/*
   Copyright (C) 2012 Ed Bueler
  
   This file is part of BETTERHYDRO.
  
   BETTERHYDRO is free software; you can redistribute it and/or modify it under the
   terms of the GNU General Public License as published by the Free Software
   Foundation; either version 2 of the License, or (at your option) any later
   version.
  
   BETTERHYDRO is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
   FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
   details.
  
   You should have received a copy of the GNU General Public License
   along with BETTERHYDRO; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

static char help[] = "Solves time-dependent subglacial hydrology model\n"
"with evolving water thickness (W) and water pressure (P) variables.\n"
"This model avoids using a function relating P and W, so P is determined\n"
"by an elliptic PDE in the W>0 region.  See phoenix.pdf for details.\n";

/* 
Example usage:

Get help:
    ./wph -fd -help | grep wph_

Finite difference evaluation of Jacobian using coloring:
    ./wph -fd
FIXME:  Jacobian evaluation routine not yet written

Options which should reproduce result of porous, the verification case:
    ./wph -fd -wph_Creep 0 -wph_Cmelt 0
    ./porous   # same computation for W

Minimal movie:
    ./wph -fd -da_grid_x 100 -da_grid_y 100 -ts_monitor_solution -draw_pause 0.5

Minimal convergence info:
    ./wph -fd -wph_steps 1 -wph_converge_check
    
More reporting:
    ./wph -fd -snes_monitor -ts_monitor -snes_vi_monitor -ksp_converged_reason
    
Matlab initial and final output:
    ./wph -fd -mfile foo.m

Relevant TS options:
  -ts_type <cn>: TS method (one of) euler beuler cn pseudo gl ssp theta alpha
  -ts_monitor_draw: <FALSE> Monitor timestep size graphically (TSMonitorLG)
  -ts_monitor_solution: <FALSE> Monitor solution graphically (TSMonitorSolution)
*/

#include <petscdmda.h>
#include <petscts.h>
#include "matlabprint.h"  /* utilities for putting petsc objects in .m file */

typedef struct {
  PetscReal W,  /* water thickness */
            P;  /* water pressure */
} WPnode;

/* we overload dof=2 Vecs to use these other node types: */

/* geometry */
typedef struct {
  PetscReal H,     /* thickness */
            b;     /* bed elevation */
} geomnode;

/* needed for water sources */
typedef struct {
  PetscReal S,      /* surface/englacial (not wall) melt; kg m-2 s-1 */
            dPsqr; /* not currently in use */
} sourcenode;


/* application context for a (W,P)-space hydrology solver */
typedef struct {
   DM          da;
   Vec         geom, sources;
   PetscReal   rhoi, rhow, g, Lfusion, Aglen, nglen,  /* standard glacier consts */
               K0,         /* constant hydraulic conductivity */
               Cmelt,      /* coefficient to turn off melt terms */
               Creep,      /* coefficient of creep; set from 1-day-to-close criterion */
               H0,         /* typical ice sheet thickness */
               dx, dy, L,  /* grid description; L=half-width of square domain */
               lastWnorm1, maxrnorm;
   PetscBool   run_silent, not_converged_warning, not_conservation_warning;
   PetscInt    fcncount,   /* count function evaluations */
               expernum;   /* various experiments can be run */
} WPCtx;

extern PetscErrorCode DefaultContext(WPCtx *user);
extern PetscErrorCode SetMarginProfile(WPCtx*,Vec,PetscReal);
extern PetscErrorCode FillSetupInitial(WPCtx*,PetscReal,PetscReal,PetscReal,Vec);
extern PetscErrorCode FormPositivityBounds(SNES,Vec,Vec);
extern PetscErrorCode checkPositivity(WPCtx *user,Vec X);
extern PetscErrorCode getPgradP(WPCtx*,Vec,Vec,Vec);
extern PetscErrorCode FormIFunction(TS,PetscReal,Vec,Vec,Vec,void*);  
extern PetscErrorCode FormIJacobian(TS,PetscReal,Vec,Vec,PetscReal,Mat*,Mat*,MatStructure*,void*);
extern PetscErrorCode getWPnorms(WPCtx*,Vec,
                                 PetscReal*,PetscReal*,PetscReal*,PetscReal*);
extern PetscErrorCode MyTSMonitor(TS,PetscInt,PetscReal,Vec,void*);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  DM             da;                   /* structured grid topology object */
  TS             ts;                   /* time-stepping object (contains snes) */
  SNES           snes;                 /* Newton solver object */
  Vec            X,residual;           /* solution, residual */
  Mat            J;                    /* Jacobian matrix */
  PetscInt       Mx,My,fsteps,steps;
  ISColoring     iscoloring;
  PetscReal      tstart,tend,ftime,secperday=3600.0*24.0,S0;
  PetscBool      fdflg = PETSC_FALSE, mfileflg = PETSC_FALSE, optflg = PETSC_FALSE;
  char           mfile[PETSC_MAX_PATH_LEN] = "out.m";
  MatFDColoring  matfdcoloring;
  WPCtx      user;                 /* user-defined work context */

  PetscInitialize(&argc,&argv,(char *)0,help);

  ierr = DMDACreate2d(PETSC_COMM_WORLD,
             DMDA_BOUNDARY_PERIODIC, DMDA_BOUNDARY_PERIODIC,
             DMDA_STENCIL_STAR, // nonlinear diffusion but diffusivity
                                //   depends on soln W not grad W
             -20,-20,           // default to 20x20 grid but override with
                                //   -da_grid_x, -da_grid_y (or -da_refine)
             PETSC_DECIDE,PETSC_DECIDE, // num of procs in each dim
             2,                 // dof = 2:  node = (W,Y)
                                //        or node = (P,dPsqr)
                                //        or node = (ddxE,ddyN)
             1,                 // s = 1 (stencil extends out one cell)
             PETSC_NULL,PETSC_NULL, // no specify proc decomposition
             &da);CHKERRQ(ierr);
  ierr = DMSetApplicationContext(da,&user);CHKERRQ(ierr);

  /* get Vecs and Mats for this grid */
  ierr = DMCreateGlobalVector(da,&X);CHKERRQ(ierr);
  ierr = VecDuplicate(X,&residual);CHKERRQ(ierr);
  ierr = VecDuplicate(X,&user.geom);CHKERRQ(ierr);
  ierr = VecDuplicate(X,&user.sources);CHKERRQ(ierr);
  ierr = DMGetMatrix(da,MATAIJ,&J);CHKERRQ(ierr);

  /* set up contexts */
  user.expernum = 1;
  tstart   = 10.0 * secperday; /* 10 days in seconds */
  tend     = 30.0 * secperday;
  steps    = 10;
  user.da  = da;
  ierr = DefaultContext(&user);CHKERRQ(ierr);
  S0       = 1.0 * user.rhow / 31556926.0;  /* = 1 m a-1 added to water layer; for exper 2 only */

  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,
           "","options to better hydrology model wph","");CHKERRQ(ierr);
  {
    ierr = PetscOptionsInt("-wph_exper",
                            "experiment 1 = (injection), 2 = (constant S)","",
                            user.expernum,&user.expernum,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsInt("-wph_steps","number of timesteps to take","",
                           steps,&steps,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-wph_exper2_S0",
                            "constant surface melt rate in experiment 2","",
                            S0,&S0,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-wph_Cmelt",
                            "additional coefficient for amount of melt","",
                            user.Cmelt,&user.Cmelt,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-wph_Creep",
                            "creep closure coefficient","",
                            user.Creep,&user.Creep,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-wph_L","half-width of square region in meters","",
                            user.L,&user.L,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-wph_tstart_days","start time in days","",
                            tstart/secperday,&tstart,&optflg);CHKERRQ(ierr);
    if (optflg) { tstart *= secperday; }
    ierr = PetscOptionsReal("-wph_tend_days","end time in days","",
                            tend/secperday,&tend,&optflg);CHKERRQ(ierr);
    if (optflg) { tend *= secperday; }
    ierr = PetscOptionsBool("-wph_converge_check",
                            "run silent and check for Newton convergence",
                            "",user.run_silent,&user.run_silent,PETSC_NULL);
                            CHKERRQ(ierr);
    ierr = PetscOptionsString("-mfile",
                            "name of Matlab file to write results","",
                            mfile,mfile,PETSC_MAX_PATH_LEN,&mfileflg);
                            CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);

  /* fix grid parameters */
  ierr = DMDASetUniformCoordinates(da,  // square domain; handles periodic appropriately
              -user.L, user.L, -user.L, user.L, 0.0, 1.0);CHKERRQ(ierr);
  ierr = DMDAGetInfo(da,PETSC_IGNORE,&Mx,&My,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);
  user.dx = 2.0 * user.L / (PetscReal)Mx;
  user.dy = 2.0 * user.L / (PetscReal)My;

  /* setup geometry (H and b), sources (S), and initial state var (W) */
  ierr = FillSetupInitial(&user,tstart,S0,X);CHKERRQ(ierr);

  /* setup TS = timestepping object */
  ierr = TSCreate(PETSC_COMM_WORLD,&ts);CHKERRQ(ierr);
  ierr = TSSetType(ts,TSCN);CHKERRQ(ierr);
  ierr = TSSetRHSFunction(ts,residual,RHSFunction,&user);CHKERRQ(ierr);

  /* use coloring to compute rhs Jacobian efficiently */
  ierr = PetscOptionsGetBool(PETSC_NULL,"-fd",&fdflg,PETSC_NULL);CHKERRQ(ierr);
  if (fdflg){
    ierr = DMGetColoring(da,IS_COLORING_GLOBAL,MATAIJ,&iscoloring);CHKERRQ(ierr);
    ierr = MatFDColoringCreate(J,iscoloring,&matfdcoloring);CHKERRQ(ierr);
    ierr = MatFDColoringSetFromOptions(matfdcoloring);CHKERRQ(ierr);
    ierr = ISColoringDestroy(&iscoloring);CHKERRQ(ierr);
    ierr = MatFDColoringSetFunction(matfdcoloring,
             (PetscErrorCode (*)(void))RHSFunction,&user);CHKERRQ(ierr);
    ierr = TSSetRHSJacobian(ts,J,J,TSDefaultComputeJacobianColor,
             matfdcoloring);CHKERRQ(ierr);
  } else { /* default case */
    ierr = TSSetRHSJacobian(ts,J,J,RHSJacobian,&user);CHKERRQ(ierr);
  }

  /* set up times for time-stepping */
  ierr = TSSetInitialTimeStep(ts,tstart,
           (tend - tstart) / (PetscReal)steps);CHKERRQ(ierr);
  ierr = TSSetDuration(ts,steps,tend);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts,PETSC_TRUE);CHKERRQ(ierr);
  ierr = TSMonitorSet(ts,MyTSMonitor,&user,PETSC_NULL);CHKERRQ(ierr);

  /* Set SNESVI type and supply upper and lower bounds. */
  ierr = TSGetSNES(ts,&snes);CHKERRQ(ierr);
  ierr = SNESSetApplicationContext(snes,&user);CHKERRQ(ierr);
  ierr = SNESVISetComputeVariableBounds(snes,FormPositivityBounds);
        CHKERRQ(ierr);

  /* ask user to finalize settings */
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);

  /* report on setup */
  if (!user.run_silent) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
      "setup done: square       side length = %.3f km\n"
      "            grid               Mx,My = %d,%d\n"
      "            spacing            dx,dy = %.3f,%.3f m\n"
      "            times     tstart:dt:tend = %.3f:%.3f:%.3f days\n",
      2.0 * user.L / 1000.0, Mx, My, user.dx, user.dy,
      tstart / secperday, (tend-tstart)/(steps*secperday), tend / secperday);
      CHKERRQ(ierr);
  }
  if (mfileflg) {
    if (!user.run_silent) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,
        "writing initial fields W,Y,P,etc. to Matlab/Octave file %s ...\n",mfile);CHKERRQ(ierr);
    }
    ierr = print2vecmatlab(da,X,"W_init","P_init",mfile,PETSC_FALSE);CHKERRQ(ierr);
    ierr = print2vecmatlab(da,user.geom,"H","b",mfile,PETSC_TRUE);CHKERRQ(ierr);
  }

  /* run time-stepping with implicit steps  */
  ierr = TSSolve(ts,X,&ftime);CHKERRQ(ierr);

  /* make a report on run and final state */
  ierr = TSGetTimeStepNumber(ts,&fsteps);CHKERRQ(ierr);
  if ((!user.run_silent) && (ftime != tend)) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
    "***WARNING3***:  reported final time wrong:  ftime(=%.12e) != tend(=%.12e) (days)\n",
    ftime / secperday, tend / secperday);CHKERRQ(ierr); }
  if ((!user.run_silent) && (fsteps != steps)) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
    "***WARNING4***:  reported number of steps wrong:  fsteps(=%D) != steps(=%D)\n",
    fsteps, steps);CHKERRQ(ierr); }

  if (mfileflg) {
    if (!user.run_silent) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,
        "writing final fields to %s ...\n",mfile);CHKERRQ(ierr);
    }
    ierr = getPgradP(&user,X,user.PdP,PETSC_NULL);CHKERRQ(ierr);
    ierr = print2vecmatlab(da,X,"W_final","Y_final",mfile,PETSC_TRUE);CHKERRQ(ierr);
    ierr = print2vecmatlab(da,user.PdP,"P_final","dPsqr_final",mfile,PETSC_TRUE);CHKERRQ(ierr);
  }

  if (user.run_silent) {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "%6d  %6d  %9.3f  %.12e\n",
                       Mx, My, (tend-tstart)/secperday, user.maxrnorm);CHKERRQ(ierr);
  }

  /* Free work space.  */
  ierr = MatDestroy(&J);CHKERRQ(ierr);
  if (fdflg) { ierr = MatFDColoringDestroy(&matfdcoloring);CHKERRQ(ierr); }
  ierr = VecDestroy(&X);CHKERRQ(ierr);
  ierr = VecDestroy(&user.geom);CHKERRQ(ierr);
  ierr = VecDestroy(&user.sources);CHKERRQ(ierr);
  ierr = VecDestroy(&residual);CHKERRQ(ierr);
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = DMDestroy(&da);CHKERRQ(ierr);

  ierr = PetscFinalize();CHKERRQ(ierr);

  PetscFunctionReturn((PetscInt)(user.not_converged_warning
                                 || user.not_conservation_warning));
}


#undef __FUNCT__
#define __FUNCT__ "DefaultContext"
PetscErrorCode DefaultContext(WPCtx *user) {
  PetscReal  Kmax  = 0.01,   /* m s-1;  maximum of hydraulic conductivity; F&C 2002 */
             Kmin  = 0.0001; /* m s-1;  minimum of hydraulic conductivity; F&C 2002 */
  /* physical constants */
  user->rhoi  = 910.0;       /* kg m-3; ice density; also PISM default */
  user->rhow  = 1000.0;      /* kg m-3; freshwater density; also PISM default */
  user->g     = 9.81;        /* m s-2;  acceleration of gravity */
  user->Lfusion = 3.34e5;    /* J kg-1; latent heat of fusion (Greve&Blatter) */
  user->Aglen = 3.1689e-24;  /* Pa-3 s-1; ice softness (EISMINT96) */
  user->nglen = 3;           /* pure; Glen exponent (EISMINT96) */
  /* parameters for this model */
  user->K0 = sqrt(Kmax * Kmin); /* geometric average */
  user->Cmelt = 1.0;         /* pure;   additional coefficient for amount of melt */
  user->Creep = 1.9014e-04;  /* pure;   creep closure coefficient; this value gives
                                1/e time of one day from pure creep closure of
                                cavity:
                                  >> AN3 = 3.1689e-24 * (910 * 9.81 * 3000)^3
                                  AN3 =  0.060870
                                  >> Creep = 1/(3600*24*AN3)
                                  Creep =  1.9014e-04 */
  user->H0     = 3000.0;     /* m;      typical ice sheet thickness */
  /* domain and grid */
  user->L  = 2.0e4;          /* m */
  user->dx = -1.0;
  user->dy = -1.0;
  /* utility */
  user->lastWnorm1 = -1.0;
  user->maxrnorm = -1.0;
  user->run_silent = PETSC_FALSE;
  user->not_converged_warning = PETSC_FALSE;
  user->not_conservation_warning = PETSC_FALSE;
  user->fcncount = 0;
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "SetMarginProfile"
/*  SetMarginProfile(): compute H(x,y) = H(x) = 'bueler profile' 
    see subsection 5.6.3 in Greve & Blatter */
PetscErrorCode  SetMarginProfile(WPCtx *user, Vec geom, PetscReal H0) {
  PetscErrorCode ierr;
  PetscInt       i,j,xs,ys,xm,ym;
  PetscReal      s, phi, profileL = user->L - 2000.0;
  DM             da = user->da, coordDA;
  Vec            coordinates;
  DMDACoor2d     **coords;
  geomnode       **hb;
  PetscFunctionBegin;

  ierr = DMDAGetCorners(da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL);CHKERRQ(ierr);

  ierr = DMDAGetCoordinateDA(da, &coordDA);CHKERRQ(ierr);
  ierr = DMDAGetCoordinates(da, &coordinates);CHKERRQ(ierr);

  ierr = DMDAVecGetArray(coordDA, coordinates, &coords);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da, geom, &hb);CHKERRQ(ierr);
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      s = PetscAbs(coords[j][i].x / profileL);  /* scaled coordinate */
      if (s <= 1.0) {
        phi = 4.0 * s - 3.0 * pow(s,4.0/3.0) + 3.0 * pow(1.0 - s,4.0/3.0) - 1.0;
        hb[j][i].H = (H0/pow(2.0,3.0/8.0)) * pow(phi,3.0/8.0);
      } else {
        hb[j][i].H = 0.0;
      }
    }
  }
  ierr = DMDAVecRestoreArray(da, geom, &hb);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(coordDA, coordinates, &coords);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "FillSetupInitial"
/*  FillSetupInitial(): fill Vecs in user context, and initial state of X */
PetscErrorCode FillSetupInitial(WPCtx *user, PetscReal tstart, PetscReal Y0, PetscReal S0,
                                Vec X) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  /* settings common to all experiments */
  ierr = VecStrideSet(user->geom,1,0.0);CHKERRQ(ierr); /* b(x,y) = 0  */
  ierr = VecStrideSet(X,1,Y0);CHKERRQ(ierr);           /* Y(x,y) = Y0 from option */
  ierr = VecStrideSet(user->S,1,-999.0);CHKERRQ(ierr);  /* unused = -999  */

  switch (user->expernum) {
    case 1:
      /* experiment 1 is injection of water;  W = barenblatt at tstart */
      ierr = BarenblattState(user,tstart,X);CHKERRQ(ierr);
      /* H(x,y) = H0; no overburden gradient */
      ierr = VecStrideSet(user->geom,0,user->H0);CHKERRQ(ierr);  
      /* S(x,y) = 0; no surface water */
      ierr = VecStrideSet(user->S,0,0.0);CHKERRQ(ierr);
      break;
    case 2:
      /* experiment 2 is steady melt under ice sheet margin;  W = 0 initially */
      ierr = VecStrideSet(X,0,0.0);CHKERRQ(ierr);  
      /* H(x,y) = bueler profile */
      ierr = SetMarginProfile(user,user->geom,user->H0);CHKERRQ(ierr);
      /* S(x,y) = S0 from option */
      ierr = VecStrideSet(user->S,0,S0);CHKERRQ(ierr);      
      break;
    default:
      SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
          "illegal experiment number user.expernum=%d\n",user->expernum);
  }
  
  /* finalize P and grad P from X */
  ierr = getPgradP(user,X,user->PdP,PETSC_NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "FormPositivityBounds"
/*  FormPositivityBounds() for call-back: tell SNESVI (variational inequality)
  that we want positive solutions for both W and Y */
PetscErrorCode FormPositivityBounds(SNES snes, Vec Xl, Vec Xu) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = VecStrideSet(Xl,0,0.0);CHKERRQ(ierr);  /* W >= 0 */
  ierr = VecStrideSet(Xl,1,0.0);CHKERRQ(ierr);  /* Y >= 0 */
  ierr = VecStrideSet(Xu,0,SNES_VI_INF);CHKERRQ(ierr);
  ierr = VecStrideSet(Xu,1,SNES_VI_INF);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "checkPositivity"
PetscErrorCode checkPositivity(WPCtx *user,Vec X) {
  PetscErrorCode ierr;
  PetscInt       i,j,xs,ys,xm,ym;
  DM             da = user->da;
  WYnode         **wy;
  PetscFunctionBegin;
  ierr = DMDAGetCorners(da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,X,&wy);CHKERRQ(ierr);
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      if (wy[j][i].W < 0) {
        SETERRQ3(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
          "negative water thickness W at (i,j)=(%d,%d) during function eval %d",
          i,j,user->fcncount);
      }
      if (wy[j][i].Y < 0) {
        SETERRQ3(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
          "negative capacity thickness Y at (i,j)=(%d,%d) during function eval %d",
          i,j,user->fcncount);
      }
    }
  }
  ierr = DMDAVecRestoreArray(da,X,&wy);CHKERRQ(ierr);
  PetscFunctionReturn(0);
} 


#undef __FUNCT__
#define __FUNCT__ "getdPsqr"
/* compute    |grad P|^2 from P, using centered finite differences */
PetscErrorCode getdPsqr(WPCtx *user, Vec P, Vec dPsqr) {
  PetscErrorCode ierr;
  PetscInt       i,j,xs,ys,xm,ym,Mx,My;
  PetscReal      Px,Py,dx=user->dx,dy=user->dy;
  DM             da = user->da;
  Vec            localP;
  WPnode         **wp;
  geomnode       **hb;
  sourcenode     **reg;

  PetscFunctionBegin;
  ierr = DMDAGetInfo(da,PETSC_IGNORE,&Mx,&My,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);
  ierr = DMDAGetCorners(da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL);CHKERRQ(ierr);

  /* we will difference P so we need a local, ghosted, copy */
  ierr = DMGetLocalVector(da,&localP);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(da,P,INSERT_VALUES,localP);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da,P,INSERT_VALUES,localP);CHKERRQ(ierr);

  ierr = DMDAVecGetArray(da,localP,&wp);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,dPsqr,&reg);CHKERRQ(ierr);
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      Px = (wp[j][i+1].P - wp[j][i-1].P) / (2.0 * dx);
      Py = (wp[j+1][i].P - wp[j-1][i].P) / (2.0 * dy);
      reg[j][i].dPsqr = Px * Px + Py * Py;
    }
  }
  ierr = DMDAVecRestoreArray(da,dPsqr,&reg);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da,localP,&wp);CHKERRQ(ierr);

  ierr = DMRestoreLocalVector(da,&localP);CHKERRQ(ierr);
  PetscFunctionReturn(0);
} 


/*
we solve the b=0 and v_b=0 case of the phoenix.pdf hydrology model:

FIXME

W_t = c1 div (W grad P) + Cmelt c2 W |grad P|^2 + S / rhow

Y_t = Cmelt c3 W |grad P|^2 - Creep A max{0,rhoi g H - P}^n Y

where

P = rhoi g H (W/(Y+Ymin))^sigma

and where

  c1 = K / (rhow g)
  c2 = c1 / (rhow L)
  c3 = c1 / (rhoi L)
*/



/*  FormIFunction:  implicit DAE is  F(t,X,Xdot) = 0 */
#undef __FUNCT__
#define __FUNCT__ "FormIFunction"
PetscErrorCode FormIFunction(TS ts, PetscReal t, Vec X, Vec Xdot, Vec F, void *ctx) {
  PetscErrorCode ierr;
  WPCtx      *user = (WPCtx*)ctx;
  DM             da = user->da;
  Vec            dPstag,localX;
  PetscInt       i,j,Mx,My,xs,ys,xm,ym;
  PetscReal      c1,c2,c3;
  PetscReal      dx = user->dx, dy = user->dy;
  PetscReal      Wij, Weast, Wwest, Wnorth, Wsouth, dQx, dQy, Nn;
  WYnode         **wy, **f;
  geomnode         **hb;
  Snode          **ss;
  Pnode          **reg;
  staggradnode   **stag;

  PetscFunctionBegin;
  ierr = DMDAGetInfo(da,PETSC_IGNORE,&Mx,&My,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);
  ierr = DMDAGetCorners(da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL);CHKERRQ(ierr);

  user->fcncount = user->fcncount + 1;
  ierr = checkPositivity(user,X);CHKERRQ(ierr);

  c1 = user->Kconst / (user->rhow * user->g);
  c2 = c1 / (user->rhow * user->Lfusion);
  c3 = c1 / (user->rhoi * user->Lfusion);

  ierr = DMGetLocalVector(da,&dPstag);CHKERRQ(ierr); /* space for staggered grad */
  ierr = getPgradP(user,X,user->PdP,dPstag);CHKERRQ(ierr); /* compute P and grad P */
  /* we will be differencing dPstag = (ddxE,ddyN) */
  ierr = DMDALocalToLocalBegin(da,dPstag,INSERT_VALUES,dPstag);CHKERRQ(ierr);
  ierr = DMDALocalToLocalEnd(da,dPstag,INSERT_VALUES,dPstag);CHKERRQ(ierr);

  ierr = DMGetLocalVector(da,&localX);CHKERRQ(ierr);
  /* we will be differencing X = (W,Y) */
  ierr = DMGlobalToLocalBegin(da,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  ierr = DMDAVecGetArray(da,localX,&wy);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,F,&f);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,user->geom,&hb);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,user->PdP,&reg);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,user->S,&ss);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,dPstag,&stag);CHKERRQ(ierr);
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
        Wij     = wy[j][i].W;
        Weast   = 0.5 * (wy[j][i+1].W + Wij);
        Wwest   = 0.5 * (wy[j][i-1].W + Wij);
        Wnorth  = 0.5 * (wy[j+1][i].W + Wij);
        Wsouth  = 0.5 * (wy[j-1][i].W + Wij);

        /* W_t = c1 div (W grad P) + Cmelt c2 W |grad P|^2 + S / rhow*/
        dQx = Weast  * stag[j][i].ddxE - Wwest  * stag[j][i-1].ddxE;
        dQy = Wnorth * stag[j][i].ddyN - Wsouth * stag[j-1][i].ddyN;
        f[j][i].W = (c1/dx) * dQx + (c1/dy) * dQy
                      + user->Cmelt * c2 * Wij * reg[j][i].dPsqr
                      + ss[j][i].S / user->rhow;

        /* Y_t = Cmelt c3 W |grad P|^2 - Creep A max{0,p_i - P}^n Y */
        Nn = PetscMax(0.0,user->rhoi * user->g * hb[j][i].H - reg[j][i].P);
        Nn = pow(Nn,user->nglen);
        f[j][i].Y = user->Cmelt * c3 * Wij * reg[j][i].dPsqr
                      - user->Creep * user->Aglen * Nn * wy[j][i].Y;
    }
  }
  ierr = DMDAVecRestoreArray(da,localX,&wy);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da,F,&f);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da,user->geom,&hb);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da,user->PdP,&reg);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da,user->S,&ss);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da,dPstag,&stag);CHKERRQ(ierr);

  ierr = DMRestoreLocalVector(da,&localX);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(da,&dPstag);CHKERRQ(ierr);
  PetscFunctionReturn(0); 
}


#undef __FUNCT__
#define __FUNCT__ "FormIJacobian"
PetscErrorCode FormIJacobian(TS ts, PetscReal t, Vec X, Vec Xdot, PetscReal a, Mat *J, Mat *Jpre, MatStructure *str, void *ctx) {
  SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER,"this Jacobian is not yet implemented");
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "getWYPnorms"
PetscErrorCode getWYPnorms(WPCtx *user,Vec X,
                           PetscReal *Wav, PetscReal *Wnorminf,
                           PetscReal *Yav, PetscReal *Ynorminf,
                           PetscReal *Pav, PetscReal *Pnorminf) {
  PetscErrorCode ierr;
  PetscInt       Mx, My;
  PetscReal      onenorms[2],infnorms[2];
  PetscFunctionBegin;  
  ierr = VecStrideNormAll(X,NORM_1,onenorms);CHKERRQ(ierr);
  ierr = VecStrideNormAll(X,NORM_INFINITY,infnorms);CHKERRQ(ierr);
  ierr = DMDAGetInfo(user->da,PETSC_IGNORE,&Mx,&My,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);
  *Wav      = onenorms[0] / ((PetscReal)Mx * (PetscReal)My);
  *Yav      = onenorms[1] / ((PetscReal)Mx * (PetscReal)My);
  *Wnorminf = infnorms[0];
  *Ynorminf = infnorms[1];
  ierr = VecStrideNormAll(user->PdP,NORM_1,onenorms);CHKERRQ(ierr);
  ierr = VecStrideNormAll(user->PdP,NORM_INFINITY,infnorms);CHKERRQ(ierr);
  *Pav      = onenorms[0] / ((PetscReal)Mx * (PetscReal)My);
  *Pnorminf = infnorms[0];
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "MyTSMonitor"
/* MyTSMonitor:  stdout report at every time step */
PetscErrorCode MyTSMonitor(TS ts,PetscInt step,PetscReal ptime,Vec X,void *ptr)
{
  PetscErrorCode ierr;
  PetscReal      twonorms[2],rnorm2,Wav,Wnorminf,Yav,Ynorminf,Pav,Pnorminf;
  PetscReal      secperday=3600.0*24.0,CD;
  MPI_Comm       comm;
  WPCtx      *user = (WPCtx*)ptr;
  SNES           snes;

  PetscFunctionBegin;  
  ierr = getWYPnorms(user,X,&Wav,&Wnorminf,&Yav,&Ynorminf,&Pav,&Pnorminf);CHKERRQ(ierr);
  CD = (user->Kconst * user->sigma) / (user->rhow * user->g);
  ierr = VecStrideNormAll(X,NORM_2,twonorms);CHKERRQ(ierr);
  ierr = PetscObjectGetComm((PetscObject)ts,&comm);CHKERRQ(ierr);
  /* summary */
  if (!user->run_silent) {
    if (user->fcncount <= 2) {
      ierr = PetscPrintf(comm,
           "  step  time(days)       av W(m)    max W(m)  max D(m2 s-1)      av Y(m)     max Y(m)    max P(bar)\n");CHKERRQ(ierr);
    }
    ierr = PetscPrintf(comm,
           "   %3d     %7G  %12.9f %11.6f    %11.6f  %11.6f  %11.6f   %11.6f\n",
           step,ptime/secperday,Wav,Wnorminf,CD*Pnorminf,
           Yav,Ynorminf,Pnorminf/1.0e5);CHKERRQ(ierr);
  }
  /* warning if solution not small */
  ierr = TSGetSNES(ts,&snes);CHKERRQ(ierr);
  ierr = SNESGetFunctionNorm(snes,&rnorm2);CHKERRQ(ierr);
  if (rnorm2 > 1.0e-10 * (twonorms[0]+twonorms[1])) {
    user->not_converged_warning = PETSC_TRUE;
    if (!user->run_silent) {
      ierr = PetscPrintf(comm,
           "***WARNING1***: residual norm not small (> 1e-10 * (|W|_2+|Y|_2)) at step %d\n",
           step);CHKERRQ(ierr); }
  }

  /* update max of rnorm (relative) so far */
  if (twonorms[0] > 0.0) user->maxrnorm = PetscMax(user->maxrnorm,rnorm2);

  PetscFunctionReturn(0);
}

