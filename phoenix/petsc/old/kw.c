/*
   Copyright (C) 2012 Ed Bueler
  
   This file is part of BETTERHYDRO.
  
   BETTERHYDRO is free software; you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation; either version 2 of the License, or (at your option) any later
   version.
  
   BETTERHYDRO is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
   FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
   details.
  
   You should have received a copy of the GNU General Public License
   along with BETTERHYDRO; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

static char help[] = "Solves time-dependent subglacial hydrology model\n"
"with evolving water thickness (W) and hydraulic conductivity (K) variables.\n"
"See vanPeltBueler.pdf for more details\n";

/* 
Example usage:

Get help:
    ./kw -fd -help | grep kw_

Finite difference evaluation of Jacobian using coloring:
    ./kw -fd
FIXME:  Jacobian evaluation routine not yet written

Options which should reproduce result of porous, the verification case:
    ./kw -fd -kw_alpha 0 -kw_Cmelt 0
    ./porous   # same computation, but for W only, without evolving K

Minimal movie:
    ./kw -fd -da_grid_x 101 -da_grid_y 101 -ts_monitor_solution -draw_pause 0.5

Minimal convergence info:
    ./kw -fd -kw_steps 1 -kw_converge_check

More reporting:
    ./kw -fd -snes_monitor -snes_vi_monitor -ksp_converged_reason

Matlab initial and final output:
    ./kw -fd -mfile foo.m

Relevant TS options:
  -ts_type <cn>: TS method (one of) euler beuler cn pseudo gl ssp theta alpha
  -ts_monitor_draw: <FALSE> Monitor timestep size graphically (TSMonitorLG)
  -ts_monitor_solution: <FALSE> Monitor solution graphically (TSMonitorSolution)
*/

#include <petscdmda.h>
#include <petscts.h>

#include "matlabprint.h"  /* utilities for putting petsc objects in .m file */
#include "context.h"      /* structs and default values for parameters */

extern PetscErrorCode InitialState(DM,PorousCtx*,PetscReal,PetscReal,Vec);
extern PetscErrorCode FormBounds(SNES,Vec,Vec);
extern PetscErrorCode checkBounds(PorousCtx *user,Vec X);
extern PetscErrorCode RHSFunction(TS,PetscReal,Vec,Vec,void*);
extern PetscErrorCode RHSJacobian(TS,PetscReal,Vec,Mat*,Mat*,
                                  MatStructure*,void*);
extern PetscErrorCode getWKnorms(PorousCtx*,Vec,
                                 PetscReal*,PetscReal*,PetscReal*,PetscReal*);
extern PetscErrorCode maxDiffusivity(DM,PorousCtx*,Vec,Vec,PetscReal*);
extern PetscErrorCode MyTSMonitor(TS,PetscInt,PetscReal,Vec,void*);


typedef struct {
  PetscReal W,  /* water thickness */
            K;  /* hydraulic conductivity */
} WKnode;

/* we overload dof=2 Vecs to use an additional node type: */
typedef struct {
  PetscReal H,     /* thickness */
            b;     /* bed elevation */
} Hbnode;


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
  PetscReal      tstart,tend,ftime,secperday=3600.0*24.0;
  PetscBool      fdflg = PETSC_FALSE, mfileflg = PETSC_FALSE, optflg = PETSC_FALSE;
  char           mfile[PETSC_MAX_PATH_LEN] = "out.m";
  MatFDColoring  matfdcoloring;
  PorousCtx      user;                 /* user-defined work context */

  PetscInitialize(&argc,&argv,(char *)0,help);

  ierr = DMDACreate2d(PETSC_COMM_WORLD,
             DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_NONE, // correct for zero Dirichlet
             DMDA_STENCIL_STAR, // nonlinear diffusion but diffusivity
                                //   depends on soln W not grad W
             -21,-21,           // default to 20x20 grid but override with
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
  ierr = DMGetMatrix(da,MATAIJ,&J);CHKERRQ(ierr);

  /* set up contexts */
  tstart   = 10.0 * secperday; /* 10 days in seconds */
  tend     = 30.0 * secperday;
  steps    = 10;  /* same default as porous */
  user.da  = da;
  ierr = DefaultContext(&user);CHKERRQ(ierr);

  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,
           "","options to (W,K)-space better-hydrology-for-PISM model kw","");
           CHKERRQ(ierr);
  {
    ierr = PetscOptionsReal("-kw_sigma","nonlinear power","",
                            user.sigma,&user.sigma,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-kw_Wcrit",
                            "critial thickness; P, dK/dt are functions of W/Wcrit","",
                            user.Wcrit,&user.Wcrit,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-kw_alpha",
                            "coefficient in dK/dt equation","",
                            user.alpha,&user.alpha,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-kw_Cmelt",
                            "coefficient (pure) for fraction of modeled melt","",
                            user.Cmelt,&user.Cmelt,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-kw_L","half-width of square region in meters","",
                            user.L,&user.L,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-kw_tstart_days","start time in days","",
                            tstart/secperday,&tstart,&optflg);CHKERRQ(ierr);
    if (optflg) { tstart *= secperday; }
    ierr = PetscOptionsReal("-kw_tend_days","end time in days","",
                            tend/secperday,&tend,&optflg);CHKERRQ(ierr);
    if (optflg) { tend *= secperday; }
    ierr = PetscOptionsInt("-kw_steps","number of timesteps to take","",
                           steps,&steps,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsBool("-kw_converge_check",
                            "run silent and check for convergence",
                            "",user.run_silent,&user.run_silent,PETSC_NULL);
                            CHKERRQ(ierr);
    ierr = PetscOptionsString("-mfile",
                            "name of Matlab file to write results","",
                            mfile,mfile,PETSC_MAX_PATH_LEN,&mfileflg);
                            CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);

  /* fix remaining parameters */
  ierr = DerivedConstants(&user);CHKERRQ(ierr);
  ierr = VecStrideSet(user.geom,0,user.H0);CHKERRQ(ierr);  /* H(x,y) = H0 */
  ierr = VecStrideSet(user.geom,1,0.0);CHKERRQ(ierr);      /* b(x,y) = 0  */
  ierr = DMDASetUniformCoordinates(da,  // square domain
              -user.L, user.L, -user.L, user.L, 0.0, 1.0);CHKERRQ(ierr);
  ierr = DMDAGetInfo(da,PETSC_IGNORE,&Mx,&My,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);
  user.dx = 2.0 * user.L / (Mx-1);
  user.dy = 2.0 * user.L / (My-1);

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

  /* set initial state:  W = barenblatt @ tstart, K = Kconst */
/* FIXME */  ierr = InitialState(da,&user,tstart,user.Kconst,X);CHKERRQ(ierr);

  /* set up times for time-stepping */
  ierr = TSSetInitialTimeStep(ts,tstart,
           (tend - tstart) / (PetscReal)steps);CHKERRQ(ierr);
  ierr = TSSetDuration(ts,steps,tend);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts,PETSC_TRUE);CHKERRQ(ierr);
  ierr = TSMonitorSet(ts,MyTSMonitor,&user,PETSC_NULL);CHKERRQ(ierr);

  /* Set SNESVI type and supply upper and lower bounds. */
  ierr = TSGetSNES(ts,&snes);CHKERRQ(ierr);
  ierr = SNESSetApplicationContext(snes,(void*)(&user));CHKERRQ(ierr);
  ierr = SNESVISetComputeVariableBounds(snes,FormBounds);CHKERRQ(ierr);

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
        "writing initial W,K and geometry H,b to Matlab file %s ...\n",
        mfile);CHKERRQ(ierr);
    }
    ierr = print2vecmatlab(da,X,"W_init","K_init",mfile,PETSC_FALSE);CHKERRQ(ierr);
    ierr = print2vecmatlab(da,user.geom,"H","b",mfile,PETSC_TRUE);CHKERRQ(ierr);
  }

  /* run time-stepping with implicit steps  */
  ierr = TSSolve(ts,X,&ftime);CHKERRQ(ierr);

  /* make a report on run and final state */
  ierr = TSGetTimeStepNumber(ts,&fsteps);CHKERRQ(ierr);
  if ((!user.run_silent) && (ftime != tend)) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
    "***WARNING3***:  reported final time wrong:  %.12e != %.12e = tend (days)\n",
    ftime / secperday, tend / secperday);CHKERRQ(ierr); }
  if ((!user.run_silent) && (fsteps != steps)) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
    "***WARNING4***:  reported number of steps wrong:  %D != steps = %D\n",
    fsteps, steps);CHKERRQ(ierr); }

  if (mfileflg) {
    if (!user.run_silent) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,
        "writing final fields to %s ...\n",mfile);CHKERRQ(ierr);
    }
    ierr = print2vecmatlab(da,X,"W_final","K_final",mfile,PETSC_TRUE);CHKERRQ(ierr);
    ierr = printfigurematlab(da,2,"W_init","W_final",mfile,PETSC_TRUE);CHKERRQ(ierr);
    ierr = printfigurematlab(da,3,"K_init","K_final",mfile,PETSC_TRUE);CHKERRQ(ierr);
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
  ierr = VecDestroy(&residual);CHKERRQ(ierr);
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = DMDestroy(&da);CHKERRQ(ierr);

  ierr = PetscFinalize();CHKERRQ(ierr);

  PetscFunctionReturn((PetscInt)(user.not_converged_warning));
}


#undef __FUNCT__
#define __FUNCT__ "InitialState"
/* InitialState:  sets initial state of (W,K):  W = barenblatt at t, K = K0 */
PetscErrorCode InitialState(DM da, PorousCtx* user, PetscReal t, PetscReal K0,
                            Vec X) {
  PetscErrorCode ierr;
  ierr = BarenblattState(user,t,X);CHKERRQ(ierr);  /* sets W */
  ierr = VecStrideSet(X,1,K0);CHKERRQ(ierr);          /* sets K = K0 */
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "FormBounds"
/*  FormBounds() answers call-back: tell SNESVI (variational inequality)
  that we want positive solutions for W and lower and upper bnds on K */
PetscErrorCode FormBounds(SNES snes, Vec Xl, Vec Xu) {
  PetscErrorCode ierr;
  PorousCtx      *user;
  ierr = SNESGetApplicationContext(snes,&user);CHKERRQ(ierr);
  ierr = VecStrideSet(Xl,0,0.0);CHKERRQ(ierr);         /* W >= 0 */
  ierr = VecStrideSet(Xl,1,user->Kmin);CHKERRQ(ierr);  /* K >= Kmin */
  ierr = VecStrideSet(Xu,0,SNES_VI_INF);CHKERRQ(ierr);
  ierr = VecStrideSet(Xu,1,user->Kmax);CHKERRQ(ierr);  /* K <= Kmax */
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "checkBounds"
PetscErrorCode checkBounds(PorousCtx *user,Vec X) {
  PetscErrorCode ierr;
  PetscInt       i,j,xs,ys,xm,ym;
  DM             da = user->da;
  WKnode         **wk;
  PetscFunctionBegin;
  ierr = DMDAGetCorners(da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,X,&wk);CHKERRQ(ierr);
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      if (wk[j][i].W < 0) {
        SETERRQ3(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
          "negative water thickness W at (i,j)=(%d,%d) during function eval %d",
          i,j,user->fcncount);
      }
      if (wk[j][i].K < user->Kmin) {
        SETERRQ3(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
          "hydraulic conductivity K at (i,j)=(%d,%d) below Kmin (fcn eval %d)",
          i,j,user->fcncount);
      }
      if (wk[j][i].K > user->Kmax) {
        SETERRQ3(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
          "hydraulic conductivity K at (i,j)=(%d,%d) above Kmax (fcn eval %d)",
          i,j,user->fcncount);
      }
    }
  }
  ierr = DMDAVecRestoreArray(da,X,&wk);CHKERRQ(ierr);
  PetscFunctionReturn(0);
} 


/*
we solve the b=0, m=0, S=0 case of the hydrology model,
in its (W,K) formulation:

  W_t = c1 div (K W grad (H W^sigma)) + Cmelt c2 K W |grad(H W^sigma)|^2

  K_t = alpha (2 Kmax - K) ((W/Wcrit)^2 - 1)

where 
  c1 = rhoi / (rhow Wcrit^sigma)
  c2 = (rhoi/rhow)^2 * g / (L Wcrit^{2 sigma})
and rhoi, rhow, Wcrit, sigma, alpha, Kmax, L are positive constants,
and Cmelt is to allow turning-off the melt term
*/

#undef __FUNCT__
#define __FUNCT__ "RHSFunction"
/* RHSFunction:  evaluates nonlinear function in ODE form X' = F(X,t).
   But our case is autonomous, so F = F(X) = F(W,K) and F has no t dependence. */
PetscErrorCode RHSFunction(TS ts,PetscReal t_unused,Vec X,Vec F,void *ptr)
{
  PetscErrorCode ierr;
  PorousCtx      *user = (PorousCtx*)ptr;
  DM             da = user->da;
  Vec            localX, localgeom;
  PetscInt       i,j,Mx,My,xs,ys,xm,ym;
  PetscReal      dx = user->dx, dy = user->dy,
                 c1, c2, sig, Wcrit,
                 Wij, Wrel, Kij,
                 KWij, KWeast, KWwest, KWnorth, KWsouth,
                 HWsig, dHWsigeast, dHWsigwest, dHWsignorth, dHWsigsouth,
                 dQx, dQy, dHWsigdx, dHWsigdy, gradsqr, Cmelt;
  WKnode         **wk, **f;
  Hbnode         **hb;

  PetscFunctionBegin;
  ierr = DMDAGetInfo(da,PETSC_IGNORE,&Mx,&My,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);
  ierr = DMDAGetCorners(da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL);CHKERRQ(ierr);

  user->fcncount = user->fcncount + 1;
  ierr = checkBounds(user,X);CHKERRQ(ierr);

  sig = user->sigma;
  Wcrit = user->Wcrit;
  Cmelt = user->Cmelt;
  c1 = user->rhoi / (user->rhow * pow(Wcrit,sig));
  c2 = user->rhoi / user->rhow;
  c2 = (c2 * c2) * user->g / (user->Lfusion * pow(Wcrit,2.0*sig));

  /* we will be differencing X = (W,K) */
  ierr = DMGetLocalVector(da,&localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(da,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  /* we will be differencing geom = (H,b) */
  ierr = DMGetLocalVector(da,&localgeom);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(da,user->geom,INSERT_VALUES,localgeom);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da,user->geom,INSERT_VALUES,localgeom);CHKERRQ(ierr);

  ierr = DMDAVecGetArray(da,localX,&wk);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,localgeom,&hb);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,F,&f);CHKERRQ(ierr);
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      if (i == 0 || j == 0 || i == Mx-1 || j == My-1) {
        f[j][i].W = 0.0;  /* no change at Dirichlet boundary conditions */
        f[j][i].K = 0.0;  /* same */
      } else {
        Wij     = wk[j][i].W;
        Kij     = wk[j][i].K;

        /* W_t = c1 div (K W grad(H W^sigma)) + c2 K W |grad(H W^sigma)|^2 */
        KWij    = Kij * Wij;
        KWeast  = 0.5 * (wk[j][i+1].K * wk[j][i+1].W + KWij);
        KWwest  = 0.5 * (wk[j][i-1].K * wk[j][i-1].W + KWij);
        KWnorth = 0.5 * (wk[j+1][i].K * wk[j+1][i].W + KWij);
        KWsouth = 0.5 * (wk[j-1][i].K * wk[j-1][i].W + KWij);
        HWsig   = hb[j][i].H * pow(Wij, sig);
        dHWsigeast  = hb[j][i+1].H * pow(wk[j][i+1].W, sig) - HWsig;
        dHWsigwest  = HWsig - hb[j][i-1].H * pow(wk[j][i-1].W, sig);
        dHWsignorth = hb[j+1][i].H * pow(wk[j+1][i].W, sig) - HWsig;
        dHWsigsouth = HWsig - hb[j-1][i].H * pow(wk[j-1][i].W, sig);
        dQx     = KWeast * (dHWsigeast / dx) - KWwest * (dHWsigwest / dx);
        dQy     = KWnorth * (dHWsignorth / dy) - KWsouth * (dHWsigsouth / dy);
        dHWsigdx = (dHWsigeast + dHWsigwest) / (2.0 * dx);
        dHWsigdy = (dHWsignorth + dHWsigsouth) / (2.0 * dy);
        gradsqr = dHWsigdx * dHWsigdx + dHWsigdy * dHWsigdy;
        f[j][i].W = c1 * (dQx / dx + dQy / dy) + Cmelt * c2 * KWij * gradsqr;

        /* K_t = alpha (2 Kmax - K) ((W/Wcrit)^2 - 1) */
        Wrel = Wij / Wcrit;
        f[j][i].K = user->alpha * (2 * user->Kmax - Kij) * (Wrel * Wrel - 1);
      }
    }
  }
  ierr = DMDAVecRestoreArray(da,localX,&wk);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da,localgeom,&hb);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da,F,&f);CHKERRQ(ierr);

  ierr = DMRestoreLocalVector(da,&localX);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(da,&localgeom);CHKERRQ(ierr);
  PetscFunctionReturn(0); 
}


#undef __FUNCT__
#define __FUNCT__ "RHSJacobian"
PetscErrorCode RHSJacobian(TS ts,PetscReal t,Vec W,Mat *J,Mat *Jpre,
                           MatStructure *str,void *ctx) {
  SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER,"this Jacobian is not yet implemented");
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "getWKnorms"
PetscErrorCode getWKnorms(PorousCtx *user, Vec X,
                          PetscReal *Wav, PetscReal *Wnorminf,
                          PetscReal *Kav, PetscReal *Knorminf) {
  PetscErrorCode ierr;
  PetscInt       Mx,My;
  PetscReal      onenorms[2],infnorms[2];
  PetscFunctionBegin;  
  ierr = VecStrideNormAll(X,NORM_1,onenorms);CHKERRQ(ierr);
  ierr = VecStrideNormAll(X,NORM_INFINITY,infnorms);CHKERRQ(ierr);
  ierr = DMDAGetInfo(user->da,PETSC_IGNORE,&Mx,&My,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);
  *Wav      = onenorms[0] / ((PetscReal)Mx * (PetscReal)My);
  *Kav      = onenorms[1] / ((PetscReal)Mx * (PetscReal)My);
  *Wnorminf = infnorms[0];
  *Knorminf = infnorms[1];
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "maxDiffusivity"
PetscErrorCode maxDiffusivity(DM da, PorousCtx *user, Vec geom, Vec X,
                              PetscReal *maxD) {
  PetscErrorCode ierr;
  PetscInt       i,j,Mx,My,xs,ys,xm,ym;
  PetscReal      sig = user->sigma, CD, locmaxD;
  WKnode         **wk;
  Hbnode         **hb;
 
  /* constant in diffusivity */
  CD = sig * user->rhoi / (user->rhow * pow(user->Wcrit,sig));

  PetscFunctionBegin;
  ierr = DMDAGetInfo(da,PETSC_IGNORE,&Mx,&My,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);
  ierr = DMDAGetCorners(da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL);CHKERRQ(ierr);

  locmaxD = -1.0;
  ierr = DMDAVecGetArray(da,X,&wk);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,user->geom,&hb);CHKERRQ(ierr);
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      if (!(i == 0 || j == 0 || i == Mx-1 || j == My-1)) {
        locmaxD = PetscMax(locmaxD,
           CD * hb[j][i].H * wk[j][i].K * pow(wk[j][i].W,sig));
      }
    }
  }
  ierr = DMDAVecRestoreArray(da,X,&wk);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da,user->geom,&hb);CHKERRQ(ierr);

  ierr = MPI_Allreduce(&locmaxD,maxD,1,MPIU_REAL,MPIU_MAX,PETSC_COMM_WORLD);
      CHKERRQ(ierr);

  PetscFunctionReturn(0);
} 


#undef __FUNCT__  
#define __FUNCT__ "MyTSMonitor"
/* MyTSMonitor:  stdout report at every time step */
PetscErrorCode MyTSMonitor(TS ts,PetscInt step,PetscReal ptime,Vec X,void *ptr)
{
  PetscErrorCode ierr;
  PetscReal      twonorms[2],rnorm2,Wav,Kav,Wnorminf,Knorminf;
  PetscReal      secperday=3600.0*24.0,maxD;
  MPI_Comm       comm;
  PorousCtx      *user = (PorousCtx*)ptr;
  SNES           snes;

  PetscFunctionBegin;  
  ierr = getWKnorms(user,X,&Wav,&Wnorminf,&Kav,&Knorminf);CHKERRQ(ierr);
  ierr = maxDiffusivity(user->da, user, user->geom, X, &maxD);CHKERRQ(ierr);
  ierr = VecStrideNormAll(X,NORM_2,twonorms);CHKERRQ(ierr);
  ierr = PetscObjectGetComm((PetscObject)ts,&comm);CHKERRQ(ierr);
  /* summary */
  if (!user->run_silent) {
    if (user->fcncount <= 2) {
      ierr = PetscPrintf(comm,
           "  step  time(days)      av W(m)   max W(m)  max D(m2 s-1)  av K(m s-1)  max K(m s-1)\n");CHKERRQ(ierr);
    }
    ierr = PetscPrintf(comm,
           "   %3d     %7G  %.9f  %9.6f    %11.6f  %11.6f      %.6f\n",
           step,ptime/secperday,Wav,Wnorminf,maxD,
           Kav,Knorminf);CHKERRQ(ierr);
  }
  /* warning if solution not small */
  ierr = TSGetSNES(ts,&snes);CHKERRQ(ierr);
  ierr = SNESGetFunctionNorm(snes,&rnorm2);CHKERRQ(ierr);
  if (rnorm2 > 1.0e-10 * (twonorms[0]+twonorms[1])) {
    user->not_converged_warning = PETSC_TRUE;
    if (!user->run_silent) {
      ierr = PetscPrintf(comm,
           "***WARNING1***: residual norm not small (> 1e-10 * (|W|_2+|K|_2)) at step %d\n",
           step);CHKERRQ(ierr); }
  }

  /* update max of rnorm (relative) so far */
  if (twonorms[0] > 0.0) user->maxrnorm = PetscMax(user->maxrnorm,rnorm2);

  PetscFunctionReturn(0);
}

