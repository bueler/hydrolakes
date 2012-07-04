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
"with evolving water thickness (W) and water pressure (P) variables.\n"
"Capacity (Y) may be computed diagnostically as needed.\n"
"See vanPeltBueler.pdf for more details\n";

/* 
Example usage:

Get help:
    ./alt -fd -help | grep alt_

Finite difference evaluation of Jacobian using coloring:
    ./alt -fd
FIXME:  Jacobian evaluation routine not yet written

Options which should (nearly:FIXME?) reproduce result of porous, the verification case:
    ./alt -fd -alt_Creep 0 -alt_Cmelt 0 -alt_Wmin 0.0001 -alt_steps 10
    ./porous   # same computation, but for W only, and without Wmin>0 regularization

Minimal movie:
    ./alt -fd -da_grid_x 101 -da_grid_y 101 -ts_monitor_solution -draw_pause 0.5

Minimal convergence info:
    ./alt -fd -alt_steps 1 -alt_converge_check
    
More reporting:
    ./alt -fd -snes_monitor -snes_vi_monitor -ksp_converged_reason
    
Matlab initial and final output:
    ./alt -fd -mfile foo.m

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
extern PetscErrorCode FormPositivityBounds(SNES,Vec,Vec);
extern PetscErrorCode checkPositivity(PorousCtx *user,Vec X);
extern PetscErrorCode RHSFunction(TS,PetscReal,Vec,Vec,void*);
extern PetscErrorCode RHSJacobian(TS,PetscReal,Vec,Mat*,Mat*,
                                  MatStructure*,void*);
extern PetscErrorCode getWPnorms(PorousCtx*,Vec,
                                 PetscReal*,PetscReal*,PetscReal*,PetscReal*);
extern PetscErrorCode MyTSMonitor(TS,PetscInt,PetscReal,Vec,void*);


typedef struct {
  PetscReal W,  /* water thickness */
            P;  /* capacity thickness */
} WPnode;

/* we overload dof=2 Vecs to use these other node types: */

typedef struct {
  PetscReal ddxE,  /* (P_i+1,j - P_i,j) / dx  at  i+1/2,j */
            ddyN;  /* (P_i,j+1 - P_i,j) / dy  at  i,j+1/2 */
} staggradnode;

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
  PetscReal      tstart,tend,ftime,secperday=3600.0*24.0,Y0;
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
  steps    = 20;
  Y0       = 1.0;              /* initial value of Y, for computing initial
                                  value of P; note Ymin = 0.1 is different */
  user.da = da;
  ierr = DefaultContext(&user);CHKERRQ(ierr);

  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,
           "","options to (W,P)-space better hydrology model alt","");CHKERRQ(ierr);
  {
    ierr = PetscOptionsReal("-alt_sigma","nonlinear power","",
                            user.sigma,&user.sigma,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-alt_Ymin",
                            "min capacity thickness (esp. in pressure computation)","",
                            user.Ymin,&user.Ymin,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-alt_Wmin",
                            "min water amount (esp. in pressure computation)","",
                            user.Wmin,&user.Wmin,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-alt_Y0",
                            "constant initial capacity thickness","",
                            Y0,&Y0,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-alt_Cmelt",
                            "additional coefficient for amount of melt","",
                            user.Cmelt,&user.Cmelt,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-alt_Creep",
                            "creep closure coefficient","",
                            user.Creep,&user.Creep,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-alt_L","half-width of square region in meters","",
                            user.L,&user.L,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-alt_tstart_days","start time in days","",
                            tstart/secperday,&tstart,&optflg);CHKERRQ(ierr);
    if (optflg) { tstart *= secperday; }
    ierr = PetscOptionsReal("-alt_tend_days","end time in days","",
                            tend/secperday,&tend,&optflg);CHKERRQ(ierr);
    if (optflg) { tend *= secperday; }
    ierr = PetscOptionsInt("-alt_steps","number of timesteps to take","",
                           steps,&steps,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsBool("-alt_converge_check",
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

  /* set initial state:  W = barenblatt, P = pi (W/Y0)^sigma */
  ierr = InitialState(da,&user,tstart,Y0,X);CHKERRQ(ierr);

  /* set up times for time-stepping */
  ierr = TSSetInitialTimeStep(ts,tstart,
           (tend - tstart) / (PetscReal)steps);CHKERRQ(ierr);
  ierr = TSSetDuration(ts,steps,tend);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts,PETSC_TRUE);CHKERRQ(ierr);
  ierr = TSMonitorSet(ts,MyTSMonitor,&user,PETSC_NULL);CHKERRQ(ierr);

  /* Set SNESVI type and supply upper and lower bounds. */
  ierr = TSGetSNES(ts,&snes);CHKERRQ(ierr);
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
        "writing initial W,P and geometry H,b to Matlab file %s ...\n",
        mfile);CHKERRQ(ierr);
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
    ierr = print2vecmatlab(da,X,"W_final","P_final",mfile,PETSC_TRUE);CHKERRQ(ierr);
    ierr = printfigurematlab(da,2,"W_init","W_final",mfile,PETSC_TRUE);CHKERRQ(ierr);
    ierr = printfigurematlab(da,3,"P_init","P_final",mfile,PETSC_TRUE);CHKERRQ(ierr);
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
/* InitialState:  sets initial state of (W,P).  specifically,
   W = barenblatt at t, P = pi (W/Y0)^sigma */
PetscErrorCode InitialState(DM da, PorousCtx* user, PetscReal t, PetscReal Y0,
                            Vec X) {
  PetscErrorCode ierr;
  PetscInt       i,j,xs,ys,xm,ym;
  PetscReal      pi;
  WPnode         **wp;
  Hbnode         **hb;
  ierr = BarenblattState(da,user,t,X);CHKERRQ(ierr);  /* sets W */
  ierr = DMDAGetCorners(da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,X,&wp);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,user->geom,&hb);CHKERRQ(ierr);
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      pi = user->rhoi * user->g * hb[j][i].H;
      wp[j][i].P = pi * pow(wp[j][i].W / Y0, user->sigma);
    }
  }
  ierr = DMDAVecRestoreArray(da,user->geom,&hb);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da,X,&wp);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "FormPositivityBounds"
/*  FormPositivityBounds() for call-back: tell SNESVI (variational inequality)
  that we want positive solutions for both W and P */
PetscErrorCode FormPositivityBounds(SNES snes, Vec Xl, Vec Xu) {
  VecStrideSet(Xl,0,0.0);  /* W >= 0 */
  VecStrideSet(Xl,1,0.0);  /* P >= 0 */
  VecStrideSet(Xu,0,SNES_VI_INF);
  VecStrideSet(Xu,1,SNES_VI_INF);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "checkPositivity"
PetscErrorCode checkPositivity(PorousCtx *user,Vec X) {
  PetscErrorCode ierr;
  PetscInt       i,j,xs,ys,xm,ym;
  DM             da = user->da;
  WPnode         **wp;
  PetscFunctionBegin;
  ierr = DMDAGetCorners(da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,X,&wp);CHKERRQ(ierr);
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      if (wp[j][i].W < 0) {
        SETERRQ3(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
          "negative water thickness W at (i,j)=(%d,%d) during function eval %d",
          i,j,user->fcncount);
      }
      if (wp[j][i].P < 0) {
        SETERRQ3(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
          "negative water pressure P at (i,j)=(%d,%d) during function eval %d",
          i,j,user->fcncount);
      }
    }
  }
  ierr = DMDAVecRestoreArray(da,X,&wp);CHKERRQ(ierr);
  PetscFunctionReturn(0);
} 


/*
we solve the S=0, b=0, and v_base=0 case of the hydrology model,
in its (W,P) alternative formulation:

W_t = c1 div (W grad P) + Cmelt c2 W |grad P|^2

P_t = sigma P/(W+Wmin) [ W_t - Cmelt c3 (P/pi)^q) W |grad P|^2
                         + Creep A max{0,(pi - P)^n} (W - (P/pi)^q Ymin) ]

where
  q = 1/sigma
  pi = rhoi g H
  c1 = K / (rhow g)
  c2 = c1 / (rhow L)
  c3 = c1 / (rhoi L)
*/

#undef __FUNCT__
#define __FUNCT__ "RHSFunction"
/* RHSFunction:  evaluates nonlinear function in ODE form X' = F(X,t).
   But our case is autonomous, so F = F(X) = F(W,P) and F has no t dependence. */
PetscErrorCode RHSFunction(TS ts,PetscReal t_unused,Vec X,Vec F,void *ptr)
{
  PetscErrorCode ierr;
  PorousCtx      *user = (PorousCtx*)ptr;
  DM             da = user->da;
  Vec            dPstag,localX;
  PetscInt       i,j,Mx,My,xs,ys,xm,ym;
  PetscReal      dx = user->dx, dy = user->dy,
                 c1, c2, c3, sig, q, Ymin, Wmin,
                 Wij, Weast, Wwest, Wnorth, Wsouth,
                 dQx, dQy, divQ, 
                 dPcentx, dPcenty, WdPsqr,
                 Pij, pi, Nn, Prelpow, zz;
  WPnode         **wp, **f;
  Hbnode         **hb;
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
  sig = user->sigma;
  q = 1.0 / sig;
  Ymin = user->Ymin;
  Wmin = user->Wmin;

  /* we will be differencing X = (W,P) */
  ierr = DMGetLocalVector(da,&localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(da,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  ierr = DMGetLocalVector(da,&dPstag);CHKERRQ(ierr); /* space for staggered grad */
  ierr = DMDAVecGetArray(da,localX,&wp);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,dPstag,&stag);CHKERRQ(ierr);
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      if (i == Mx-1) {
        stag[j][i].ddxE = 0.0;  /* value will not be referenced */
      } else {
        stag[j][i].ddxE = (wp[j][i+1].P - wp[j][i].P) / dx;
      }
      if (j == My-1) {
        stag[j][i].ddyN = 0.0;  /* value will not be referenced */
      } else {
        stag[j][i].ddyN = (wp[j+1][i].P - wp[j][i].P) / dy;
      }
    }
  }
  ierr = DMDAVecRestoreArray(da,dPstag,&stag);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da,localX,&wp);CHKERRQ(ierr);

  /* we will be differencing dPstag = (ddxE,ddyN) */
  ierr = DMDALocalToLocalBegin(da,dPstag,INSERT_VALUES,dPstag);CHKERRQ(ierr);
  ierr = DMDALocalToLocalEnd(da,dPstag,INSERT_VALUES,dPstag);CHKERRQ(ierr);

  ierr = DMDAVecGetArray(da,localX,&wp);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,F,&f);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,user->geom,&hb);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,dPstag,&stag);CHKERRQ(ierr);
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      if (i == 0 || j == 0 || i == Mx-1 || j == My-1) {
        f[j][i].W = 0.0;  /* no change at Dirichlet boundary conditions */
        f[j][i].P = 0.0;
      } else {
        Wij     = wp[j][i].W;
        Weast   = 0.5 * (wp[j][i+1].W + Wij);
        Wwest   = 0.5 * (wp[j][i-1].W + Wij);
        Wnorth  = 0.5 * (wp[j+1][i].W + Wij);
        Wsouth  = 0.5 * (wp[j-1][i].W + Wij);

        dPcentx = (wp[j][i+1].P - wp[j][i-1].P) / dx;
        dPcenty = (wp[j+1][i].P - wp[j-1][i].P) / dy;
        WdPsqr  = Wij * ( dPcentx * dPcentx + dPcenty * dPcenty );

        /* W_t = c1 div (W grad P) + c2 W |grad P|^2 */
        dQx       = Weast  * stag[j][i].ddxE - Wwest  * stag[j][i-1].ddxE;
        dQy       = Wnorth * stag[j][i].ddyN - Wsouth * stag[j-1][i].ddyN;
        divQ      = dQx / dx + dQy / dy;
        f[j][i].W = c1 * divQ + user->Cmelt * c2 * WdPsqr;

        /* P_t = sigma P/(W+Wmin) [ W_t - Cmelt c3 (P/pi)^q) W |grad P|^2
                                    + Creep A max{0,N}^n (W - (P/pi)^q Ymin) ] */
        pi        = user->rhoi * user->g * hb[j][i].H;
        Pij       = wp[j][i].P;
        Nn        = pow(PetscMax(0.0,pi - Pij),user->nglen);
        Prelpow   = pow(Pij / pi,q);
        zz        = f[j][i].W - user->Cmelt * c3 * Prelpow * WdPsqr
                      + user->Creep * user->Aglen * Nn * (Wij - Prelpow * Ymin);
        f[j][i].P = sig * (Pij / (Wij + Wmin)) * zz;
      }
    }
  }
  ierr = DMDAVecRestoreArray(da,localX,&wp);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da,F,&f);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da,user->geom,&hb);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da,dPstag,&stag);CHKERRQ(ierr);

  ierr = DMRestoreLocalVector(da,&localX);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(da,&dPstag);CHKERRQ(ierr);
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
#define __FUNCT__ "getWPnorms"
PetscErrorCode getWPnorms(PorousCtx *user, Vec X,
                          PetscReal *Wnorm1, PetscReal *Wnorminf,
                          PetscReal *Pnorm1, PetscReal *Pnorminf) {
  PetscErrorCode ierr;
  PetscInt       j;
  PetscReal      onenorms[2],infnorms[2];
  PetscFunctionBegin;  
  ierr = VecStrideNormAll(X,NORM_1,onenorms);CHKERRQ(ierr);
  for (j=0; j<2; j++)  onenorms[j] *= user->dx * user->dy;
  ierr = VecStrideNormAll(X,NORM_INFINITY,infnorms);CHKERRQ(ierr);
  *Wnorm1   = onenorms[0];  *Pnorm1   = onenorms[1];
  *Wnorminf = infnorms[0];  *Pnorminf = infnorms[1];
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "MyTSMonitor"
/* MyTSMonitor:  stdout report at every time step */
PetscErrorCode MyTSMonitor(TS ts,PetscInt step,PetscReal ptime,Vec X,void *ptr)
{
  PetscErrorCode ierr;
  PetscReal      twonorms[2],rnorm2,Wnorm1,Wnorminf,Pnorm1,Pnorminf;
  PetscReal      secperday=3600.0*24.0,CD;
  MPI_Comm       comm;
  PorousCtx      *user = (PorousCtx*)ptr;
  SNES           snes;

  PetscFunctionBegin;  
  ierr = getWPnorms(user,X,&Wnorm1,&Wnorminf,&Pnorm1,&Pnorminf);CHKERRQ(ierr);
  CD = (user->Kconst * user->sigma) / (user->rhow * user->g);
  ierr = VecStrideNormAll(X,NORM_2,twonorms);CHKERRQ(ierr);
  ierr = PetscObjectGetComm((PetscObject)ts,&comm);CHKERRQ(ierr);
  /* summary */
  if (!user->run_silent) {
/*
    ierr = PetscPrintf(comm,
           "step %3d at %7G days:  |W|_1=%.9e m3,  max W=%.3f m,\n"
           "                           max D=%.3f m2s-1,  max P=%.3f bar\n",
           step,ptime/secperday,Wnorm1,Wnorminf,
           CD*Pnorminf,Pnorminf/1.0e5);CHKERRQ(ierr);
*/
    if (user->fcncount <= 2) {
      ierr = PetscPrintf(comm,
           "  step  time(days)     |W|_1(m3)    max W(m)  max D(m2 s-1)    max P(bar)\n");CHKERRQ(ierr);
    }
    ierr = PetscPrintf(comm,
           "   %3d     %7G  %11.6e %11.6f    %11.6f   %11.6f\n",
           step,ptime/secperday,Wnorm1,Wnorminf,CD*Pnorminf,
           Pnorminf/1.0e5);CHKERRQ(ierr);
  }
  /* warning if solution not small */
  ierr = TSGetSNES(ts,&snes);CHKERRQ(ierr);
  ierr = SNESGetFunctionNorm(snes,&rnorm2);CHKERRQ(ierr);
  if (rnorm2 > 1.0e-10 * (twonorms[0]+twonorms[1])) {
    user->not_converged_warning = PETSC_TRUE;
    if (!user->run_silent) {
      ierr = PetscPrintf(comm,
           "***WARNING1***: residual norm not small (> 1e-10 * (|W|_2+|P|_2)) at step %d\n",
           step);CHKERRQ(ierr); }
  }

  /* update max of rnorm (relative) so far */
  if (twonorms[0] > 0.0) user->maxrnorm = PetscMax(user->maxrnorm,rnorm2);

  PetscFunctionReturn(0);
}

