/*
   Copyright (C) 2011,2012 Ed Bueler
  
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

static char help[] = "Solves time-dependent, vertically-integrated porous\n"
"medium equation in 2D for subglacial water thickness.\n"
"Uses PETSc TS and SNES to do implicit time steps.  Defaults to Crank-Nicolson.\n"
"Uses SNESVI to constrain for nonnegative values of the water thickness.\n";

/* 
We solve this vertically-integrated porous medium equation (VPME)

    W_t = gamma div ( W grad (W^sigma) )

where  W = W(t,x,y)  is a water depth and gamma is a constant:

    gamma = (K H0 rhoi) / (rhow Ybaren^sigma).  

All of the values on the right side of this equation for gamma are
constant in this verification-focused case.

Let (x_i,y_j) be a grid in a rectangle, with uniform spacing dx,dy.
Let W(t) = W[j][i](t) be an approximation to W(t,x_i,y_j).
Semi-discretization by finite differences in space gives:

         gamma
(1) W' = ----- ( Weast (W[j][i+1]^s - W^s) - Wwest (W^s - W[j][i-1]^s) )
          dx^2
           gamma
         + ----- ( Wnorth (W[j+1][i]^s - W^s) - Wsouth (W^s - W[j-1][i]^s) )
            dy^2

where  s = sigma  and  W = W[j][i]  (notation for simplicity), and where
Weast,Wwest, Wnorth,Wsouth are values of W average onto the staggered grid (e.g.
Weast = (W[j][i] + W[j][i+1])/2).

Equation (1) is computed by RHSFunction().  The derivatives of (1) are 
computed by RHSJacobian().

By default we use Crank-Nicolson time-stepping.  By default we use a
SNESVI type with bounds on W to avoid negative values.

The exit status of the program is 0 if the SNES reports at least 10^-10
reduction in residual norm at every time step (converged) and 1 otherwise
(not converged).

Example usage follows.

Get help:
  ./porous -help | grep por_

Minimal movie:
  ./porous -da_grid_x 101 -da_grid_y 101 \
       -ts_monitor_solution -draw_pause 0.5 -por_steps 20

Parallel runs, spatial refinement only:
  for M in 21 41 81 161 321; do
    echo "case M=$M:"
    mpiexec -n 4 ./porous -da_grid_x $M -da_grid_y $M |grep -A 1 errors
  done

Parallel runs, temporal and spatial refinement:
  for N in 2 4 8 16 32 64 128; do
    (( M = 10*N+1 ))
    echo "case N=$N, M=$M:"
    mpiexec -n 4 ./porous -por_steps $N -da_grid_x $M -da_grid_y $M |grep -A 1 errors
  done

Minimal convergence info:
  ./porous -por_steps 1 -por_converge_check
For a plot of Newton convergence, see newtonconverge.sh.

Fails to converge:
  mpiexec -n 4 ./porous -da_grid_x 551 -da_grid_y 551 \
      -snes_monitor -ksp_converged_reason
But this succeeds:
  mpiexec -n 4 ./porous -da_grid_x 551 -da_grid_y 551 \
      -snes_monitor -ksp_converged_reason -pc_type asm -sub_pc_type lu
Just doubling the number of steps also succeeds (default is -por_steps 10):
  mpiexec -n 4 ./porous -da_grid_x 551 -da_grid_y 551 \
      -snes_monitor -ksp_converged_reason -por_steps 20
    
With finite difference evaluation of Jacobian using coloring:
  ./porous -fd

More reporting:
  ./porous -snes_monitor -ts_monitor_solution -draw_pause 0.5 -mfile foo.m

Here are some relevant TS options:

  -ts_type <cn>: TS method (one of) euler beuler cn pseudo gl ssp theta alpha
  -ts_monitor_draw: <FALSE> Monitor timestep size graphically (TSMonitorLG)
  -ts_monitor_solution: <FALSE> Monitor solution graphically (TSMonitorSolution)
*/

#include <petscdmda.h>
#include <petscts.h>

#include "matlabprint.h"  /* utilities for putting petsc objects in .m file */
#include "context.h"      /* default values for parameters, and Barenblatt soln */

extern PetscErrorCode FormPositivityBounds(SNES,Vec,Vec);
extern PetscErrorCode RHSFunction(TS,PetscReal,Vec,Vec,void*);
extern PetscErrorCode RHSJacobian(TS,PetscReal,Vec,Mat*,Mat*,MatStructure*,void*);
extern PetscErrorCode maxDiffusivity(DM da,PorousCtx *user,Vec X, PetscReal *maxD);
extern PetscErrorCode MyTSMonitor(TS,PetscInt,PetscReal,Vec,void*);


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  TS             ts;                   /* time-stepping object (contains snes) */
  SNES snes;
  Vec            W,Wexact,r;           /* solution, exact solution, residual vector */
  Mat            J;                    /* Jacobian matrix */
  PetscInt       Mx,My,fsteps;
  DM             da;
  ISColoring     iscoloring;
  PetscReal      dx,dy,ftime,secperday=3600.0*24.0,Wnorm1,error1,errorinf;
  PetscBool      fdflg = PETSC_FALSE, mfileflg = PETSC_FALSE, optflg = PETSC_FALSE;
  char           mfile[PETSC_MAX_PATH_LEN] = "porout.m";
  MatFDColoring  matfdcoloring;
  PorousCtx      user;                 /* user-defined work context */
  PetscReal      tstart, tend;
  PetscInt       steps;

  PetscInitialize(&argc,&argv,(char *)0,help);

  ierr = DMDACreate2d(PETSC_COMM_WORLD,
             DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_NONE, // correct for zero Dirichlet
             DMDA_STENCIL_STAR, // nonlinear diffusion but diffusivity
                                //   depends on soln W not grad W
             -21,-21,           // default to 20x20 grid but override with
                                //   -da_grid_x, -da_grid_y (or -da_refine)
             PETSC_DECIDE,PETSC_DECIDE, // num of procs in each dim
             1,1,               // dof = 1, s = 1 (stencil extends out one cell)
             PETSC_NULL,PETSC_NULL, // no specify proc decomposition
             &da);CHKERRQ(ierr);
  ierr = DMSetApplicationContext(da,&user);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(da,&W);CHKERRQ(ierr);
  ierr = VecDuplicate(W,&Wexact);CHKERRQ(ierr);
  ierr = VecDuplicate(W,&r);CHKERRQ(ierr);

  /* default constants & setup context */
  tstart   = 10.0 * secperday; /* 10 days in seconds */
  tend     = 30.0 * secperday;
  steps    = 10;
  user.da = da;
  ierr = DefaultContext(&user);CHKERRQ(ierr);

  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,
           "","options to porous medium equation","");
  {
    ierr = PetscOptionsReal("-por_sigma","nonlinear power","",
                            user.sigma,&user.sigma,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-por_Ybaren",
                            "critical water thickness in verification case","",
                            user.Ybaren,&user.Ybaren,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-por_L","half-width of square region in meters","",
                            user.L,&user.L,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-por_tstart_days","start time in days","",
                            tstart/secperday,&tstart,&optflg);CHKERRQ(ierr);
    if (optflg) { tstart *= secperday; }
    ierr = PetscOptionsReal("-por_tend_days","end time in days","",
                            tend/secperday,&tend,&optflg);CHKERRQ(ierr);
    if (optflg) { tend *= secperday; }
    ierr = PetscOptionsInt("-por_steps","number of timesteps to take","",
                           steps,&steps,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsBool("-por_converge_check",
                            "run silent and check for Newton convergence",
                            "",user.run_silent,&user.run_silent,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsString("-mfile",
                            "name of Matlab file to write results","",
                            mfile,mfile,PETSC_MAX_PATH_LEN,&mfileflg);CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);

  ierr = DerivedConstants(&user);CHKERRQ(ierr);
  ierr = DMDASetUniformCoordinates(da,  // square domain
              -user.L, user.L, -user.L, user.L, 0.0, 1.0);CHKERRQ(ierr);

  ierr = TSCreate(PETSC_COMM_WORLD,&ts);CHKERRQ(ierr);
  ierr = TSSetType(ts,TSCN);CHKERRQ(ierr);
  ierr = TSSetRHSFunction(ts,r,RHSFunction,&user);CHKERRQ(ierr);

  /* let TS know about Jacobian */
  ierr = DMGetMatrix(da,MATAIJ,&J);CHKERRQ(ierr);

  /* use coloring to compute rhs Jacobian efficiently; does SNESVI need this? */
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
  } else {
    /* default case */
    ierr = TSSetRHSJacobian(ts,J,J,RHSJacobian,&user);CHKERRQ(ierr);
  }

  /* set initial state */
  ierr = BarenblattState(&user,tstart,W);CHKERRQ(ierr);
  if (mfileflg) {
    ierr = printvecmatlab(da,W,"Winitial",mfile,PETSC_FALSE);CHKERRQ(ierr);
  }

  /* set up time-stepping */
  ierr = TSSetInitialTimeStep(ts,tstart,
           (tend - tstart) / (PetscReal)steps);CHKERRQ(ierr);
  ierr = TSSetDuration(ts,steps,tend);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts,PETSC_TRUE);CHKERRQ(ierr);
  ierr = TSMonitorSet(ts,MyTSMonitor,&user,PETSC_NULL);CHKERRQ(ierr);

  /* Set SNESVI type and supply upper and lower bounds.
  The following is a workaround to fix a bug in petsc3.2; see
     http://petsc.cs.iit.edu/petsc/petsc-dev/rev/2d5e73b75a6c#l1.2
  We would like to use SNESVISetVariableBounds, but use
  a call-back instead to avoid some error checking issue. */
  ierr = TSGetSNES(ts,&snes);CHKERRQ(ierr);
  ierr = SNESVISetComputeVariableBounds(snes,FormPositivityBounds);
        CHKERRQ(ierr);

  /* ask user to finalize settings */
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);

  /* report on setup */
  ierr = DMDAGetInfo(da,PETSC_IGNORE,&Mx,&My,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);
  dx = 2.0 * user.L / (Mx-1);
  dy = 2.0 * user.L / (My-1);
  if (!user.run_silent) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
    "setup done: square       side length = %.3f km\n"
    "            grid               Mx,My = %d,%d\n"
    "            spacing            dx,dy = %.3f,%.3f m\n"
    "            times     tstart:dt:tend = %.3f:%.3f:%.3f days\n",
    2.0 * user.L / 1000.0,
    Mx, My,
    dx, dy,
    tstart / secperday, (tend-tstart)/(steps*secperday), tend / secperday);
    CHKERRQ(ierr); }

  /* run time-stepping with implicit steps  */
  ierr = TSSolve(ts,W,&ftime);CHKERRQ(ierr);

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
    ierr = printvecmatlab(da,W,"Wfinal",mfile,PETSC_TRUE);CHKERRQ(ierr);
    ierr = printfigurematlab(da,1,"Winitial","Wfinal",mfile,PETSC_TRUE);CHKERRQ(ierr);
  }

  /* compare to exact */
  ierr = VecNorm(W,NORM_1,&Wnorm1);CHKERRQ(ierr);  Wnorm1 *= dx*dy;
  ierr = BarenblattState(&user,ftime,Wexact);CHKERRQ(ierr);
  if (mfileflg) {
    ierr = printvecmatlab(da,Wexact,"Wfinalexact",mfile,PETSC_TRUE);CHKERRQ(ierr);
    /* uncomment next line if Mat wanted;
       note "-mat_ascii_output_large" option probably needed */
    /* ierr = printmatmatlab(da,J,"J",mfile,PETSC_TRUE);CHKERRQ(ierr); */  
  }
  ierr = VecWAXPY(r,-1.0,Wexact,W);CHKERRQ(ierr);  /* r = W - Wexact */
  ierr = VecNorm(r,NORM_1,&error1);CHKERRQ(ierr);
  error1 /= (PetscReal)Mx * (PetscReal)My;
  ierr = VecNorm(r,NORM_INFINITY,&errorinf);CHKERRQ(ierr);
  if (!user.run_silent) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
    "errors:    av |W-Wexact|  = %.3e m\n"
    "           |W-Wexact|_inf = %.3e m\n",
    error1,errorinf);CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "%6d  %6d  %9.3f  %.12e\n",
                       Mx, My, (tend-tstart)/secperday, user.maxrnorm);CHKERRQ(ierr);
  }

  /* Free work space.  */
  ierr = MatDestroy(&J);CHKERRQ(ierr);
  if (fdflg) { ierr = MatFDColoringDestroy(&matfdcoloring);CHKERRQ(ierr); }
  ierr = VecDestroy(&W);CHKERRQ(ierr);
  ierr = VecDestroy(&Wexact);CHKERRQ(ierr);
  ierr = VecDestroy(&r);CHKERRQ(ierr);
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = DMDestroy(&da);CHKERRQ(ierr);

  ierr = PetscFinalize();CHKERRQ(ierr);

  PetscFunctionReturn((PetscInt)(user.not_converged_warning
                                 || user.not_conservation_warning));
}


#undef __FUNCT__
#define __FUNCT__ "FormPositivityBounds"
/*  FormPositivityBounds:  for call-back: tell SNESVI (variational inequality),
  that we want positive solutions */
PetscErrorCode FormPositivityBounds(SNES snes, Vec xl, Vec xu) {
  VecSet(xl,0.0);
  VecSet(xu,SNES_VI_INF);
  return(0);
}


#undef __FUNCT__
#define __FUNCT__ "RHSFunction"
/* RHSFunction - Evaluates nonlinear function in ODE form W' = F(W,t).
   But our case is autonomous, so F = F(W). */
PetscErrorCode RHSFunction(TS ts,PetscReal t,Vec W,Vec F,void *ptr)
{
  PetscErrorCode ierr;
  PorousCtx      *user = (PorousCtx*)ptr;
  DM             da = (DM)user->da;
  Vec            localW;
  PetscInt       i,j,Mx,My,xs,ys,xm,ym;
  PetscScalar    hx,hy,Cx,Cy,sig;
  PetscScalar    Wij, Wpow, Weast, Wwest, Wnorth, Wsouth, Qx, Qy;
  PetscScalar    **w,**f;

  PetscFunctionBegin;
  ierr = DMDAGetInfo(da,PETSC_IGNORE,&Mx,&My,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);
  ierr = DMDAGetCorners(da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL);CHKERRQ(ierr);

  hx = 2.0 * user->L / (PetscReal)(Mx-1);
  Cx = user->gamma / (hx*hx);
  hy = 2.0 * user->L / (PetscReal)(My-1);
  Cy = user->gamma / (hy*hy);
  sig = user->sigma;

  /* We will difference current state W so we get a local vector which can
     store ghost points, and then we scatter values from W to the local 
     vector.  */
  ierr = DMGetLocalVector(da,&localW);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(da,W,INSERT_VALUES,localW);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da,W,INSERT_VALUES,localW);CHKERRQ(ierr);

  ierr = DMDAVecGetArray(da,localW,&w);CHKERRQ(ierr);
  /* check for nonnegativity */
  user->fcncount = user->fcncount + 1;
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      if (w[j][i] < 0) {
        SETERRQ3(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
        "negative water thickness W value at (i,j)=(%d,%d) during function eval %d",
        i,j,user->fcncount);
      }
    }
  }

  ierr = DMDAVecGetArray(da,F,&f);CHKERRQ(ierr);
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      if (i == 0 || j == 0 || i == Mx-1 || j == My-1) {
        f[j][i] = 0.0;  /* no change at Dirichlet boundary conditions */
      } else {
        Wij     = w[j][i];
        Weast   = 0.5 * (w[j][i+1] + Wij);
        Wwest   = 0.5 * (w[j][i-1] + Wij);
        Wnorth  = 0.5 * (w[j+1][i] + Wij);
        Wsouth  = 0.5 * (w[j-1][i] + Wij);
        Wpow    = pow(Wij,sig);
        Qx      = Weast  * (pow(w[j][i+1],sig) - Wpow)
                      - Wwest  * (Wpow - pow(w[j][i-1],sig));
        Qy      = Wnorth * (pow(w[j+1][i],sig) - Wpow)
                      - Wsouth * (Wpow - pow(w[j-1][i],sig));
        f[j][i] = Cx * Qx + Cy * Qy;
      }
    }
  }
  ierr = DMDAVecRestoreArray(da,localW,&w);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da,F,&f);CHKERRQ(ierr);

  ierr = DMRestoreLocalVector(da,&localW);CHKERRQ(ierr);

  PetscFunctionReturn(0); 
}


#define G(a,b,s) ( ((s+1)/2) * pow(a,s) + (b/2) * (s * pow(a,s-1) - pow(b,s-1)) )

#undef __FUNCT__
#define __FUNCT__ "RHSJacobian"
PetscErrorCode RHSJacobian(TS ts,PetscReal t,Vec W,Mat *J,Mat *Jpre,
                           MatStructure *str,void *ctx) {
  PetscErrorCode ierr;
  PorousCtx      *user = (PorousCtx*)ctx;
  DM             da = (DM)user->da;
  Mat            jac = *Jpre;
  Vec            localW;
  PetscInt       Mx,My,xs,ys,xm,ym,i,j;
  MatStencil     col[5],row;
  PetscScalar    v[5], **w;
  PetscReal      hx,hy,Cx,Cy;
  PetscReal      sig = user->sigma;
  PetscScalar    Wij, Wpow, Wwest,  Weast,  Wsouth, Wnorth, tmp1, tmp2;

  PetscFunctionBegin;
  ierr = DMDAGetInfo(da,PETSC_IGNORE,&Mx,&My,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);
  ierr = DMDAGetCorners(da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL);CHKERRQ(ierr);

  hx = 2.0 * user->L / (PetscReal)(Mx-1);
  Cx = user->gamma / (hx*hx);
  hy = 2.0 * user->L / (PetscReal)(My-1);
  Cy = user->gamma / (hy*hy);

  /* We will difference current state W so we get a local vector which can
     store ghost points, and then we scatter values from W to it.  */
  ierr = DMGetLocalVector(da,&localW);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(da,W,INSERT_VALUES,localW);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da,W,INSERT_VALUES,localW);CHKERRQ(ierr);

  /* Compute entries for the locally owned part of the Jacobian.
      - Currently, all PETSc parallel matrix formats are partitioned by
        contiguous chunks of rows across the processors. 
      - Each processor needs to insert only elements that it owns
        locally (but any non-local elements will be sent to the
        appropriate processor during matrix assembly). 
      - Here, we set all entries for a particular row at once.  */
  ierr = DMDAVecGetArray(da,localW,&w);CHKERRQ(ierr);
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      row.j = j; row.i = i;  /* identifies the equation:  k = k(i,j) */
      if (i == 0 || j == 0 || i == Mx-1 || j == My-1) {
        /* boundary points */
        v[0] = 0.0;
        ierr = MatSetValuesStencil(jac,1,&row,1,&row,v,INSERT_VALUES);
            CHKERRQ(ierr);
      } else {
        /* interior grid points; col[.] identifies the variable:  l = l(i,j)
           note w[j][i] = W_m(x_i,y_j)   (a current iterate value) */
        Wij     = w[j][i];
        Weast   = 0.5 * (w[j][i+1] + Wij);
        Wwest   = 0.5 * (w[j][i-1] + Wij);
        Wnorth  = 0.5 * (w[j+1][i] + Wij);
        Wsouth  = 0.5 * (w[j-1][i] + Wij);
        Wpow    = pow(Wij,sig);
        /* for \partial F_k / \partial W : */
        col[0].j = j;     col[0].i = i;
        tmp1 = 0.5 * (pow(w[j][i+1],sig) - Wpow) - 0.5 * (Wpow - pow(w[j][i-1],sig))
                 - sig * pow(Wij,sig-1) * (Weast + Wwest);
        tmp2 = 0.5 * (pow(w[j+1][i],sig) - Wpow) - 0.5 * (Wpow - pow(w[j-1][i],sig))
                 - sig * pow(Wij,sig-1) * (Wnorth + Wsouth);
        v[0] = Cx * tmp1 + Cy * tmp2;
        /* for \partial F_k / \partial W[j][i-1] : */
        col[1].j = j;     col[1].i = i-1;
        v[1] = Cx * G(w[j][i-1],Wij,sig);
        /* for \partial F_k / \partial W[j][i+1] : */
        col[2].j = j;     col[2].i = i+1;
        v[2] = Cx * G(w[j][i+1],Wij,sig);
        /* for \partial F_k / \partial W[j-1][i] : */
        col[3].j = j - 1; col[3].i = i;
        v[3] = Cy * G(w[j-1][i],Wij,sig);
        /* for \partial F_k / \partial W[j+1][i] : */
        col[4].j = j + 1; col[4].i = i;
        v[4] = Cy * G(w[j+1][i],Wij,sig);
        /* done with row; insert it */
        ierr = MatSetValuesStencil(jac,1,&row,5,col,v,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
  }
  ierr = DMDAVecRestoreArray(da,localW,&w);CHKERRQ(ierr);

  ierr = DMRestoreLocalVector(da,&localW);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* Tell the matrix we will never add a new nonzero location to the
     matrix. If we do, it will generate an error.                    */
  ierr = MatSetOption(jac,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "maxDiffusivity"
PetscErrorCode maxDiffusivity(DM da,PorousCtx *user,Vec X, PetscReal *maxD) {
  PetscErrorCode ierr;
  PetscInt       i,j,Mx,My,xs,ys,xm,ym;
  PetscScalar    sig = user->sigma, CD, locmaxD, **x;

  /* constant in diffusivity */
  CD = sig * user->H0 * user->Kconst * user->rhoi
          / (user->rhow * pow(user->Ybaren,sig));

  PetscFunctionBegin;
  ierr = DMDAGetInfo(da,PETSC_IGNORE,&Mx,&My,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);
  ierr = DMDAGetCorners(da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL);CHKERRQ(ierr);
  locmaxD = -1.0;
  ierr = DMDAVecGetArray(da,X,&x);CHKERRQ(ierr);
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      if (i == 0 || j == 0 || i == Mx-1 || j == My-1) {
        /* no action */
      } else {
        locmaxD = PetscMax(locmaxD, CD * pow(PetscAbs(x[j][i]),sig));
      }
    }
  }
  ierr = DMDAVecRestoreArray(da,X,&x);CHKERRQ(ierr);
  ierr = MPI_Allreduce(&locmaxD,maxD,1,MPIU_REAL,MPIU_MAX,PETSC_COMM_WORLD);
      CHKERRQ(ierr);

  PetscFunctionReturn(0);
} 


#undef __FUNCT__  
#define __FUNCT__ "MyTSMonitor"
PetscErrorCode MyTSMonitor(TS ts,PetscInt step,PetscReal ptime,Vec W,void *ptr)
{
  PetscErrorCode ierr;
  PetscReal      dx,dy,norm[2],maxW,maxD,rnorm2,secperday=3600.0*24.0;
  PetscInt       Mx,My;
  MPI_Comm       comm;
  PorousCtx      *user = (PorousCtx*)ptr;
  DM             da = (DM)user->da;
  SNES           snes;

  PetscFunctionBegin;  
  /* stdout report at every step */
  ierr = DMDAGetInfo(da,PETSC_IGNORE,&Mx,&My,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);
  dx = 2.0 * user->L / (Mx-1);
  dy = 2.0 * user->L / (My-1);
  ierr = VecNorm(W,NORM_1_AND_2,norm);CHKERRQ(ierr);
  norm[0] /= (PetscReal)Mx * (PetscReal)My;
  ierr = VecNorm(W,NORM_INFINITY,&maxW);CHKERRQ(ierr);
  ierr = maxDiffusivity(da,user,W,&maxD);CHKERRQ(ierr);  
  ierr = PetscObjectGetComm((PetscObject)ts,&comm);CHKERRQ(ierr);
  if (!user->run_silent) {
    if (user->fcncount <= 2) {
      ierr = PetscPrintf(comm,
        "  step  time(days)      av W(m)   max W(m)  max D(m2 s-1)  tau_expl(days)\n");
        CHKERRQ(ierr);
    }
    ierr = PetscPrintf(comm,
        "   %3d     %7G  %.9f  %9.6f    %11.6f     %11.6f\n",
        step,ptime/secperday,norm[0],maxW,maxD,dx*dy/(maxD*secperday));
        CHKERRQ(ierr);
  }
  /* warning if solution not small */
  ierr = TSGetSNES(ts,&snes);CHKERRQ(ierr);
  ierr = SNESGetFunctionNorm(snes,&rnorm2);CHKERRQ(ierr);
  if (rnorm2 > 1.0e-10 * norm[1]) {
    user->not_converged_warning = PETSC_TRUE;
    if (!user->run_silent) {
      ierr = PetscPrintf(comm,
           "***WARNING1***: residual norm not small (> 1e-10 * |W|_2) at step %d\n",
           step);CHKERRQ(ierr); }
  }
  /* update max of rnorm (relative) so far */
  if (norm[1] > 0.0) user->maxrnorm = PetscMax(user->maxrnorm,rnorm2);
  /* check for conservation error */
  if (user->lastWnorm1 > 0) {
    if (PetscAbs(norm[0] - user->lastWnorm1)/user->lastWnorm1 > 1.0e-8) {
      user->not_conservation_warning = PETSC_TRUE;
      if (!user->run_silent) {
        ierr = PetscPrintf(comm,
           "***WARNING2***: conservation error detected (d|W|_1 > 1e-8 * |W|_1) at step %d\n",
           step);CHKERRQ(ierr); }
    }
  }
  user->lastWnorm1 = norm[0];
  PetscFunctionReturn(0);
}

