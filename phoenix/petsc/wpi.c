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

static char help[] =
"Solves implicit time-step equations for W = (water thickness) and\n"
"psi = (hydraulic potential) form of the hydrology model.  Uses SNESVI.\n"
"Includes physically-based opening and closing processes based on\n"
"cavity-opening sliding and creep closure of cavities.  Includes a wall\n"
"melt model including turbulent dissipation heating, friction heating,\n"
"and conductive heating.  Allows general bed elevations and ice thickness.\n"
"Applies Dirichlet boundary conditions for psi and inflow boundary condition\n"
"for W.  Compares to exact solution in verification cases 1 and 2.\n"
"See notes phoenix.pdf.\n";

/* Example usage.

Get help:
    ./wpi -help | less

Verification cases, with finite difference evaluation of Jacobian using coloring:
    ./wpi -fd -case 1   # default case
    ./wpi -fd -case 2
(Direct Jacobian evaluation not yet implemented.)

Parallel runs, spatial refinement:
  case=1; for M in 21 41 81 161; do
    echo "case $case with M=$M points:"
    mpiexec -n 4 ./wpi -fd -case $case -da_grid_x $M -da_grid_y $M |grep errors
  done

Minimal movie:
    ./wpi -fd -da_grid_x 101 -da_grid_y 101 -snes_monitor_solution -draw_pause 0.5
    
More reporting:
    ./wpi -fd -snes_monitor -snes_vi_monitor -ksp_converged_reason
    
Matlab initial and final output:
    ./wpi -fd -mfile foo.m
*/

#include <petscdmda.h>
#include <petscsnes.h>

#include "matlabprint.h"  /* utilities for putting petsc objects in .m file */

#include "context.h"
#include "verifcase1.h"
#include "verifcase2.h"
#include "expcase3.h"


/* for both psi and W^{l+1} equations */
extern PetscErrorCode checkBounds(WPCtx*,Vec,Vec);
/* for psi equation */
extern PetscErrorCode FormBounds(SNES,Vec,Vec);
extern PetscErrorCode FormFunction(SNES,Vec,Vec,void*);
extern PetscErrorCode FormJacobian(SNES,Vec,Mat*,Mat*,MatStructure*,void*);
/* for W^{l+1} equation */
extern PetscErrorCode FormBoundsUpdateW(SNES,Vec,Vec);
extern PetscErrorCode FormFunctionUpdateW(SNES,Vec,Vec,void*);
extern PetscErrorCode FormJacobianUpdateW(SNES,Vec,Mat*,Mat*,MatStructure*,void*);


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  SNES           snes, snesUpdateW;
  Vec            x,   /* solution (P or Wnew, resp.) */
                 r,   /* residual; same type */
                 initialWnew, /* put Wnew guess here */
                 psiexact, Wnewexact;  /* available in verification cases */
  Mat            J;   /* space for Jacobian */
  PetscInt       N=1, /* number of time steps */
                 Mx, My, /* spatial grid points: set with -da_grid_{x,y} */
                 its, /* SNES iteration reporting */
                 n,   /* timestep counter */
                 verifcase=1;

  SNESConvergedReason reason;
  ISColoring          iscoloring;
  MatFDColoring       matfdcoloring;

  WPCtx          user;
  PetscReal      error1, errorinf, exactinf;
  PetscBool      fdflg = PETSC_FALSE, mfileflg = PETSC_FALSE,
                 Lsetflg = PETSC_FALSE;
  char           mfile[PETSC_MAX_PATH_LEN] = "wpiout.m";

  PetscInitialize(&argc,&argv,(char *)0,help);

  /* set up contexts */
  ierr = DefaultContext(&user);CHKERRQ(ierr);

  ierr = DMDACreate2d(PETSC_COMM_WORLD,
             DMDA_BOUNDARY_NONE, /* Dirichlet conditions, so no ghosts out there */
             DMDA_BOUNDARY_NONE,
             DMDA_STENCIL_STAR, /* divergence form 2nd order */
             -11,-11,           /* default to 10x10 grid but override with
                                     -da_grid_x, -da_grid_y 
                                   (*but*, because periodic, must be multiples of 5) */
             PETSC_DECIDE,PETSC_DECIDE, /* num of procs in each dim */
             1,1,               /* dof = 2, s = 1 (stencil extends out one cell) */
             PETSC_NULL,PETSC_NULL, /* no specify proc decomposition */
             &user.da);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(user.da,&x);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&r);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&initialWnew);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&psiexact);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&Wnewexact);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&(user.psiforW));CHKERRQ(ierr);
  ierr = VecDuplicate(x,&(user.Wl));CHKERRQ(ierr);
  ierr = VecDuplicate(x,&(user.S));CHKERRQ(ierr);
  ierr = VecDuplicate(x,&(user.vbspeed));CHKERRQ(ierr);
  ierr = VecDuplicate(x,&(user.fric));CHKERRQ(ierr);
  ierr = VecDuplicate(x,&(user.Ggeo));CHKERRQ(ierr);
  ierr = VecDuplicate(x,&(user.H));CHKERRQ(ierr);
  ierr = VecDuplicate(x,&(user.b));CHKERRQ(ierr);
  ierr = VecDuplicate(x,&(user.psiBC));CHKERRQ(ierr);
  ierr = VecDuplicate(x,&(user.WnewBC));CHKERRQ(ierr);

  ierr = DMGetMatrix(user.da,MATAIJ,&J);CHKERRQ(ierr);

  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,
           "","options to implicit-timestepping (W,psi) space hydrology solver","");
  {
    ierr = PetscOptionsString("-mfile",
               "name of Matlab file to write results","",
               mfile,mfile,PETSC_MAX_PATH_LEN,&mfileflg);CHKERRQ(ierr);
    ierr = PetscOptionsBool("-fd",
               "use coloring to compute Jacobian by finite differences",
               PETSC_NULL,fdflg,&fdflg,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsInt("-N",
               "number of time steps to take","",
               N,&N,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsInt("-case",
               "1 or 2 for verification case; 3 for experiment","",
               verifcase,&verifcase,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-Creep",
               "creep closure coefficient","",
               user.Creep,&user.Creep,PETSC_NULL);CHKERRQ(ierr);
    /* FIXME  also allow control of Cmelt, omega0 */
    ierr = PetscOptionsReal("-L","half-width of square region in meters","",
               user.L,&user.L,&Lsetflg);CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);

  /* there should be an API call for this but I do not see it: */
  ierr = PetscOptionsSetValue("-snes_vi_ignore_function_sign","");CHKERRQ(ierr);

  /* configure DA */
  ierr = DMSetApplicationContext(user.da,&user);CHKERRQ(ierr);
  ierr = DMDAGetInfo(user.da,PETSC_IGNORE,&Mx,&My,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);

  if ((verifcase < 1) & (verifcase > 3)) {
    SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
             "%d is invalid verification or experiment case; nothing else is implemented\n",
             verifcase);
  }

  /* DA is compatible with geometry */
  if (!Lsetflg) {
    if (verifcase == 1) {
      user.L = 20000.0;
    } else if (verifcase == 2) {
      user.L = 250000.0;
    } else if (verifcase == 3) {
      user.L = 250000.0;
    }
  }
  ierr = DMDASetUniformCoordinates(user.da,  /* square domain */
              -user.L, user.L, -user.L, user.L, 0.0, 1.0);CHKERRQ(ierr);
  user.dx = 2.0 * user.L / (PetscReal)(Mx-1);  /* non-periodic grid */
  user.dy = 2.0 * user.L / (PetscReal)(My-1);

  /* FULLY initialize verification cases */
  if (verifcase == 1) {
    VC1Ctx vc1;
    ierr = VerifCase1Context(&user,&vc1); CHKERRQ(ierr);
    ierr = VerifCase1SetGeometry(&user,&vc1,user.b,user.H); CHKERRQ(ierr);
    ierr = VerifCase1SetHeating(&user,&vc1,user.fric,user.Ggeo); CHKERRQ(ierr);
    ierr = VerifCase1SetCurrentW(&user,&vc1,user.Wl); CHKERRQ(ierr);
    ierr = VerifCase1SetInitialpsi(&user,&vc1,x); CHKERRQ(ierr);
    ierr = VerifCase1SetInitialWnew(&user,&vc1,initialWnew); CHKERRQ(ierr);
    ierr = VerifCase1Manufacture(&user,&vc1,
               psiexact,Wnewexact,user.S,user.vbspeed); CHKERRQ(ierr);
    ierr = checkBounds(&user,psiexact,Wnewexact);CHKERRQ(ierr);
    ierr = VerifCase1SetpsiBC(&user,&vc1,user.psiBC); CHKERRQ(ierr);
    ierr = VerifCase1SetWnewBC(&user,&vc1,user.WnewBC); CHKERRQ(ierr);
  } else if (verifcase == 2) {
    VC2Ctx vc2;
    ierr = VerifCase2Context(&user,&vc2); CHKERRQ(ierr);
    ierr = VerifCase2SetGeometry(&user,&vc2,user.b,user.H); CHKERRQ(ierr);
    ierr = VerifCase2SetHeating(&user,&vc2,user.fric,user.Ggeo); CHKERRQ(ierr);
    ierr = VerifCase2SetCurrentW(&user,&vc2,user.Wl); CHKERRQ(ierr);
    ierr = VerifCase2SetInitialpsi(&user,&vc2,x); CHKERRQ(ierr);
    ierr = VerifCase2SetInitialWnew(&user,&vc2,initialWnew); CHKERRQ(ierr);
    ierr = VerifCase2Manufacture(&user,&vc2,
               psiexact,Wnewexact,user.S,user.vbspeed); CHKERRQ(ierr);
    ierr = checkBounds(&user,psiexact,Wnewexact);CHKERRQ(ierr);
    ierr = VerifCase2SetpsiBC(&user,&vc2,user.psiBC); CHKERRQ(ierr);
    ierr = VerifCase2SetWnewBC(&user,&vc2,user.WnewBC); CHKERRQ(ierr);
  } else if (verifcase == 3) {
    EC3Ctx ec3;
    ierr = VecSet(psiexact,-1.0); CHKERRQ(ierr);
    ierr = VecSet(Wnewexact,-1.0); CHKERRQ(ierr);
    ierr = ExpCase3Context(&user,&ec3); CHKERRQ(ierr);
    ierr = ExpCase3SetGeometry(&user,&ec3,user.b,user.H); CHKERRQ(ierr);
    ierr = ExpCase3SetHeating(&user,&ec3,user.fric,user.Ggeo); CHKERRQ(ierr);
    ierr = ExpCase3SetCurrentW(&user,&ec3,user.Wl); CHKERRQ(ierr);
    ierr = ExpCase3SetInitialpsi(&user,&ec3,x); CHKERRQ(ierr);
    ierr = ExpCase3SetInitialWnew(&user,&ec3,initialWnew); CHKERRQ(ierr);
    ierr = ExpCase3SetDrainageSliding(&user,&ec3,user.S,user.vbspeed); CHKERRQ(ierr);
    ierr = ExpCase3SetpsiBC(&user,&ec3,user.psiBC); CHKERRQ(ierr);
    ierr = ExpCase3SetWnewBC(&user,&ec3,user.WnewBC); CHKERRQ(ierr);
  }

  ierr = checkBounds(&user,x,user.Wl);CHKERRQ(ierr);

  /* set up the solver and Jacobian for the psi equation */
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
  ierr = SNESSetDM(snes,user.da);CHKERRQ(ierr);
  ierr = SNESSetApplicationContext(snes,&user);CHKERRQ(ierr);
  ierr = SNESSetType(snes,SNESVI);CHKERRQ(ierr);
  ierr = SNESVISetComputeVariableBounds(snes,&FormBounds);CHKERRQ(ierr);
  ierr = SNESSetFunction(snes,r,FormFunction,(void*)&user);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(PETSC_NULL,"-fd",&fdflg,PETSC_NULL);CHKERRQ(ierr);
  if (fdflg){  /* use coloring to compute finite difference Jacobian efficiently */
    ierr = DMGetColoring(user.da,IS_COLORING_GLOBAL,MATAIJ,&iscoloring);CHKERRQ(ierr);
    ierr = MatFDColoringCreate(J,iscoloring,&matfdcoloring);CHKERRQ(ierr);
    ierr = MatFDColoringSetFromOptions(matfdcoloring);CHKERRQ(ierr);
    ierr = ISColoringDestroy(&iscoloring);CHKERRQ(ierr);
    ierr = MatFDColoringSetFunction(matfdcoloring,
             (PetscErrorCode (*)(void))FormFunction,&user);CHKERRQ(ierr);
    ierr = SNESSetJacobian(snes,J,J,SNESDefaultComputeJacobianColor,
             matfdcoloring);CHKERRQ(ierr);
  } else { /* default case */
    ierr = SNESSetJacobian(snes,J,J,FormJacobian,(void*)&user); CHKERRQ(ierr);
  }
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

  /* set up the solver and Jacobian for the Wnew equation */
  ierr = SNESCreate(PETSC_COMM_WORLD,&snesUpdateW);CHKERRQ(ierr);
  ierr = SNESSetDM(snesUpdateW,user.da);CHKERRQ(ierr);
  ierr = SNESSetApplicationContext(snesUpdateW,&user);CHKERRQ(ierr);
  ierr = SNESSetType(snesUpdateW,SNESVI);CHKERRQ(ierr);
  ierr = SNESVISetComputeVariableBounds(snesUpdateW,&FormBoundsUpdateW);CHKERRQ(ierr);
  ierr = SNESSetFunction(snesUpdateW,r,FormFunctionUpdateW,(void*)&user);CHKERRQ(ierr);
  if (fdflg){  /* use coloring to compute finite difference Jacobian efficiently */
    ierr = DMGetColoring(user.da,IS_COLORING_GLOBAL,MATAIJ,&iscoloring);CHKERRQ(ierr);
    ierr = MatFDColoringCreate(J,iscoloring,&matfdcoloring);CHKERRQ(ierr);
    ierr = MatFDColoringSetFromOptions(matfdcoloring);CHKERRQ(ierr);
    ierr = ISColoringDestroy(&iscoloring);CHKERRQ(ierr);
    ierr = MatFDColoringSetFunction(matfdcoloring,
             (PetscErrorCode (*)(void))FormFunctionUpdateW,&user);CHKERRQ(ierr);
    ierr = SNESSetJacobian(snesUpdateW,J,J,SNESDefaultComputeJacobianColor,
             matfdcoloring);CHKERRQ(ierr);
  } else { /* default case */
    ierr = SNESSetJacobian(snesUpdateW,J,J,FormJacobianUpdateW,(void*)&user); CHKERRQ(ierr);
  }
  ierr = SNESSetFromOptions(snesUpdateW);CHKERRQ(ierr);

  /* report on setup */
  ierr = PetscPrintf(PETSC_COMM_WORLD,
      "setup done for verification case %d\n"
      "    domain: square of side length = %.3f km:\n"
      "    grid:   Mx,My = %d,%d;    spacing:  dx,dy = %.3f,%.3f m\n",
      verifcase, 2.0 * user.L / 1000.0, Mx, My, user.dx, user.dy); CHKERRQ(ierr);
  if (mfileflg) {
    ierr = printvecmatlab(user.da,x,"psi_init",mfile,PETSC_FALSE);CHKERRQ(ierr);
    ierr = printvecmatlab(user.da,initialWnew,"Wnew_init",mfile,PETSC_TRUE);CHKERRQ(ierr);
    ierr = printvecmatlab(user.da,user.b,"b",mfile,PETSC_TRUE);CHKERRQ(ierr);
    ierr = printvecmatlab(user.da,user.H,"H",mfile,PETSC_TRUE);CHKERRQ(ierr);
    ierr = printvecmatlab(user.da,user.S,"S",mfile,PETSC_TRUE);CHKERRQ(ierr);
    ierr = printvecmatlab(user.da,user.vbspeed,"vbspeed",mfile,PETSC_TRUE);CHKERRQ(ierr);
    if (verifcase < 3) {
      ierr = printvecmatlab(user.da,psiexact,"Pexact",mfile,PETSC_TRUE);CHKERRQ(ierr);
      ierr = printvecmatlab(user.da,Wnewexact,"Wnewexact",mfile,PETSC_TRUE);CHKERRQ(ierr);
    }
  }

  /* main time-stepping loop */
  for (n=0; n<N; n++) {
    if (N>1) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,
          "step %3d;  updating to t_{l+1} = %.3f days\n",
          n+1, (n+1) * user.dt / (3600.0*24.0));CHKERRQ(ierr);
    }

    /* solve nonlinear system for psi */
    ierr = PetscPrintf(PETSC_COMM_WORLD,"  solving for psi ... ");CHKERRQ(ierr);
    user.fcncount = 0;
    ierr = SNESSolve(snes,PETSC_NULL,x);CHKERRQ(ierr); 
    ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr);
    ierr = SNESGetConvergedReason(snes,&reason);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,
            "done: '%s'\n"
            "    number of Newton iterations = %D; function evaluations = %D\n",
            SNESConvergedReasons[reason],its,user.fcncount);CHKERRQ(ierr);

    if (n == N-1) { /* if last step */
      if (mfileflg) {
        ierr = printvecmatlab(user.da,x,"psi",mfile,PETSC_TRUE);CHKERRQ(ierr);
        /* uncomment to save ascii matrix:
        ierr = printmatmatlab(user.da,J,"J_psi",mfile,PETSC_TRUE);CHKERRQ(ierr); */
      }
      if (verifcase < 3) { /* if we know exact */
        /* compare P to exact */
        ierr = VecWAXPY(r,-1.0,psiexact,x);CHKERRQ(ierr);  /* r = P - Pexact */
        ierr = VecNorm(r,NORM_1,&error1);CHKERRQ(ierr);
        error1 /= (PetscReal)(Mx * My);
        ierr = VecNorm(r,NORM_INFINITY,&errorinf);CHKERRQ(ierr);
        ierr = VecNorm(psiexact,NORM_INFINITY,&exactinf);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,
          "errors:    av |psi-exact| = %.3e (Pa),   max |psi-exact| / max |exact| = %.3e\n",
          error1,errorinf/exactinf);CHKERRQ(ierr);
        if (mfileflg) {
          ierr = printerrorfigurematlab(user.da,1,"psi","Pexact",
                                        mfile,PETSC_TRUE);CHKERRQ(ierr);
        }
      }
    }

    /* now put x in app context as psi AND set initial condition for W = W^{l+1} */
    ierr = VecCopy(x,user.psiforW);CHKERRQ(ierr);
    if (n == 0) {
      ierr = VecCopy(initialWnew,x);CHKERRQ(ierr);
    } else {
      ierr = VecCopy(user.Wl,x);CHKERRQ(ierr);
    }

    /* solve nonlinear system for Wnew */
    ierr = PetscPrintf(PETSC_COMM_WORLD,"  solving for W_{l+1} ... ");CHKERRQ(ierr);
    user.fcncount = 0;
    ierr = SNESSolve(snesUpdateW,PETSC_NULL,x);CHKERRQ(ierr); 
    ierr = SNESGetIterationNumber(snesUpdateW,&its);CHKERRQ(ierr);
    ierr = SNESGetConvergedReason(snesUpdateW,&reason);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,
            "done: '%s'\n"
            "    number of Newton iterations = %D; function evaluations = %D\n",
            SNESConvergedReasons[reason],its,user.fcncount);CHKERRQ(ierr);

    if (n == N-1) { /* if last step */
      if (mfileflg) {
        ierr = printvecmatlab(user.da,user.Wl,"Wl",mfile,PETSC_TRUE);CHKERRQ(ierr);
        ierr = printvecmatlab(user.da,x,"Wnew",mfile,PETSC_TRUE);CHKERRQ(ierr);
        /* uncomment to save ascii matrix:
        ierr = printmatmatlab(user.da,J,"J_Wnew",mfile,PETSC_TRUE);CHKERRQ(ierr); */
      }
      if (verifcase < 3) { /* if we know exact */
        /* compare W_{l+1} to exact */
        ierr = VecWAXPY(r,-1.0,Wnewexact,x);CHKERRQ(ierr);  /* r = Wnew - Wnewexact */
        ierr = VecNorm(r,NORM_1,&error1);CHKERRQ(ierr);
        error1 /= (PetscReal)(Mx * My);
        ierr = VecNorm(r,NORM_INFINITY,&errorinf);CHKERRQ(ierr);
        ierr = VecNorm(Wnewexact,NORM_INFINITY,&exactinf);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,
          "errors:    av |Wnew-exact| = %.3e (m),   max |Wnew-exact| / max |exact| = %.3e\n",
          error1,errorinf,exactinf);CHKERRQ(ierr);
        if (mfileflg) {
          ierr = printerrorfigurematlab(user.da,2,"Wnew","Wnewexact",
                                        mfile,PETSC_TRUE);CHKERRQ(ierr);
        }
      }
    }
    
    /* Wnew --> W^l */
    ierr = VecCopy(x,user.Wl); CHKERRQ(ierr);  
    /* stored psi is best guess for next solve */
    ierr = VecCopy(user.psiforW,x); CHKERRQ(ierr);
  } /* end time stepping */

  /* free work space */
  ierr = SNESDestroy(&snes);CHKERRQ(ierr);
  ierr = SNESDestroy(&snesUpdateW);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&r);CHKERRQ(ierr);
  ierr = VecDestroy(&initialWnew);CHKERRQ(ierr);
  ierr = VecDestroy(&psiexact);CHKERRQ(ierr);
  ierr = VecDestroy(&Wnewexact);CHKERRQ(ierr);
  ierr = VecDestroy(&(user.psiforW));CHKERRQ(ierr);
  ierr = VecDestroy(&(user.Wl));CHKERRQ(ierr);
  ierr = VecDestroy(&(user.S));CHKERRQ(ierr);
  ierr = VecDestroy(&(user.vbspeed));CHKERRQ(ierr);
  ierr = VecDestroy(&(user.fric));CHKERRQ(ierr);
  ierr = VecDestroy(&(user.Ggeo));CHKERRQ(ierr);
  ierr = VecDestroy(&(user.H));CHKERRQ(ierr);
  ierr = VecDestroy(&(user.b));CHKERRQ(ierr);
  ierr = VecDestroy(&(user.psiBC));CHKERRQ(ierr);
  ierr = VecDestroy(&(user.WnewBC));CHKERRQ(ierr);

  if (fdflg) {
    ierr = MatFDColoringDestroy(&matfdcoloring);CHKERRQ(ierr);
  }
  ierr = MatDestroy(&J);CHKERRQ(ierr);
  ierr = DMDestroy(&user.da);CHKERRQ(ierr);

  ierr = PetscFinalize();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "checkBounds"
PetscErrorCode checkBounds(WPCtx *user, Vec Psi, Vec W) {
  PetscErrorCode ierr;
  PetscInt       i,j,xs,ys,xm,ym;
  PetscReal      **thk, **bed, **psi, **w, psilower, psiupper;

  PetscFunctionBegin;
  ierr = DMDAGetCorners(user->da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, Psi, &psi);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, W, &w);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, user->H, &thk);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, user->b, &bed);CHKERRQ(ierr);
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      if (w[j][i] < 0.0) {
        SETERRQ4(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
          "negative water thickness W = %.6e at (i,j)=(%d,%d) during function eval %d",
          w[j][i],i,j,user->fcncount);
      }
      psilower = user->rhow * user->g * bed[j][i];
      psiupper = user->rhoi * user->g * thk[j][i] + psilower;
      if (psi[j][i] < psilower) {
        SETERRQ4(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
          "potential psi = %.6e too low (negative pressure) at (i,j)=(%d,%d)\n"
          "during function eval %d",
          psi[j][i],i,j,user->fcncount);
      }
      if (psi[j][i] > psiupper) {
        SETERRQ4(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
          "potential psi = %.6e too high (pressure exceeds overburden) at (i,j)=(%d,%d)\n"
          "during function eval %d",
          psi[j][i],i,j,user->fcncount);
      }
    }
  }
  ierr = DMDAVecRestoreArray(user->da, user->H, &thk);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, user->b, &bed);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, Psi, &psi);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, W, &w);CHKERRQ(ierr);
  PetscFunctionReturn(0);
} 


#undef __FUNCT__
#define __FUNCT__ "FormBounds"
/*  FormBounds() for call-back: tell SNESVI we want
   rhow g b  <=  psi  <=  p_i + rhow g b
where p_i = rhoi g H is overburden pressure */
PetscErrorCode FormBounds(SNES snes, Vec Xl, Vec Xu) {
  PetscErrorCode ierr;
  WPCtx          *user;
  PetscFunctionBegin;
  ierr = SNESGetApplicationContext(snes,&user);CHKERRQ(ierr);
  /* Xl = rhow g b */
  ierr = VecCopy(user->b,Xl);CHKERRQ(ierr);
  ierr = VecScale(Xl, user->rhow * user->g);CHKERRQ(ierr);
  /* Xu = rhoi g H + rhow g b */
  ierr = VecCopy(user->H,Xu);CHKERRQ(ierr);
  ierr = VecScale(Xu, user->rhoi * user->g);CHKERRQ(ierr);
  ierr = VecAXPY(Xu,1.0,Xl);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "FormFunction"
/* FormFunction - Evaluates nonlinear function, F(X), on local process patch,
for psi equation */
PetscErrorCode FormFunction(SNES snes, Vec X, Vec F, void *ptr) {
  PetscErrorCode ierr;
  WPCtx          *user = (WPCtx*)ptr;
  PetscInt       i,j,xs,ys,xm,ym,Mx,My;
  PetscReal      dx = user->dx, dy = user->dy, Pover, psib,
                 c0, c1, c2, c3, tom0, scale,
                 psieast, psiwest, psinorth, psisouth,
                 divQ, dpsidx, dpsidy, dpsisqr, Melt, OC,
                 Wij, psiij, east, west, north, south,
                 **w, **psi, **f, **h, **bed,
                 **vb, **s, **fricheat, **geoheat, **bc;
  Vec            localpsi, localW, localH;

  PetscFunctionBegin;
  user->fcncount = user->fcncount + 1;
/*  ierr = checkBounds(user,X,user->Wl); CHKERRQ(ierr);  FIXME: fails with -fd; restore? */

  ierr = DMDAGetCorners(user->da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetInfo(user->da,PETSC_IGNORE,&Mx,&My,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);

  ierr = DMGetLocalVector(user->da,&localpsi); CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(user->da,X,INSERT_VALUES,localpsi); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(user->da,X,INSERT_VALUES,localpsi); CHKERRQ(ierr);

  ierr = DMGetLocalVector(user->da,&localW); CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(user->da,user->Wl,INSERT_VALUES,localW); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(user->da,user->Wl,INSERT_VALUES,localW); CHKERRQ(ierr);

  ierr = DMGetLocalVector(user->da,&localH); CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(user->da,user->H,INSERT_VALUES,localH); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(user->da,user->H,INSERT_VALUES,localH); CHKERRQ(ierr);

  ierr = DMDAVecGetArray(user->da, user->b, &bed);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, user->vbspeed, &vb);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, user->fric, &fricheat);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, user->Ggeo, &geoheat);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, user->S, &s);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, user->psiBC, &bc);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, localW, &w);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, localpsi, &psi);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, localH, &h);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, F, &f);CHKERRQ(ierr);
  c0    = user->K0 / (user->rhow * user->g);
  c1    = (1.0 / user->rhoi) - (1.0 / user->rhow);
  c2    = user->Creep * user->Aglen;
  c3    = user->Cmelt * c0 / user->Lfusion;
  tom0  = 2.0 * user->omega0;
  scale = dx * dy / c0;  /* make entries in Jacobian close to O(1) */
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      if (i == 0 || j == 0 || i == Mx-1 || j == My-1) {
        /* along bdry of computational domain set Dirichlet b.c. */
        f[j][i] = psi[j][i] - bc[j][i];
      } else {
        if (h[j][i] < 1.0) {
          /* pressure is atmospheric if ice is absent (even if marine ... the ice is absent!) */
          psib = user->rhow * user->g * bed[j][i];
          f[j][i] = psi[j][i] - psib;
        } else {
          Wij = w[j][i];
          psiij = psi[j][i];
          psib = user->rhow * user->g * bed[j][i];
          Pover = user->rhoi * user->g * h[j][i];          
          psieast  = (h[j][i+1] < 1.0) ? (Pover + psib) : psi[j][i+1];
          psiwest  = (h[j][i-1] < 1.0) ? (Pover + psib) : psi[j][i-1];
          psinorth = (h[j+1][i] < 1.0) ? (Pover + psib) : psi[j+1][i];
          psisouth = (h[j-1][i] < 1.0) ? (Pover + psib) : psi[j-1][i];
          /* Melt(x,y) = (1/L) (-tau_b . v_b + G) + (Cmelt c0 / L) x y
             Melt = Melt(Wl, dpsisqr)                                      */
          dpsidx = (psieast - psiwest) / (2.0 * dx);
          dpsidy = (psinorth - psisouth) / (2.0 * dy);
          dpsisqr = dpsidx * dpsidx + dpsidy * dpsidy;
          Melt = (fricheat[j][i] + geoheat[j][i]) / user->Lfusion + c3 * Wij * dpsisqr;
          /* OC(x,y) = Cavit |v_b| - Creep A x^n y
             OC = OC(psi_i - psi, Wl)                                   */
          OC = user->Cavit * vb[j][i] - c2 * pow(Pover + psib - psiij,user->nglen) * Wij;
          /* c0 div(Wl grad psi) - c1 Melt - OC + S/rhow = 0 */
          east  = 0.5 * (w[j][i+1] + Wij + tom0) * (psieast - psiij) / dx;
          west  = 0.5 * (Wij + w[j][i-1] + tom0) * (psiij - psiwest) / dx;
          north = 0.5 * (w[j+1][i] + Wij + tom0) * (psinorth - psiij) / dy;
          south = 0.5 * (Wij + w[j-1][i] + tom0) * (psiij - psisouth) / dy;
          divQ  = c0 * ((east - west) / dx + (north - south) / dy);
          f[j][i] = scale * (divQ - c1 * Melt - OC + s[j][i] / user->rhow);
        }
      }
    }
  }
  ierr = DMDAVecRestoreArray(user->da, user->b, &bed);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, user->vbspeed, &vb);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, user->fric, &fricheat);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, user->Ggeo, &geoheat);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, user->S, &s);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, user->psiBC, &bc);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, localW, &w);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, localpsi, &psi);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, localH, &h);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, F, &f);CHKERRQ(ierr);

  ierr = DMRestoreLocalVector(user->da, &localpsi);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(user->da, &localW);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(user->da, &localH);CHKERRQ(ierr);

  /*ierr = PetscLogFlops(XX * info->ym * info->xm);CHKERRQ(ierr);*/
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "FormJacobian"
/* FormJacobian - Evaluates Jacobian matrix on local process patch */
PetscErrorCode FormJacobian(SNES snes,Vec X,Mat *J,Mat *B,MatStructure *flag,
                            void *ptr) {
  SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER,"FormJacobian() is not yet implemented");
  PetscFunctionReturn(0); 
} 


#undef __FUNCT__
#define __FUNCT__ "FormBoundsUpdateW"
/*  FormBoundsUpdateW() for call-back: tell SNESVI we want 0 <= W_{l+1} */
PetscErrorCode FormBoundsUpdateW(SNES snes, Vec Xl, Vec Xu) {
  PetscErrorCode ierr;
  WPCtx          *user;
  PetscFunctionBegin;
  ierr = SNESGetApplicationContext(snes,&user);CHKERRQ(ierr);
  ierr = VecSet(Xl,0.0);CHKERRQ(ierr);         /* W >= 0 */
  ierr = VecSet(Xu,SNES_VI_INF);CHKERRQ(ierr); /* no upper bound */
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "FormFunctionUpdateW"
/* FormFunctionUpdateW() computes the residual F(X) for the equation we solve
   to update W, computing X = W_{l+1};  needs psi, Wl, S */
PetscErrorCode FormFunctionUpdateW(SNES snes, Vec X, Vec F, void *ptr) {
  PetscErrorCode ierr;
  WPCtx          *user = (WPCtx*)ptr;
  DM             da = user->da;
  PetscInt       i,j,xs,ys,xm,ym,Mx,My;
  PetscReal      dx = user->dx, dy = user->dy, dt = user->dt, nux, nuy,
                 c0, c3, divQ, dpsidx, dpsidy, dpsisqr, Melt,
                 Wij, psiij, east, west, north, south,
                 vWeast, vWwest, vWnorth, vWsouth,
                 **w, **psi, **fricheat, **geoheat, **s, **wold, **h, **f, **bc;
  Vec            localpsi, localW, localH;

  PetscFunctionBegin;
  user->fcncount = user->fcncount + 1;
  ierr = DMDAGetCorners(da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetInfo(da,PETSC_IGNORE,&Mx,&My,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);

  ierr = DMGetLocalVector(da,&localpsi); CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(da,user->psiforW,INSERT_VALUES,localpsi); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da,user->psiforW,INSERT_VALUES,localpsi); CHKERRQ(ierr);

  ierr = DMGetLocalVector(da,&localW); CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(da,X,INSERT_VALUES,localW); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da,X,INSERT_VALUES,localW); CHKERRQ(ierr);

  ierr = DMGetLocalVector(da,&localH); CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(da,user->H,INSERT_VALUES,localH); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da,user->H,INSERT_VALUES,localH); CHKERRQ(ierr);

  nux = dt / dx;
  nuy = dt / dy;
  c0  = user->K0 / (user->rhow * user->g);
  c3  = user->Cmelt * c0 / user->Lfusion;

  ierr = DMDAVecGetArray(da, localpsi, &psi);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da, localW, &w);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da, localH, &h);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da, user->fric, &fricheat);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da, user->Ggeo, &geoheat);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da, user->S, &s);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da, user->Wl, &wold);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da, user->WnewBC, &bc);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da, F, &f);CHKERRQ(ierr);
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      if (i == 0 || j == 0 || i == Mx-1 || j == My-1) {
        f[j][i] = w[j][i] - bc[j][i];
      } else {
        /* FIXME:  more thought needed on inflow b.c. */
        if (h[j][i] < 1.0) {
          f[j][i] = w[j][i] - 0.0;  /* no water if ice is absent */
        /* ALTERNATE VERSION: generates NaN???
          if ( (h[j][i] < 1.0) || (h[j][i+1] < 1.0) || (h[j][i-1] < 1.0)
               || (h[j+1][i] < 1.0) || (h[j-1][i] < 1.0) )  {
            f[j][i] = w[j][i] - 0.0;   no water if next to ice-free
        */
        } else {
          Wij = w[j][i];
          psiij = psi[j][i];
          /* Melt(a,b) = (1/L) (-tau_b . v_b + G) + (Cmelt c0 / L) a b
             Melt = Melt(W, dpsisqr)                                      */
          dpsidx = (psi[j][i+1] - psi[j][i-1]) / (2.0 * dx);
          dpsidy = (psi[j+1][i] - psi[j-1][i]) / (2.0 * dy);
          dpsisqr = dpsidx * dpsidx + dpsidy * dpsidy;
          Melt = (fricheat[j][i] + geoheat[j][i]) / user->Lfusion + c3 * Wij * dpsisqr;
          /* think finite volume: east,...,south are velocities at
               midpoints of boundaries of (i,j)-centered-cell */
          east  = - c0 * (psi[j][i+1] - psiij) / dx;
          west  = - c0 * (psiij - psi[j][i-1]) / dx;
          north = - c0 * (psi[j+1][i] - psiij) / dy;
          south = - c0 * (psiij - psi[j-1][i]) / dy;
          vWeast  = (east  >= 0.0) ? east  * Wij       : east  * w[j][i+1];
          vWwest  = (west  >= 0.0) ? west  * w[j][i-1] : west  * Wij;
          vWnorth = (north >= 0.0) ? north * Wij       : north * w[j+1][i];
          vWsouth = (south >= 0.0) ? south * w[j-1][i] : south * Wij;
          divQ  = nux * (vWeast - vWwest) + nuy * (vWnorth - vWsouth);
          /* no scaling required because Jacobain already close to identity when dt small */
          f[j][i] = Wij - wold[j][i] + divQ - (dt / user->rhow) * (Melt + s[j][i]);
        }
      }
    }
  }
  ierr = DMDAVecRestoreArray(da, localpsi, &psi);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da, localW, &w);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da, localH, &h);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da, user->fric, &fricheat);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da, user->Ggeo, &geoheat);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da, user->S, &s);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da, user->Wl, &wold);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da, user->WnewBC, &bc);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da, F, &f);CHKERRQ(ierr);

  ierr = DMRestoreLocalVector(da, &localpsi);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(da, &localW);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(da, &localH);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "FormJacobianUpdateW"
/* FormJacobianUpdateW() - Evaluates Jacobian matrix on local process patch */
PetscErrorCode FormJacobianUpdateW(SNES snes,Vec X,Mat *J,Mat *B,MatStructure *flag,
                            void *ptr) {
  SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER,
          "FormJacobianUpdateW() is not yet implemented");
  PetscFunctionReturn(0); 
} 

