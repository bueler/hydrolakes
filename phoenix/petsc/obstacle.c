
static char help[] = "Solves obstacle problem in 2D to demonstrate this\n"
"functionality of SNESVI.  An elliptic problem with solution  u  constrained\n"
"to be above a given function  psi.  Exact solution is known.\n";

/* 
Solve on a square  R = {-2<x<2,-2<y<2}:

    0 = u_xx + u_yy   on the set {(x,y) | u(x,y) > psi(x,y)},

where psi is the upper hemisphere of the unit ball, and we require
u >= psi on all of R.  On the boundary of R we have nonhomogenous
Dirichlet boundary conditions coming from the exact solution.

The exact solution is known for the given psi and boundary values in
question.  See http://www.dms.uaf.edu/~bueler/obstacleDOC.pdf

Example usage follows.

Get help:
  ./obstacle -help

Parallel runs, spatial refinement only:
  for M in 21 41 81 161 321; do
    echo "case M=$M:"
    mpiexec -n 4 ./obstacle -da_grid_x $M -da_grid_y $M -snes_monitor
  done

With finite difference evaluation of Jacobian using coloring:
  ./obstacle -fd

Ascii Matlab-format output:
  ./obstacle -mfile foo.m
*/

#include <petscdmda.h>
#include <petscsnes.h>

#include "matlabprint.h"  /* utilities for putting petsc objects in .m file */


/* application context for obstacle problem solver */
typedef struct {
   DM          da;
   Vec         psi, uexact;
} ObsCtx;


extern PetscErrorCode FormPsiAndInitialGuess(ObsCtx*,Vec);
extern PetscErrorCode FormBounds(SNES,Vec,Vec);
extern PetscErrorCode FormFunctionLocal(DMDALocalInfo*,
                        PetscScalar**,PetscScalar**,ObsCtx*);
extern PetscErrorCode FormJacobianLocal(DMDALocalInfo*,PetscScalar**,Mat,ObsCtx*);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  SNES snes;
  Vec            u, r;   /* solution, residual vector */
  PetscInt       Mx,My,its;
  SNESConvergedReason reason;
  DM             da;
  ObsCtx         user;
  PetscReal      dx,dy,error1,errorinf;
  PetscBool      fdflg = PETSC_FALSE, mfileflg = PETSC_FALSE;
  char           mfile[PETSC_MAX_PATH_LEN] = "obout.m";

  PetscInitialize(&argc,&argv,(char *)0,help);

  ierr = DMDACreate2d(PETSC_COMM_WORLD,
             DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_NONE,
             DMDA_STENCIL_STAR, // nonlinear diffusion but diffusivity
                                //   depends on soln W not grad W
             -11,-11,           // default to 10x10 grid but override with
                                //   -da_grid_x, -da_grid_y (or -da_refine)
             PETSC_DECIDE,PETSC_DECIDE, // num of procs in each dim
             1,1,               // dof = 1, s = 1 (stencil extends out one cell)
             PETSC_NULL,PETSC_NULL, // no specify proc decomposition
             &da);CHKERRQ(ierr);
  user.da = da;

  ierr = DMCreateGlobalVector(da,&u);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&r);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&(user.uexact));CHKERRQ(ierr);
  ierr = VecDuplicate(u,&(user.psi));CHKERRQ(ierr);

  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,
           "","options to obstacle problem","");
  {
    ierr = PetscOptionsString("-mfile",
               "name of Matlab file to write results","",
               mfile,mfile,PETSC_MAX_PATH_LEN,&mfileflg);CHKERRQ(ierr);
    ierr = PetscOptionsBool("-fd",
               "use coloring to compute Jacobian by finite differences",
               PETSC_NULL,fdflg,&fdflg,PETSC_NULL);CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  /* there should be an API call for this but I do not see it: */
  ierr = PetscOptionsSetValue("-snes_vi_ignore_function_sign","");CHKERRQ(ierr);

  ierr = DMDASetUniformCoordinates(da,  // square domain
              -2.0, 2.0, -2.0, 2.0, 0.0, 1.0);CHKERRQ(ierr);
  ierr = DMSetApplicationContext(da,&user);CHKERRQ(ierr);

  ierr = FormPsiAndInitialGuess(&user,u);CHKERRQ(ierr);

  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
  ierr = SNESSetDM(snes,da);CHKERRQ(ierr);
  ierr = SNESSetApplicationContext(snes,&user);CHKERRQ(ierr);
  ierr = SNESSetType(snes,SNESVI);CHKERRQ(ierr);
  ierr = SNESVISetComputeVariableBounds(snes,&FormBounds);CHKERRQ(ierr);

  ierr = DMDASetLocalFunction(da,(DMDALocalFunction1)FormFunctionLocal);
             CHKERRQ(ierr);
  if (!fdflg) {
    ierr = DMDASetLocalJacobian(da,(DMDALocalFunction1)FormJacobianLocal);
               CHKERRQ(ierr); 
  }

  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

  /* report on setup */
  ierr = DMDAGetInfo(da,PETSC_IGNORE,&Mx,&My,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);
  dx = 4.0 / (PetscReal)(Mx-1);
  dy = 4.0 / (PetscReal)(My-1);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
    "setup done: square       side length = %.3f\n"
    "            grid               Mx,My = %d,%d\n"
    "            spacing            dx,dy = %.3f,%.3f\n",
    4.0, Mx, My, dx, dy); CHKERRQ(ierr);

  if (mfileflg) {
    ierr = printvecmatlab(da,user.psi,"psi",mfile,PETSC_FALSE);CHKERRQ(ierr);
    ierr = printvecmatlab(da,u,"u0",mfile,PETSC_TRUE);CHKERRQ(ierr);
    ierr = printvecmatlab(da,user.uexact,"uexact",mfile,PETSC_TRUE);CHKERRQ(ierr);
  }

  /* solve nonlinear system */
  ierr = SNESSolve(snes,PETSC_NULL,u);CHKERRQ(ierr); 
  ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr);
  ierr = SNESGetConvergedReason(snes,&reason);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Newton iterations = %D; result = %s\n",
            its,SNESConvergedReasons[reason]);CHKERRQ(ierr);

  if (mfileflg) {
    ierr = printvecmatlab(da,u,"u",mfile,PETSC_TRUE);CHKERRQ(ierr);
  }

  /* compare to exact */
  ierr = VecWAXPY(r,-1.0,user.uexact,u);CHKERRQ(ierr);  /* r = W - Wexact */
  ierr = VecNorm(r,NORM_1,&error1);CHKERRQ(ierr);
  error1 /= (PetscReal)Mx * (PetscReal)My;
  ierr = VecNorm(r,NORM_INFINITY,&errorinf);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
    "errors:    av |u-uexact|  = %.3e\n"
    "           |u-uexact|_inf = %.3e\n",
    error1,errorinf);CHKERRQ(ierr);

  /* Free work space.  */
  ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&r);CHKERRQ(ierr);
  ierr = VecDestroy(&(user.psi));CHKERRQ(ierr);
  ierr = VecDestroy(&(user.uexact));CHKERRQ(ierr);

  ierr = SNESDestroy(&snes);CHKERRQ(ierr);
  ierr = DMDestroy(&da);CHKERRQ(ierr);

  ierr = PetscFinalize();CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "FormPsiAndInitialGuess"
PetscErrorCode FormPsiAndInitialGuess(ObsCtx *user, Vec U0) {
  PetscErrorCode ierr;
  PetscInt       i,j,Mx,My,xs,ys,xm,ym;
  DM             coordDA;
  Vec            coordinates;
  DMDACoor2d     **coords;
  PetscReal      **psi, **u0, **uexact,
                 x, y, r,
                 afree = 0.69797, A = 0.68026, B = 0.47152, pi = 3.1415926536;

  PetscFunctionBegin;
  ierr = DMDAGetInfo(user->da,PETSC_IGNORE,&Mx,&My,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);
  ierr = DMDAGetCorners(user->da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL);CHKERRQ(ierr);

  ierr = DMDAGetCoordinateDA(user->da, &coordDA);CHKERRQ(ierr);
  ierr = DMDAGetCoordinates(user->da, &coordinates);CHKERRQ(ierr);

  ierr = DMDAVecGetArray(coordDA, coordinates, &coords);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, user->psi, &psi);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, U0, &u0);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, user->uexact, &uexact);CHKERRQ(ierr);
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      x = coords[j][i].x;
      y = coords[j][i].y;
      r = sqrt(x * x + y * y);
      if (r <= 1.0) {
        psi[j][i] = sqrt(1.0 - r * r);
      } else {
        psi[j][i] = -1.0;        
      }
      if (r <= afree) {
        uexact[j][i] = psi[j][i];  /* on the obstacle */
      } else {
        uexact[j][i] = - A * log(r) + B;   /* solves the laplace eqn */
      }
      if (i == 0 || j == 0 || i == Mx-1 || j == My-1) {
        u0[j][i] = uexact[j][i];
      } else {
        /* initial guess is admissible: it is above the obstacle */
        u0[j][i] = uexact[j][i] + cos(pi * x / 4) * cos(pi * y / 4);
      }
    }
  }
  ierr = DMDAVecRestoreArray(user->da, user->psi, &psi);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, U0, &u0);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, user->uexact, &uexact);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(coordDA, coordinates, &coords);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "FormBounds"
/*  FormBounds() for call-back: tell SNESVI (variational inequality)
  that we want u >= psi */
PetscErrorCode FormBounds(SNES snes, Vec Xl, Vec Xu) {
  PetscErrorCode ierr;
  ObsCtx         *user;
  PetscFunctionBegin;
  ierr = SNESGetApplicationContext(snes,&user);CHKERRQ(ierr);
  ierr = VecCopy(user->psi,Xl);CHKERRQ(ierr);  /* u >= psi */
  ierr = VecSet(Xu,SNES_VI_INF);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "FormFunctionLocal"
/* FormFunctionLocal - Evaluates nonlinear function, F(x) on local process patch */
PetscErrorCode FormFunctionLocal(DMDALocalInfo *info,
                                 PetscScalar **x,PetscScalar **f,ObsCtx *user) {
  PetscErrorCode ierr;
  PetscInt       i,j;
  PetscReal      dx,dy,uxx,uyy,
                 **uexact;  /* for boundary values only */

  PetscFunctionBegin;
  dx = 4.0 / (PetscReal)(info->mx-1);
  dy = 4.0 / (PetscReal)(info->my-1);

  ierr = DMDAVecGetArray(user->da, user->uexact, &uexact);CHKERRQ(ierr);
  for (j=info->ys; j<info->ys+info->ym; j++) {
    for (i=info->xs; i<info->xs+info->xm; i++) {
      if (i == 0 || j == 0 || i == info->mx-1 || j == info->my-1) {
        f[j][i] = x[j][i] - uexact[j][i];
      } else {
        uxx     = (x[j][i-1] - 2.0 * x[j][i] + x[j][i+1]) / (dx*dx);
        uyy     = (x[j-1][i] - 2.0 * x[j][i] + x[j+1][i]) / (dy*dy);
        f[j][i] = uxx + uyy;
      }
    }
  }
  ierr = DMDAVecRestoreArray(user->da, user->uexact, &uexact);CHKERRQ(ierr);

  ierr = PetscLogFlops(10.0*info->ym*info->xm);CHKERRQ(ierr);
  PetscFunctionReturn(0); 
} 


#undef __FUNCT__
#define __FUNCT__ "FormJacobianLocal"
/* FormJacobianLocal - Evaluates Jacobian matrix on local process patch */
PetscErrorCode FormJacobianLocal(DMDALocalInfo *info,PetscScalar **x,Mat jac,
                                 ObsCtx *user) {
  PetscErrorCode ierr;
  PetscInt       i,j;
  MatStencil     col[5],row;
  PetscReal      v[5],dx,dy,oxx,oyy;

  PetscFunctionBegin;
  dx = 4.0 / (PetscReal)(info->mx-1);
  dy = 4.0 / (PetscReal)(info->my-1);
  oxx = 1.0 / (dx * dx);
  oyy = 1.0 / (dy * dy);

  for (j=info->ys; j<info->ys+info->ym; j++) {
    for (i=info->xs; i<info->xs+info->xm; i++) {
      row.j = j; row.i = i;
      if (i == 0 || j == 0 || i == info->mx-1 || j == info->my-1) { /* boundary */
        v[0] = 1.0;
        ierr = MatSetValuesStencil(jac,1,&row,1,&row,v,INSERT_VALUES);CHKERRQ(ierr);
      } else { /* interior grid points */
        v[0] = oyy;                 col[0].j = j - 1;  col[0].i = i;
        v[1] = oxx;                 col[1].j = j;      col[1].i = i - 1;
        v[2] = -2.0 * (oxx + oyy);  col[2].j = j;      col[2].i = i;
        v[3] = oxx;                 col[3].j = j;      col[3].i = i + 1;
        v[4] = oyy;                 col[4].j = j + 1;  col[4].i = i;
        ierr = MatSetValuesStencil(jac,1,&row,5,col,v,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
  }

  /* Assemble matrix, using the 2-step process: */
  ierr = MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  /* Tell the matrix we will never add a new nonzero location to the
     matrix. If we do, it will generate an error. */
  ierr = MatSetOption(jac,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);

  ierr = PetscLogFlops(2.0*info->ym*info->xm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

