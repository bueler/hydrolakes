static char help[] = "Solves a single Crank-Nicolson timestep\n\
for the vertically-integrated porous medium\n\
equation in 2d, namely the degenerate diffusion equation\n\
    W_t = gamma div ( W grad (W^sigma) )\n\
where gamma > 0 and sigma > 0.\n\
TODO: Verification using the Barenblatt similarity solution.\n\n";


/* try this kind of thing:
   mpiexec -n 4 ./onestep -vi -snes_monitor -ksp_monitor

situation:
  ./onestep -vi -da_grid_x 31 -da_grid_y 31   // default; works
  ./onestep -vi -da_grid_x 51 -da_grid_y 51   // works; W has tiny neg entries
  ./onestep -vi -da_grid_x 101 -da_grid_y 101 // solve fails (DIVERGED_LINE_SEARCH)
  ./onestep -vi -da_grid_x 101 -da_grid_y 101 -por_dtdays 10 // succeeds with shorter time-step

view sparsity:
  ./onestep -vi -mat_view_draw -draw_pause 0.5 -da_grid_x 5 -da_grid_y 5
*/

/* ------------------------------------------------------------------------

Consider this vertically-integrated porous medium equation (VPME)

    W_t = gamma div ( W grad (W^sigma) )

where  W = W(t,x,y)  is a water depth and

    gamma = (K0 H0 rhoi) / (rhow Wcrit^sigma).  

Let  t = t_n  be some time and suppose  W_n = W(t,x,y)  is known.  Let

    W(x,y) be an approximation to W_{n+1} = W(t+dt,x,y)

One implicit time-step by the Crank-Nicolson (trapezoid) rule requires
us to solve this equation for  W:

    W - W_n   gamma                              gamma
    ------- = ----- div ( W grad ( W^sigma ) ) + ----- div ( W_n grad ( W_n^sigma ) )
      dt        2                                  2

Let  beta = gamma dt / 2.  Rearranging as an elliptic equation for  W,

(*)    W - beta div (W grad ( W^sigma ) ) = f

where

    f = W_n + beta div (W_n grad ( W_n^sigma ) )

This code solves equation (*) and it gets  W_n,  to compute  f,  from the
Barenblatt solution.  We check our result our result by comparing the
computed  W  to the Barenblatt solution  W(t+dt).

FIXME:  measure difference with Barenblatt and see/measure O(dt^2)?, O(dx^z)?
convergence
 
------------------------------------------------------------------------- */

#include <petscdmda.h>
#include <petscsnes.h>

#include "../matlabprint.h"  /* utilities for putting petsc objects in .m file */
#include "../context.h"      /* define and set defaults for parameters */
#include "../barenblatt.h"

/* Application-provided routines (call-back) */
extern PetscErrorCode FormSource(DM,PorousCtx*,Vec,Vec);
extern PetscErrorCode FormFunctionLocal(DMDALocalInfo*,
                        PetscScalar**,PetscScalar**,PorousCtx*);
extern PetscErrorCode FormJacobianLocal(DMDALocalInfo*,
                        PetscScalar**,Mat,PorousCtx*);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  SNES                   snes;                 /* nonlinear solver */
  Vec                    x;                    /* solution vector */
  PorousCtx              user;                 /* user-defined work context */
  PetscInt               its;                  /* iterations for convergence */
  SNESConvergedReason    reason;               /* reason for (con/di)vergence */
  PetscErrorCode         ierr;
  PetscBool              fdflg = PETSC_FALSE, viflg = PETSC_FALSE,
                         mfileflg = PETSC_FALSE, optflg = PETSC_FALSE;
  char                   mfile[PETSC_MAX_PATH_LEN] = "porout.m";
  DM                     da;
  PetscInt    Mx, My;
  PetscReal   dx, maxD, stabilitydt, secperday = 3600.0 * 24.0;

  PetscInitialize(&argc,&argv,(char *)0,help);

  user.da = da;
  ierr = DefaultContext(&user);CHKERRQ(ierr);
  
  /* initialize problem parameters which correspond to options */
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,
           "","options to porous medium equation","");
  {
    ierr = PetscOptionsReal("-por_sigma","nonlinear power","",
                            user.sigma,&user.sigma,PETSC_NULL);CHKERRQ(ierr);
    /* default region is 40km x 40km: */
    ierr = PetscOptionsReal("-por_L","half-width of square region","",
                            user.L,&user.L,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-por_tsdays","start time in days","",
                            user.tstart/secperday,&user.tstart,&optflg);CHKERRQ(ierr);
    if (optflg) { user.tstart *= secperday; }
    ierr = PetscOptionsReal("-por_dtdays","duration of time step in days","",
                            user.dt/secperday,&user.dt,&optflg);CHKERRQ(ierr);
    if (optflg) { user.dt *= secperday; }
    ierr = PetscOptionsString("-mfile",
                            "name of Matlab file to write results","",
                            mfile,mfile,PETSC_MAX_PATH_LEN,&mfileflg);
                            CHKERRQ(ierr);
    ierr = PetscOptionsBool("-fd","use finite difference evaluation of jacobian","",
                            fdflg,&fdflg,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsBool("-vi","use SNESVI for variational inequality","",
                            viflg,&viflg,PETSC_NULL);CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);

  ierr = DerivedConstants(&user);CHKERRQ(ierr);

  /* now do PETSc setup */
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);

  /* Create distributed array (DMDA) to manage parallel grid and vectors */
  ierr = DMDACreate2d(PETSC_COMM_WORLD,
             DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_NONE, // correct for zero Dirichlet
             DMDA_STENCIL_STAR, // nonlinear diffusion but diffusivity
                                //   depends on soln W not grad W
             -31,-31,           // default to 30x30 grid but override with
                                //   -da_grid_x, -da_grid_y (or -da_refine)
             PETSC_DECIDE,PETSC_DECIDE, // num of procs in each dim
             1,1,               // dof = 1, s = 1 (stencil extends out one cell)
             PETSC_NULL,PETSC_NULL, // no specify proc decomposition
             &da);CHKERRQ(ierr);
  ierr = DMDASetUniformCoordinates(da,  // square domain
              -user.L, user.L, -user.L, user.L, 0.0, 1.0);CHKERRQ(ierr);
  ierr = DMSetApplicationContext(da,&user);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(da,&x);CHKERRQ(ierr);  // x = solution
  ierr = VecDuplicate(x,&(user.f));CHKERRQ(ierr);

  ierr = DMDASetLocalFunction(da,(DMDALocalFunction1)FormFunctionLocal);
           CHKERRQ(ierr);
  if (!fdflg) {
    ierr = DMDASetLocalJacobian(da,(DMDALocalFunction1)FormJacobianLocal);
        CHKERRQ(ierr); 
  }

  ierr = SNESSetDM(snes,da);CHKERRQ(ierr);

  /* supply upper and lower bounds to SNESVI type */
  /* The following is a workaround to fix a bug in petsc3.2; see
     http://petsc.cs.iit.edu/petsc/petsc-dev/rev/2d5e73b75a6c#l1.2
  We would like to use SNESVISetVariableBounds, but use
  a call-back instead to avoid some error checking issue. */
  if (viflg) {
    ierr = SNESVISetComputeVariableBounds(snes,FormPositivityBounds);
        CHKERRQ(ierr);
  }

  /* Customize nonlinear solver; set runtime options */
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

  ierr = DMDAGetInfo(da,PETSC_IGNORE,&Mx,&My,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);
  dx = 2.0 * user.L / (Mx-1);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
    "  setup done: square of side length %.3f km\n"
    "              spacing  dx = %.3f m\n"
    "              initial time  tstart = %.3f days\n"
    "              time step  dt = %.3f days\n",
    2.0 * user.L / 1000.0, dx,
    user.tstart / secperday, (user.dt) / secperday);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,
    "  forming initial state Wn from Barenblatt soln ...\n");CHKERRQ(ierr);
  ierr = BarenblattState(da,&user,user.tstart,x);CHKERRQ(ierr);
  ierr = FormSource(da,&user,x,user.f);CHKERRQ(ierr);

  ierr = maxDiffusivity(da,&user,x,&maxD);CHKERRQ(ierr);
  stabilitydt = dx*dx/maxD;
  ierr = PetscPrintf(PETSC_COMM_WORLD,
    "  state Wn:   max D = %f m2s-1;\n"
    "              explicit stability dt = %.5f days\n"
    "              time-step is %.2f times initial stability\n",
    maxD, stabilitydt / secperday,
    user.dt / stabilitydt);CHKERRQ(ierr);

  if (mfileflg) {
    ierr = printvecmatlab(da,x,"Wn",mfile,PETSC_FALSE);CHKERRQ(ierr);
    ierr = printvecmatlab(da,user.f,"f",mfile,PETSC_TRUE);CHKERRQ(ierr);
  }

  ierr = PetscPrintf(PETSC_COMM_WORLD,
    "  solving for final state W at time t = %.3f days ...\n",
    (user.tstart+user.dt) / secperday);CHKERRQ(ierr);

  /* Solve nonlinear system */
  ierr = SNESSolve(snes,PETSC_NULL,x);CHKERRQ(ierr);  
  ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr);
  ierr = SNESGetConvergedReason(snes,&reason);CHKERRQ(ierr);
  if (reason > 0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
      "  solve SUCCEEDED!\n");CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
      "ERROR:  solve FAILED with reason = %s\n",
      SNESConvergedReasons[reason]);CHKERRQ(ierr);
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,
      "  did %d Newton iterations and %d function evaluations\n",
      its,user.fcncount);CHKERRQ(ierr);
  if (mfileflg) {
    ierr = printvecmatlab(da,x,"W",mfile,PETSC_TRUE);CHKERRQ(ierr);
    ierr = printfigurematlab(da,"Wn","W",mfile);CHKERRQ(ierr);
  }

  ierr = maxDiffusivity(da,&user,x,&maxD);CHKERRQ(ierr);
  stabilitydt = dx*dx/maxD;
  ierr = PetscPrintf(PETSC_COMM_WORLD,
    "  state W:    max D = %f m2s-1;\n"
    "              explicit stability dt = %.5f days\n"
    "              time-step is %.2f times final stability\n",
    maxD, stabilitydt/secperday,
    user.dt / stabilitydt);CHKERRQ(ierr);

  if (mfileflg) {
    KSP         ksp;
    Mat         J;
    ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
    ierr = KSPGetOperators(ksp,&J,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
    ierr = printmatmatlab(da,J,"J",mfile,PETSC_TRUE);CHKERRQ(ierr);
  }

  /* All PETSc objects should be destroyed when they are no longer needed. */
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&(user.f));CHKERRQ(ierr);
  ierr = SNESDestroy(&snes);CHKERRQ(ierr);
  ierr = DMDestroy(&da);CHKERRQ(ierr);
  ierr = PetscFinalize();

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "FormSource"
/*  FormSource:   Fill in source term f(x,y) in PDE from state X.
Does finite-differencing on X, so care is required with
differencing when working in parallel. */
PetscErrorCode FormSource(DM da,PorousCtx *user,Vec X,Vec F)
{
  PetscInt       i,j,Mx,My,xs,ys,xm,ym;
  PetscErrorCode ierr;
  PetscReal      hx,hy;
  PetscReal      L, sig, beta;
  PetscScalar    **x, **f;
  Vec            Xlocal;
  PetscScalar    W, Wpow, Weast, Wwest, Wnorth, Wsouth,
                 Qx, Qy;

  PetscFunctionBegin;

  /* Get global and then local grid boundaries (for 2-dimensional DMDA):
       Mx, My   - total number of grid points in each dimension
       xs, ys   - starting grid indices (no ghost points)
       xm, ym   - widths of local grid (no ghost points)   */
  ierr = DMDAGetInfo(da,PETSC_IGNORE,&Mx,&My,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);
  ierr = DMDAGetCorners(da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL);CHKERRQ(ierr);

  /* set local constants */
  L      = user->L;
  sig    = user->sigma;
  beta   = user->beta;
  hx     = (2.0 * L) / (PetscReal)(Mx-1);
  hy     = (2.0 * L) / (PetscReal)(My-1);

  /* we are going to difference X, so get a local vec and
     communicate X into it */
  ierr = DMCreateLocalVector(da,&Xlocal);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(da,X,INSERT_VALUES,Xlocal);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da,X,INSERT_VALUES,Xlocal);CHKERRQ(ierr);

  /* Compute source term in PDE over the locally owned part of
     the grid by finite-differencing initial guess */
  ierr = DMDAVecGetArray(da,F,&f);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,Xlocal,&x);CHKERRQ(ierr);
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      if (i == 0 || j == 0 || i == Mx-1 || j == My-1) {
        /* no other value makes sense; can't difference */
        f[j][i] = 0.0;
      } else {
        W       = x[j][i];
        Weast   = 0.5 * (x[j][i+1] + W);
        Wwest   = 0.5 * (x[j][i-1] + W);
        Wnorth  = 0.5 * (x[j+1][i] + W);
        Wsouth  = 0.5 * (x[j-1][i] + W);
        Wpow    = pow(W,sig);
        Qx      = ( Weast  * (pow(x[j][i+1],sig) - Wpow)
                      - Wwest  * (Wpow - pow(x[j][i-1],sig)) );
        Qy      = ( Wnorth * (pow(x[j+1][i],sig) - Wpow)
                      - Wsouth * (Wpow - pow(x[j-1][i],sig)) );
        f[j][i] = W + beta * (Qx / (hx*hx) + Qy / (hy*hy));
      }
    }
  }
  ierr = DMDAVecRestoreArray(da,F,&f);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da,Xlocal,&x);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "FormFunctionLocal"
/*  FormFunctionLocal - Evaluates nonlinear residual function r(x) on
    local process patch.  */
PetscErrorCode FormFunctionLocal(DMDALocalInfo *info,
                                 PetscScalar **x,PetscScalar **r,PorousCtx *user) {
  PetscErrorCode ierr;
  PetscInt       i,j;
  PetscReal      hx,hy,hxohy,hyohx,carea;
  PetscReal      beta = user->beta, sig = user->sigma;
  PetscScalar    W, Wpow, Weast, Wwest, Wnorth, Wsouth,
                 Qx, Qy;
  PetscScalar    **f;

  PetscFunctionBegin;
  
  user->fcncount = user->fcncount + 1;
  for (j=info->ys; j<info->ys+info->ym; j++) {
    for (i=info->xs; i<info->xs+info->xm; i++) {
      if (x[j][i] < 0) {
        SETERRQ3(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
        "negative solution value at (i,j)=(%d,%d) during function eval %d",
        i,j,user->fcncount);
      }
    }
  }

  hx     = (2.0 * user->L) / (PetscReal)(info->mx-1);
  hy     = (2.0 * user->L) / (PetscReal)(info->my-1);
  hxohy  = hx / hy;
  hyohx  = hy / hx;
  carea  = hx * hy;

  /* Compute function over the locally owned part of the grid */
  ierr = DMDAVecGetArray(info->da,user->f,&f);CHKERRQ(ierr);
  for (j=info->ys; j<info->ys+info->ym; j++) {
    for (i=info->xs; i<info->xs+info->xm; i++) {
      if (i == 0 || j == 0 || i == info->mx-1 || j == info->my-1) {
        /* since Dirichlet condition, residual at boundary is current value */
        r[j][i] = carea * x[j][i];
      } else {
        W       = x[j][i];
        Weast   = 0.5 * (x[j][i+1] + W);
        Wwest   = 0.5 * (x[j][i-1] + W);
        Wnorth  = 0.5 * (x[j+1][i] + W);
        Wsouth  = 0.5 * (x[j-1][i] + W);
        Wpow    = pow(W,sig);
        Qx      = ( Weast  * (pow(x[j][i+1],sig) - Wpow)
                      - Wwest  * (Wpow - pow(x[j][i-1],sig)) );
        Qy      = ( Wnorth * (pow(x[j+1][i],sig) - Wpow)
                      - Wsouth * (Wpow - pow(x[j-1][i],sig)) );
        r[j][i] = carea * (W - f[j][i]) - beta * (hyohx * Qx + hxohy * Qy);
      }
    }
  }
  ierr = DMDAVecRestoreArray(info->da,user->f,&f);CHKERRQ(ierr);

  PetscFunctionReturn(0); 
} 


#undef __FUNCT__
#define __FUNCT__ "FormJacobianLocal"
/* FormJacobianLocal - Evaluates Jacobian matrix on local process patch
Evaluates partial derivatives \partial F_k / \partial W_l =0 on
local process patch. */
PetscErrorCode FormJacobianLocal(DMDALocalInfo *info,
                                 PetscScalar **x,Mat jac,PorousCtx *user) {
  PetscErrorCode ierr;
  PetscInt       i,j;
  MatStencil     col[5],row;
  PetscScalar    v[5];
  PetscReal      hx,hy,hxohy,hyohx,carea;
  PetscReal      beta = user->beta, sig = user->sigma;
  PetscScalar    W, Wpow, Wpowm1, tmp1, tmp2,
                 Wleft,  Wright, Wdown,  Wup,
                 Wwest,  Weast,  Wsouth, Wnorth;

  PetscFunctionBegin;
  hx     = (2.0 * user->L) / (PetscReal)(info->mx-1);
  hy     = (2.0 * user->L) / (PetscReal)(info->my-1);
  hxohy  = hx / hy;
  hyohx  = hy / hx;
  carea  = hx * hy;

  /* Compute entries for the locally owned part of the Jacobian.
      - Currently, all PETSc parallel matrix formats are partitioned by
        contiguous chunks of rows across the processors. 
      - Each processor needs to insert only elements that it owns
        locally (but any non-local elements will be sent to the
        appropriate processor during matrix assembly). 
      - Here, we set all entries for a particular row at once.  */
  for (j=info->ys; j<info->ys+info->ym; j++) {
    for (i=info->xs; i<info->xs+info->xm; i++) {
      row.j = j; row.i = i;  /* identifies the equation:  k = k(i,j) */
      if (i == 0 || j == 0 || i == info->mx-1 || j == info->my-1) {
        /* boundary points */
        v[0] = carea;
        ierr = MatSetValuesStencil(jac,1,&row,1,&row,v,INSERT_VALUES);
            CHKERRQ(ierr);
      } else {
        /* interior grid points; col[.] identifies the variable:  l = l(i,j)
           note x[j][i] = W_m(x_i,y_j)   (a current iterate value) */
        W       = x[j][i];
        Wleft   = x[j][i-1];
        Wright  = x[j][i+1];
        Wdown   = x[j-1][i];
        Wup     = x[j+1][i];
        Wwest   = 0.5 * (Wleft  + W);
        Weast   = 0.5 * (Wright + W);
        Wsouth  = 0.5 * (Wdown  + W);
        Wnorth  = 0.5 * (Wup    + W);
        Wpow    = pow(W,sig);
        Wpowm1  = pow(W,sig-1);
        /* for \partial F_k / \partial W : */
        col[0].j = j;     col[0].i = i;
        tmp1 = 0.5 * ( pow(Wright,sig) - Wpow ) - sig * Weast * Wpowm1
                 - 0.5 * ( Wpow - pow(Wleft,sig) ) - sig * Wwest * Wpowm1;
        tmp2 = 0.5 * ( pow(Wup,sig) - Wpow ) - sig * Wnorth * Wpowm1
                 - 0.5 * ( Wpow - pow(Wdown,sig) ) - sig * Wsouth * Wpowm1;
        v[0] = carea - beta * hyohx * tmp1 - beta * hxohy * tmp2;
        /* for \partial F_k / \partial Wleft : */
        col[1].j = j;     col[1].i = i-1;
        v[1] = -0.5 * ( Wpow - pow(Wleft,sig) ) + sig * Wwest * pow(Wleft,sig-1);
        v[1] = - beta * hyohx * v[1];
        /* for \partial F_k / \partial Wright : */
        col[2].j = j;     col[2].i = i+1;
        v[2] = 0.5 * ( pow(Wright,sig) - Wpow ) + sig * Weast * pow(Wright,sig-1);
        v[2] = - beta * hyohx * v[2];
        /* for \partial F_k / \partial Wdown : */
        col[3].j = j - 1; col[3].i = i;
        v[3] = -0.5 * ( Wpow - pow(Wdown,sig) ) + sig * Wsouth * pow(Wdown,sig-1);
        v[3] = - beta * hxohy * v[3];
        /* for \partial F_k / \partial Wup : */
        col[4].j = j + 1; col[4].i = i;
        v[4] = 0.5 * ( pow(Wup,sig) - Wpow ) + sig * Wnorth * pow(Wup,sig-1);
        v[4] = - beta * hxohy * v[4];
        /* done with row; insert it */
        ierr = MatSetValuesStencil(jac,1,&row,5,col,v,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
  }

  ierr = MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* Tell the matrix we will never add a new nonzero location to the
     matrix. If we do, it will generate an error.                    */
  ierr = MatSetOption(jac,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

