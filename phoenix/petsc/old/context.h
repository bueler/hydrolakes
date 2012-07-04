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

#ifndef __context_h
#define __context_h

#include <petscdmda.h>
#include <math.h>

/* application context for basal hydrology solvers */
typedef struct {
   DM          da;
   Vec         geom, PdP, S;
   PetscReal   rhoi, rhow, g, Lfusion, Aglen, nglen,  /* standard glacier consts */
               sigma,  /* power to get P:  P = p_i (W/(Y+Ymin))^sigma */
               Kmax,   /* upper bound on hydraulic conductivity (some models) */
               Kmin,   /* lower bound on hydraulic conductivity (some models) */
   /* (W,K) model: */
               alpha,  /* coefficient in dK/dt equation */
               Wcrit,  /* dK/dt and P are functions of W/Wcrit */
   /* (W,Y) and (W,P) models: */
               Kconst, /* constant hydraulic conductivity used in some models */
               Ymin,   /* used in P and dP/dt computations */
               Wmin,   /* in dP/dt computation, avoid div by W=0 */
               Cmelt,  /* coefficient to turn off melt terms; = zero for verif case */
               Creep,  /* coefficient of creep; set from time to close empty cavity */
   /* verification case and Barenblatt: */
               gamma, Ybaren, H0, /* constants in verification case */
               W0, t0, R0, /* Barenblatt solution parameters */
               dx, dy, L,  /* grid description; L=half-width of square domain */
               lastWnorm1, maxrnorm;
   PetscBool   run_silent, not_converged_warning, not_conservation_warning;
   PetscInt    fcncount,   /* count function evaluations */
               expernum;   /* various experiments can be run */
} PorousCtx;


extern PetscErrorCode DefaultContext(PorousCtx *user);
extern PetscErrorCode BarenblattConstants(PorousCtx *user);
extern PetscErrorCode DerivedConstants(PorousCtx *user);
extern PetscErrorCode BarenblattState(PorousCtx *user,PetscReal t,Vec X);


#undef __FUNCT__
#define __FUNCT__ "DefaultContext"
PetscErrorCode DefaultContext(PorousCtx *user) {
  PetscErrorCode ierr;

  /* physical constants */
  user->rhoi  = 910.0;       /* kg m-3; ice density; also PISM default */
  user->rhow  = 1000.0;      /* kg m-3; freshwater density; also PISM default */
  user->g     = 9.81;        /* m s-2;  acceleration of gravity */
  user->Lfusion = 3.34e5;    /* J kg-1; latent heat of fusion (Greve&Blatter) */
  user->Aglen = 3.1689e-24;  /* Pa-3 s-1; ice softness (EISMINT96) */
  user->nglen = 3;           /* pure; Glen exponent (EISMINT96) */

  /* reference values from Flowers & Clarke (2002) = F&C */
  user->sigma = 7.0/2.0;     /* pure; */
  user->Kmax  = 0.01,        /* m s-1;  maximum of hydraulic conductivity */
  user->Kmin  = 0.0001;      /* m s-1;  minimum of hydraulic conductivity */

  /* controls in (W,K) model */
  user->Wcrit = 1.0;         /* m */
  user->alpha = (0.1)/(3600.0*24.0); /* s-1;  = (1/10) per day */

  /* regularizations and controls in (W,Y) and (W,P) models */
  user->Kconst = sqrt(user->Kmax * user->Kmin); /* geometric average */
  user->Ymin  = 0.1;
  user->Wmin  = 0.1;
  user->Cmelt = 1.0;         /* pure;   additional coefficient for amount of melt */
  user->Creep = 1.9014e-04;  /* pure;   creep closure coefficient; this value gives
                                1/e time of one day from pure creep closure of
                                cavity:
                                  >> AN3 = 3.1689e-24 * (910 * 9.81 * 3000)^3
                                  AN3 =  0.060870
                                  >> Creep = 1/(3600*24*AN3)
                                  Creep =  1.9014e-04 */

  /* verification case */
  user->H0     = 3000.0;      /* m;      typical ice sheet thickness */
  user->Ybaren = 1.0;        /* m;      verification critical water thickness */

  /* space-time domain and grid */
  user->L  = 2.0e4;            /* m */
  user->dx = -1.0;  user->dy = -1.0;

  /* utility */
  user->lastWnorm1 = -1.0;
  user->maxrnorm = -1.0;
  user->run_silent = PETSC_FALSE;
  user->not_converged_warning = PETSC_FALSE;
  user->not_conservation_warning = PETSC_FALSE;
  user->fcncount = 0;

  /* fill the rest */
  ierr = DerivedConstants(user);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DerivedConstants"
PetscErrorCode DerivedConstants(PorousCtx *user) {
  PetscErrorCode ierr;
  user->gamma  = user->Kconst * user->H0 * user->rhoi
                   / (user->rhow * pow(user->Ybaren,user->sigma));
  ierr = BarenblattConstants(user);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "BarenblattConstants"
PetscErrorCode BarenblattConstants(PorousCtx *user) {
  PetscReal sig = user->sigma;
  user->W0     = 200.0;
  user->R0     = 2.0 * (sig + 1) * user->W0 / sqrt(sig);
  user->t0     = ( (sig + 1) * user->rhow * pow(user->Ybaren,sig) ) 
                   / ( sig * user->rhoi * user->Kconst * user->H0
                         * pow(user->W0,sig-2) );
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "BarenblattState"
/*  BarenblattState:  Compute t value of Barenblatt solution and
  put it in global Vec X (if dof=1) or in first component of X (if dof>1).
  Used for initial guess and verification. */
PetscErrorCode BarenblattState(PorousCtx *user,PetscReal t,Vec X)
{
  PetscErrorCode ierr;
  DM             da = user->da, coordDA;
  Vec            coordinates;
  DMDACoor2d     **coords;
  PetscInt       i,j,xs,ys,xm,ym,dof;
  PetscReal      trp, A, psi, rsqr, W0, sig;
  PetscScalar    **xx, W;
  typedef struct { PetscReal W, Y; } WYnode;
  WYnode**       xxx;

  PetscFunctionBegin;

  ierr = DMDAGetCorners(da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetInfo(da,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     &dof,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);

  ierr = DMDAGetCoordinateDA(da, &coordDA);CHKERRQ(ierr);
  ierr = DMDAGetCoordinates(da, &coordinates);CHKERRQ(ierr);

  /* set local constants */
  W0     = user->W0;
  sig    = user->sigma;

  /* compute Barenblatt at t */
  trp = pow( user->t0 / t, 1.0 / (sig+1.0) );
  A   = 4.0 * (sig+1.0) * (sig+1.0) * W0 * W0;
  A   = (sig * trp) / A;

  /* Compute initial guess over the locally owned part of the grid */
  if (dof == 1) { ierr = DMDAVecGetArray(da,X,&xx);CHKERRQ(ierr); }
  else          { ierr = DMDAVecGetArray(da,X,&xxx);CHKERRQ(ierr); }
  ierr = DMDAVecGetArray(coordDA, coordinates, &coords);CHKERRQ(ierr);
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      /* Matlab:   psi = 1 - A * rsqr;  psi = max(psi,0.0);
                   W   = W0 * trp * psi.^(1/sigma); */
      rsqr = coords[j][i].x * coords[j][i].x + coords[j][i].y * coords[j][i].y;
      psi  = 1.0 - A * rsqr;
      psi  = PetscMax(psi,0.0);
      W = W0 * trp * pow(psi,(1.0/sig));
      if (dof == 1)  xx[j][i] = W;
      else           xxx[j][i].W = W;
    }
  }
  ierr = DMDAVecRestoreArray(coordDA, coordinates, &coords);CHKERRQ(ierr);
  if (dof == 1) { ierr = DMDAVecRestoreArray(da,X,&xx);CHKERRQ(ierr); }
  else          { ierr = DMDAVecRestoreArray(da,X,&xxx);CHKERRQ(ierr); }

  PetscFunctionReturn(0);
}

#endif

