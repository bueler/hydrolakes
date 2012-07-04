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

#ifndef __context_h
#define __context_h

#include <petscdmda.h>

/* application context for a (W,psi)-space hydrology solver */
typedef struct {
   DM          da;
   Vec         Wl, psiforW,
               H, b,
               fric, Ggeo,
               S, vbspeed,
               psiBC, WnewBC;
   PetscReal   g, Lfusion, rhoi, rhow, /* physical consts */
               Aglen, nglen,  /* Glen power law for flow */
               K0,         /* constant hydraulic conductivity */
               Cavit,      /* cavitation coefficient */
               Cmelt,      /* coefficient of dissipation heating in wall melt */
               Creep,      /* coefficient of creep */
               omega0,     /* regularization water thickness for pressure problem */
               dx, dy,     /* grid cell dimensions */
               L,          /* half-width of square domain */
               dt;         /* time step in seconds */
   PetscBool   not_converged_warning, not_conservation_warning;
   PetscInt    fcncount;   /* count function evaluations */
} WPCtx;


#undef __FUNCT__
#define __FUNCT__ "DefaultContext"
PetscErrorCode DefaultContext(WPCtx *user) {
  PetscReal      Kmax, Kmin;
  /* physical constants */
  user->rhoi  = 910.0;       /* kg m-3; ice density; also PISM default */
  user->rhow  = 1000.0;      /* kg m-3; freshwater density; also PISM default */
  user->g     = 9.81;        /* m s-2;  acceleration of gravity */
  user->Lfusion = 3.34e5;    /* J kg-1; latent heat of fusion (Greve&Blatter) */
  user->Aglen = 3.1689e-24;  /* Pa-3 s-1; ice softness (EISMINT96) */
  user->nglen = 3.0;         /* pure; Glen exponent (EISMINT96) */
  /* reference values from Flowers & Clarke (2002) = F&C */
  Kmax  = 0.01,        /* m s-1;  maximum of hydraulic conductivity */
  Kmin  = 0.0001;      /* m s-1;  minimum of hydraulic conductivity */
  user->K0    = sqrt(Kmax * Kmin); /* geometric average */
  /* controls */
  user->Cavit = 0.001;       /* pure */
  user->Cmelt = 1.0;         /* pure */
  user->Creep = 1.9014e-04;  /* pure;   creep closure coefficient; this value gives
                                1/e time of one day from pure creep closure of
                                cavity:
                                  >> AN3 = 3.1689e-24 * (910 * 9.81 * 3000)^3
                                  AN3 =  0.060870
                                  >> Creep = 1/(3600*24*AN3)
                                  Creep =  1.9014e-04 */
  user->omega0 = 0.1;        /* m */
  /* utility */
  user->dx = -1.0;  user->dy = -1.0;
  user->dt = -1.0;
  user->L = -1.0;
  user->not_converged_warning = PETSC_FALSE;
  user->not_conservation_warning = PETSC_FALSE;
  user->fcncount = 0;

  PetscFunctionReturn(0);
}

#endif	/* __context_h */

