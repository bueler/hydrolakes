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

#ifndef __verifcase2_h
#define __verifcase2_h

#include <petscdmda.h>
#include "context.h"


/* context for verification case 2 */
typedef struct {
   PetscReal   H0,     /* typical ice sheet thickness */
               alpha0, /* spatial rate for current water thickness */
               alpha1, /* spatial rate for new water thickness */
               beta,   /* spatial rate for pressure */
               Ls,     /* radius of ice cap profile */
               P0,     /* constant overburden pressure */
               delta0; /* coefficient for water thicknesses */
} VC2Ctx;


#undef __FUNCT__
#define __FUNCT__ "VerifCase2Context"
PetscErrorCode VerifCase2Context(WPCtx *user, VC2Ctx *v) {
  PetscReal secperday=3600.0*24.0, pi=3.1415926536;

  user->Cmelt  = 0.0;
  user->dt     = 100.0 * secperday;
  user->L      = 250000.0;
  user->omega0 = 0.0;

  v->alpha0   = 45000.0;
  v->alpha1   = 50000.0;
  v->Ls       = 200000.0;
  v->beta     = pi / (3.0 * v->Ls);
  v->H0       = 2000.0;
  v->P0       = 0.9 * (user->rhoi * user->g * v->H0);
  v->delta0   = 0.0001;
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VerifCase2SetGeometry"
/* sets  geometry = b, H
   where  H(r) = (Bueler profile) */
PetscErrorCode VerifCase2SetGeometry(WPCtx *user, VC2Ctx *v, Vec b, Vec H) {
  PetscErrorCode ierr;
  PetscInt       i,j,xs,ys,xm,ym;
  PetscReal      r, n, m, q, s, CC, phi, **h;
  DM             coordDA;
  Vec            coordinates;
  DMDACoor2d     **coords;

  PetscFunctionBegin;
  ierr = VecSet(b,0.0);CHKERRQ(ierr);   /* b = 0 */

  ierr = DMDAGetCorners(user->da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetCoordinateDA(user->da, &coordDA);CHKERRQ(ierr);
  ierr = DMDAGetCoordinates(user->da, &coordinates);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(coordDA, coordinates, &coords);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, H, &h);CHKERRQ(ierr);
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      r = sqrt(coords[j][i].x * coords[j][i].x + coords[j][i].y * coords[j][i].y);
      s = r / v->Ls;
      if (s < 1.0) {
        n = user->nglen;
        m = n / (2.0 * n + 2.0);
        q = (n + 1.0) / n;
        CC = v->H0 / pow(1.0 - (1.0 / n),m);
        phi = q * s - (1.0 / n) + pow(1.0 - s,q) - pow(s,q);
        h[j][i] = CC * pow(phi,m);
      } else {
        h[j][i] = 0.0;
      }
    }
  }
  ierr = DMDAVecRestoreArray(user->da, user->H, &h);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(coordDA, coordinates, &coords);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VerifCase2SetHeating"
/* sets  heating  = fric, Ggeo */
PetscErrorCode VerifCase2SetHeating(WPCtx *user, VC2Ctx *v, Vec fric, Vec Ggeo) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = VecSet(fric,0.0);CHKERRQ(ierr); /* - tau_b . v_b = 0 */
  ierr = VecSet(Ggeo,0.0);CHKERRQ(ierr); /* G = 0 (= conductive heating) */
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VerifCase2SetCurrentW"
/* sets Wl */
PetscErrorCode VerifCase2SetCurrentW(WPCtx *user, VC2Ctx *v, Vec Wl) {
  PetscErrorCode ierr;
  PetscInt       i,j,xs,ys,xm,ym;
  DM             coordDA;
  Vec            coordinates;
  DMDACoor2d     **coords;
  PetscReal      **w, r;

  PetscFunctionBegin;
  ierr = DMDAGetCorners(user->da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetCoordinateDA(user->da, &coordDA);CHKERRQ(ierr);
  ierr = DMDAGetCoordinates(user->da, &coordinates);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(coordDA, coordinates, &coords);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, Wl, &w);CHKERRQ(ierr);
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      r = sqrt(coords[j][i].x * coords[j][i].x + coords[j][i].y * coords[j][i].y);
      if (r < v->Ls) {
        w[j][i] = v->delta0 * r * exp(-r / v->alpha0);
      } else {
        w[j][i] = 0.0;
      }
    }
  }
  ierr = DMDAVecRestoreArray(user->da, Wl, &w);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(coordDA, coordinates, &coords);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VerifCase2SetInitialpsi"
/* sets psi */
PetscErrorCode VerifCase2SetInitialpsi(WPCtx *user, VC2Ctx *v, Vec Psi) {
  PetscErrorCode ierr;
  PetscInt       i,j,xs,ys,xm,ym;
  DM             coordDA;
  Vec            coordinates;
  DMDACoor2d     **coords;
  PetscReal      **psi, **h;

  PetscFunctionBegin;
  ierr = DMDAGetCorners(user->da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetCoordinateDA(user->da, &coordDA);CHKERRQ(ierr);
  ierr = DMDAGetCoordinates(user->da, &coordinates);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(coordDA, coordinates, &coords);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, user->H, &h);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, Psi, &psi);CHKERRQ(ierr);
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      /* P = (80% of overburden) */
      psi[j][i] = 0.8 * user->rhoi * user->g * h[j][i];
    }
  }
  ierr = DMDAVecRestoreArray(user->da, user->H, &h);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, Psi, &psi);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(coordDA, coordinates, &coords);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VerifCase2SetInitialWnew"
/* sets Wnew */
PetscErrorCode VerifCase2SetInitialWnew(WPCtx *user, VC2Ctx *v, Vec Wnew) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  /* use previous as guess */
  ierr = VecCopy(user->Wl,Wnew);CHKERRQ(ierr);  
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VerifCase2Manufacture"
/* sets Pexact, Wnewexact, S, vbspeed */
PetscErrorCode VerifCase2Manufacture(WPCtx *user, VC2Ctx *v,
                                     Vec Pexact, Vec Wnewexact,
                                     Vec S, Vec vbspeed) {
  PetscErrorCode ierr;
  PetscInt       i,j,xs,ys,xm,ym;
  DM             coordDA;
  Vec            coordinates;
  DMDACoor2d     **coords;
  PetscReal      pi=3.1415926536, r, Y, Z0, Z1, Wl, c0, CC1, CC2, CC3, pover,
                 **h, **pexact, **wexact, **s, **vb;

  PetscFunctionBegin;
  ierr = DMDAGetCorners(user->da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL);CHKERRQ(ierr);  
  ierr = DMDAGetCoordinateDA(user->da, &coordDA);CHKERRQ(ierr);
  ierr = DMDAGetCoordinates(user->da, &coordinates);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(coordDA, coordinates, &coords);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, Pexact, &pexact);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, Wnewexact, &wexact);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, S, &s);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, vbspeed, &vb);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, user->H, &h);CHKERRQ(ierr);
  c0  = user->K0 / (user->rhow * user->g);
  CC1 = user->Creep * user->Aglen / user->Cavit;
  CC2 = (4.0 / 3.0) * user->rhow * c0 * v->P0 * v->delta0 * v->beta;
  CC3 = 4.0 * pi * c0 * v->P0 * v->delta0 / (9.0 * v->Ls * user->Cavit);
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      r = sqrt(coords[j][i].x * coords[j][i].x + coords[j][i].y * coords[j][i].y);
      if (r < v->Ls) {
        Y  = 2.0 * v->beta * r;
        Z0 = r / v->alpha0;
        Z1 = r / v->alpha1;
        Wl = v->delta0 * r * exp(-Z0);
        pexact[j][i] = (v->P0 / 3.0) * (2.0 * cos(Y) + 1.0);
        wexact[j][i] = v->delta0 * r * exp(-Z1);
        pover        = user->rhoi * user->g * h[j][i];
        s[j][i]      = (user->rhow / user->dt) * (wexact[j][i] - Wl)
                       + CC2 * exp(-Z1) * ( (2.0 - Z1) * sin(Y) + Y * cos(Y) );
        vb[j][i]     = - CC3 * exp(-Z0) * ( (2.0 - Z0) * sin(Y) + Y * cos(Y) )
                       + CC1 * pow(pover - pexact[j][i],user->nglen) * Wl
                       + s[j][i] / (user->rhow * user->Cavit);
      } else {
        pexact[j][i] = 0.0;
        wexact[j][i] = 0.0;
        s[j][i]      = 0.0;
        vb[j][i]     = 0.0;
      }
    }
  }
  ierr = DMDAVecRestoreArray(user->da, user->H, &h);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, Pexact, &pexact);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, Wnewexact, &wexact);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, S, &s);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, vbspeed, &vb);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(coordDA, coordinates, &coords);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VerifCase2SetpsiBC"
/* sets psiBC:  along bdry of computational domain, exact solution is
                Dirichlet b.c. for psi */
PetscErrorCode VerifCase2SetpsiBC(WPCtx *user, VC2Ctx *v, Vec psiBC) {
  PetscErrorCode ierr;
  PetscInt       i,j,xs,ys,xm,ym,Mx,My;
  DM             coordDA;
  Vec            coordinates;
  DMDACoor2d     **coords;
  PetscReal      r, Y, **psi;

  PetscFunctionBegin;
  ierr = DMDAGetCorners(user->da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetInfo(user->da,PETSC_IGNORE,&Mx,&My,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);
  ierr = DMDAGetCoordinateDA(user->da, &coordDA);CHKERRQ(ierr);
  ierr = DMDAGetCoordinates(user->da, &coordinates);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(coordDA, coordinates, &coords);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, psiBC, &psi);CHKERRQ(ierr);
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      if (i == 0 || j == 0 || i == Mx-1 || j == My-1) {
        r = sqrt(coords[j][i].x * coords[j][i].x + coords[j][i].y * coords[j][i].y);
        if (r < v->Ls) {
          Y = 2.0 * v->beta * r;
          psi[j][i] = (v->P0 / 3.0) * (2.0 * cos(Y) + 1.0);
        } else {
          psi[j][i] = 0.0;
        }
      } else {
        psi[j][i] = NAN;
      }
    }
  }
  ierr = DMDAVecRestoreArray(user->da, psiBC, &psi);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(coordDA, coordinates, &coords);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VerifCase2SetWnewBC"
/* sets WnewBC: along bdry of computational domain, exact solution is
                Dirichlet b.c. for Wnew (FIXME: more thought needed on inflow) */
PetscErrorCode VerifCase2SetWnewBC(WPCtx *user, VC2Ctx *v, Vec WnewBC) {
  PetscErrorCode ierr;
  PetscInt       i,j,xs,ys,xm,ym,Mx,My;
  DM             coordDA;
  Vec            coordinates;
  DMDACoor2d     **coords;
  PetscReal      r, **w;

  PetscFunctionBegin;
  ierr = DMDAGetCorners(user->da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetInfo(user->da,PETSC_IGNORE,&Mx,&My,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);
  ierr = DMDAGetCoordinateDA(user->da, &coordDA);CHKERRQ(ierr);
  ierr = DMDAGetCoordinates(user->da, &coordinates);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(coordDA, coordinates, &coords);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, WnewBC, &w);CHKERRQ(ierr);
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      if (i == 0 || j == 0 || i == Mx-1 || j == My-1) {
        r = sqrt(coords[j][i].x * coords[j][i].x + coords[j][i].y * coords[j][i].y);
        if (r < v->Ls) {
          w[j][i] = v->delta0 * r * exp(- r / v->alpha1);
        } else {
          w[j][i] = 0.0;
        }
      } else {
        w[j][i] = NAN;
      }
    }
  }
  ierr = DMDAVecRestoreArray(user->da, WnewBC, &w);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(coordDA, coordinates, &coords);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#endif	/* __verifcase2_h */

