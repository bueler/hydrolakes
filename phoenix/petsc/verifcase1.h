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

#ifndef __verifcase1_h
#define __verifcase1_h

#include <petscdmda.h>
#include "context.h"


/* context for verification case 1 */
typedef struct {
   PetscReal   H0,     /* typical ice sheet thickness */
               alpha,  /* spatial rate for water thickness */
               beta,   /* spatial rate for pressure */
               gamma,  /* amount of depression of pressure */
               P0;     /* constant overburden pressure */
} VC1Ctx;


#undef __FUNCT__
#define __FUNCT__ "VerifCase1Context"
PetscErrorCode VerifCase1Context(WPCtx *user, VC1Ctx *v) {
  PetscReal secperday=3600.0*24.0;

  user->Cmelt  = 0.0;
  user->dt     = 100.0 * secperday;
  user->L      = 20000.0;
  user->omega0 = 0.0;

  v->alpha  = 3.0e-4;
  v->beta   = 0.0002;
  v->gamma  = 0.15;
  v->H0     = 1000.0;
  v->P0     = user->rhoi * user->g * v->H0;
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VerifCase1SetGeometry"
/* sets  geometry = b, H */
PetscErrorCode VerifCase1SetGeometry(WPCtx *user, VC1Ctx *v, Vec b, Vec H) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = VecSet(b,0.0);CHKERRQ(ierr);   /* b = 0 */
  ierr = VecSet(H,v->H0);CHKERRQ(ierr); /* H = H_0 */
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VerifCase1SetHeating"
/* sets  heating  = fric, Ggeo */
PetscErrorCode VerifCase1SetHeating(WPCtx *user, VC1Ctx *v, Vec fric, Vec Ggeo) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = VecSet(fric,0.0);CHKERRQ(ierr); /* - tau_b . v_b = 0 */
  ierr = VecSet(Ggeo,0.0);CHKERRQ(ierr); /* G = 0 (= conductive heating) */
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VerifCase1SetCurrentW"
/* sets Wl */
PetscErrorCode VerifCase1SetCurrentW(WPCtx *user, VC1Ctx *v, Vec Wl) {
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
      w[j][i] = 1.0 + 0.5 * cos(v->alpha * r);
    }
  }
  ierr = DMDAVecRestoreArray(user->da, Wl, &w);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(coordDA, coordinates, &coords);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VerifCase1SetInitialpsi"
/* sets psi */
PetscErrorCode VerifCase1SetInitialpsi(WPCtx *user, VC1Ctx *v, Vec Psi) {
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
#define __FUNCT__ "VerifCase1SetInitialWnew"
/* sets Wnew */
PetscErrorCode VerifCase1SetInitialWnew(WPCtx *user, VC1Ctx *v, Vec Wnew) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  /* should not use W^l because that is exact; note x=0 does not work because
       it is against the VI bound */
  ierr = VecSet(Wnew,1.0);CHKERRQ(ierr);  
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VerifCase1Manufacture"
/* sets Pexact, Wnewexact, S, vbspeed */
PetscErrorCode VerifCase1Manufacture(WPCtx *user, VC1Ctx *v,
                                     Vec Pexact, Vec Wnewexact,
                                     Vec S, Vec vbspeed) {
  PetscErrorCode ierr;
  PetscInt       i,j,xs,ys,xm,ym;
  DM             coordDA;
  Vec            coordinates;
  DMDACoor2d     **coords;
  PetscReal      r, Y, Z, Wl, cp, c0, CC0, CC1,
                 **pexact, **wexact, **s, **vb;

  PetscFunctionBegin;
  ierr = DMDAGetCorners(user->da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL);CHKERRQ(ierr);  
  ierr = DMDAGetCoordinateDA(user->da, &coordDA);CHKERRQ(ierr);
  ierr = DMDAGetCoordinates(user->da, &coordinates);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(coordDA, coordinates, &coords);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, Pexact, &pexact);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, Wnewexact, &wexact);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, S, &s);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, vbspeed, &vb);CHKERRQ(ierr);
  c0  = user->K0 / (user->rhow * user->g);
  CC0 = c0 * user->rhow * v->P0 * v->gamma * v->beta;
  CC1 = user->Creep * user->Aglen / user->Cavit;
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      r  = sqrt(coords[j][i].x * coords[j][i].x + coords[j][i].y * coords[j][i].y);
      Y  = 2.0 * v->beta * r;
      Z = (r > 0.0) ? sin(Y) / r : 2.0 * v->beta;
      cp = cos(v->beta * r);
      pexact[j][i] = v->P0 * (1.0 - v->gamma * cp * cp);
      Wl = 1.0 + 0.5 * cos(v->alpha * r);
      wexact[j][i] = Wl;
      s[j][i]  = - CC0 * ( Wl * (Z + 2.0 * v->beta * cos(Y))
                           - 0.5 * v->alpha * sin(v->alpha * r) * sin(Y) );
      vb[j][i] = CC1 * pow(v->P0 * v->gamma * cp * cp, user->nglen) * Wl;
    }
  }
  ierr = DMDAVecRestoreArray(user->da, Pexact, &pexact);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, Wnewexact, &wexact);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, S, &s);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, vbspeed, &vb);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(coordDA, coordinates, &coords);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VerifCase1SetpsiBC"
/* sets psiBC:  along bdry of computational domain, exact solution is
                Dirichlet b.c. for psi */
PetscErrorCode VerifCase1SetpsiBC(WPCtx *user, VC1Ctx *v, Vec psiBC) {
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
        Y = cos(v->beta * r);
        psi[j][i] = v->P0 * (1.0 - v->gamma * Y * Y);
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
#define __FUNCT__ "VerifCase1SetWnewBC"
/* sets WnewBC: along bdry of computational domain, exact solution is
                Dirichlet b.c. for Wnew (FIXME: more thought needed on inflow) */
PetscErrorCode VerifCase1SetWnewBC(WPCtx *user, VC1Ctx *v, Vec WnewBC) {
  PetscErrorCode ierr;
  PetscInt       i,j,xs,ys,xm,ym,Mx,My;
  PetscReal      **w, **wold;

  PetscFunctionBegin;
  ierr = DMDAGetCorners(user->da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetInfo(user->da,PETSC_IGNORE,&Mx,&My,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, WnewBC, &w);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, user->Wl, &wold);CHKERRQ(ierr);
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      if (i == 0 || j == 0 || i == Mx-1 || j == My-1) {
        w[j][i] = wold[j][i];
      } else {
        w[j][i] = NAN;
      }
    }
  }
  ierr = DMDAVecRestoreArray(user->da, WnewBC, &w);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, user->Wl, &wold);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#endif	/* __verifcase1_h */

