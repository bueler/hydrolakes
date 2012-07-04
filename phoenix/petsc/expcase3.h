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

#ifndef __expcase3_h
#define __expcase3_h

#include <petscdmda.h>
#include "context.h"


/* in numerical experiment case 3:
   -- we don't know exact solution
   -- we use the same surface elevation as in verif case 2
   -- we use a bed with radial troughs
   -- we use a sliding distribution  vbspeed  that follows the troughs
   -- we set friction heating (fric) and geothermal heat (Ggeo) to zero
   -- we use an axially-symmetric S
   -- we start with zero water
   -- we intend to run to steady state
 */

/* context for numerical experiment case 3 */
typedef struct {
   PetscReal   H0,     /* typical ice sheet thickness */
               Ls,     /* radius of ice cap profile */
               b0,     /* magnitude of bed troughs at r=Ls */
               n0,     /* number of bed troughs */
               w0,     /* drainage rate scale as velocity (m/s) */
               v0;     /* sliding speed scale */
} EC3Ctx;


#undef __FUNCT__
#define __FUNCT__ "ExpCase3Context"
PetscErrorCode ExpCase3Context(WPCtx *user, EC3Ctx *ectx) {
  PetscReal secperday=3600.0*24.0, secpera = 31556926.0;

  user->dt    = 100.0 * secperday;
  user->L     = 250000.0;

  ectx->H0       = 2000.0;
  ectx->Ls       = 200000.0;
  ectx->b0       = 500.0;
  ectx->n0       = 5.0;
  ectx->w0       = 1.0 / secpera;     /* 1 m/a  so  max s(r) ~~ 0.15 m/a */
  ectx->v0       = 1000.0 / secpera;  /* 1000 m/a */
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ExpCase3SetGeometry"
/* sets  geometry = b, H
   where  H(r) = (Bueler profile)   and  b(r,theta) = [radiating troughs]  */
PetscErrorCode ExpCase3SetGeometry(WPCtx *user, EC3Ctx *ectx, Vec b, Vec H) {
  PetscErrorCode ierr;
  PetscInt       i,j,xs,ys,xm,ym;
  PetscReal      r, theta, s, n, m, q, CC, phi, usurf, **thk, **bed;
  DM             coordDA;
  Vec            coordinates;
  DMDACoor2d     **coords;

  PetscFunctionBegin;
  ierr = DMDAGetCorners(user->da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetCoordinateDA(user->da, &coordDA);CHKERRQ(ierr);
  ierr = DMDAGetCoordinates(user->da, &coordinates);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(coordDA, coordinates, &coords);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, b, &bed);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, H, &thk);CHKERRQ(ierr);
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      r = sqrt(coords[j][i].x * coords[j][i].x + coords[j][i].y * coords[j][i].y);
      if (r > 0.0)  theta = atan2(coords[j][i].y, coords[j][i].x);
      else          theta = 0.0;
      s = r / ectx->Ls;
      bed[j][i] = ectx->b0 * s * s * cos(ectx->n0 * theta);
      if (s < 1.0) {
        n = user->nglen;
        m = n / (2.0 * n + 2.0);
        q = (n + 1.0) / n;
        CC = ectx->H0 / pow(1.0 - (1.0 / n),m);
        phi = q * s - (1.0 / n) + pow(1.0 - s,q) - pow(s,q);
        usurf = CC * pow(phi,m);
        if (usurf > bed[j][i]) {
          thk[j][i] = usurf - bed[j][i];
        } else {
          thk[j][i] = 0.0;
        }
      } else {
        thk[j][i] = 0.0;
      }
    }
  }
  ierr = DMDAVecRestoreArray(user->da, user->b, &bed);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, user->H, &thk);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(coordDA, coordinates, &coords);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ExpCase3SetDrainageSliding"
/* sets S and vbspeed; requires geometry to already be set */
PetscErrorCode ExpCase3SetDrainageSliding(WPCtx *user, EC3Ctx *ectx, Vec S, Vec vbspeed) {
  PetscErrorCode ierr;
  PetscInt       i,j,xs,ys,xm,ym;
  PetscReal      r, s, bb, **bed, **drain, **vb;
  DM             coordDA;
  Vec            coordinates;
  DMDACoor2d     **coords;

  PetscFunctionBegin;
  ierr = DMDAGetCorners(user->da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetCoordinateDA(user->da, &coordDA);CHKERRQ(ierr);
  ierr = DMDAGetCoordinates(user->da, &coordinates);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(coordDA, coordinates, &coords);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, user->b, &bed);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, S, &drain);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, vbspeed, &vb);CHKERRQ(ierr);
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      r = sqrt(coords[j][i].x * coords[j][i].x + coords[j][i].y * coords[j][i].y);
      s = r / ectx->Ls;
      if (s < 1.0) {
        drain[j][i] = ectx->w0 * user->rhow * s * s * (1.0 - s);
        if (bed[j][i] < 0) {
          bb = bed[j][i] / ectx->b0;
          vb[j][i] = ectx->v0 * bb * bb;
        } else {
          vb[j][i] = 0.0;
        }
      } else {
        drain[j][i] = 0.0;
        vb[j][i] = 0.0;
      }
    }
  }
  ierr = DMDAVecRestoreArray(user->da, user->b, &bed);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, S, &drain);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, vbspeed, &vb);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(coordDA, coordinates, &coords);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ExpCase3SetHeating"
/* sets  heating  = fric, Ggeo */
PetscErrorCode ExpCase3SetHeating(WPCtx *user, EC3Ctx *ectx, Vec fric, Vec Ggeo) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = VecSet(fric,0.0);CHKERRQ(ierr); /* - tau_b . v_b = 0 */
  ierr = VecSet(Ggeo,0.0);CHKERRQ(ierr); /* G = 0 (= conductive heating) */
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ExpCase3SetCurrentW"
/* sets initial value for Wl to zero */
PetscErrorCode ExpCase3SetCurrentW(WPCtx *user, EC3Ctx *ectx, Vec Wl) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = VecSet(Wl,0.0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ExpCase3SetInitialpsi"
/* sets psi initial guess to 0.8 * overburden; requires geometry to already be set
   (recall if psi_init == (overburden) then SNESVI thinks its done)  */
PetscErrorCode ExpCase3SetInitialpsi(WPCtx *user, EC3Ctx *ectx, Vec Psi) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  /*  psi = psi_i = (0.8 rhoi g H) + (rhow g b)  */
  ierr = VecCopy(user->b, Psi);CHKERRQ(ierr);  
  ierr = VecScale(Psi, user->rhow * user->g);CHKERRQ(ierr);  
  ierr = VecAXPY(Psi, 0.8 * user->rhoi * user->g, user->H);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ExpCase3SetInitialWnew"
/* sets Wnew initial guess to 0.5;
   (recall if Wnew_init == 0.0 then SNESVI thinks its done)  */
PetscErrorCode ExpCase3SetInitialWnew(WPCtx *user, EC3Ctx *ectx, Vec Wnew) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = VecSet(Wnew,0.5);CHKERRQ(ierr);  
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ExpCase3SetpsiBC"
/* sets psiBC:  along bdry of computational domain, Dirichlet b.c. for
  pressure is zero, so psi = psi_b */
PetscErrorCode ExpCase3SetpsiBC(WPCtx *user, EC3Ctx *ectx, Vec psiBC) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  /* along boundaries,  psi = psi_b = rhow g b */
  ierr = VecCopy(user->b, psiBC);CHKERRQ(ierr);  
  ierr = VecScale(psiBC, user->rhow * user->g);CHKERRQ(ierr);  
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ExpCase3SetWnewBC"
/* sets WnewBC: along bdry of computational domain, Dirichlet b.c. for Wnew is zero */
PetscErrorCode ExpCase3SetWnewBC(WPCtx *user, EC3Ctx *ectx, Vec WnewBC) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = VecSet(WnewBC,0.0);CHKERRQ(ierr);  
  PetscFunctionReturn(0);
}

#endif	/* __expcase3_h */

