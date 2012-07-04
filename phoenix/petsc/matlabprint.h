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

#ifndef __matlabprint_h
#define __matlabprint_h

/* some utilities for viewing Vec and Mat objects into a common
   matlab file, and drawing figures */

#include <petscdmda.h>
#include <petscsnes.h>

PetscErrorCode printvecmatlab(DM da, Vec x,
                              const char* vecname, const char* filename,
                              PetscBool append) {
  PetscErrorCode ierr;
  PetscViewer    viewer;
  PetscInt       Mx,My;

  PetscFunctionBegin;
  ierr = DMDAGetInfo(da,PETSC_IGNORE,&Mx,&My,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);
  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer, PETSCVIEWERASCII);CHKERRQ(ierr);
  if (append) {
    ierr = PetscViewerFileSetMode(viewer,FILE_MODE_APPEND); CHKERRQ(ierr);
  }
  ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, filename);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) x, vecname); CHKERRQ(ierr);
  ierr = VecView(x, viewer);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"%s = reshape(%s,%d,%d);\n\n",
           vecname,vecname,Mx,My);  CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


PetscErrorCode print2vecmatlab(DM da, Vec x,
                               const char* vecname1, const char* vecname2,
                               const char* filename,
                               PetscBool append) {
  PetscErrorCode ierr;
  PetscViewer    viewer;
  PetscInt       Mx,My;

  PetscFunctionBegin;
  ierr = DMDAGetInfo(da,PETSC_IGNORE,&Mx,&My,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);
  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer, PETSCVIEWERASCII);CHKERRQ(ierr);
  if (append) {
    ierr = PetscViewerFileSetMode(viewer,FILE_MODE_APPEND); CHKERRQ(ierr);
  }
  ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, filename);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"clear tmpvecdof2\n");  CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) x, "tmpvecdof2"); CHKERRQ(ierr);
  ierr = VecView(x, viewer);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"tmpvecdof2 = reshape(tmpvecdof2,2,%d);\n",
           Mx*My);  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"%s = tmpvecdof2(1,:);\n",
           vecname1);  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"%s = reshape(%s,%d,%d);\n\n",
           vecname1,vecname1,Mx,My);  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"%s = tmpvecdof2(2,:);\n",
           vecname2);  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"%s = reshape(%s,%d,%d);\n\n",
           vecname2,vecname2,Mx,My);  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"clear tmpvecdof2\n");  CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


PetscErrorCode printmatmatlab(DM da, Mat A,
                              const char* matname, const char* filename,
                              PetscBool append) {
  PetscErrorCode ierr;
  PetscViewer    viewer;

  PetscFunctionBegin;
  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer, PETSCVIEWERASCII);CHKERRQ(ierr);
  if (append) {
    ierr = PetscViewerFileSetMode(viewer,FILE_MODE_APPEND); CHKERRQ(ierr);
  }
  ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, filename);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) A, matname); CHKERRQ(ierr);
  ierr = MatView(A, viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


PetscErrorCode printerrorfigurematlab(DM da, PetscInt figurenum,
                                      const char* vecname, const char* exactname,
                                      const char* filename,
                                      PetscBool append) {
  PetscErrorCode ierr;
  PetscViewer    viewer;
  PetscInt       Mx,My;

  PetscFunctionBegin;
  ierr = DMDAGetInfo(da,PETSC_IGNORE,&Mx,&My,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);
  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer, PETSCVIEWERASCII);CHKERRQ(ierr);
  if (append) {
    ierr = PetscViewerFileSetMode(viewer,FILE_MODE_APPEND); CHKERRQ(ierr);
  }
  ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, filename);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
           "err = abs(%s-%s) / max(abs(%s(:)));\n"
           "figure(%d), imagesc(err), axis equal, colorbar\n"
           "title('abs(%s-%s) / |%s|_{inf}')\n",
           vecname,exactname,exactname,
           figurenum,
           vecname,exactname,exactname);  CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


PetscErrorCode printfigurematlab(DM da, PetscInt figurenum,
                                 const char* vecnameleft, const char* vecnameright,
                                 const char* filename,
                                 PetscBool append) {
  PetscErrorCode ierr;
  PetscViewer    viewer;
  PetscInt       Mx,My;

  PetscFunctionBegin;
  ierr = DMDAGetInfo(da,PETSC_IGNORE,&Mx,&My,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);
  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer, PETSCVIEWERASCII);CHKERRQ(ierr);
  if (append) {
    ierr = PetscViewerFileSetMode(viewer,FILE_MODE_APPEND); CHKERRQ(ierr);
  }
  ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, filename);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
           "figure(%d)\n"
           "subplot(1,2,1), imagesc(%s), axis equal, title('%s'), colorbar\n"
           "subplot(1,2,2), imagesc(%s), axis equal, title('%s'), colorbar\n",
           figurenum,
           vecnameleft,vecnameleft,
           vecnameright,vecnameright);  CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#endif

