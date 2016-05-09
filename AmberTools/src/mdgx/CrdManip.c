#include <math.h>
#include <stdlib.h>
#include "Matrix.h"
#include "pmeRecip.h"
#include "CrdManip.h"
#include "mdgxVector.h"

#include "pmeDirectDS.h"
#include "TopologyDS.h"

/***=======================================================================***/
/*** CreateCoord: create a coordinate structure.                           ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   natom:  the number of atoms to allocate for                         ***/
/***=======================================================================***/
coord CreateCoord(int natom)
{
  coord crd;

  crd.natom = natom;
  crd.atmid = (int*)calloc(natom, sizeof(int));
  crd.loc = (double*)calloc(3*natom, sizeof(double));
  crd.prvloc = (double*)calloc(3*natom, sizeof(double));
  crd.scrloc = (double*)calloc(3*natom, sizeof(double));
  crd.vel = (double*)calloc(3*natom, sizeof(double));
  crd.prvvel = (double*)calloc(3*natom, sizeof(double));
  crd.frc = (double*)calloc(3*natom, sizeof(double));
  crd.prvfrc = (double*)calloc(3*natom, sizeof(double));
  crd.scrfrc = (double*)calloc(3*natom, sizeof(double));
  crd.U = CreateDmat(3, 3, 0);
  crd.invU = CreateDmat(3, 3, 0);
  crd.fcnorm = CreateDmat(26, 3, 0);

  return crd;
}

/***=======================================================================***/
/*** CopyCoord: copy coord struct crd into Xcrd.                           ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   crd:  the original coord struct                                     ***/
/***=======================================================================***/
coord CopyCoord(coord *crd)
{
  int i;
  coord Xcrd;

  Xcrd.natom  = crd->natom;
  Xcrd.isortho = crd->isortho;
  Xcrd.atmid  = CpyIVec(crd->atmid, crd->natom);
  Xcrd.loc    = CpyDVec(crd->loc, 3*crd->natom);
  Xcrd.prvloc = CpyDVec(crd->prvloc, 3*crd->natom);
  Xcrd.scrloc = CpyDVec(crd->scrloc, 3*crd->natom);
  Xcrd.vel    = CpyDVec(crd->vel, 3*crd->natom);
  Xcrd.prvvel = CpyDVec(crd->prvvel, 3*crd->natom);
  Xcrd.frc    = CpyDVec(crd->frc, 3*crd->natom);
  Xcrd.prvfrc = CpyDVec(crd->prvfrc, 3*crd->natom);
  Xcrd.scrfrc = CpyDVec(crd->scrfrc, 3*crd->natom);
  for (i = 0; i < 6; i++) {
    Xcrd.gdim[i] = crd->gdim[i];
  }
  for (i = 0; i < 3; i++) {
    Xcrd.hgdim[i] = crd->hgdim[i];
  }
  Xcrd.U = CreateDmat(3, 3, 0);
  Xcrd.invU = CreateDmat(3, 3, 0);
  CompXfrm(Xcrd.gdim, Xcrd.U, Xcrd.invU);
  CopyDmat(&Xcrd.fcnorm, &crd->fcnorm, 0);

  return Xcrd;
}

/***=======================================================================***/
/*** DestroyCoord: free all memory associated with a coord struct.         ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   crd: the coord struct to destroy                                    ***/
/***=======================================================================***/
void DestroyCoord(coord *crd)
{
  free(crd->atmid);
  free(crd->loc);
  free(crd->prvloc);
  free(crd->scrloc);
  free(crd->vel);
  free(crd->prvvel);
  free(crd->frc);
  free(crd->prvfrc);
  free(crd->scrfrc);
  DestroyDmat(&crd->U);
  DestroyDmat(&crd->invU);
  DestroyDmat(&crd->fcnorm);
}

/***=======================================================================***/
/*** TransCrd: translates a set of coordinates by step*tvec[].             ***/
/***=======================================================================***/
void TransCrd(double* crds, int natom, double* tvec, double step)
{
  int i;
  double nvec[3];

  for (i = 0; i < 3; i++) {
    nvec[i] = tvec[i]*step;
  }
  for (i = 0; i < natom; i++) {
    crds[3*i] += nvec[0];
    crds[3*i+1] += nvec[1];
    crds[3*i+2] += nvec[2];
  }
}

/***=======================================================================***/
/*** RotateCrd: rotates a set of coordinates using matrix U.               ***/
/***=======================================================================***/
void RotateCrd(double* crds, int natom, dmat U)
{
  int i;
  double tmp_crds[3];
  double* dtmp;

  dtmp = &crds[0];
  for (i = 0; i < natom; i++) {
    tmp_crds[0] = U.data[0]*dtmp[0] + U.data[1]*dtmp[1] + U.data[2]*dtmp[2];
    tmp_crds[1] = U.data[3]*dtmp[0] + U.data[4]*dtmp[1] + U.data[5]*dtmp[2];
    tmp_crds[2] = U.data[6]*dtmp[0] + U.data[7]*dtmp[1] + U.data[8]*dtmp[2];
    dtmp[0] = tmp_crds[0];
    dtmp[1] = tmp_crds[1];
    dtmp[2] = tmp_crds[2];
    dtmp = &dtmp[3];
  }
}

/***=======================================================================***/
/*** CompXfrm: compute the transformation matrices U and invU that take a  ***/
/***           set of coordinates into and out of box space.  The matrices ***/
/***           U and invU must be pre-allocated.                           ***/
/***=======================================================================***/
void CompXfrm(double* cd, dmat U, dmat invU)
{
  double dx, dy;

  dx = (cos(cd[4])*cos(cd[5]) - cos(cd[3])) / (sin(cd[4]) * sin(cd[5]));
  dy = sqrt(1.0 - dx*dx);
  invU.data[0] = cd[0];
  invU.data[1] = cd[1]*cos(cd[5]);
  invU.data[2] = cd[2]*cos(cd[4]);
  invU.data[4] = cd[1]*sin(cd[5]);
  invU.data[5] = -cd[2]*sin(cd[4])*dx;
  invU.data[8] = cd[2]*sin(cd[4])*dy;
  ttInv(invU, U);
}

/***=======================================================================***/
/*** OrthoReim: re-image a single set of three coordinates (such as a      ***/
/***            displacement) given an orthorhombic box and the knowledge  ***/
/***            that nothing is more than one box length from the primary  ***/
/***            image.                                                     ***/
/***=======================================================================***/
void OrthoReim(double *dx, double *dy, double *dz, double* hgdim)
{
  if (*dx >= hgdim[0]) {
    *dx -= 2.0*hgdim[0];
  }
  else if (*dx < -hgdim[0]) {
    *dx += 2.0*hgdim[0];
  }
  if (*dy >= hgdim[1]) {
    *dy -= 2.0*hgdim[1];
  }
  else if (*dy < -hgdim[1]) {
    *dy += 2.0*hgdim[1];
  }
  if (*dz >= hgdim[2]) {
    *dz -= 2.0*hgdim[2];
  }
  else if (*dz < -hgdim[2]) {
    *dz += 2.0*hgdim[2];
  }
}

/***=======================================================================***/
/*** NonOrthoReim: re-image a single set of three coordinates (such as a   ***/
/***               displacement) given an non-orthorhombic box and the     ***/
/***               appropriate transformation matrices.  No assumptions    ***/
/***               are made about the initial distances between            ***/
/***               coordinates.                                            ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   dx:                                                                 ***/
/***   dy:     Cartesian x, y, and z coordinates                           ***/
/***   dz:                                                                 ***/
/***   [inv]U: transformation matrices [into] and out of unit cell space   ***/
/***=======================================================================***/
void NonOrthoReim(double *dx, double *dy, double *dz, dmat *U, dmat *invU)
{
  double ndx, ndy, ndz;
  double *dtmp;

  dtmp = U->data;
  ndx = dtmp[0]*(*dx) + dtmp[1]*(*dy) + dtmp[2]*(*dz);
  ndy = dtmp[3]*(*dx) + dtmp[4]*(*dy) + dtmp[5]*(*dz);
  ndz = dtmp[6]*(*dx) + dtmp[7]*(*dy) + dtmp[8]*(*dz);
  if (ndx >= 0.5) {
    ndx -= floor(ndx+0.5);
  }
  else if (ndx < -0.5) {
    ndx -= ceil(ndx-0.5);
  }
  if (ndy >= 0.5) {
    ndy -= floor(ndy+0.5);
  }
  else if (ndy < -0.5) {
    ndy -= ceil(ndy-0.5);
  }
  if (ndz >= 0.5) {
    ndz -= floor(ndz+0.5);
  }
  else if (ndz < -0.5) {
    ndz -= ceil(ndz-0.5);
  }
  dtmp = invU->data;
  *dx = dtmp[0]*ndx + dtmp[1]*ndy + dtmp[2]*ndz;
  *dy = dtmp[3]*ndx + dtmp[4]*ndy + dtmp[5]*ndz;
  *dz = dtmp[6]*ndx + dtmp[7]*ndy + dtmp[8]*ndz;
}

/***=======================================================================***/
/*** ImageBondedGroups: align atoms in each bonded group to be in the same ***/
/***                    periodic image.  This routine loops over all atoms ***/
/***                    in the group, starting with the first, and         ***/
/***                    re-images other atoms to be in the minimum image.  ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   crd:      the coordinates (this routine operates on a unified list  ***/
/***             of coordinates, not a cell grid)                          ***/
/***   tp:       the topology (contains bonded lists of atoms)             ***/
/***=======================================================================***/
void ImageBondedGroups(coord *crd, prmtop *tp)
{
  int h, i, tgai3;
  double cenx, ceny, cenz, dx, dy, dz;
  lgrp *tg;

  /*** Take all atoms into box space ***/
  RotateCrd(crd->loc, crd->natom, crd->U);

  /*** Loop over all groups ***/
  for (h = 0; h < tp->ngrp; h++) {
    tg = &tp->lgrps[h];
    cenx = crd->loc[3*tg->atoms[0]];
    ceny = crd->loc[3*tg->atoms[0]+1];
    cenz = crd->loc[3*tg->atoms[0]+2];
    for (i = 1; i < tg->natom; i++) {
      tgai3 = 3*tg->atoms[i];
      dx = crd->loc[tgai3] - cenx;
      dy = crd->loc[tgai3+1] - ceny;
      dz = crd->loc[tgai3+2] - cenz;
      dx -= floor(dx+0.5);
      dy -= floor(dy+0.5);
      dz -= floor(dz+0.5);
      crd->loc[tgai3] = cenx + dx;
      crd->loc[tgai3+1] = ceny + dy;
      crd->loc[tgai3+2] = cenz + dz;
    }
  }

  /*** Take all atoms back to real space ***/
  RotateCrd(crd->loc, crd->natom, crd->invU);
}
