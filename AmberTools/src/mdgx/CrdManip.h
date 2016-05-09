#ifndef CrdManipHeadings
#define CrdManipHeadings

#include "CrdManipDS.h"
#include "MatrixDS.h"
#include "TopologyDS.h"

coord CreateCoord(int natm);

coord CopyCoord(coord *crd);

void DestroyCoord(coord *crd);

void TransCrd(double* crds, int natom, double* tvec, double step);

void RotateCrd(double* crds, int natom, dmat U);

void CompXfrm(double* cd, dmat U, dmat invU);

void OrthoReim(double *dx, double *dy, double *dz, double* hgdim);

void NonOrthoReim(double *dx, double *dy, double *dz, dmat *U, dmat *invU);

void ImageBondedGroups(coord *crd, prmtop *tp);

#endif
