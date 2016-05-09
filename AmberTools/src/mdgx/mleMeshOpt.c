#include <math.h>
#include <stdlib.h>
#include "Grid.h"
#include "BSpline.h"
#include "mdgxVector.h"
#include "pmeRecip.h"
#include "mleRecip.h"
#include "ChargeMap.h"
#include "fftw3.h"

#include "mleOptDS.h"

/***=======================================================================***/
/*** MeshBasisFunctions: this function takes a reciprocal space pair       ***/
/***                     potential and numbers all the unique points in it ***/
/***                     using the criteria of 1.) location relative to    ***/
/***                     the origin and 2.) value of the potential at the  ***/
/***                     point.                                            ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   Urec: the reciprocal space pair potential (in real space)           ***/
/***   crd:  the system coordinates (needed for unit cell dimensions)      ***/
/***=======================================================================***/
bssf MeshBasisFunctions(dbook *Urec, coord *crd)
{
  int i, j, k, ii, jj, kk, nivec, njvec, nkvec, iloc, jloc, kloc;
  int ivec[2], jvec[2], kvec[2], ng[3];
  int *itmp;
  int* bcon;
  double dx, dy, dz, mmx, mmy, mmz, r2ijk;
  double *dtmp, *ivd;
  double* mx;
  double* my;
  double* mz;
  dbook r2;
  ibook bn;
  bssf basis;

  /*** Determine each point's distance from the origin ***/
  ng[0] = Urec->pag;
  ng[1] = Urec->row;
  ng[2] = Urec->col;
  r2 = CreateDbook(ng[0], ng[1], ng[2], 0);
  mx = LoadMVec(ng[0]);
  my = LoadMVec(ng[1]);
  mz = LoadMVec(ng[2]);
  ivd = crd->invU.data;
  for (i = 0; i < ng[0]; i++) {
    dx = mx[i];
    for (j = 0; j < ng[1]; j++) {
      dy = my[j];
      dtmp = r2.map[i][j];
      for (k = 0; k < ng[2]; k++) {
        dz = mz[k];
        mmx = ivd[0]*dx + ivd[1]*dy + ivd[2]*dz;
        mmy = ivd[3]*dx + ivd[4]*dy + ivd[5]*dz;
        mmz = ivd[6]*dx + ivd[7]*dy + ivd[8]*dz;
        dtmp[k] = mmx*mmx + mmy*mmy + mmz*mmz;
      }
    }
  }

  /*** Count the number of basis functions ***/
  bn = CreateIbook(ng[0], ng[1], ng[2]);
  SetIVec(bn.data, ng[0]*ng[1]*ng[2], -1);

  basis.nbss = 0;
  for (i = 0; i < ng[0]; i++) {
    ivec[0] = i;
    if (i > 0) {
      ivec[1] = ng[0]-i;
      nivec = 2;
    }
    else {
      nivec = 1;
    }
    for (j = 0; j < ng[1]; j++) {
      jvec[0] = j;
      if (j > 0) {
        jvec[1] = ng[1]-j;
        njvec = 2;
      }
      else {
        njvec = 1;
      }
      itmp = bn.map[i][j];
      for (k = 0; k < ng[2]; k++) {
        kvec[0] = k;
        if (k > 0) {
          kvec[1] = ng[2]-k;
          nkvec = 2;
        }
        else {
          nkvec = 1;
        }

        /*** Skip if this element is already part of a basis function ***/
        if (itmp[k] != -1) {
          continue;
        }

        /*** Mark this spot as part of a basis function ***/
        itmp[k] = basis.nbss;
        r2ijk = r2.map[i][j][k];

        /*** Look for symmetry-related points ***/
        for (ii = 0; ii < nivec; ii++) {
          iloc = ivec[ii];
          for (jj = 0; jj < njvec; jj++) {
            jloc = jvec[jj];
            for (kk = 0; kk < nkvec; kk++) {
              kloc = kvec[kk];
	      if (bn.map[iloc][jloc][kloc] >= 0) {
		continue;
	      }
              if (fabs(r2.map[iloc][jloc][kloc] - r2ijk) < 1.0e-9) {
                bn.map[iloc][jloc][kloc] = basis.nbss;
              }
              if (ng[0] == ng[1] && 
                  fabs(r2.map[jloc][iloc][kloc] - r2ijk) < 1.0e-9) {
                bn.map[jloc][iloc][kloc] = basis.nbss;
              }
              if (ng[0] == ng[2] && 
                  fabs(r2.map[kloc][jloc][iloc] - r2ijk) < 1.0e-9) {
                bn.map[kloc][jloc][iloc] = basis.nbss;
              }
              if (ng[1] == ng[2] && 
                  fabs(r2.map[iloc][kloc][jloc] - r2ijk) < 1.0e-9) {
                bn.map[iloc][kloc][jloc] = basis.nbss;
              }
              if (ng[0] == ng[2] && ng[1] == ng[2]) {
                if (fabs(r2.map[kloc][iloc][jloc] - r2ijk) < 1.0e-9) {
                  bn.map[kloc][iloc][jloc] = basis.nbss;
                }
                if (fabs(r2.map[jloc][kloc][iloc] - r2ijk) < 1.0e-9) {
                  bn.map[jloc][kloc][iloc] = basis.nbss;
                }
              }
            }
          }
        }
        basis.nbss += 1;
      }
    }
  }

  /*** Record coordinates for all basis functions ***/
  basis.bsslim = (int*)calloc(basis.nbss+1, sizeof(int));
  basis.bsscrd = (int*)malloc(3*ng[0]*ng[1]*ng[2]*sizeof(int));
  bcon = (int*)calloc(basis.nbss, sizeof(int));
  for (i = 0; i < ng[0]; i++) {
    for (j = 0; j < ng[1]; j++) {
      itmp = bn.map[i][j];
      for (k = 0; k < ng[2]; k++) {
        basis.bsslim[itmp[k]+1] += 1;
      }
    }
  }
  j = 0;
  for (i = 0; i <= basis.nbss; i++) {
    j += basis.bsslim[i];
    basis.bsslim[i] = j;
  }
  for (i = 0; i < ng[0]; i++) {
    for (j = 0; j < ng[1]; j++) {
      itmp = bn.map[i][j];
      for (k = 0; k < ng[2]; k++) {
        ii = basis.bsslim[itmp[k]] + bcon[itmp[k]];
        basis.bsscrd[3*ii] = i;
        basis.bsscrd[3*ii+1] = j;
        basis.bsscrd[3*ii+2] = k;
        bcon[itmp[k]] += 1;
      }
    }
  }

  /*** Free allocated memory ***/
  DestroyIbook(&bn);
  DestroyDbook(&r2);
  free(bcon);
  free(mx);
  free(my);
  free(mz);

  return basis;
}
