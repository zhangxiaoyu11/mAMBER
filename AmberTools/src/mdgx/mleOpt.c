#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "Grid.h"
#include "BSpline.h"
#include "mdgxVector.h"
#include "pmeRecip.h"
#include "mleOpt.h"
#include "mleRecip.h"
#include "ChargeMap.h"
#include "fftw3.h"

#include "MatrixDS.h"

/***=======================================================================***/
/*** PrepareLayerTargets: this function selects a number of layers from a  ***/
/***                      charge grid to use in generating targets for     ***/
/***                      optimizing a coarsened MLE reciprocal space pair ***/
/***                      potential.                                       ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   Uact:   the fine-resolution, "actual" reciprocal space pair         ***/
/***             potential                                                 ***/
/***   Qact:   the fine-resolution, "actual" charge mesh                   ***/
/***   Qcrs:   the coarsened charge mesh (pre-computed for efficiency)     ***/
/***   nlyr:   the number of the fine-resolution mesh layer to be          ***/
/***             approximated; nlyr has an acceptable range [0, Uact.pag)  ***/
/***   Utrg:   the electrostatic potential targets (pre-allocated array)   ***/
/***             obtained by convoluting the charge mesh Qact with the     ***/
/***             potential Uact                                            ***/
/***   Qtrg:   the charge mesh layers to be used to obtain electrostatic   ***/
/***             potential targets (pre-allocated array), obtained by      ***/
/***             extracting a layer from Qact and then interpolating it    ***/
/***             to a coarsened mesh (note that Uact and Qact are inputs,  ***/
/***             and both have the finest possible resolution, whereas     ***/
/***             Utrg and Qcrs are outputs, with different resolutions     ***/
/***             because Qcrs is pre-interpolated for speed)               ***/
/***   rcinp:  reciprocal space control data                               ***/
/***=======================================================================***/
void PrepareLayerTargets(reccon *rcinp, dbook *Qact, dbook *Qcrs, int nlyr,
			 dmat* Utrg, dmat* Qtrg)
{
  int i, j, k, ntrg;
  int* UseAsTarget;
  double invpag, invi, qureal, quimag, normfac, bfac, byfac;
  double* By;
  double* Bz;
  dmat pUact;
  dbook *Uact;
  fftw_plan forwp, backp;
  fftw_complex *cmplx, *cmplx2;

  /*** At this stage, the complete, fine resolution reciprocal space pair ***/
  /*** potential is available, but the other elements of the Urec array   ***/
  /*** have not been computed.                                            ***/
  Uact = &rcinp->Urec[rcinp->nlev];

  /*** Make a copy of the layer of interest ***/
  pUact = CreateDmat(Uact->row, Uact->col, 1);
  cmplx = pUact.fdata;
  forwp = fftw_plan_dft_r2c_2d(pUact.row, pUact.col, pUact.data, pUact.fdata,
			       FFTW_ESTIMATE);
  ExtractDpage(Uact, &pUact, nlyr, 1);
  fftw_execute(forwp);
  fftw_destroy_plan(forwp);

  /*** Select a series of targets ***/
  UseAsTarget = (int*)calloc(Uact->pag, sizeof(int));
  i = 0;
  ntrg = rcinp->ntrg;
  if (ntrg > Uact->pag) {
    ntrg = Uact->pag;
  }
  invpag = 1.0/Uact->pag;
  j = 0;
  for (k = 0; k < 2; k++) {
    for (i = 0; i < Uact->pag; i++) {
      invi = (i == 0) ? 0.0 : 1.0/i;
      if (j*invi < ntrg*invpag) {
        UseAsTarget[i] = 1;
        j++;
      }
    }
  }
  j = 0;
  for (i = 0; i < Uact->pag; i++) {
    if (UseAsTarget[i] == 1) {
      UseAsTarget[j] = i;
      j++;
    }
  }

  /*** Normalization for the convolutions ***/
  normfac = 1.0 / (pUact.row * pUact.col);

  /*** Create a series of electrostatic potential targets ***/
  for (i = 0; i < ntrg; i++) {
    Utrg[i] = CreateDmat(Uact->row, Uact->col, 1);
    Utrg[i].col = Uact->col;
    forwp = fftw_plan_dft_r2c_2d(Utrg[i].row, Utrg[i].col, Utrg[i].data,
				 Utrg[i].fdata, FFTW_ESTIMATE);
    backp = fftw_plan_dft_c2r_2d(Utrg[i].row, Utrg[i].col, Utrg[i].fdata,
				 Utrg[i].data, FFTW_ESTIMATE);
    cmplx2 = Utrg[i].fdata;
    ExtractDpage(Qact, &Utrg[i], UseAsTarget[i], 1);
    fftw_execute(forwp);
    for (j = 0; j < pUact.row*(pUact.col/2+1); j++) {
      qureal = cmplx[j][0]*cmplx2[j][0] - cmplx[j][1]*cmplx2[j][1];
      quimag = cmplx[j][0]*cmplx2[j][1] + cmplx[j][1]*cmplx2[j][0];
      cmplx2[j][0] = qureal*normfac;
      cmplx2[j][1] = quimag*normfac;
    }
    fftw_execute(backp);
    fftw_destroy_plan(forwp);
    fftw_destroy_plan(backp);
  }

  /*** Extract the layers of the charge mesh corresponding to the    ***/
  /*** electrostatic potential targets, then interpolate them to the ***/
  /*** higher level mesh, apply the forward FFT, and fold in the     ***/
  /*** Euler exponential splines and normalization factor in         ***/
  /*** preparation for fitting.                                      ***/
  By = LoadPrefac(rcinp->ggordr, Qcrs->row);
  Bz = LoadPrefac(rcinp->ggordr, Qcrs->col);
  normfac = 1.0 / (Qcrs->row * Qcrs->col);
  for (i = 0; i < ntrg; i++) {
    Qtrg[i] = CreateDmat(Qcrs->row, Qcrs->col, 1);
    forwp = fftw_plan_dft_r2c_2d(Qtrg[i].row, Qtrg[i].col, Qtrg[i].data,
				 Qtrg[i].fdata, FFTW_ESTIMATE);
    ExtractDpage(Qcrs, &Qtrg[i], UseAsTarget[i], 1);
    fftw_execute(forwp);
    fftw_destroy_plan(forwp);
    for (j = 0; j < Qtrg[i].row; j++) {
      cmplx = Qtrg[i].fmap[j];
      byfac = By[j]*normfac;
      for (k = 0; k < Qtrg[i].col/2 + 1; k++) {
	bfac = byfac*Bz[k];
        cmplx[k][0] *= bfac;
        cmplx[k][1] *= bfac;
      }
    }
  }
  free(By);
  free(Bz);
}

/***=======================================================================***/
/*** LayerBasisFunctions: this function works much like MeshBasisFunctions ***/
/***                      (see above), but on a specific layer.            ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   pLrec:   the electrostatic potential layer                          ***/
/***   nlyr:    the number of this layer; range [0, maxlyr)                ***/
/***   maxlyr:  the number of this layer; range [0, maxlyr)                ***/
/***   crd:     the coordinates (used to get unit cell dimensions)         ***/
/***=======================================================================***/
bssf LayerBasisFunctions(dmat *pLrec, int nlyr, int maxlyr, coord *crd)
{
  int i, j, nivec, njvec, ii, jj, iloc, jloc;
  int *itmp;
  int ng[2], ivec[2], jvec[2];
  int* bcon;
  double dx, dy, dz, mmx, mmy, mmz, r2ij;
  double *dtmp, *ivd;
  double* my;
  double* mz;
  imat bn;
  dmat r2;
  bssf basis;

  /*** First, compute the distance of each point from the origin ***/
  ng[0] = pLrec->row;
  ng[1] = pLrec->col;
  dx = (nlyr <= maxlyr/2) ? nlyr : nlyr-maxlyr; 
  my = LoadMVec(ng[0]);
  mz = LoadMVec(ng[1]);
  r2 = CreateDmat(ng[0], ng[1], 0);
  ivd = crd->invU.data;
  for (i = 0; i < ng[0]; i++) {
    dtmp = r2.map[i];
    dy = my[i];
    for (j = 0; j < ng[1]; j++) {
      dz = mz[j];
      mmx = ivd[0]*dx + ivd[1]*dy + ivd[2]*dz;
      mmy = ivd[3]*dx + ivd[4]*dy + ivd[5]*dz;
      mmz = ivd[6]*dx + ivd[7]*dy + ivd[8]*dz;
      dtmp[j] = mmx*mmx + mmy*mmy + mmz*mmz;
    }
  }

  /*** As before (in the mesh basis function computation), ***/
  /*** count the number of basis functions                 ***/
  bn = CreateImat(ng[0], ng[1]);
  SetIVec(bn.data, ng[0]*ng[1], -1);
  basis.nbss = 0;
  for (i = 0; i < ng[0]; i++) {
    itmp = bn.map[i];
    ivec[0] = i;
    if (i == 0) {
      nivec = 1;
    }
    else {
      ivec[1] = ng[0]-i;
      nivec = 2;
    }
    for (j = 0; j < ng[1]; j++) {
      if (bn.map[i][j] != -1) {
	continue;
      }
      jvec[0] = j;
      if (j == 0) {
	njvec = 1;
      }
      else {
	jvec[1] = ng[1]-j;
	njvec = 2;
      }
      itmp[j] = basis.nbss;
      r2ij = r2.map[i][j];
      for (ii = 0; ii < nivec; ii++) {
	iloc = ivec[ii];
	for (jj = 0; jj < njvec; jj++) {
	  jloc = jvec[jj];
	  if (bn.map[iloc][jloc] < 0 &&
	      fabs(r2.map[iloc][jloc] - r2ij) < 1.0e-9) {
	    bn.map[iloc][jloc] = basis.nbss;
	  }
	  if (ng[0] == ng[1] && bn.map[jloc][iloc] < 0 &&
	      fabs(r2.map[jloc][iloc] - r2ij) < 1.0e-9) {
	    bn.map[jloc][iloc] = basis.nbss;
	  }
	}
      }
      basis.nbss += 1;
    }
  }

  /*** Record coordinates for all basis functions ***/
  basis.bsslim = (int*)calloc(basis.nbss+1, sizeof(int));
  basis.bsscrd = (int*)malloc(2*ng[0]*ng[1]*sizeof(int));
  bcon = (int*)calloc(basis.nbss, sizeof(int));
  for (i = 0; i < ng[0]; i++) {
    itmp = bn.map[i];
    for (j = 0; j < ng[1]; j++) {
      basis.bsslim[itmp[j]+1] += 1;
    }
  }
  j = 0;
  for (i = 0; i <= basis.nbss; i++) {
    j += basis.bsslim[i];
    basis.bsslim[i] = j;
  }
  for (i = 0; i < ng[0]; i++) {
    itmp = bn.map[i];
    for (j = 0; j < ng[1]; j++) {
      ii = basis.bsslim[itmp[j]] + bcon[itmp[j]];
      basis.bsscrd[2*ii] = i;
      basis.bsscrd[2*ii+1] = j;
      bcon[itmp[j]] += 1;
    }
  }

  /*** Free allocated memory ***/
  DestroyImat(&bn);
  DestroyDmat(&r2);
  free(bcon);
  free(my);
  free(mz);

  return basis;
}

/***=======================================================================***/
/*** DestroyBasis: free allocated memory for a basis function structure;   ***/
/***               works for both and layer basis functions.               ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   basis:  the basis function set to free                              ***/
/***=======================================================================***/
static void DestroyBasis(bssf *basis)
{
  free(basis->bsslim);
  free(basis->bsscrd);
}

/***=======================================================================***/
/*** ConvLayer: convolute a series of layers of the charge mesh with a     ***/
/***            putative reciprocal space pair potential (i.e. one of its  ***/
/***            basis functions).                                          ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***                                                                       ***/
/***=======================================================================***/
static void ConvLayer(dmat *ThtBSS, dmat* Qtrg, dmat* UCguess, dmat* Uguess,
		      fftw_plan *forwp, fftw_plan* backp, reccon *rcinp,
		      int lcon)
{
  int h, i, j;
  fftw_complex *qtmp, *utmp, *qutmp;
  dbook UCbook, Ubook;

  /*** Tranform the putative pair potential forward ***/
  fftw_execute(*forwp);

  /*** Convolute each of the charge mesh layers ***/
  const int jlim = ThtBSS->col/2 + 1;
  for (h = 0; h < rcinp->ntrg; h++) {
    for (i = 0; i < ThtBSS->row; i++) {
      qtmp = Qtrg[h].fmap[i];
      utmp = ThtBSS->fmap[i];
      qutmp = UCguess[h].fmap[i];
      for (j = 0; j < jlim; j++) {
	qutmp[j][0] = qtmp[j][0]*utmp[j][0] - qtmp[j][1]*utmp[j][1];
	qutmp[j][1] = qtmp[j][0]*utmp[j][1] + qtmp[j][1]*utmp[j][0];
      }
    }
    fftw_execute(backp[h]);

    /*** Fake a dbook structure out of the array of dmats.    ***/
    /*** This is done to save memory and time that would      ***/
    /*** otherwise be needed to shuffle around all that data. ***/
    UCbook = Dmat2Dbook(&UCguess[h]);
    Ubook = Dmat2Dbook(&Uguess[h]);
    IntrpBook(Ubook, UCbook, rcinp->SPrv[lcon], rcinp->SPcv[lcon]);

    /*** Free memory allocated for the fake book structures ***/
    free(UCbook.map);
    free(UCbook.fmap);
    free(Ubook.map);
    free(Ubook.fmap);
  }
}

/***=======================================================================***/
/*** LayerRMSE: compute the RMS error between a guess of the electrostatic ***/
/***            potential and the actual potential.  The guess is composed ***/
/***            of two parts, a base and a delta, with a scaling factor    ***/
/***            for the delta.                                             ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   base:    the base guesses to the electrostatic potential            ***/
/***   delta:   the delta to the electrostatic potential guesses           ***/
/***   scl:     the scaling factor for the delta; the guesses as to the    ***/
/***              electrostatic potential are, in all:                     ***/
/***                                                                       ***/
/***                GUESS = BASE + SCL * DELTA                             ***/
/***                                                                       ***/
/***   actual:  the actual electrostatic potentials                        ***/
/***   ntrg:    the number of targets                                      ***/
/***=======================================================================***/
static double LayerRMSE(dmat* base, dmat* delta, double scl, dmat* actual,
			int ntrg)
{
  int h, i, j;
  double *btmp, *dtmp, *atmp;
  double dx, rmserr;

  const int ilim = base[0].row;
  const int jlim = base[0].col;
  rmserr = 0.0;
  for (h = 0; h < ntrg; h++) {
    for (i = 0; i < ilim; i++) {
      btmp = base[h].map[i];
      dtmp = delta[h].map[i];
      atmp = actual[h].map[i];
      for (j = 0; j < jlim; j++) {
	dx = atmp[j] - (btmp[j] + scl*dtmp[j]);
	rmserr += dx*dx;
      }
    }
  }

  return sqrt(rmserr/(base->row*base->col*ntrg));
}

/***=======================================================================***/
/*** LayerGradient: compute the gradient in the approximate reciprocal     ***/
/***                space pair potential score based on the accuracy of    ***/
/***                the electrostatic potential produced for several       ***/
/***                target layers of the charge mesh.                      ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***                                                                       ***/
/***=======================================================================***/
static double LayerGradient(dmat* Uguess, dmat* dUguess, dmat* Utrg, int ntrg)
{
  double gpscr, gmscr;

  gpscr = LayerRMSE(Uguess, dUguess, 1.0, Utrg, ntrg);
  gmscr = LayerRMSE(Uguess, dUguess, -1.0, Utrg, ntrg);
 
  return gpscr - gmscr;
}

/***=======================================================================***/
/*** EvalLayerBasis: compute the gradient of the reciprocal space          ***/
/***                 electrostatic potential score function obtained by    ***/
/***                 convoluting layers of the charge mesh with layers of  ***/
/***                 the coarsened electrostatic potential.                ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   ThtI:   a guess for the reciprocal space pair potential layer       ***/
/***   Qtrg:   the coarsened charge mesh layers, represented in Fourier    ***/
/***             space with Euler exponential spline coefficients folded   ***/
/***             in for convenience and efficiency                         ***/
/***   Utrg:   the fine electrostatic potential computed using the true    ***/
/***             pair potential, represented in real space                 ***/
/***   lyrbss:  the basis functions that can be used to optimize ThtI      ***/
/***   grad:   the gradients of the scoring function for all basis         ***/
/***             functions (pre-allocated array with one more spaces than  ***/
/***             there are basis functions; the gradient in ThtI itself is ***/
/***             stored in the last element)                               ***/
/***   rcinp: reciprocal space control data                                ***/
/***   lcon:  the level of the reciprocal space mesh                       ***/
/***   dg:    testing increment for finite difference gradient computation ***/
/***=======================================================================***/
static double EvalLayerBasis(dmat *ThtI, dmat* Qtrg, dmat* Utrg, bssf *lyrbss,
			     double* grad, reccon *rcinp, int lcon, double dg)
{
  int i, j;
  int *tlim, *tcrd;
  double rmse0, rmseN, rmseB, tstep, dstep;
  dmat ThtBSS;
  dmat* UCguess;
  dmat* Uguess;
  dmat* dUguess;
  fftw_plan forwp;
  fftw_plan* backp;

  /*** Unpack structures ***/
  tlim = lyrbss->bsslim;
  tcrd = lyrbss->bsscrd;
  const int nbss = lyrbss->nbss;

  /*** Prepare scratch space and plans ***/
  ThtBSS = CreateDmat(ThtI->row, ThtI->col, 1);
  forwp = fftw_plan_dft_r2c_2d(ThtBSS.row, ThtBSS.col, ThtBSS.data,
			       ThtBSS.fdata, FFTW_MEASURE);
  UCguess = (dmat*)malloc(rcinp->ntrg*sizeof(dmat));
  Uguess = (dmat*)malloc(rcinp->ntrg*sizeof(dmat));
  dUguess = (dmat*)malloc(rcinp->ntrg*sizeof(dmat));
  backp = (fftw_plan*)malloc(rcinp->ntrg*sizeof(fftw_plan));
  for (i = 0; i < rcinp->ntrg; i++) {
    UCguess[i] = CreateDmat(ThtI->row, ThtI->col, 1);
    Uguess[i] = CreateDmat(Utrg[i].row, Utrg[i].col, 1);
    dUguess[i] = CreateDmat(Utrg[i].row, Utrg[i].col, 1);
    backp[i] = fftw_plan_dft_c2r_2d(UCguess[i].row, UCguess[i].col,
				    UCguess[i].fdata, UCguess[i].data,
				    FFTW_MEASURE);
  }

  /*** Compute the potential with ThtI ***/
  CopyDmat(&ThtBSS, ThtI, 1);
  ConvLayer(&ThtBSS, Qtrg, UCguess, Uguess, &forwp, backp, rcinp, lcon);
  rmse0 = LayerRMSE(Uguess, Uguess, 0.0, Utrg, rcinp->ntrg);

  /*** Loop over all basis functions ***/
  for (i = 0; i < nbss; i++) {

    /*** Make the basis function ***/
    SetDVec(ThtBSS.data, 2*ThtBSS.row*(ThtBSS.col/2+1), 0.0);
    for (j = tlim[i]; j < tlim[i+1]; j++) {
      ThtBSS.map[tcrd[2*j]][tcrd[2*j+1]] = 0.1*dg;
    }

    /*** Compute the gradient in this basis function ***/
    ConvLayer(&ThtBSS, Qtrg, UCguess, dUguess, &forwp, backp, rcinp, lcon);
    grad[i] = LayerGradient(Uguess, dUguess, Utrg, rcinp->ntrg);
  }
  Normalize(grad, nbss);
  for (i = 0; i < nbss; i++) {
    grad[i] *= -dg;
  }

  /*** Compute the new guess based on      ***/
  /*** the gradient in all basis functions ***/
  SetDVec(ThtBSS.data, 2*ThtBSS.row*(ThtBSS.col/2+1), 0.0);
  for (i = 0; i < nbss; i++) {
    for (j = tlim[i]; j < tlim[i+1]; j++) {
      ThtBSS.map[tcrd[2*j]][tcrd[2*j+1]] += grad[i];
    }
  }
  ConvLayer(&ThtBSS, Qtrg, UCguess, dUguess, &forwp, backp, rcinp, lcon);

  /*** Line minimization ***/
  dstep = 1.0;
  tstep = 0.0;
  rmseB = rmse0;
  while (dstep > 0.0625) {
    rmseN = LayerRMSE(Uguess, dUguess, tstep+dstep, Utrg, rcinp->ntrg);
    if (rmseN < rmseB) {
      rmseB = rmseN;
      tstep += dstep;
      dstep *= 1.1;
    }
    else {
      dstep *= 0.5;
    }
  }

  /*** Store the best solution found thus far.  Note that ThtBSS ***/
  /*** must be reset as it has been hit with a forward FFT.      ***/
  SetDVec(ThtBSS.data, 2*ThtBSS.row*(ThtBSS.col/2+1), 0.0);
  for (i = 0; i < nbss; i++) {
    for (j = tlim[i]; j < tlim[i+1]; j++) {
      ThtBSS.map[tcrd[2*j]][tcrd[2*j+1]] += grad[i];
    }
  }
  AddDmat(ThtI, &ThtBSS, 1.0, tstep);

  /*** Free allocated memory ***/
  fftw_destroy_plan(forwp);
  for (i = 0; i < rcinp->ntrg; i++) {
    fftw_destroy_plan(backp[i]);
    DestroyDmat(&UCguess[i]);
    DestroyDmat(&Uguess[i]);
    DestroyDmat(&dUguess[i]);
  }
  free(backp);
  free(UCguess);
  free(Uguess);
  free(dUguess);

  return rmse0 - rmseB;
}

/***=======================================================================***/
/*** OptimizeLayer: this function optimizes one layer of a coarsened       ***/
/***                reciprocal space pair potential.                       ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***                                                                       ***/
/***=======================================================================***/
void OptimizeLayer(reccon *rcinp, int lcon, dbook *Lrec, dbook *Qact,
		   dbook *Qcrs, int nlyr, coord *crd)
{
  int i, niter;
  double derr, dg;
  double* grad;
  dmat pLrec;
  dmat* Utrg;
  dmat* Qtrg;
  bssf lyrbss;

  /*** Make a copy of the layer from Lrec;   ***/
  /*** the matrix is allocated automatically ***/
  ExtractDpage(Lrec, &pLrec, nlyr, 0);

  /*** Prepare a series of targets ("given this layer of the finest      ***/
  /*** resolution charge mesh, produce this electrostatic potential...") ***/
  Utrg = (dmat*)malloc(rcinp->ntrg*sizeof(dmat));
  Qtrg = (dmat*)malloc(rcinp->ntrg*sizeof(dmat));
  PrepareLayerTargets(rcinp, Qact, Qcrs, nlyr, Utrg, Qtrg);

  /*** Decide on a set of basis functions ***/
  lyrbss = LayerBasisFunctions(&pLrec, nlyr, Qcrs->pag, crd);

  /*** Iteratively refine the potential ***/
  niter = 0;
  grad = (double*)calloc(lyrbss.nbss, sizeof(double));
  dg = rcinp->dg0;
  while (niter < rcinp->maxiter && dg > rcinp->mindg) {

    /*** Loop over all basis functions, determine the gradient, perform ***/
    /*** line minimization, and return the improvement in the error.    ***/
    derr = EvalLayerBasis(&pLrec, Qtrg, Utrg, &lyrbss, grad, rcinp, lcon, dg);

    /*** Modify the test increment based on the success of this iteration ***/
    dg = (derr < rcinp->MeshIterTol) ? dg*0.5 : dg*1.1;
    niter++;
  }

  /*** Copy the final result of the optimization back into Lrec ***/
  ScribeDpage(&pLrec, Lrec, nlyr);

  /*** Free allocated memory ***/
  DestroyDmat(&pLrec);
  DestroyBasis(&lyrbss);
  for (i = 0; i < rcinp->ntrg; i++) {
    DestroyDmat(&Utrg[i]);
    DestroyDmat(&Qtrg[i]);
  }
  free(Utrg);
  free(Qtrg);
  free(grad);
}

/***=======================================================================***/
/*** OptimizeMesh: this function optimizes a coarsened reciprocal space    ***/
/***               pair potential, one layer at a time.                    ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   Uact:   the actual reciprocal space pair potential                  ***/
/***   Lrec:   on input, the initial guess as to the coarsened reciprocal  ***/
/***             space pair potential; on output, the optimized coarsened  ***/
/***             reciprocal space pair potential                           ***/
/***   rcinp:  the reciprocal space control data                           ***/
/***   crd:    the system coordinates                                      ***/
/***   tp:     the system topology                                         ***/
/***   lcon:   the level index (some information, such as the size of the  ***/
/***             higher level mesh Lrec, can be obtained from Lrec itself, ***/
/***             but lcon is authoritative on the actual mesh level)       ***/
/***=======================================================================***/
void OptimizeMesh(dbook *Uact, dbook *Lrec, reccon *rcinp, coord *crd,
		  prmtop *tp, int lcon)
{
  int i, llim1, llim2, hlim1, hlim2;
  bmap pmmap;
  dbook Qact, Qcrs;

  /*** Compute the charge mesh and interpolate ***/
  /*** it to the desired (higher) level        ***/
  pmmap = CreateBmap(rcinp, crd->natom);
  SplCoeff(crd, &pmmap, rcinp);
  Qact = CreateDbook(Uact->pag, Uact->row, Uact->col, 1);
  Qcrs = CreateDbook(Lrec->pag, Lrec->row, Lrec->col, 1);
  Charges2Book(&pmmap, rcinp, tp, &Qact);
  SpreadBook(Qcrs, Qact, rcinp->SPrv[lcon], rcinp->SPcv[lcon]);

  /*** Optimize each layer, one by one ***/
  for (i = 0; i < Lrec->pag; i++) {
    MeshLevelLimits(rcinp, lcon, &llim1, &hlim1, &llim2, &hlim2);
    if ((i >= llim1 && i < hlim1) || (i >= llim2 && i < hlim2)) {
      OptimizeLayer(rcinp, lcon, Lrec, &Qact, &Qcrs, i, crd);
    }
  }

  /*** Free allocated memory ***/
  DestroyBmap(&pmmap);
  DestroyDbook(&Qact);
  DestroyDbook(&Qcrs);
}
