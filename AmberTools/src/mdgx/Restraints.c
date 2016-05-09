#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/***=======================================================================***/
/*** GridRestraint: compute restraint forces on atoms due to a grid-based  ***/
/***                potential.                                             ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***                                                                       ***/
/***                                                                       ***/
/***=======================================================================***/
void GridRestraints(cellgrid *CG, prmtop *tp, dbook *kgrd)
{

  /*** Loop over all cells and interpolate ***/
  /*** the restraint force on each atom    ***/
  for (i = 0; i < CG->MyCellCount; i++) {

  }
}

/***=======================================================================***/
/*** SymmetryOnGrid: for one specific grid point, compute the set of       ***/
/***                 symmetry-related points and their weights by a cubic  ***/
/***                 interpolation scheme.                                 ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   kgrd:     the grid                                                  ***/
/***   kspc:     the spacing of the grid, distance between neighboring     ***/
/***             points                                                    ***/
/***   p[x,y,z]: the grid location of the initial point                    ***/
/***   xfrm:     the set of symmetry operations                            ***/
/***   crd:      the coordinates (for box transformation matrices)         ***/
/***=======================================================================***/
void SymmetryOnGrid(fbook *kgrd, dmat *kspc, double* korig, int px, int py,
                    int pz, symop *xfrm, coord *crd)
{
  int i, j, k;
  double prx, pry, prz;
  double sr[3];
  double *Utmp, *spctmp;
  dbook A;

  /*** Compute the symmetry-related points ***/
  spctmp = &kspc->data;
  ptwt = 1.0/xfrm->nt;
  for (i = 0; i < xfrm->nt; i++) {

    /*** Check for the identity transform and skip it ***/
    if (fabs(spctmp[0] - 1.0) < 1.0e-8 && fabs(spctmp[1]) < 1.0e-8 &&
	fabs(spctmp[2]) < 1.0e-8 && fabs(spctmp[3]) < 1.0e-8 &&
	fabs(spctmp[4] - 1.0) < 1.0e-8 && fabs(spctmp[5]) < 1.0e-8 &&
	fabs(spctmp[6]) < 1.0e-8 && fabs(spctmp[7]) < 1.0e-8 &&
	fabs(spctmp[8] - 1.0) < 1.0e-8) {
      continue;
    }

    /*** Compute the symmetry-related point location ***/
    prx = px*spctmp[0] + py*spctmp[1] + pz*spctmp[2] + korig[0];
    pry = px*spctmp[3] + py*spctmp[4] + pz*spctmp[5] + korig[1];
    prz = px*spctmp[6] + py*spctmp[7] + pz*spctmp[8] + korig[2];
    Utmp = &xfrm->rmat[i].data;
    sr[0] = prx*Utmp[0] + pry*Utmp[1] + prz*Utmp[2] + xfrm->tvec[3*i];
    sr[1] = prx*Utmp[3] + pry*Utmp[4] + prz*Utmp[5] + xfrm->tvec[3*i+1];
    sr[2] = prx*Utmp[6] + pry*Utmp[7] + prz*Utmp[8] + xfrm->tvec[3*i+2];

    /*** For each symmetry-related point, ***/
    /*** compute the interpolation domain ***/
    Interp2Grid(sr, ptwt, A, crd->U, crd->invU);
  }

}

/***=======================================================================***/
/*** SymmetricPotential: with the potential (or restraint potential) due   ***/
/***                     to the asymmetric unit plotted, spread the        ***/
/***                     potential throughout the unit cell according to   ***/
/***                     symmetry operations.                              ***/
/***=======================================================================***/
void SymmetricPotential(fbook *kgrd, dmat *kspc, double* korig, symop *xfrm)
{
  fbook ngrd;

  /*** Create a new grid to accumulate the symmetry-related potential ***/
  ngrd = CreateFBook(kgrd->npag, kgrd->nrow, kgrd->ncol);

  for (i = 0; i < kgrd->npag; i++) {
    for (j = 0; j < kgrd->nrow; j++) {
      for (k = 0; k < kgrd->ncol; k++) {
	SymmetryOnGrid(kgrd, kspc, korig, i, j, k, xfrm);
      }
    }
  }
}
