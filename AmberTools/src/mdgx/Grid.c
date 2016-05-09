#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Matrix.h"
#include "Grid.h"

/***=======================================================================***/
/*** CreateIbook: create an integer book of dimensions M x N x P.          ***/
/***=======================================================================***/
ibook CreateIbook(int M, int N, int P)
{
  int i, j;
  ibook A;

  /*** Allocate memory ***/
  A.data = (int*)calloc(M*N*P, sizeof(int));

  /*** Make the map ***/
  A.map = (int***)malloc(M*sizeof(int**));
  for (i = 0; i < M; i++) {
    A.map[i] = (int**)malloc(N*sizeof(int*));
    for (j = 0; j < N; j++) {
      A.map[i][j] = &A.data[i*N*P + j*P];
    }
  }

  /*** Store the dimensionality ***/
  A.pag = M;
  A.row = N;
  A.col = P;

  return A;
}

/***=======================================================================***/
/*** CreateDbook: create a double-precision real book with M x N x P       ***/
/***              indices.  There is added functionality for preparing the ***/
/***              array for 1, 2, and 3-dimensional real-to-complex (R2C)  ***/
/***              Discrete Fourier Transforms (DFTs).                      ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   M:       the number of pages (the most slowly varying index)        ***/
/***   N:       the number of rows (the second most slowly varying index)  ***/
/***   P:       the number of columns (the most rapidly varying index)     ***/
/***   prepFFT: flag to request preparation for R2C transforms--setting    ***/
/***            prepFFT=1 pads along the 3rd dimension.  This will prepare ***/
/***            the array for a single 3D R2C DFT, as well as a series of  ***/
/***            2D or 1D R2C DFTs.                                         ***/
/***=======================================================================***/
dbook CreateDbook(int M, int N, int P, int prepFFT)
{
  int i, j, Pp;
  dbook A;

  /*** Prepare this book for 3D FFTs ***/
  A.pfft = prepFFT;
  Pp = (prepFFT == 1) ? 2*(P/2+1) : P;

  /*** Allocate memory ***/
  A.data = (double*)calloc(M*N*Pp, sizeof(double));

  /*** This is pointer-array abuse, but I know of no other way to do it ***/
  if (prepFFT > 0) {
    A.fdata = (fftw_complex*)A.data;
  }

  /*** Make the map ***/
  A.map = (double***)malloc(M*sizeof(double**));
  for (i = 0; i < M; i++) {
    A.map[i] = (double**)malloc(N*sizeof(double*));
    for (j = 0; j < N; j++) {
      A.map[i][j] = &A.data[(i*N+j)*Pp];
    }
  }

  /*** Set up the complex maps ***/
  if (prepFFT == 1) {
    A.fmap = (fftw_complex***)malloc(M*sizeof(fftw_complex**));
    for (i = 0; i < M; i++) {
      A.fmap[i] = (fftw_complex**)malloc(N*sizeof(fftw_complex*));
      for (j = 0; j < N; j++) {
	A.fmap[i][j] = &A.fdata[(i*N+j)*Pp/2];
      }
    }
  }

  /*** Store the dimensionality ***/
  A.pag = M;
  A.row = N;
  A.col = P;

  return A;
}

/***=======================================================================***/
/*** CreateFbook: create a single-precision real book with M x N x P       ***/
/***              indices.  This is the format for grids read from files,  ***/
/***              or grids that represent potential fields for which there ***/
/***              is not a ton of memory.                                  ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   M:       the number of pages (the most slowly varying index)        ***/
/***   N:       the number of rows (the second most slowly varying index)  ***/
/***   P:       the number of columns (the most rapidly varying index)     ***/
/***=======================================================================***/
fbook CreateFbook(int M, int N, int P)
{
  int i, j;
  fbook A;

  /*** Allocate memory ***/
  A.data = (float*)calloc(M*N*P, sizeof(float));

  /*** Make the map ***/
  A.map = (float***)malloc(M*sizeof(float**));
  for (i = 0; i < M; i++) {
    A.map[i] = (float**)malloc(N*sizeof(float*));
    for (j = 0; j < N; j++) {
      A.map[i][j] = &A.data[(i*N+j)*P];
    }
  }

  /*** Store the dimensionality ***/
  A.pag = M;
  A.row = N;
  A.col = P;

  /*** Allocate the spacing ***/
  A.lvec = CreateDmat(3, 3, 0);

  return A;
}

/***=======================================================================***/
/*** CreateCbook: create a character (or unsigned 8-bit integer) book of   ***/
/***              dimensions M x N x P.                                    ***/
/***=======================================================================***/
cbook CreateCbook(int M, int N, int P)
{
  int i, j;
  cbook A;

  /*** Allocate memory ***/
  A.data = (char*)calloc(M*N*P, sizeof(char));

  /*** Make the map ***/
  A.map = (char***)malloc(M*sizeof(char**));
  for (i = 0; i < M; i++) {
    A.map[i] = (char**)malloc(N*sizeof(char*));
    for (j = 0; j < N; j++) {
      A.map[i][j] = &A.data[i*N*P + j*P];
    }
  }

  /*** Store the dimensionality ***/
  A.pag = M;
  A.row = N;
  A.col = P;

  return A;
}

/***=======================================================================***/
/*** DestroyIbook: destroy an integer book.                                ***/
/***=======================================================================***/
void DestroyIbook(ibook *A)
{
  int i;

  /*** Free the map ***/
  for (i = 0; i < A->pag; i++) {
    free(A->map[i]);
  }
  free(A->map);

  /*** Free the data ***/
  free(A->data);
}

/***=======================================================================***/
/*** DestroyFbook: destroy a single-precision real book.                   ***/
/***=======================================================================***/
void DestroyFbook(fbook *A)
{
  int i;

  /*** Free the map ***/
  for (i = 0; i < A->pag; i++) {
    free(A->map[i]);
  }
  free(A->map);

  /*** Free the grid spacing matrix ***/
  DestroyDmat(&A->lvec);

  /*** Free the data ***/
  free(A->data);
}

/***=======================================================================***/
/*** DestroyDbook: destroy a double-precision real book.                   ***/
/***=======================================================================***/
void DestroyDbook(dbook *A)
{
  int i;

  /*** Free the map ***/
  for (i = 0; i < A->pag; i++) {
    free(A->map[i]);
  }
  free(A->map);
  if (A->pfft == 1) {
    for (i = 0; i < A->pag; i++) {
      free(A->fmap[i]);
    }
    free(A->fmap);
  }

  /*** Free the data ***/
  free(A->data);
}

/***=======================================================================***/
/*** DestroyCbook: destroy a character book.                               ***/
/***=======================================================================***/
void DestroyCbook(cbook *A)
{
  int i;

  /*** Free the map ***/
  for (i = 0; i < A->pag; i++) {
    free(A->map[i]);
  }
  free(A->map);

  /*** Free the data ***/
  free(A->data);
}

/***=======================================================================***/
/*** CompressQ: convert a 64-bit double-precision real book to a 32-bit    ***/
/***            integer book, with a list of exceptions needing special    ***/
/***            treatment.  The integers range from roughly -2.1e9 to      ***/
/***            +2.1e9, so charges are stored as integers in the range     ***/
/***            -2.1 to +2.1, with nine decimal places of precision.  The  ***/
/***            vast majority grid points will fall in this category.  A   ***/
/***            separate list enumerates the few exceptional cases.        ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   cQ: the compressed charge book                                      ***/
/***   Q:  the uncompressed charge book                                    ***/
/***                                                                       ***/
/*** Note that the compressed charge book is assumed to be pre-allocated,  ***/
/*** as this routine is expected to be used many times to convert the same ***/
/*** charge book into a compressed version.  See also DecompressQ in this  ***/
/*** library.                                                              ***/
/***=======================================================================***/
void CompressQ(qbook *cQ, dbook *Q)
{
  int i, j, M, N, P;
  int *Qcomp;
  double *Qorig;

  /*** Set the exception counter to zero ***/
  cQ->nexcp = 0;

  /*** Set pointers and abbreviate constants ***/
  M = cQ->qdata.pag;
  N = cQ->qdata.row;
  P = cQ->qdata.col;
  Qcomp = cQ->qdata.data;
  Qorig = Q->data;

  /*** Checkbook dimensions ***/
  if (M != Q->pag || N != Q->row || P != Q->col) {
    printf("CompressQ >> Error.  Book dimensions are inconsistent.\n"
	   "CompressQ >> Original book   = [ %4d pages of %4d x %4d ]\n"
	   "CompressQ >> Compressed book = [ %4d pages of %4d x %4d ]\n",
	   M, N, P, Q->pag, Q->row, Q->col);
    exit(1);
  }

  /*** Loop over all data ***/
  for (i = 0; i < M*N*P; i++) {
    if (Qorig[i] >= -2.14748364 && Qorig[i] <= 2.14748364) {
      Qcomp[i] = Qorig[i]*1.0e9;
    }
    else {
      Qcomp[i] = 0.0;
      j = cQ->nexcp;
      cQ->eidx[j] = i;
      cQ->eval[j] = Qorig[i];
      cQ->nexcp = ++j;
      if (j > cQ->maxexcp) {
	cQ->maxexcp += 1024;
	cQ->eidx = (int*)realloc(cQ->eidx, cQ->maxexcp*sizeof(int));
	cQ->eval = (double*)realloc(cQ->eval, cQ->maxexcp*sizeof(double));
      }
    }
  }
}

/***=======================================================================***/
/*** DecompressQ: convert a 32-bit integer book back into a 64-bit         ***/
/***              double-precision real book, abiding a list of exceptions ***/
/***              showing grid points that contain lots of charge.         ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   Q:  the uncompressed charge book                                    ***/
/***   cQ: the compressed charge book                                      ***/
/***                                                                       ***/
/*** As in CompressQ, the uncompressed charge book must be pre-allocated;  ***/
/*** this routine is designed to be used many times with the same data.    ***/
/*** structures.                                                           ***/
/***=======================================================================***/
void DecompressQ(dbook *Q, qbook *cQ)
{
  int i, M, N, P;
  int *Qcomp;
  double *Qnew;

  /*** Set pointers and abbreviate constants ***/
  M = cQ->qdata.row;
  N = cQ->qdata.col;
  P = cQ->qdata.pag;
  Qcomp = cQ->qdata.data;
  Qnew = Q->data;

  /*** Loop over all data ***/
  for (i = 0; i < M*N*P; i++) {
    Qnew[i] = Qcomp[i];
  }

  /*** Loop over exceptions ***/
  for (i = 0; i < cQ->nexcp; i++) {
    Qnew[cQ->eidx[i]] = cQ->eval[i];
  }
}

/***=======================================================================***/
/*** ExtractDpage: extract a page from a book and copy it into a dmat      ***/
/***               struct.                                                 ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   B:   the book of double-precision real numbers                      ***/
/***   P:   the matrix of double-precision real numbers                    ***/
/***   n:   the page number to extract                                     ***/
/***   ll:  flag to tell whether the page is pre-allocated (1) or must be  ***/
/***          allocated (0)                                                ***/
/***=======================================================================***/
void ExtractDpage(dbook *B, dmat *P, int n, int ll)
{
  int i, j;
  double *btmp, *ptmp;

  /*** Allocate the matrix (page) if necessary.  Note that the page will ***/
  /*** be prepared for FFTs if the original book is set up that way.     ***/
  if (ll == 0) {
    *P = CreateDmat(B->row, B->col, B->pfft);
  }

  /*** Check ***/
  if (B->row > P->row) {
    printf("ExtractDpage >> Error.  %d rows in the book page, %d rows "
	   "available in the matrix.\n", B->row, P->row);
    exit(1);
  }
  if (B->col > P->col) {
    printf("ExtractDpage >> Error.  %d columns in the book page, %d columns "
	   "available in the matrix.\n", B->col, P->col);
    exit(1);
  }

  /*** Copy values ***/
  for (i = 0; i < B->row; i++) {
    btmp = B->map[n][i];
    ptmp = P->map[i];
    for (j = 0; j < B->col; j++) {
      ptmp[j] = btmp[j];
    }
  }
}

/***=======================================================================***/
/*** ScribeDpage: copy a double-precision matrix into a dbook struct.      ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   P:   the matrix of double-precision real numbers                    ***/
/***   B:   the book of double-precision real numbers                      ***/
/***   n:   the page number to which the matrix values should be written   ***/
/***=======================================================================***/
void ScribeDpage(dmat *P, dbook *B, int n)
{
  int i, j;
  double *btmp, *ptmp;

  /*** Check ***/
  if (P->row > B->row) {
    printf("ScribeDpage >> Error.  %d rows in the book page, %d rows "
           "available in the matrix.\n", P->row, B->row);
    exit(1);
  }
  if (P->col > B->col) {
    printf("ScribeDpage >> Error.  %d columns in the book page, %d columns "
           "available in the matrix.\n", P->col, B->col);
    exit(1);
  }

  /*** Copy values ***/
  for (i = 0; i < P->row; i++) {
    ptmp = P->map[i];
    btmp = B->map[n][i];
    for (j = 0; j < P->col; j++) {
      btmp[j] = ptmp[j];
    }
  }
}

/***=======================================================================***/
/*** Dmat2Dbook: promote a double-precision M x N real matrix into a       ***/
/***             1 x M x N double precision real book.  FFT preparations   ***/
/***             are preserved.  No additional data is actually allocated  ***/
/***             for the book--its single page simply points directly to   ***/
/***             the matrix data.  Therefore, the returned values of this  ***/
/***             function should not be freed with the usual destructor.   ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   A:  the matrix to promote                                           ***/
/***=======================================================================***/
dbook Dmat2Dbook(dmat *A)
{
  dbook Ab;

  Ab.pag = 1;
  Ab.row = A->row;
  Ab.col = A->col;
  Ab.data = A->data;
  Ab.pfft = A->pfft;
  Ab.map = (double***)malloc(sizeof(double**));
  Ab.map[0] = A->map;
  if (Ab.pfft == 1) {
    Ab.fdata = A->fdata;
    Ab.fmap = (fftw_complex***)malloc(sizeof(fftw_complex**));
    Ab.fmap[0] = A->fmap;
  }

  return Ab;
}

/***=======================================================================***/
/*** SumDbook: compute the sum of a double-precision real book, accounting ***/
/***           for possible 3D-FFT padding.                                ***/
/***=======================================================================***/
double SumDbook(dbook *A)
{
  int i, j, k;
  double qs;
  double *dtmp;

  qs = 0.0;
  const int klim = A->col;
  for (i = 0; i < A->pag; i++) {
    for (j = 0; j < A->row; j++) {
      dtmp = A->map[i][j];
      for (k = 0; k < klim; k++) {
	qs += dtmp[k];
      }
    }
  }

  return qs;
}

#if 0
/***=======================================================================***/
/*** TriLinInterp: trilinear interpolation of a single point to a grid.    ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   crd:    the coordinates of the points                               ***/
/***   val:    the value of the function at the points                     ***/
/***   npts:   the number of points                                        ***/
/***   A:      the three-dimensional grid, complete with origin and        ***/
/***           spacing information                                         ***/
/***   U:      the simulation cell tranformation matrix                    ***/
/***   invU:   the inverse simulation cell tranformation matrix            ***/
/***=======================================================================***/
void Interp2Grid(double* crd, double* val, int npts, fbook *A, dmat *U,
                 dmat *invU)
{
  int i, j, k;
  double dx, dy, dz, Aox, Aoy, Aoz, Agx, Agy, Agz, ldx, ldy, ldz;

  /*** Determine the grid coordinates ***/
  dx = (crd[0] - A->ox);
  dy = (crd[1] - A->oy);
  dz = (crd[2] - A->oz);
  NonOrthoReim(&dx, &dy, &dz, U, invU);
  dx *= A->invgx;
  dy *= A->invgy;
  dz *= A->invgz;
  i = dx;
  dx -= i;
  j = dy;
  dy -= j;
  k = dz;
  dz -= k;
  ldx = val*(1.0-dx);
  ldy = 1.0-dy;
  ldz = 1.0-dz;
  dx *= val;
  A->map[i][j][k] += ldx*ldy*ldz;
  A->map[i+1][j][k] += dx*ldy*ldz;
  A->map[i][j+1][k] += ldx*dy*ldz;
  A->map[i+1][j+1][k] += dx*dy*ldz;
  A->map[i][j][k+1] += ldx*ldy*dz;
  A->map[i+1][j][k+1] += dx*ldy*dz;
  A->map[i][j+1][k+1] += ldx*dy*dz;
  A->map[i+1][j+1][k+1] += dx*dy*dz;
}
#endif
