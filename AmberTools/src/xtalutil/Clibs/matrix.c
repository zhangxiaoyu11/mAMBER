#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include "vector.h"
#include "macros.h"

/***=======================================================================***/
/*** CreateImat: create an N x M integer matrix, initialized to zero.      ***/
/***=======================================================================***/
imat CreateImat(int N, int M)
{
  int i;
  imat A;

  A.row = N;
  A.col = M;
  A.map = (int**)malloc(N*sizeof(int*));
  A.data = (int*)calloc(N*M, sizeof(int));
  for (i = 0; i < N; i++) {
    A.map[i] = &A.data[M*i];
  }

  return A;
}

/***=======================================================================***/
/*** DestroyImat: destroy an integer matrix.                               ***/
/***=======================================================================***/
void DestroyImat(imat *A)
{
  free(A->data);
  free(A->map);
}

/***=======================================================================***/
/*** CreateDmat: create an N x M double-precision real matrix, initialized ***/
/***             to zero.                                                  ***/
/***=======================================================================***/
dmat CreateDmat(int N, int M)
{
  int i;
  dmat A;

  A.row = N;
  A.col = M;
  A.map = (double**)malloc(N*sizeof(double*));
  A.data = (double*)calloc(N*M, sizeof(double));
  for (i = 0; i < N; i++) {
    A.map[i] = &A.data[M*i];
  }

  return A;
}

/***=======================================================================***/
/*** DestroyDmat: destroy a double-precision real matrix.                  ***/
/***=======================================================================***/
void DestroyDmat(dmat *A)
{
  free(A->data);
  free(A->map);
}

/***=======================================================================***/
/*** CreateCmat: create an N x M character matrix, initialized to null.    ***/
/***=======================================================================***/
cmat CreateCmat(int N, int M)
{
  int i;
  cmat A;

  A.row = N;
  A.col = M;
  A.map = (char**)malloc(N*sizeof(char*));
  A.data = (char*)calloc(N*M, sizeof(char));
  for (i = 0; i < N; i++) {
    A.map[i] = &A.data[M*i];
  }

  return A;
}

/***=======================================================================***/
/*** DestroyCmat: destroy a character matrix.                              ***/
/***=======================================================================***/
void DestroyCmat(cmat *A)
{
  free(A->data);
  free(A->map);
}

/***=======================================================================***/
/*** DMatMult: multiply two double-precision real matrices A and B to form ***/
/***           another matrix C.                                           ***/
/***=======================================================================***/
void DMatMult(dmat A, dmat B, dmat C)
{
  int i, j, k;
  double tempval;
  double* tmp;
  double* tm2p;

  if (A.col != B.row) {
    printf("DMatMult >> Error: internal dimensions do not agree.\n"
           "DMatMult >> Taking minimum value.\n\n");
    if (A.col > B.row) {
      A.col = B.row;
    }
    else {
      B.row = A.col;
    }
  }

  for (i = 0; i < A.row; i++) {
    tmp = A.map[i];
    tm2p = C.map[i];
    for (j = 0; j < B.col; j++) {
      tempval = 0.0;
      for (k = 0; k < A.col; k++) {
	tempval += tmp[k]*B.map[k][j];
      }
      tm2p[j] = tempval;
    }
  }
}

/***=======================================================================***/
/*** IMatMult: multiply two integer matrices A and B to form matrix C.     ***/
/***=======================================================================***/
void IMatMult(imat A, imat B, imat C)
{
  int i, j, k;
  int tempval;
  int* tmp;
  int* tm2p;

  if (A.col != B.row) {
    printf("IMatMult >> Error: internal dimensions do not agree.\n"
           "IMatMult >> Taking minimum value.\n\n");
    if (A.col > B.row) {
      A.col = B.row;
    }
    else {
      B.row = A.col;
    }
  }

  for (i = 0; i < A.row; i++) {
    tmp = A.map[i];
    tm2p = C.map[i];
    for (j = 0; j < B.col; j++) {
      tempval = 0.0;
      for (k = 0; k < A.col; k++) {
	tempval += tmp[k]*B.map[k][j];
      }
      tm2p[j] = tempval;
    }
  }
}

/***=======================================================================***/
/*** ttInv: inverts a 3 x 3 matrix a and returns the result as inva.  The  ***/
/***        algorithm is HIGHLY INEFFICIENT (scales as N!) but for the     ***/
/***        3 x 3 case (3! = 6) we don't care.                             ***/
/***=======================================================================***/
void ttInv(dmat A, dmat invA)
{
  double deta;
  double* pq;
  double* rs;

  /*** Compute the determinant of a ***/
  deta = A.map[0][0]*(A.map[1][1]*A.map[2][2] - A.map[1][2]*A.map[2][1]) +
         A.map[0][1]*(A.map[1][2]*A.map[2][0] - A.map[1][0]*A.map[2][2]) +
         A.map[0][2]*(A.map[1][0]*A.map[2][1] - A.map[1][1]*A.map[2][0]);
  deta = 1.0/deta;

  pq = A.map[1];
  rs = A.map[2];
  ttInvRow(invA, pq, rs, 0, deta);
  pq = A.map[0];
  ttInvRow(invA, pq, rs, 1, deta);
  rs = A.map[1];
  ttInvRow(invA, pq, rs, 2, deta);
  invA.map[0][1] *= -1.0;
  invA.map[1][0] *= -1.0;
  invA.map[1][2] *= -1.0;
  invA.map[2][1] *= -1.0;
}

/***=======================================================================***/
/*** ttInvRow: accessory to ttInv above, makes the rows of the inverse.    ***/
/***=======================================================================***/
void ttInvRow(dmat invA, double* pq, double* rs, int j, double deta)
{
  invA.map[0][j] = deta*(pq[1]*rs[2]-rs[1]*pq[2]);
  invA.map[1][j] = deta*(pq[0]*rs[2]-rs[0]*pq[2]);
  invA.map[2][j] = deta*(pq[0]*rs[1]-rs[0]*pq[1]);
}

/***=======================================================================***/
/*** DMatAdd: computes the expression da*A + db*B = C, where da and db     ***/
/***          are scalars and A, B, and C are matrices.  Note that C may   ***/
/***          point to addresses of either A or B so that the result will  ***/
/***          overwrite A or B.                                            ***/
/***=======================================================================***/
void DMatAdd(dmat A, dmat B, double da, double db, dmat C)
{
  int i;

  /*** Check ***/
  if (A.row != B.row || A.col != B.col || A.row != C.row || A.col != C.col) {
    printf("DMatAdd >> Error.  Matrix dimensions do not agree.\n"
	   "DMatAdd >> A = [ %6d x %6d ]\nDMatAdd >> B = [ %6d x %6d ]\n"
	   "DMatAdd >> C = [ %6d x %6d ]\n", A.row, A.col, B.row, B.col,
	   C.row, C.col);
    exit(1);
  }

  for (i = 0; i < A.row*A.col; i++) {
    C.data[i] = da*A.data[i] + db*B.data[i];
  }
}

/***=======================================================================***/
/*** ROTATION_MATRIX: constructs a rotation matrix for rotating about a    ***/
/***                  particular vector.                                   ***/
/***=======================================================================***/
void RotationMatrix(dmat *mat, double* vec, double angle)
{
  double cT = cos(angle), cT1 = 1-cT;
  double sT = sin(angle);
  double x = vec[0], y = vec[1], z = vec[2];

  /*** Generate the rotation matrix ***/
  mat->map[0][0] = 1 + cT1 * (x*x - 1);
  mat->map[0][1] = -z*sT + cT1*x*y;
  mat->map[0][2] = y*sT + cT1*x*z;
  mat->map[1][0] = z*sT + cT1*x*y;
  mat->map[1][1] = 1 + cT1*(y*y - 1);
  mat->map[1][2] = -x*sT + cT1*y*z;
  mat->map[2][0] = -y*sT + cT1*x*z;
  mat->map[2][1] = x*sT + cT1*y*z;
  mat->map[2][2] = 1 + cT1*(z*z - 1);
}

/***=======================================================================***/
/*** MatVecMult: multiply matrix A times vector B and put it in vector C.  ***/
/***=======================================================================***/
void MatVecMult(dmat *A, double* B, double* C)
{
  int i;

  for (i = 0; i < A->row; i++) {
    C[i] = DotP(A->map[i], B, A->col);
  }
}

/***=======================================================================***/
/*** AxbQRRxc: function for solving linear least-squares problems.  This   ***/
/***           performs the first stage of the solution by taking a        ***/
/***           matrix problem Ax = b, where A { m by n, m >= n, and        ***/
/***           decomposes A by the "modern classical" QR algorithm to      ***/
/***           recast the problem as Rx = c, where R is the the upper-     ***/
/***           triangular component of A and the vector c is implicitly    ***/
/***           computed as (Q*b).  This algorithm may be found in [REF]:   ***/
/***                                                                       ***/
/***           Trefethen, Lloyd N. and Bau, David III. "Numerical Linear   ***/
/***           Algebra." pp.73.  Society for Industrial and Applied        ***/
/***           Mathematics, Philadelphia.  1997.                           ***/
/***=======================================================================***/
void AxbQRRxc(dmat A, double* b, int update_user)
{
  int i, j, k, iks;
  double tnm_v, tempval, sign_v;
  double* v;
  double* vprime;
  double* tmp;

  if (A.row >= A.col) {
    v = (double*)malloc(A.row*sizeof(double));
    vprime = (double*)malloc(A.row*sizeof(double));
    for (k = 0; k < A.col; k++) {

      /*** Update the user ***/
      if (update_user == 1) {
        fprintf(stderr, "\rAxbQRRxc >> Computing for column %5d", k);
        fflush(stderr);
      }

      /*** Compute the kth column of Q* ***/
      for (i = 0; i < A.row-k; i++) {
        v[i] = A.map[i+k][k];
      }
      sign_v = SIGN(v[0]);
      tnm_v = DMag(v, A.row-k);
      v[0] += sign_v*tnm_v;
      tnm_v = 1.0/DMag(v, A.row-k);
      for (i = 0; i < A.row-k; i++) {
        v[i] = v[i]*tnm_v;
      }

      /*** Update A as R evolves ***/
      for (i = 0; i < A.col-k; i++) {
        tempval = 0.0;
        iks = i+k;
        for (j = 0; j < A.row-k; j++) {
          tempval += v[j]*A.map[j+k][iks];
        }
        vprime[i] = tempval;
      }
      for (i = 0; i < A.row-k; i++) {
        tmp = A.map[i+k];
        tempval = 2.0*v[i];
        for (j = 0; j < A.col-k; j++) {
          tmp[j+k] -= tempval*vprime[j];
        }
      }

      /*** Update b as Q* evolves ***/
      for (i = 0; i < A.row-k; i++) {
        vprime[i] = b[i+k];
      }
      tempval = 2.0*DotP(v, vprime, A.row-k);
      for (i = 0; i < A.row-k; i++) {
        b[i+k] -= tempval*v[i];
      }
    }

    /*** Free Allocated Memory ***/
    free(v);
    free(vprime);

    if (update_user == 1) {
      printf("\n\n");
    }
  }
  else {
    printf("AxbQRRxc >> Error: Matrix A is rank deficient.\n");
  }
}

/***=======================================================================***/
/*** BackSub: solve the equation Rx = b, where R is an upper triangular    ***/
/***          matrix of dimension n, and b is a vector of dimension n.     ***/
/***          Results are returned in the vector b.                        ***/
/***=======================================================================***/
void BackSub(dmat R, double* b)
{
  int i, j;
  double multval, temp_b, pivot;

  for (i = R.col-1; i > 0; i--) {
    pivot = 1.0/R.map[i][i];
    temp_b = b[i];
    for (j = i-1; j >= 0; j--) {
      multval = R.map[j][i]*pivot;
      b[j] -= multval*temp_b;
    }
    b[i] *= pivot;
  }
  b[0] /= R.map[0][0];
}

/***=======================================================================***/
/*** TRED2: reduces a symmetric matrix to tridiagonal form.  It is a       ***/
/***        modification more suited for C programs from [REF]:            ***/
/***                                                                       ***/
/***        William H. Press, Saul A. Teukolsky, William T. Vetterling,    ***/
/***        and Brian P. Flannery.  Numerical Recipes in C, Second         ***/
/***        Edition.  Cambridge University Press, 1992.                    ***/
/***=======================================================================***/
void TRed2(dmat A, double* d, double* e)
{
  int i, j, k, l, n;
  double scale, hh, h, g, f, invscale, invh;
  double* tmp;
  double* tm2p;

  n = A.row;
  for (i = n-1; i >= 1; i--) {
    l = i - 1;
    h = 0.0;
    scale = 0.0;
    tmp = A.map[i];
    if (l > 0) {
      for (k = 0; k <= l; k++) {
        scale += fabs(tmp[k]);
      }
      if (scale == 0.0) {
        e[i] = tmp[l];
      }
      else {
        invscale = 1.0/scale;
        for (k = 0; k <= l; k++) {
          tmp[k] *= invscale;
          h += tmp[k]*tmp[k];
        }
        f = tmp[l];
        g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
        e[i] = scale*g;
        h -= f*g;
        tmp[l] = f - g;
        f = 0.0;
        invh = 1.0/h;
        for (j = 0; j <= l; j++) {
          tm2p = A.map[j];
          tm2p[i] = tmp[j]*invh;
          g = 0.0;
          for (k = 0; k <= j; k++) {
            g += tm2p[k]*tmp[k];
          }
          for (k = j+1; k <=l; k++) {
            g += A.map[k][j]*tmp[k];
          }
          e[j] = g*invh;
          f += e[j]*tmp[j];
        }
        hh = 0.5*f*invh;
        for (j = 0; j <= l; j++) {
          f = tmp[j];
          e[j] = g = e[j] - hh*f;
          tm2p = A.map[j];
          for (k = 0; k <= j; k++) {
            tm2p[k] -= (f*e[k] + g*tmp[k]);
          }
        }
      }
    }
    else {
      e[i] = tmp[l];
    }
    d[i] = h;
  }
  d[0] = 0.0;
  e[0] = 0.0;
  for (i = 0; i < n; i++) {
    tmp = A.map[i];
    l = i - 1;
    if (d[i]) {
      for (j = 0; j <= l; j++) {
        g = 0.0;
        for (k = 0; k <= l; k++) {
          g += tmp[k]*A.map[k][j];
        }
        for (k = 0; k <= l; k++) {
          A.map[k][j] -= g*A.map[k][i];
        }
      }
    }
    d[i] = tmp[i];
    tmp[i] = 1.0;
    for (j = 0; j <= l; j++) {
      A.map[j][i] = tmp[j] = 0.0;
    }
  }
}

/***=======================================================================***/
/*** TQLI: this function is also modified from [REF]:                      ***/
/***                                                                       ***/
/***        William H. Press, Saul A. Teukolsky, William T. Vetterling,    ***/
/***        and Brian P. Flannery.  Numerical Recipes in C, Second         ***/
/***        Edition.  Cambridge University Press, 1992.                    ***/
/***                                                                       ***/
/***        The input is the super- and sub- diagonal vectors d and e.     ***/
/***        The eigenvalues are returned in d and the eigenvectors are     ***/
/***        returned in z.                                                 ***/
/***=======================================================================***/
void tqli(double* d, double* e, dmat z)
{
  int m, l, n, iter, i, k;
  double s, r, p, g, f, dd, c, b, invr;
  double* tmp;

  n = z.row;
  for (i = 1; i < n; i++) {
    e[i-1] = e[i];
  }
  e[n-1] = 0.0;
  for (l = 0; l < n; l++) {
    iter = 0;

    /*** Added to initialize m ***/
    m = l - 1;

    while (m != l) {
      for (m = l; m < n - 1; m++) {
        dd = fabs(d[m]) + fabs(d[m+1]);
        if (fabs(e[m]+dd) == dd) break;
      }
      if (m != l) {
        if (iter++ == 30) {
          printf("tqli >> Error: Too many iterations in tqli\n");
	  exit(1);
        }
        g = (d[l+1]-d[l])/(2.0*e[l]);
        r = pythag(g,1.0);
        g = d[m]-d[l]+e[l]/(g+SIGN2(r,g));
        c = 1.0;
        s = 1.0;
        p = 0.0;
        for (i = m-1; i >= l; i--) {
          f = s*e[i];
          b = c*e[i];
          e[i+1] = (r=pythag(f,g));
          if (r == 0.0) {
            d[i+1] -= p;
            e[m] = 0.0;
            break;
          }
          invr = 1.0/r;
          s = f*invr;
          c = g*invr;
          g = d[i+1]-p;
          r = (d[i] - g)*s + 2.0*c*b;
          d[i+1] = g + (p = s*r);
          g = c*r - b;
          for (k = 0; k < n; k++) {
            tmp = z.map[k];
            f = tmp[i+1];
            tmp[i+1] = s*tmp[i] + c*f;
            tmp[i] = c*tmp[i] - s*f;
          }
        }
        if (r == 0.0 && i >= l) continue;
        d[l] -= p;
        e[l] = g;
        e[m] = 0.0;
      }
    }
  }
}
