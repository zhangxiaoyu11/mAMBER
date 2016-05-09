#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_dge_sum_mv_s_d(enum blas_order_type order, int m, int n,
			 double alpha, const float *a, int lda,
			 const double *x, int incx,
			 double beta, const float *b, int ldb,
			 double *y, int incy)

/*
 * Purpose
 * =======
 *
 * Computes y = alpha * A * x + beta * B * y, 
 *     where A, B are general matricies.
 *
 * Arguments
 * =========
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * m            (input) int
 *              Row Dimension of A, B, length of output vector y
 *
 * n            (input) int
 *              Column Dimension of A, B and the length of vector x
 *
 * alpha        (input) double
 *              
 * A            (input) const float*
 *
 * lda          (input) int 
 *              Leading dimension of A
 *
 * x            (input) const double*
 * 
 * incx         (input) int
 *              The stride for vector x.
 *
 * beta         (input) double
 *
 * b            (input) const float*
 *
 * ldb          (input) int 
 *              Leading dimension of B
 *
 * y            (input/output) double*
 *
 * incy         (input) int
 *              The stride for vector y.
 * 
 */
{
  /* Routine name */
  static const char routine_name[] = "BLAS_dge_sum_mv_s_d";
  int i, j;
  int xi, yi;
  int x_starti, y_starti, incxi, incyi;
  int lda_min;
  int ai;
  int incai;
  int aij;
  int incaij;
  int bi;
  int incbi;
  int bij;
  int incbij;

  const float *a_i = a;
  const float *b_i = b;
  const double *x_i = x;
  double *y_i = y;
  double alpha_i = alpha;
  double beta_i = beta;
  float a_elem;
  float b_elem;
  double x_elem;
  double prod;
  double sumA;
  double sumB;
  double tmp1;
  double tmp2;



  /* m is number of rows */
  /* n is number of columns */

  if (m == 0 || n == 0)
    return;


  /* all error calls */
  if (order == blas_rowmajor) {
    lda_min = n;
    incai = lda;		/* row stride */
    incbi = ldb;
    incbij = incaij = 1;	/* column stride */
  } else if (order == blas_colmajor) {
    lda_min = m;
    incai = incbi = 1;		/*row stride */
    incaij = lda;		/* column stride */
    incbij = ldb;
  } else {
    /* error, order not blas_colmajor not blas_rowmajor */
    BLAS_error(routine_name, -1, order, 0);
    return;
  }

  if (m < 0)
    BLAS_error(routine_name, -2, m, 0);
  else if (n < 0)
    BLAS_error(routine_name, -3, n, 0);
  if (lda < lda_min)
    BLAS_error(routine_name, -6, lda, 0);
  else if (ldb < lda_min)
    BLAS_error(routine_name, -11, ldb, 0);
  else if (incx == 0)
    BLAS_error(routine_name, -8, incx, 0);
  else if (incy == 0)
    BLAS_error(routine_name, -13, incy, 0);

  incxi = incx;
  incyi = incy;







  if (incxi > 0)
    x_starti = 0;
  else
    x_starti = (1 - n) * incxi;

  if (incyi > 0)
    y_starti = 0;
  else
    y_starti = (1 - m) * incyi;



  if (alpha_i == 0.0) {
    if (beta_i == 0.0) {
      /* alpha, beta are 0.0 */
      for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	y_i[yi] = 0.0;
      }
    } else if (beta_i == 1.0) {
      /* alpha is 0.0, beta is 1.0 */


      bi = 0;
      for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {

	sumB = 0.0;
	bij = bi;
	for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	  x_elem = x_i[xi];

	  b_elem = b_i[bij];
	  prod = b_elem * x_elem;
	  sumB = sumB + prod;
	  bij += incbij;
	}
	/* now put the result into y_i */
	y_i[yi] = sumB;

	bi += incbi;
      }
    } else {
      /* alpha is 0.0, beta not 1.0 nor 0.0 */


      bi = 0;
      for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {

	sumB = 0.0;
	bij = bi;
	for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	  x_elem = x_i[xi];

	  b_elem = b_i[bij];
	  prod = b_elem * x_elem;
	  sumB = sumB + prod;
	  bij += incbij;
	}
	/* now put the result into y_i */
	tmp1 = sumB * beta_i;
	y_i[yi] = tmp1;

	bi += incbi;
      }
    }
  } else if (alpha_i == 1.0) {
    if (beta_i == 0.0) {
      /* alpha is 1.0, beta is 0.0 */

      ai = 0;

      for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	sumA = 0.0;
	aij = ai;

	for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	  x_elem = x_i[xi];
	  a_elem = a_i[aij];
	  prod = a_elem * x_elem;
	  sumA = sumA + prod;
	  aij += incaij;

	}
	/* now put the result into y_i */
	y_i[yi] = sumA;
	ai += incai;

      }
    } else if (beta_i == 1.0) {
      /* alpha is 1.0, beta is 1.0 */

      ai = 0;
      bi = 0;
      for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	sumA = 0.0;
	aij = ai;
	sumB = 0.0;
	bij = bi;
	for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	  x_elem = x_i[xi];
	  a_elem = a_i[aij];
	  prod = a_elem * x_elem;
	  sumA = sumA + prod;
	  aij += incaij;
	  b_elem = b_i[bij];
	  prod = b_elem * x_elem;
	  sumB = sumB + prod;
	  bij += incbij;
	}
	/* now put the result into y_i */
	tmp1 = sumA;
	tmp2 = sumB;
	tmp1 = tmp1 + tmp2;
	y_i[yi] = tmp1;
	ai += incai;
	bi += incbi;
      }
    } else {
      /* alpha is 1.0, beta is other */

      ai = 0;
      bi = 0;
      for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	sumA = 0.0;
	aij = ai;
	sumB = 0.0;
	bij = bi;
	for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	  x_elem = x_i[xi];
	  a_elem = a_i[aij];
	  prod = a_elem * x_elem;
	  sumA = sumA + prod;
	  aij += incaij;
	  b_elem = b_i[bij];
	  prod = b_elem * x_elem;
	  sumB = sumB + prod;
	  bij += incbij;
	}
	/* now put the result into y_i */
	tmp1 = sumA;
	tmp2 = sumB * beta_i;
	tmp1 = tmp1 + tmp2;
	y_i[yi] = tmp1;
	ai += incai;
	bi += incbi;
      }
    }
  } else {
    if (beta_i == 0.0) {
      /* alpha is other, beta is 0.0 */

      ai = 0;

      for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	sumA = 0.0;
	aij = ai;

	for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	  x_elem = x_i[xi];
	  a_elem = a_i[aij];
	  prod = a_elem * x_elem;
	  sumA = sumA + prod;
	  aij += incaij;

	}
	/* now put the result into y_i */
	tmp1 = sumA * alpha_i;
	y_i[yi] = tmp1;
	ai += incai;

      }
    } else if (beta_i == 1.0) {
      /* alpha is other, beta is 1.0 */

      ai = 0;
      bi = 0;
      for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	sumA = 0.0;
	aij = ai;
	sumB = 0.0;
	bij = bi;
	for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	  x_elem = x_i[xi];
	  a_elem = a_i[aij];
	  prod = a_elem * x_elem;
	  sumA = sumA + prod;
	  aij += incaij;
	  b_elem = b_i[bij];
	  prod = b_elem * x_elem;
	  sumB = sumB + prod;
	  bij += incbij;
	}
	/* now put the result into y_i */
	tmp1 = sumA * alpha_i;
	tmp2 = sumB;
	tmp1 = tmp1 + tmp2;
	y_i[yi] = tmp1;
	ai += incai;
	bi += incbi;
      }
    } else {
      /* most general form, alpha, beta are other */

      ai = 0;
      bi = 0;
      for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	sumA = 0.0;
	aij = ai;
	sumB = 0.0;
	bij = bi;
	for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	  x_elem = x_i[xi];
	  a_elem = a_i[aij];
	  prod = a_elem * x_elem;
	  sumA = sumA + prod;
	  aij += incaij;
	  b_elem = b_i[bij];
	  prod = b_elem * x_elem;
	  sumB = sumB + prod;
	  bij += incbij;
	}
	/* now put the result into y_i */
	tmp1 = sumA * alpha_i;
	tmp2 = sumB * beta_i;
	tmp1 = tmp1 + tmp2;
	y_i[yi] = tmp1;
	ai += incai;
	bi += incbi;
      }
    }
  }


}				/* end BLAS_dge_sum_mv_s_d */
