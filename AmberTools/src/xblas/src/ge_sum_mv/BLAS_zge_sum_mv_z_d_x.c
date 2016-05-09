#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_zge_sum_mv_z_d_x(enum blas_order_type order, int m, int n,
			   const void *alpha, const void *a, int lda,
			   const double *x, int incx,
			   const void *beta, const void *b, int ldb,
			   void *y, int incy, enum blas_prec_type prec)

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
 * alpha        (input) const void*
 *              
 * A            (input) const void*
 *
 * lda          (input) int 
 *              Leading dimension of A
 *
 * x            (input) const double*
 * 
 * incx         (input) int
 *              The stride for vector x.
 *
 * beta         (input) const void*
 *
 * b            (input) const void*
 *
 * ldb          (input) int 
 *              Leading dimension of B
 *
 * y            (input/output) void*
 *
 * incy         (input) int
 *              The stride for vector y.
 * 
 * prec   (input) enum blas_prec_type
 *        Specifies the internal precision to be used.
 *        = blas_prec_single: single precision.
 *        = blas_prec_double: double precision.
 *        = blas_prec_extra : anything at least 1.5 times as accurate
 *                            than double, and wider than 80-bits.
 *                            We use double-double in our implementation.
 *
 */
{
  /* Routine name */
  static const char routine_name[] = "BLAS_zge_sum_mv_z_d";
  switch (prec) {
  case blas_prec_single:
  case blas_prec_indigenous:
  case blas_prec_double:
    {
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

      const double *a_i = (double *) a;
      const double *b_i = (double *) b;
      const double *x_i = x;
      double *y_i = (double *) y;
      double *alpha_i = (double *) alpha;
      double *beta_i = (double *) beta;
      double a_elem[2];
      double b_elem[2];
      double x_elem;
      double prod[2];
      double sumA[2];
      double sumB[2];
      double tmp1[2];
      double tmp2[2];



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
	incai = incbi = 1;	/*row stride */
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

      incyi *= 2;
      incai *= 2;
      incaij *= 2;
      incbi *= 2;
      incbij *= 2;

      if (incxi > 0)
	x_starti = 0;
      else
	x_starti = (1 - n) * incxi;

      if (incyi > 0)
	y_starti = 0;
      else
	y_starti = (1 - m) * incyi;



      if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
	if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	  /* alpha, beta are 0.0 */
	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	    y_i[yi] = 0.0;
	    y_i[yi + 1] = 0.0;
	  }
	} else if ((beta_i[0] == 1.0 && beta_i[1] == 0.0)) {
	  /* alpha is 0.0, beta is 1.0 */


	  bi = 0;
	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {

	    sumB[0] = sumB[1] = 0.0;
	    bij = bi;
	    for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	      x_elem = x_i[xi];

	      b_elem[0] = b_i[bij];
	      b_elem[1] = b_i[bij + 1];
	      {
		prod[0] = b_elem[0] * x_elem;
		prod[1] = b_elem[1] * x_elem;
	      }
	      sumB[0] = sumB[0] + prod[0];
	      sumB[1] = sumB[1] + prod[1];
	      bij += incbij;
	    }
	    /* now put the result into y_i */
	    y_i[yi] = sumB[0];
	    y_i[yi + 1] = sumB[1];

	    bi += incbi;
	  }
	} else {
	  /* alpha is 0.0, beta not 1.0 nor 0.0 */


	  bi = 0;
	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {

	    sumB[0] = sumB[1] = 0.0;
	    bij = bi;
	    for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	      x_elem = x_i[xi];

	      b_elem[0] = b_i[bij];
	      b_elem[1] = b_i[bij + 1];
	      {
		prod[0] = b_elem[0] * x_elem;
		prod[1] = b_elem[1] * x_elem;
	      }
	      sumB[0] = sumB[0] + prod[0];
	      sumB[1] = sumB[1] + prod[1];
	      bij += incbij;
	    }
	    /* now put the result into y_i */
	    {
	      tmp1[0] =
		(double) sumB[0] * beta_i[0] - (double) sumB[1] * beta_i[1];
	      tmp1[1] =
		(double) sumB[0] * beta_i[1] + (double) sumB[1] * beta_i[0];
	    }
	    y_i[yi] = tmp1[0];
	    y_i[yi + 1] = tmp1[1];

	    bi += incbi;
	  }
	}
      } else if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {
	if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	  /* alpha is 1.0, beta is 0.0 */

	  ai = 0;

	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	    sumA[0] = sumA[1] = 0.0;
	    aij = ai;

	    for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	      x_elem = x_i[xi];
	      a_elem[0] = a_i[aij];
	      a_elem[1] = a_i[aij + 1];
	      {
		prod[0] = a_elem[0] * x_elem;
		prod[1] = a_elem[1] * x_elem;
	      }
	      sumA[0] = sumA[0] + prod[0];
	      sumA[1] = sumA[1] + prod[1];
	      aij += incaij;

	    }
	    /* now put the result into y_i */
	    y_i[yi] = sumA[0];
	    y_i[yi + 1] = sumA[1];
	    ai += incai;

	  }
	} else if ((beta_i[0] == 1.0 && beta_i[1] == 0.0)) {
	  /* alpha is 1.0, beta is 1.0 */

	  ai = 0;
	  bi = 0;
	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	    sumA[0] = sumA[1] = 0.0;
	    aij = ai;
	    sumB[0] = sumB[1] = 0.0;
	    bij = bi;
	    for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	      x_elem = x_i[xi];
	      a_elem[0] = a_i[aij];
	      a_elem[1] = a_i[aij + 1];
	      {
		prod[0] = a_elem[0] * x_elem;
		prod[1] = a_elem[1] * x_elem;
	      }
	      sumA[0] = sumA[0] + prod[0];
	      sumA[1] = sumA[1] + prod[1];
	      aij += incaij;
	      b_elem[0] = b_i[bij];
	      b_elem[1] = b_i[bij + 1];
	      {
		prod[0] = b_elem[0] * x_elem;
		prod[1] = b_elem[1] * x_elem;
	      }
	      sumB[0] = sumB[0] + prod[0];
	      sumB[1] = sumB[1] + prod[1];
	      bij += incbij;
	    }
	    /* now put the result into y_i */
	    tmp1[0] = sumA[0];
	    tmp1[1] = sumA[1];
	    tmp2[0] = sumB[0];
	    tmp2[1] = sumB[1];
	    tmp1[0] = tmp1[0] + tmp2[0];
	    tmp1[1] = tmp1[1] + tmp2[1];
	    y_i[yi] = tmp1[0];
	    y_i[yi + 1] = tmp1[1];
	    ai += incai;
	    bi += incbi;
	  }
	} else {
	  /* alpha is 1.0, beta is other */

	  ai = 0;
	  bi = 0;
	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	    sumA[0] = sumA[1] = 0.0;
	    aij = ai;
	    sumB[0] = sumB[1] = 0.0;
	    bij = bi;
	    for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	      x_elem = x_i[xi];
	      a_elem[0] = a_i[aij];
	      a_elem[1] = a_i[aij + 1];
	      {
		prod[0] = a_elem[0] * x_elem;
		prod[1] = a_elem[1] * x_elem;
	      }
	      sumA[0] = sumA[0] + prod[0];
	      sumA[1] = sumA[1] + prod[1];
	      aij += incaij;
	      b_elem[0] = b_i[bij];
	      b_elem[1] = b_i[bij + 1];
	      {
		prod[0] = b_elem[0] * x_elem;
		prod[1] = b_elem[1] * x_elem;
	      }
	      sumB[0] = sumB[0] + prod[0];
	      sumB[1] = sumB[1] + prod[1];
	      bij += incbij;
	    }
	    /* now put the result into y_i */
	    tmp1[0] = sumA[0];
	    tmp1[1] = sumA[1];
	    {
	      tmp2[0] =
		(double) sumB[0] * beta_i[0] - (double) sumB[1] * beta_i[1];
	      tmp2[1] =
		(double) sumB[0] * beta_i[1] + (double) sumB[1] * beta_i[0];
	    }
	    tmp1[0] = tmp1[0] + tmp2[0];
	    tmp1[1] = tmp1[1] + tmp2[1];
	    y_i[yi] = tmp1[0];
	    y_i[yi + 1] = tmp1[1];
	    ai += incai;
	    bi += incbi;
	  }
	}
      } else {
	if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	  /* alpha is other, beta is 0.0 */

	  ai = 0;

	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	    sumA[0] = sumA[1] = 0.0;
	    aij = ai;

	    for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	      x_elem = x_i[xi];
	      a_elem[0] = a_i[aij];
	      a_elem[1] = a_i[aij + 1];
	      {
		prod[0] = a_elem[0] * x_elem;
		prod[1] = a_elem[1] * x_elem;
	      }
	      sumA[0] = sumA[0] + prod[0];
	      sumA[1] = sumA[1] + prod[1];
	      aij += incaij;

	    }
	    /* now put the result into y_i */
	    {
	      tmp1[0] =
		(double) sumA[0] * alpha_i[0] - (double) sumA[1] * alpha_i[1];
	      tmp1[1] =
		(double) sumA[0] * alpha_i[1] + (double) sumA[1] * alpha_i[0];
	    }
	    y_i[yi] = tmp1[0];
	    y_i[yi + 1] = tmp1[1];
	    ai += incai;

	  }
	} else if ((beta_i[0] == 1.0 && beta_i[1] == 0.0)) {
	  /* alpha is other, beta is 1.0 */

	  ai = 0;
	  bi = 0;
	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	    sumA[0] = sumA[1] = 0.0;
	    aij = ai;
	    sumB[0] = sumB[1] = 0.0;
	    bij = bi;
	    for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	      x_elem = x_i[xi];
	      a_elem[0] = a_i[aij];
	      a_elem[1] = a_i[aij + 1];
	      {
		prod[0] = a_elem[0] * x_elem;
		prod[1] = a_elem[1] * x_elem;
	      }
	      sumA[0] = sumA[0] + prod[0];
	      sumA[1] = sumA[1] + prod[1];
	      aij += incaij;
	      b_elem[0] = b_i[bij];
	      b_elem[1] = b_i[bij + 1];
	      {
		prod[0] = b_elem[0] * x_elem;
		prod[1] = b_elem[1] * x_elem;
	      }
	      sumB[0] = sumB[0] + prod[0];
	      sumB[1] = sumB[1] + prod[1];
	      bij += incbij;
	    }
	    /* now put the result into y_i */
	    {
	      tmp1[0] =
		(double) sumA[0] * alpha_i[0] - (double) sumA[1] * alpha_i[1];
	      tmp1[1] =
		(double) sumA[0] * alpha_i[1] + (double) sumA[1] * alpha_i[0];
	    }
	    tmp2[0] = sumB[0];
	    tmp2[1] = sumB[1];
	    tmp1[0] = tmp1[0] + tmp2[0];
	    tmp1[1] = tmp1[1] + tmp2[1];
	    y_i[yi] = tmp1[0];
	    y_i[yi + 1] = tmp1[1];
	    ai += incai;
	    bi += incbi;
	  }
	} else {
	  /* most general form, alpha, beta are other */

	  ai = 0;
	  bi = 0;
	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	    sumA[0] = sumA[1] = 0.0;
	    aij = ai;
	    sumB[0] = sumB[1] = 0.0;
	    bij = bi;
	    for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	      x_elem = x_i[xi];
	      a_elem[0] = a_i[aij];
	      a_elem[1] = a_i[aij + 1];
	      {
		prod[0] = a_elem[0] * x_elem;
		prod[1] = a_elem[1] * x_elem;
	      }
	      sumA[0] = sumA[0] + prod[0];
	      sumA[1] = sumA[1] + prod[1];
	      aij += incaij;
	      b_elem[0] = b_i[bij];
	      b_elem[1] = b_i[bij + 1];
	      {
		prod[0] = b_elem[0] * x_elem;
		prod[1] = b_elem[1] * x_elem;
	      }
	      sumB[0] = sumB[0] + prod[0];
	      sumB[1] = sumB[1] + prod[1];
	      bij += incbij;
	    }
	    /* now put the result into y_i */
	    {
	      tmp1[0] =
		(double) sumA[0] * alpha_i[0] - (double) sumA[1] * alpha_i[1];
	      tmp1[1] =
		(double) sumA[0] * alpha_i[1] + (double) sumA[1] * alpha_i[0];
	    }
	    {
	      tmp2[0] =
		(double) sumB[0] * beta_i[0] - (double) sumB[1] * beta_i[1];
	      tmp2[1] =
		(double) sumB[0] * beta_i[1] + (double) sumB[1] * beta_i[0];
	    }
	    tmp1[0] = tmp1[0] + tmp2[0];
	    tmp1[1] = tmp1[1] + tmp2[1];
	    y_i[yi] = tmp1[0];
	    y_i[yi + 1] = tmp1[1];
	    ai += incai;
	    bi += incbi;
	  }
	}
      }

    }
    break;

  case blas_prec_extra:
    {
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

      const double *a_i = (double *) a;
      const double *b_i = (double *) b;
      const double *x_i = x;
      double *y_i = (double *) y;
      double *alpha_i = (double *) alpha;
      double *beta_i = (double *) beta;
      double a_elem[2];
      double b_elem[2];
      double x_elem;
      double head_prod[2], tail_prod[2];
      double head_sumA[2], tail_sumA[2];
      double head_sumB[2], tail_sumB[2];
      double head_tmp1[2], tail_tmp1[2];
      double head_tmp2[2], tail_tmp2[2];

      FPU_FIX_DECL;

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
	incai = incbi = 1;	/*row stride */
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

      incyi *= 2;
      incai *= 2;
      incaij *= 2;
      incbi *= 2;
      incbij *= 2;

      if (incxi > 0)
	x_starti = 0;
      else
	x_starti = (1 - n) * incxi;

      if (incyi > 0)
	y_starti = 0;
      else
	y_starti = (1 - m) * incyi;

      FPU_FIX_START;

      if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
	if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	  /* alpha, beta are 0.0 */
	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	    y_i[yi] = 0.0;
	    y_i[yi + 1] = 0.0;
	  }
	} else if ((beta_i[0] == 1.0 && beta_i[1] == 0.0)) {
	  /* alpha is 0.0, beta is 1.0 */


	  bi = 0;
	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {

	    head_sumB[0] = head_sumB[1] = tail_sumB[0] = tail_sumB[1] = 0.0;
	    bij = bi;
	    for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	      x_elem = x_i[xi];

	      b_elem[0] = b_i[bij];
	      b_elem[1] = b_i[bij + 1];
	      {
		/* Compute complex-extra = complex-double * real. */
		double head_t, tail_t;
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = x_elem * split;
		  a1 = con - x_elem;
		  a1 = con - a1;
		  a2 = x_elem - a1;
		  con = b_elem[0] * split;
		  b1 = con - b_elem[0];
		  b1 = con - b1;
		  b2 = b_elem[0] - b1;

		  head_t = x_elem * b_elem[0];
		  tail_t =
		    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		head_prod[0] = head_t;
		tail_prod[0] = tail_t;
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = x_elem * split;
		  a1 = con - x_elem;
		  a1 = con - a1;
		  a2 = x_elem - a1;
		  con = b_elem[1] * split;
		  b1 = con - b_elem[1];
		  b1 = con - b1;
		  b2 = b_elem[1] - b1;

		  head_t = x_elem * b_elem[1];
		  tail_t =
		    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		head_prod[1] = head_t;
		tail_prod[1] = tail_t;
	      }
	      {
		double head_t, tail_t;
		double head_a, tail_a;
		double head_b, tail_b;
		/* Real part */
		head_a = head_sumB[0];
		tail_a = tail_sumB[0];
		head_b = head_prod[0];
		tail_b = tail_prod[0];
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_a + head_b;
		  bv = s1 - head_a;
		  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_a + tail_b;
		  bv = t1 - tail_a;
		  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_sumB[0] = head_t;
		tail_sumB[0] = tail_t;
		/* Imaginary part */
		head_a = head_sumB[1];
		tail_a = tail_sumB[1];
		head_b = head_prod[1];
		tail_b = tail_prod[1];
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_a + head_b;
		  bv = s1 - head_a;
		  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_a + tail_b;
		  bv = t1 - tail_a;
		  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_sumB[1] = head_t;
		tail_sumB[1] = tail_t;
	      }
	      bij += incbij;
	    }
	    /* now put the result into y_i */
	    y_i[yi] = head_sumB[0];
	    y_i[yi + 1] = head_sumB[1];

	    bi += incbi;
	  }
	} else {
	  /* alpha is 0.0, beta not 1.0 nor 0.0 */


	  bi = 0;
	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {

	    head_sumB[0] = head_sumB[1] = tail_sumB[0] = tail_sumB[1] = 0.0;
	    bij = bi;
	    for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	      x_elem = x_i[xi];

	      b_elem[0] = b_i[bij];
	      b_elem[1] = b_i[bij + 1];
	      {
		/* Compute complex-extra = complex-double * real. */
		double head_t, tail_t;
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = x_elem * split;
		  a1 = con - x_elem;
		  a1 = con - a1;
		  a2 = x_elem - a1;
		  con = b_elem[0] * split;
		  b1 = con - b_elem[0];
		  b1 = con - b1;
		  b2 = b_elem[0] - b1;

		  head_t = x_elem * b_elem[0];
		  tail_t =
		    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		head_prod[0] = head_t;
		tail_prod[0] = tail_t;
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = x_elem * split;
		  a1 = con - x_elem;
		  a1 = con - a1;
		  a2 = x_elem - a1;
		  con = b_elem[1] * split;
		  b1 = con - b_elem[1];
		  b1 = con - b1;
		  b2 = b_elem[1] - b1;

		  head_t = x_elem * b_elem[1];
		  tail_t =
		    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		head_prod[1] = head_t;
		tail_prod[1] = tail_t;
	      }
	      {
		double head_t, tail_t;
		double head_a, tail_a;
		double head_b, tail_b;
		/* Real part */
		head_a = head_sumB[0];
		tail_a = tail_sumB[0];
		head_b = head_prod[0];
		tail_b = tail_prod[0];
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_a + head_b;
		  bv = s1 - head_a;
		  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_a + tail_b;
		  bv = t1 - tail_a;
		  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_sumB[0] = head_t;
		tail_sumB[0] = tail_t;
		/* Imaginary part */
		head_a = head_sumB[1];
		tail_a = tail_sumB[1];
		head_b = head_prod[1];
		tail_b = tail_prod[1];
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_a + head_b;
		  bv = s1 - head_a;
		  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_a + tail_b;
		  bv = t1 - tail_a;
		  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_sumB[1] = head_t;
		tail_sumB[1] = tail_t;
	      }
	      bij += incbij;
	    }
	    /* now put the result into y_i */
	    {
	      /* Compute complex-extra = complex-extra * complex-double. */
	      double head_a0, tail_a0;
	      double head_a1, tail_a1;
	      double head_t1, tail_t1;
	      double head_t2, tail_t2;
	      head_a0 = head_sumB[0];
	      tail_a0 = tail_sumB[0];
	      head_a1 = head_sumB[1];
	      tail_a1 = tail_sumB[1];
	      /* real part */
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a0 * split;
		a11 = con - head_a0;
		a11 = con - a11;
		a21 = head_a0 - a11;
		con = beta_i[0] * split;
		b1 = con - beta_i[0];
		b1 = con - b1;
		b2 = beta_i[0] - b1;

		c11 = head_a0 * beta_i[0];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a0 * beta_i[0];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a1 * split;
		a11 = con - head_a1;
		a11 = con - a11;
		a21 = head_a1 - a11;
		con = beta_i[1] * split;
		b1 = con - beta_i[1];
		b1 = con - b1;
		b2 = beta_i[1] - b1;

		c11 = head_a1 * beta_i[1];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a1 * beta_i[1];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t2 = t1 + t2;
		tail_t2 = t2 - (head_t2 - t1);
	      }
	      head_t2 = -head_t2;
	      tail_t2 = -tail_t2;
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_t1 + head_t2;
		bv = s1 - head_t1;
		s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_t1 + tail_t2;
		bv = t1 - tail_t1;
		t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      head_tmp1[0] = head_t1;
	      tail_tmp1[0] = tail_t1;
	      /* imaginary part */
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a1 * split;
		a11 = con - head_a1;
		a11 = con - a11;
		a21 = head_a1 - a11;
		con = beta_i[0] * split;
		b1 = con - beta_i[0];
		b1 = con - b1;
		b2 = beta_i[0] - b1;

		c11 = head_a1 * beta_i[0];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a1 * beta_i[0];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a0 * split;
		a11 = con - head_a0;
		a11 = con - a11;
		a21 = head_a0 - a11;
		con = beta_i[1] * split;
		b1 = con - beta_i[1];
		b1 = con - b1;
		b2 = beta_i[1] - b1;

		c11 = head_a0 * beta_i[1];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a0 * beta_i[1];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t2 = t1 + t2;
		tail_t2 = t2 - (head_t2 - t1);
	      }
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_t1 + head_t2;
		bv = s1 - head_t1;
		s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_t1 + tail_t2;
		bv = t1 - tail_t1;
		t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      head_tmp1[1] = head_t1;
	      tail_tmp1[1] = tail_t1;
	    }

	    y_i[yi] = head_tmp1[0];
	    y_i[yi + 1] = head_tmp1[1];

	    bi += incbi;
	  }
	}
      } else if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {
	if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	  /* alpha is 1.0, beta is 0.0 */

	  ai = 0;

	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	    head_sumA[0] = head_sumA[1] = tail_sumA[0] = tail_sumA[1] = 0.0;
	    aij = ai;

	    for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	      x_elem = x_i[xi];
	      a_elem[0] = a_i[aij];
	      a_elem[1] = a_i[aij + 1];
	      {
		/* Compute complex-extra = complex-double * real. */
		double head_t, tail_t;
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = x_elem * split;
		  a1 = con - x_elem;
		  a1 = con - a1;
		  a2 = x_elem - a1;
		  con = a_elem[0] * split;
		  b1 = con - a_elem[0];
		  b1 = con - b1;
		  b2 = a_elem[0] - b1;

		  head_t = x_elem * a_elem[0];
		  tail_t =
		    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		head_prod[0] = head_t;
		tail_prod[0] = tail_t;
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = x_elem * split;
		  a1 = con - x_elem;
		  a1 = con - a1;
		  a2 = x_elem - a1;
		  con = a_elem[1] * split;
		  b1 = con - a_elem[1];
		  b1 = con - b1;
		  b2 = a_elem[1] - b1;

		  head_t = x_elem * a_elem[1];
		  tail_t =
		    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		head_prod[1] = head_t;
		tail_prod[1] = tail_t;
	      }
	      {
		double head_t, tail_t;
		double head_a, tail_a;
		double head_b, tail_b;
		/* Real part */
		head_a = head_sumA[0];
		tail_a = tail_sumA[0];
		head_b = head_prod[0];
		tail_b = tail_prod[0];
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_a + head_b;
		  bv = s1 - head_a;
		  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_a + tail_b;
		  bv = t1 - tail_a;
		  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_sumA[0] = head_t;
		tail_sumA[0] = tail_t;
		/* Imaginary part */
		head_a = head_sumA[1];
		tail_a = tail_sumA[1];
		head_b = head_prod[1];
		tail_b = tail_prod[1];
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_a + head_b;
		  bv = s1 - head_a;
		  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_a + tail_b;
		  bv = t1 - tail_a;
		  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_sumA[1] = head_t;
		tail_sumA[1] = tail_t;
	      }
	      aij += incaij;

	    }
	    /* now put the result into y_i */
	    y_i[yi] = head_sumA[0];
	    y_i[yi + 1] = head_sumA[1];
	    ai += incai;

	  }
	} else if ((beta_i[0] == 1.0 && beta_i[1] == 0.0)) {
	  /* alpha is 1.0, beta is 1.0 */

	  ai = 0;
	  bi = 0;
	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	    head_sumA[0] = head_sumA[1] = tail_sumA[0] = tail_sumA[1] = 0.0;
	    aij = ai;
	    head_sumB[0] = head_sumB[1] = tail_sumB[0] = tail_sumB[1] = 0.0;
	    bij = bi;
	    for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	      x_elem = x_i[xi];
	      a_elem[0] = a_i[aij];
	      a_elem[1] = a_i[aij + 1];
	      {
		/* Compute complex-extra = complex-double * real. */
		double head_t, tail_t;
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = x_elem * split;
		  a1 = con - x_elem;
		  a1 = con - a1;
		  a2 = x_elem - a1;
		  con = a_elem[0] * split;
		  b1 = con - a_elem[0];
		  b1 = con - b1;
		  b2 = a_elem[0] - b1;

		  head_t = x_elem * a_elem[0];
		  tail_t =
		    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		head_prod[0] = head_t;
		tail_prod[0] = tail_t;
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = x_elem * split;
		  a1 = con - x_elem;
		  a1 = con - a1;
		  a2 = x_elem - a1;
		  con = a_elem[1] * split;
		  b1 = con - a_elem[1];
		  b1 = con - b1;
		  b2 = a_elem[1] - b1;

		  head_t = x_elem * a_elem[1];
		  tail_t =
		    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		head_prod[1] = head_t;
		tail_prod[1] = tail_t;
	      }
	      {
		double head_t, tail_t;
		double head_a, tail_a;
		double head_b, tail_b;
		/* Real part */
		head_a = head_sumA[0];
		tail_a = tail_sumA[0];
		head_b = head_prod[0];
		tail_b = tail_prod[0];
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_a + head_b;
		  bv = s1 - head_a;
		  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_a + tail_b;
		  bv = t1 - tail_a;
		  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_sumA[0] = head_t;
		tail_sumA[0] = tail_t;
		/* Imaginary part */
		head_a = head_sumA[1];
		tail_a = tail_sumA[1];
		head_b = head_prod[1];
		tail_b = tail_prod[1];
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_a + head_b;
		  bv = s1 - head_a;
		  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_a + tail_b;
		  bv = t1 - tail_a;
		  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_sumA[1] = head_t;
		tail_sumA[1] = tail_t;
	      }
	      aij += incaij;
	      b_elem[0] = b_i[bij];
	      b_elem[1] = b_i[bij + 1];
	      {
		/* Compute complex-extra = complex-double * real. */
		double head_t, tail_t;
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = x_elem * split;
		  a1 = con - x_elem;
		  a1 = con - a1;
		  a2 = x_elem - a1;
		  con = b_elem[0] * split;
		  b1 = con - b_elem[0];
		  b1 = con - b1;
		  b2 = b_elem[0] - b1;

		  head_t = x_elem * b_elem[0];
		  tail_t =
		    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		head_prod[0] = head_t;
		tail_prod[0] = tail_t;
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = x_elem * split;
		  a1 = con - x_elem;
		  a1 = con - a1;
		  a2 = x_elem - a1;
		  con = b_elem[1] * split;
		  b1 = con - b_elem[1];
		  b1 = con - b1;
		  b2 = b_elem[1] - b1;

		  head_t = x_elem * b_elem[1];
		  tail_t =
		    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		head_prod[1] = head_t;
		tail_prod[1] = tail_t;
	      }
	      {
		double head_t, tail_t;
		double head_a, tail_a;
		double head_b, tail_b;
		/* Real part */
		head_a = head_sumB[0];
		tail_a = tail_sumB[0];
		head_b = head_prod[0];
		tail_b = tail_prod[0];
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_a + head_b;
		  bv = s1 - head_a;
		  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_a + tail_b;
		  bv = t1 - tail_a;
		  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_sumB[0] = head_t;
		tail_sumB[0] = tail_t;
		/* Imaginary part */
		head_a = head_sumB[1];
		tail_a = tail_sumB[1];
		head_b = head_prod[1];
		tail_b = tail_prod[1];
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_a + head_b;
		  bv = s1 - head_a;
		  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_a + tail_b;
		  bv = t1 - tail_a;
		  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_sumB[1] = head_t;
		tail_sumB[1] = tail_t;
	      }
	      bij += incbij;
	    }
	    /* now put the result into y_i */
	    head_tmp1[0] = head_sumA[0];
	    tail_tmp1[0] = tail_sumA[0];
	    head_tmp1[1] = head_sumA[1];
	    tail_tmp1[1] = tail_sumA[1];
	    head_tmp2[0] = head_sumB[0];
	    tail_tmp2[0] = tail_sumB[0];
	    head_tmp2[1] = head_sumB[1];
	    tail_tmp2[1] = tail_sumB[1];
	    {
	      double head_t, tail_t;
	      double head_a, tail_a;
	      double head_b, tail_b;
	      /* Real part */
	      head_a = head_tmp1[0];
	      tail_a = tail_tmp1[0];
	      head_b = head_tmp2[0];
	      tail_b = tail_tmp2[0];
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_a + head_b;
		bv = s1 - head_a;
		s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_a + tail_b;
		bv = t1 - tail_a;
		t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t = t1 + t2;
		tail_t = t2 - (head_t - t1);
	      }
	      head_tmp1[0] = head_t;
	      tail_tmp1[0] = tail_t;
	      /* Imaginary part */
	      head_a = head_tmp1[1];
	      tail_a = tail_tmp1[1];
	      head_b = head_tmp2[1];
	      tail_b = tail_tmp2[1];
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_a + head_b;
		bv = s1 - head_a;
		s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_a + tail_b;
		bv = t1 - tail_a;
		t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t = t1 + t2;
		tail_t = t2 - (head_t - t1);
	      }
	      head_tmp1[1] = head_t;
	      tail_tmp1[1] = tail_t;
	    }
	    y_i[yi] = head_tmp1[0];
	    y_i[yi + 1] = head_tmp1[1];
	    ai += incai;
	    bi += incbi;
	  }
	} else {
	  /* alpha is 1.0, beta is other */

	  ai = 0;
	  bi = 0;
	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	    head_sumA[0] = head_sumA[1] = tail_sumA[0] = tail_sumA[1] = 0.0;
	    aij = ai;
	    head_sumB[0] = head_sumB[1] = tail_sumB[0] = tail_sumB[1] = 0.0;
	    bij = bi;
	    for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	      x_elem = x_i[xi];
	      a_elem[0] = a_i[aij];
	      a_elem[1] = a_i[aij + 1];
	      {
		/* Compute complex-extra = complex-double * real. */
		double head_t, tail_t;
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = x_elem * split;
		  a1 = con - x_elem;
		  a1 = con - a1;
		  a2 = x_elem - a1;
		  con = a_elem[0] * split;
		  b1 = con - a_elem[0];
		  b1 = con - b1;
		  b2 = a_elem[0] - b1;

		  head_t = x_elem * a_elem[0];
		  tail_t =
		    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		head_prod[0] = head_t;
		tail_prod[0] = tail_t;
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = x_elem * split;
		  a1 = con - x_elem;
		  a1 = con - a1;
		  a2 = x_elem - a1;
		  con = a_elem[1] * split;
		  b1 = con - a_elem[1];
		  b1 = con - b1;
		  b2 = a_elem[1] - b1;

		  head_t = x_elem * a_elem[1];
		  tail_t =
		    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		head_prod[1] = head_t;
		tail_prod[1] = tail_t;
	      }
	      {
		double head_t, tail_t;
		double head_a, tail_a;
		double head_b, tail_b;
		/* Real part */
		head_a = head_sumA[0];
		tail_a = tail_sumA[0];
		head_b = head_prod[0];
		tail_b = tail_prod[0];
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_a + head_b;
		  bv = s1 - head_a;
		  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_a + tail_b;
		  bv = t1 - tail_a;
		  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_sumA[0] = head_t;
		tail_sumA[0] = tail_t;
		/* Imaginary part */
		head_a = head_sumA[1];
		tail_a = tail_sumA[1];
		head_b = head_prod[1];
		tail_b = tail_prod[1];
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_a + head_b;
		  bv = s1 - head_a;
		  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_a + tail_b;
		  bv = t1 - tail_a;
		  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_sumA[1] = head_t;
		tail_sumA[1] = tail_t;
	      }
	      aij += incaij;
	      b_elem[0] = b_i[bij];
	      b_elem[1] = b_i[bij + 1];
	      {
		/* Compute complex-extra = complex-double * real. */
		double head_t, tail_t;
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = x_elem * split;
		  a1 = con - x_elem;
		  a1 = con - a1;
		  a2 = x_elem - a1;
		  con = b_elem[0] * split;
		  b1 = con - b_elem[0];
		  b1 = con - b1;
		  b2 = b_elem[0] - b1;

		  head_t = x_elem * b_elem[0];
		  tail_t =
		    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		head_prod[0] = head_t;
		tail_prod[0] = tail_t;
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = x_elem * split;
		  a1 = con - x_elem;
		  a1 = con - a1;
		  a2 = x_elem - a1;
		  con = b_elem[1] * split;
		  b1 = con - b_elem[1];
		  b1 = con - b1;
		  b2 = b_elem[1] - b1;

		  head_t = x_elem * b_elem[1];
		  tail_t =
		    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		head_prod[1] = head_t;
		tail_prod[1] = tail_t;
	      }
	      {
		double head_t, tail_t;
		double head_a, tail_a;
		double head_b, tail_b;
		/* Real part */
		head_a = head_sumB[0];
		tail_a = tail_sumB[0];
		head_b = head_prod[0];
		tail_b = tail_prod[0];
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_a + head_b;
		  bv = s1 - head_a;
		  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_a + tail_b;
		  bv = t1 - tail_a;
		  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_sumB[0] = head_t;
		tail_sumB[0] = tail_t;
		/* Imaginary part */
		head_a = head_sumB[1];
		tail_a = tail_sumB[1];
		head_b = head_prod[1];
		tail_b = tail_prod[1];
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_a + head_b;
		  bv = s1 - head_a;
		  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_a + tail_b;
		  bv = t1 - tail_a;
		  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_sumB[1] = head_t;
		tail_sumB[1] = tail_t;
	      }
	      bij += incbij;
	    }
	    /* now put the result into y_i */
	    head_tmp1[0] = head_sumA[0];
	    tail_tmp1[0] = tail_sumA[0];
	    head_tmp1[1] = head_sumA[1];
	    tail_tmp1[1] = tail_sumA[1];
	    {
	      /* Compute complex-extra = complex-extra * complex-double. */
	      double head_a0, tail_a0;
	      double head_a1, tail_a1;
	      double head_t1, tail_t1;
	      double head_t2, tail_t2;
	      head_a0 = head_sumB[0];
	      tail_a0 = tail_sumB[0];
	      head_a1 = head_sumB[1];
	      tail_a1 = tail_sumB[1];
	      /* real part */
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a0 * split;
		a11 = con - head_a0;
		a11 = con - a11;
		a21 = head_a0 - a11;
		con = beta_i[0] * split;
		b1 = con - beta_i[0];
		b1 = con - b1;
		b2 = beta_i[0] - b1;

		c11 = head_a0 * beta_i[0];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a0 * beta_i[0];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a1 * split;
		a11 = con - head_a1;
		a11 = con - a11;
		a21 = head_a1 - a11;
		con = beta_i[1] * split;
		b1 = con - beta_i[1];
		b1 = con - b1;
		b2 = beta_i[1] - b1;

		c11 = head_a1 * beta_i[1];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a1 * beta_i[1];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t2 = t1 + t2;
		tail_t2 = t2 - (head_t2 - t1);
	      }
	      head_t2 = -head_t2;
	      tail_t2 = -tail_t2;
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_t1 + head_t2;
		bv = s1 - head_t1;
		s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_t1 + tail_t2;
		bv = t1 - tail_t1;
		t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      head_tmp2[0] = head_t1;
	      tail_tmp2[0] = tail_t1;
	      /* imaginary part */
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a1 * split;
		a11 = con - head_a1;
		a11 = con - a11;
		a21 = head_a1 - a11;
		con = beta_i[0] * split;
		b1 = con - beta_i[0];
		b1 = con - b1;
		b2 = beta_i[0] - b1;

		c11 = head_a1 * beta_i[0];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a1 * beta_i[0];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a0 * split;
		a11 = con - head_a0;
		a11 = con - a11;
		a21 = head_a0 - a11;
		con = beta_i[1] * split;
		b1 = con - beta_i[1];
		b1 = con - b1;
		b2 = beta_i[1] - b1;

		c11 = head_a0 * beta_i[1];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a0 * beta_i[1];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t2 = t1 + t2;
		tail_t2 = t2 - (head_t2 - t1);
	      }
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_t1 + head_t2;
		bv = s1 - head_t1;
		s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_t1 + tail_t2;
		bv = t1 - tail_t1;
		t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      head_tmp2[1] = head_t1;
	      tail_tmp2[1] = tail_t1;
	    }

	    {
	      double head_t, tail_t;
	      double head_a, tail_a;
	      double head_b, tail_b;
	      /* Real part */
	      head_a = head_tmp1[0];
	      tail_a = tail_tmp1[0];
	      head_b = head_tmp2[0];
	      tail_b = tail_tmp2[0];
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_a + head_b;
		bv = s1 - head_a;
		s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_a + tail_b;
		bv = t1 - tail_a;
		t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t = t1 + t2;
		tail_t = t2 - (head_t - t1);
	      }
	      head_tmp1[0] = head_t;
	      tail_tmp1[0] = tail_t;
	      /* Imaginary part */
	      head_a = head_tmp1[1];
	      tail_a = tail_tmp1[1];
	      head_b = head_tmp2[1];
	      tail_b = tail_tmp2[1];
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_a + head_b;
		bv = s1 - head_a;
		s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_a + tail_b;
		bv = t1 - tail_a;
		t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t = t1 + t2;
		tail_t = t2 - (head_t - t1);
	      }
	      head_tmp1[1] = head_t;
	      tail_tmp1[1] = tail_t;
	    }
	    y_i[yi] = head_tmp1[0];
	    y_i[yi + 1] = head_tmp1[1];
	    ai += incai;
	    bi += incbi;
	  }
	}
      } else {
	if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	  /* alpha is other, beta is 0.0 */

	  ai = 0;

	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	    head_sumA[0] = head_sumA[1] = tail_sumA[0] = tail_sumA[1] = 0.0;
	    aij = ai;

	    for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	      x_elem = x_i[xi];
	      a_elem[0] = a_i[aij];
	      a_elem[1] = a_i[aij + 1];
	      {
		/* Compute complex-extra = complex-double * real. */
		double head_t, tail_t;
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = x_elem * split;
		  a1 = con - x_elem;
		  a1 = con - a1;
		  a2 = x_elem - a1;
		  con = a_elem[0] * split;
		  b1 = con - a_elem[0];
		  b1 = con - b1;
		  b2 = a_elem[0] - b1;

		  head_t = x_elem * a_elem[0];
		  tail_t =
		    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		head_prod[0] = head_t;
		tail_prod[0] = tail_t;
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = x_elem * split;
		  a1 = con - x_elem;
		  a1 = con - a1;
		  a2 = x_elem - a1;
		  con = a_elem[1] * split;
		  b1 = con - a_elem[1];
		  b1 = con - b1;
		  b2 = a_elem[1] - b1;

		  head_t = x_elem * a_elem[1];
		  tail_t =
		    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		head_prod[1] = head_t;
		tail_prod[1] = tail_t;
	      }
	      {
		double head_t, tail_t;
		double head_a, tail_a;
		double head_b, tail_b;
		/* Real part */
		head_a = head_sumA[0];
		tail_a = tail_sumA[0];
		head_b = head_prod[0];
		tail_b = tail_prod[0];
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_a + head_b;
		  bv = s1 - head_a;
		  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_a + tail_b;
		  bv = t1 - tail_a;
		  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_sumA[0] = head_t;
		tail_sumA[0] = tail_t;
		/* Imaginary part */
		head_a = head_sumA[1];
		tail_a = tail_sumA[1];
		head_b = head_prod[1];
		tail_b = tail_prod[1];
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_a + head_b;
		  bv = s1 - head_a;
		  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_a + tail_b;
		  bv = t1 - tail_a;
		  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_sumA[1] = head_t;
		tail_sumA[1] = tail_t;
	      }
	      aij += incaij;

	    }
	    /* now put the result into y_i */
	    {
	      /* Compute complex-extra = complex-extra * complex-double. */
	      double head_a0, tail_a0;
	      double head_a1, tail_a1;
	      double head_t1, tail_t1;
	      double head_t2, tail_t2;
	      head_a0 = head_sumA[0];
	      tail_a0 = tail_sumA[0];
	      head_a1 = head_sumA[1];
	      tail_a1 = tail_sumA[1];
	      /* real part */
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a0 * split;
		a11 = con - head_a0;
		a11 = con - a11;
		a21 = head_a0 - a11;
		con = alpha_i[0] * split;
		b1 = con - alpha_i[0];
		b1 = con - b1;
		b2 = alpha_i[0] - b1;

		c11 = head_a0 * alpha_i[0];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a0 * alpha_i[0];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a1 * split;
		a11 = con - head_a1;
		a11 = con - a11;
		a21 = head_a1 - a11;
		con = alpha_i[1] * split;
		b1 = con - alpha_i[1];
		b1 = con - b1;
		b2 = alpha_i[1] - b1;

		c11 = head_a1 * alpha_i[1];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a1 * alpha_i[1];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t2 = t1 + t2;
		tail_t2 = t2 - (head_t2 - t1);
	      }
	      head_t2 = -head_t2;
	      tail_t2 = -tail_t2;
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_t1 + head_t2;
		bv = s1 - head_t1;
		s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_t1 + tail_t2;
		bv = t1 - tail_t1;
		t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      head_tmp1[0] = head_t1;
	      tail_tmp1[0] = tail_t1;
	      /* imaginary part */
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a1 * split;
		a11 = con - head_a1;
		a11 = con - a11;
		a21 = head_a1 - a11;
		con = alpha_i[0] * split;
		b1 = con - alpha_i[0];
		b1 = con - b1;
		b2 = alpha_i[0] - b1;

		c11 = head_a1 * alpha_i[0];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a1 * alpha_i[0];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a0 * split;
		a11 = con - head_a0;
		a11 = con - a11;
		a21 = head_a0 - a11;
		con = alpha_i[1] * split;
		b1 = con - alpha_i[1];
		b1 = con - b1;
		b2 = alpha_i[1] - b1;

		c11 = head_a0 * alpha_i[1];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a0 * alpha_i[1];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t2 = t1 + t2;
		tail_t2 = t2 - (head_t2 - t1);
	      }
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_t1 + head_t2;
		bv = s1 - head_t1;
		s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_t1 + tail_t2;
		bv = t1 - tail_t1;
		t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      head_tmp1[1] = head_t1;
	      tail_tmp1[1] = tail_t1;
	    }

	    y_i[yi] = head_tmp1[0];
	    y_i[yi + 1] = head_tmp1[1];
	    ai += incai;

	  }
	} else if ((beta_i[0] == 1.0 && beta_i[1] == 0.0)) {
	  /* alpha is other, beta is 1.0 */

	  ai = 0;
	  bi = 0;
	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	    head_sumA[0] = head_sumA[1] = tail_sumA[0] = tail_sumA[1] = 0.0;
	    aij = ai;
	    head_sumB[0] = head_sumB[1] = tail_sumB[0] = tail_sumB[1] = 0.0;
	    bij = bi;
	    for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	      x_elem = x_i[xi];
	      a_elem[0] = a_i[aij];
	      a_elem[1] = a_i[aij + 1];
	      {
		/* Compute complex-extra = complex-double * real. */
		double head_t, tail_t;
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = x_elem * split;
		  a1 = con - x_elem;
		  a1 = con - a1;
		  a2 = x_elem - a1;
		  con = a_elem[0] * split;
		  b1 = con - a_elem[0];
		  b1 = con - b1;
		  b2 = a_elem[0] - b1;

		  head_t = x_elem * a_elem[0];
		  tail_t =
		    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		head_prod[0] = head_t;
		tail_prod[0] = tail_t;
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = x_elem * split;
		  a1 = con - x_elem;
		  a1 = con - a1;
		  a2 = x_elem - a1;
		  con = a_elem[1] * split;
		  b1 = con - a_elem[1];
		  b1 = con - b1;
		  b2 = a_elem[1] - b1;

		  head_t = x_elem * a_elem[1];
		  tail_t =
		    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		head_prod[1] = head_t;
		tail_prod[1] = tail_t;
	      }
	      {
		double head_t, tail_t;
		double head_a, tail_a;
		double head_b, tail_b;
		/* Real part */
		head_a = head_sumA[0];
		tail_a = tail_sumA[0];
		head_b = head_prod[0];
		tail_b = tail_prod[0];
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_a + head_b;
		  bv = s1 - head_a;
		  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_a + tail_b;
		  bv = t1 - tail_a;
		  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_sumA[0] = head_t;
		tail_sumA[0] = tail_t;
		/* Imaginary part */
		head_a = head_sumA[1];
		tail_a = tail_sumA[1];
		head_b = head_prod[1];
		tail_b = tail_prod[1];
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_a + head_b;
		  bv = s1 - head_a;
		  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_a + tail_b;
		  bv = t1 - tail_a;
		  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_sumA[1] = head_t;
		tail_sumA[1] = tail_t;
	      }
	      aij += incaij;
	      b_elem[0] = b_i[bij];
	      b_elem[1] = b_i[bij + 1];
	      {
		/* Compute complex-extra = complex-double * real. */
		double head_t, tail_t;
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = x_elem * split;
		  a1 = con - x_elem;
		  a1 = con - a1;
		  a2 = x_elem - a1;
		  con = b_elem[0] * split;
		  b1 = con - b_elem[0];
		  b1 = con - b1;
		  b2 = b_elem[0] - b1;

		  head_t = x_elem * b_elem[0];
		  tail_t =
		    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		head_prod[0] = head_t;
		tail_prod[0] = tail_t;
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = x_elem * split;
		  a1 = con - x_elem;
		  a1 = con - a1;
		  a2 = x_elem - a1;
		  con = b_elem[1] * split;
		  b1 = con - b_elem[1];
		  b1 = con - b1;
		  b2 = b_elem[1] - b1;

		  head_t = x_elem * b_elem[1];
		  tail_t =
		    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		head_prod[1] = head_t;
		tail_prod[1] = tail_t;
	      }
	      {
		double head_t, tail_t;
		double head_a, tail_a;
		double head_b, tail_b;
		/* Real part */
		head_a = head_sumB[0];
		tail_a = tail_sumB[0];
		head_b = head_prod[0];
		tail_b = tail_prod[0];
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_a + head_b;
		  bv = s1 - head_a;
		  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_a + tail_b;
		  bv = t1 - tail_a;
		  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_sumB[0] = head_t;
		tail_sumB[0] = tail_t;
		/* Imaginary part */
		head_a = head_sumB[1];
		tail_a = tail_sumB[1];
		head_b = head_prod[1];
		tail_b = tail_prod[1];
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_a + head_b;
		  bv = s1 - head_a;
		  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_a + tail_b;
		  bv = t1 - tail_a;
		  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_sumB[1] = head_t;
		tail_sumB[1] = tail_t;
	      }
	      bij += incbij;
	    }
	    /* now put the result into y_i */
	    {
	      /* Compute complex-extra = complex-extra * complex-double. */
	      double head_a0, tail_a0;
	      double head_a1, tail_a1;
	      double head_t1, tail_t1;
	      double head_t2, tail_t2;
	      head_a0 = head_sumA[0];
	      tail_a0 = tail_sumA[0];
	      head_a1 = head_sumA[1];
	      tail_a1 = tail_sumA[1];
	      /* real part */
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a0 * split;
		a11 = con - head_a0;
		a11 = con - a11;
		a21 = head_a0 - a11;
		con = alpha_i[0] * split;
		b1 = con - alpha_i[0];
		b1 = con - b1;
		b2 = alpha_i[0] - b1;

		c11 = head_a0 * alpha_i[0];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a0 * alpha_i[0];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a1 * split;
		a11 = con - head_a1;
		a11 = con - a11;
		a21 = head_a1 - a11;
		con = alpha_i[1] * split;
		b1 = con - alpha_i[1];
		b1 = con - b1;
		b2 = alpha_i[1] - b1;

		c11 = head_a1 * alpha_i[1];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a1 * alpha_i[1];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t2 = t1 + t2;
		tail_t2 = t2 - (head_t2 - t1);
	      }
	      head_t2 = -head_t2;
	      tail_t2 = -tail_t2;
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_t1 + head_t2;
		bv = s1 - head_t1;
		s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_t1 + tail_t2;
		bv = t1 - tail_t1;
		t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      head_tmp1[0] = head_t1;
	      tail_tmp1[0] = tail_t1;
	      /* imaginary part */
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a1 * split;
		a11 = con - head_a1;
		a11 = con - a11;
		a21 = head_a1 - a11;
		con = alpha_i[0] * split;
		b1 = con - alpha_i[0];
		b1 = con - b1;
		b2 = alpha_i[0] - b1;

		c11 = head_a1 * alpha_i[0];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a1 * alpha_i[0];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a0 * split;
		a11 = con - head_a0;
		a11 = con - a11;
		a21 = head_a0 - a11;
		con = alpha_i[1] * split;
		b1 = con - alpha_i[1];
		b1 = con - b1;
		b2 = alpha_i[1] - b1;

		c11 = head_a0 * alpha_i[1];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a0 * alpha_i[1];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t2 = t1 + t2;
		tail_t2 = t2 - (head_t2 - t1);
	      }
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_t1 + head_t2;
		bv = s1 - head_t1;
		s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_t1 + tail_t2;
		bv = t1 - tail_t1;
		t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      head_tmp1[1] = head_t1;
	      tail_tmp1[1] = tail_t1;
	    }

	    head_tmp2[0] = head_sumB[0];
	    tail_tmp2[0] = tail_sumB[0];
	    head_tmp2[1] = head_sumB[1];
	    tail_tmp2[1] = tail_sumB[1];
	    {
	      double head_t, tail_t;
	      double head_a, tail_a;
	      double head_b, tail_b;
	      /* Real part */
	      head_a = head_tmp1[0];
	      tail_a = tail_tmp1[0];
	      head_b = head_tmp2[0];
	      tail_b = tail_tmp2[0];
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_a + head_b;
		bv = s1 - head_a;
		s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_a + tail_b;
		bv = t1 - tail_a;
		t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t = t1 + t2;
		tail_t = t2 - (head_t - t1);
	      }
	      head_tmp1[0] = head_t;
	      tail_tmp1[0] = tail_t;
	      /* Imaginary part */
	      head_a = head_tmp1[1];
	      tail_a = tail_tmp1[1];
	      head_b = head_tmp2[1];
	      tail_b = tail_tmp2[1];
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_a + head_b;
		bv = s1 - head_a;
		s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_a + tail_b;
		bv = t1 - tail_a;
		t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t = t1 + t2;
		tail_t = t2 - (head_t - t1);
	      }
	      head_tmp1[1] = head_t;
	      tail_tmp1[1] = tail_t;
	    }
	    y_i[yi] = head_tmp1[0];
	    y_i[yi + 1] = head_tmp1[1];
	    ai += incai;
	    bi += incbi;
	  }
	} else {
	  /* most general form, alpha, beta are other */

	  ai = 0;
	  bi = 0;
	  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
	    head_sumA[0] = head_sumA[1] = tail_sumA[0] = tail_sumA[1] = 0.0;
	    aij = ai;
	    head_sumB[0] = head_sumB[1] = tail_sumB[0] = tail_sumB[1] = 0.0;
	    bij = bi;
	    for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
	      x_elem = x_i[xi];
	      a_elem[0] = a_i[aij];
	      a_elem[1] = a_i[aij + 1];
	      {
		/* Compute complex-extra = complex-double * real. */
		double head_t, tail_t;
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = x_elem * split;
		  a1 = con - x_elem;
		  a1 = con - a1;
		  a2 = x_elem - a1;
		  con = a_elem[0] * split;
		  b1 = con - a_elem[0];
		  b1 = con - b1;
		  b2 = a_elem[0] - b1;

		  head_t = x_elem * a_elem[0];
		  tail_t =
		    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		head_prod[0] = head_t;
		tail_prod[0] = tail_t;
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = x_elem * split;
		  a1 = con - x_elem;
		  a1 = con - a1;
		  a2 = x_elem - a1;
		  con = a_elem[1] * split;
		  b1 = con - a_elem[1];
		  b1 = con - b1;
		  b2 = a_elem[1] - b1;

		  head_t = x_elem * a_elem[1];
		  tail_t =
		    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		head_prod[1] = head_t;
		tail_prod[1] = tail_t;
	      }
	      {
		double head_t, tail_t;
		double head_a, tail_a;
		double head_b, tail_b;
		/* Real part */
		head_a = head_sumA[0];
		tail_a = tail_sumA[0];
		head_b = head_prod[0];
		tail_b = tail_prod[0];
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_a + head_b;
		  bv = s1 - head_a;
		  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_a + tail_b;
		  bv = t1 - tail_a;
		  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_sumA[0] = head_t;
		tail_sumA[0] = tail_t;
		/* Imaginary part */
		head_a = head_sumA[1];
		tail_a = tail_sumA[1];
		head_b = head_prod[1];
		tail_b = tail_prod[1];
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_a + head_b;
		  bv = s1 - head_a;
		  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_a + tail_b;
		  bv = t1 - tail_a;
		  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_sumA[1] = head_t;
		tail_sumA[1] = tail_t;
	      }
	      aij += incaij;
	      b_elem[0] = b_i[bij];
	      b_elem[1] = b_i[bij + 1];
	      {
		/* Compute complex-extra = complex-double * real. */
		double head_t, tail_t;
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = x_elem * split;
		  a1 = con - x_elem;
		  a1 = con - a1;
		  a2 = x_elem - a1;
		  con = b_elem[0] * split;
		  b1 = con - b_elem[0];
		  b1 = con - b1;
		  b2 = b_elem[0] - b1;

		  head_t = x_elem * b_elem[0];
		  tail_t =
		    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		head_prod[0] = head_t;
		tail_prod[0] = tail_t;
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = x_elem * split;
		  a1 = con - x_elem;
		  a1 = con - a1;
		  a2 = x_elem - a1;
		  con = b_elem[1] * split;
		  b1 = con - b_elem[1];
		  b1 = con - b1;
		  b2 = b_elem[1] - b1;

		  head_t = x_elem * b_elem[1];
		  tail_t =
		    (((a1 * b1 - head_t) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		head_prod[1] = head_t;
		tail_prod[1] = tail_t;
	      }
	      {
		double head_t, tail_t;
		double head_a, tail_a;
		double head_b, tail_b;
		/* Real part */
		head_a = head_sumB[0];
		tail_a = tail_sumB[0];
		head_b = head_prod[0];
		tail_b = tail_prod[0];
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_a + head_b;
		  bv = s1 - head_a;
		  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_a + tail_b;
		  bv = t1 - tail_a;
		  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_sumB[0] = head_t;
		tail_sumB[0] = tail_t;
		/* Imaginary part */
		head_a = head_sumB[1];
		tail_a = tail_sumB[1];
		head_b = head_prod[1];
		tail_b = tail_prod[1];
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_a + head_b;
		  bv = s1 - head_a;
		  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_a + tail_b;
		  bv = t1 - tail_a;
		  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_sumB[1] = head_t;
		tail_sumB[1] = tail_t;
	      }
	      bij += incbij;
	    }
	    /* now put the result into y_i */
	    {
	      /* Compute complex-extra = complex-extra * complex-double. */
	      double head_a0, tail_a0;
	      double head_a1, tail_a1;
	      double head_t1, tail_t1;
	      double head_t2, tail_t2;
	      head_a0 = head_sumA[0];
	      tail_a0 = tail_sumA[0];
	      head_a1 = head_sumA[1];
	      tail_a1 = tail_sumA[1];
	      /* real part */
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a0 * split;
		a11 = con - head_a0;
		a11 = con - a11;
		a21 = head_a0 - a11;
		con = alpha_i[0] * split;
		b1 = con - alpha_i[0];
		b1 = con - b1;
		b2 = alpha_i[0] - b1;

		c11 = head_a0 * alpha_i[0];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a0 * alpha_i[0];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a1 * split;
		a11 = con - head_a1;
		a11 = con - a11;
		a21 = head_a1 - a11;
		con = alpha_i[1] * split;
		b1 = con - alpha_i[1];
		b1 = con - b1;
		b2 = alpha_i[1] - b1;

		c11 = head_a1 * alpha_i[1];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a1 * alpha_i[1];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t2 = t1 + t2;
		tail_t2 = t2 - (head_t2 - t1);
	      }
	      head_t2 = -head_t2;
	      tail_t2 = -tail_t2;
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_t1 + head_t2;
		bv = s1 - head_t1;
		s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_t1 + tail_t2;
		bv = t1 - tail_t1;
		t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      head_tmp1[0] = head_t1;
	      tail_tmp1[0] = tail_t1;
	      /* imaginary part */
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a1 * split;
		a11 = con - head_a1;
		a11 = con - a11;
		a21 = head_a1 - a11;
		con = alpha_i[0] * split;
		b1 = con - alpha_i[0];
		b1 = con - b1;
		b2 = alpha_i[0] - b1;

		c11 = head_a1 * alpha_i[0];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a1 * alpha_i[0];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a0 * split;
		a11 = con - head_a0;
		a11 = con - a11;
		a21 = head_a0 - a11;
		con = alpha_i[1] * split;
		b1 = con - alpha_i[1];
		b1 = con - b1;
		b2 = alpha_i[1] - b1;

		c11 = head_a0 * alpha_i[1];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a0 * alpha_i[1];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t2 = t1 + t2;
		tail_t2 = t2 - (head_t2 - t1);
	      }
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_t1 + head_t2;
		bv = s1 - head_t1;
		s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_t1 + tail_t2;
		bv = t1 - tail_t1;
		t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      head_tmp1[1] = head_t1;
	      tail_tmp1[1] = tail_t1;
	    }

	    {
	      /* Compute complex-extra = complex-extra * complex-double. */
	      double head_a0, tail_a0;
	      double head_a1, tail_a1;
	      double head_t1, tail_t1;
	      double head_t2, tail_t2;
	      head_a0 = head_sumB[0];
	      tail_a0 = tail_sumB[0];
	      head_a1 = head_sumB[1];
	      tail_a1 = tail_sumB[1];
	      /* real part */
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a0 * split;
		a11 = con - head_a0;
		a11 = con - a11;
		a21 = head_a0 - a11;
		con = beta_i[0] * split;
		b1 = con - beta_i[0];
		b1 = con - b1;
		b2 = beta_i[0] - b1;

		c11 = head_a0 * beta_i[0];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a0 * beta_i[0];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a1 * split;
		a11 = con - head_a1;
		a11 = con - a11;
		a21 = head_a1 - a11;
		con = beta_i[1] * split;
		b1 = con - beta_i[1];
		b1 = con - b1;
		b2 = beta_i[1] - b1;

		c11 = head_a1 * beta_i[1];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a1 * beta_i[1];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t2 = t1 + t2;
		tail_t2 = t2 - (head_t2 - t1);
	      }
	      head_t2 = -head_t2;
	      tail_t2 = -tail_t2;
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_t1 + head_t2;
		bv = s1 - head_t1;
		s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_t1 + tail_t2;
		bv = t1 - tail_t1;
		t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      head_tmp2[0] = head_t1;
	      tail_tmp2[0] = tail_t1;
	      /* imaginary part */
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a1 * split;
		a11 = con - head_a1;
		a11 = con - a11;
		a21 = head_a1 - a11;
		con = beta_i[0] * split;
		b1 = con - beta_i[0];
		b1 = con - b1;
		b2 = beta_i[0] - b1;

		c11 = head_a1 * beta_i[0];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a1 * beta_i[0];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_a0 * split;
		a11 = con - head_a0;
		a11 = con - a11;
		a21 = head_a0 - a11;
		con = beta_i[1] * split;
		b1 = con - beta_i[1];
		b1 = con - b1;
		b2 = beta_i[1] - b1;

		c11 = head_a0 * beta_i[1];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_a0 * beta_i[1];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t2 = t1 + t2;
		tail_t2 = t2 - (head_t2 - t1);
	      }
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_t1 + head_t2;
		bv = s1 - head_t1;
		s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_t1 + tail_t2;
		bv = t1 - tail_t1;
		t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      head_tmp2[1] = head_t1;
	      tail_tmp2[1] = tail_t1;
	    }

	    {
	      double head_t, tail_t;
	      double head_a, tail_a;
	      double head_b, tail_b;
	      /* Real part */
	      head_a = head_tmp1[0];
	      tail_a = tail_tmp1[0];
	      head_b = head_tmp2[0];
	      tail_b = tail_tmp2[0];
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_a + head_b;
		bv = s1 - head_a;
		s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_a + tail_b;
		bv = t1 - tail_a;
		t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t = t1 + t2;
		tail_t = t2 - (head_t - t1);
	      }
	      head_tmp1[0] = head_t;
	      tail_tmp1[0] = tail_t;
	      /* Imaginary part */
	      head_a = head_tmp1[1];
	      tail_a = tail_tmp1[1];
	      head_b = head_tmp2[1];
	      tail_b = tail_tmp2[1];
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_a + head_b;
		bv = s1 - head_a;
		s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_a + tail_b;
		bv = t1 - tail_a;
		t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t = t1 + t2;
		tail_t = t2 - (head_t - t1);
	      }
	      head_tmp1[1] = head_t;
	      tail_tmp1[1] = tail_t;
	    }
	    y_i[yi] = head_tmp1[0];
	    y_i[yi + 1] = head_tmp1[1];
	    ai += incai;
	    bi += incbi;
	  }
	}
      }
      FPU_FIX_STOP;
    }
    break;

  default:
    {
      BLAS_error(routine_name, -14, prec, 0);
    }
    break;
  }

}				/* end BLAS_zge_sum_mv_z_d */
