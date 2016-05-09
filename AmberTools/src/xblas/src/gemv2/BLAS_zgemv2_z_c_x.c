#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_zgemv2_z_c_x(enum blas_order_type order, enum blas_trans_type trans,
		       int m, int n, const void *alpha, const void *a,
		       int lda, const void *head_x, const void *tail_x,
		       int incx, const void *beta, void *y, int incy,
		       enum blas_prec_type prec)

/*
 * Purpose
 * =======
 *
 * Computes y = alpha * op(A) * head_x + alpha * op(A) * tail_x + beta * y,
 * where A is a general matrix.
 *
 * Arguments
 * =========
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans        (input) blas_trans_type
 *              Transpose of A: no trans, trans, or conjugate trans
 *
 * m            (input) int
 *              Dimension of A
 *
 * n            (input) int
 *              Dimension of A and the length of vector x and z
 *
 * alpha        (input) const void*
 *              
 * A            (input) const void*
 *
 * lda          (input) int 
 *              Leading dimension of A
 *
 * head_x
 * tail_x       (input) const void*
 * 
 * incx         (input) int
 *              The stride for vector x.
 *
 * beta         (input) const void*
 *
 * y            (input) const void*
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
  static const char routine_name[] = "BLAS_zgemv2_z_c_x";
  switch (prec) {
  case blas_prec_single:
  case blas_prec_double:
  case blas_prec_indigenous:{

      int i, j;
      int iy, jx, kx, ky;
      int lenx, leny;
      int ai, aij;
      int incai, incaij;

      const double *a_i = (double *) a;
      const float *head_x_i = (float *) head_x;
      const float *tail_x_i = (float *) tail_x;
      double *y_i = (double *) y;
      double *alpha_i = (double *) alpha;
      double *beta_i = (double *) beta;
      double a_elem[2];
      float x_elem[2];
      double y_elem[2];
      double prod[2];
      double sum[2];
      double sum2[2];
      double tmp1[2];
      double tmp2[2];


      /* all error calls */
      if (m < 0)
	BLAS_error(routine_name, -3, m, 0);
      else if (n <= 0)
	BLAS_error(routine_name, -4, n, 0);
      else if (incx == 0)
	BLAS_error(routine_name, -10, incx, 0);
      else if (incy == 0)
	BLAS_error(routine_name, -13, incy, 0);

      if ((order == blas_rowmajor) && (trans == blas_no_trans)) {
	lenx = n;
	leny = m;
	incai = lda;
	incaij = 1;
      } else if ((order == blas_rowmajor) && (trans != blas_no_trans)) {
	lenx = m;
	leny = n;
	incai = 1;
	incaij = lda;
      } else if ((order == blas_colmajor) && (trans == blas_no_trans)) {
	lenx = n;
	leny = m;
	incai = 1;
	incaij = lda;
      } else {			/* colmajor and blas_trans */
	lenx = m;
	leny = n;
	incai = lda;
	incaij = 1;
      }

      if (lda < leny)
	BLAS_error(routine_name, -7, lda, NULL);



      incx *= 2;
      incy *= 2;
      incai *= 2;
      incaij *= 2;

      if (incx > 0)
	kx = 0;
      else
	kx = (1 - lenx) * incx;
      if (incy > 0)
	ky = 0;
      else
	ky = (1 - leny) * incy;

      /* No extra-precision needed for alpha = 0 */
      if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
	if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	  iy = ky;
	  for (i = 0; i < leny; i++) {
	    y_i[iy] = 0.0;
	    y_i[iy + 1] = 0.0;
	    iy += incy;
	  }
	} else if (!(beta_i[0] == 0.0 && beta_i[1] == 0.0)) {
	  iy = ky;
	  for (i = 0; i < leny; i++) {
	    y_elem[0] = y_i[iy];
	    y_elem[1] = y_i[iy + 1];
	    {
	      tmp1[0] =
		(double) y_elem[0] * beta_i[0] -
		(double) y_elem[1] * beta_i[1];
	      tmp1[1] =
		(double) y_elem[0] * beta_i[1] +
		(double) y_elem[1] * beta_i[0];
	    }
	    y_i[iy] = tmp1[0];
	    y_i[iy + 1] = tmp1[1];
	    iy += incy;
	  }
	}
      } else {			/* alpha != 0 */
	if (trans == blas_conj_trans) {

	  /* if beta = 0, we can save m multiplies:
	     y = alpha*A*head_x + alpha*A*tail_x  */
	  if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	    if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {
	      /* save m more multiplies if alpha = 1 */
	      ai = 0;
	      iy = ky;
	      for (i = 0; i < leny; i++) {
		sum[0] = sum[1] = 0.0;
		sum2[0] = sum2[1] = 0.0;
		aij = ai;
		jx = kx;
		for (j = 0; j < lenx; j++) {
		  a_elem[0] = a_i[aij];
		  a_elem[1] = a_i[aij + 1];
		  a_elem[1] = -a_elem[1];
		  x_elem[0] = head_x_i[jx];
		  x_elem[1] = head_x_i[jx + 1];
		  {
		    prod[0] =
		      (double) a_elem[0] * x_elem[0] -
		      (double) a_elem[1] * x_elem[1];
		    prod[1] =
		      (double) a_elem[0] * x_elem[1] +
		      (double) a_elem[1] * x_elem[0];
		  }
		  sum[0] = sum[0] + prod[0];
		  sum[1] = sum[1] + prod[1];
		  x_elem[0] = tail_x_i[jx];
		  x_elem[1] = tail_x_i[jx + 1];
		  {
		    prod[0] =
		      (double) a_elem[0] * x_elem[0] -
		      (double) a_elem[1] * x_elem[1];
		    prod[1] =
		      (double) a_elem[0] * x_elem[1] +
		      (double) a_elem[1] * x_elem[0];
		  }
		  sum2[0] = sum2[0] + prod[0];
		  sum2[1] = sum2[1] + prod[1];
		  aij += incaij;
		  jx += incx;
		}
		sum[0] = sum[0] + sum2[0];
		sum[1] = sum[1] + sum2[1];
		y_i[iy] = sum[0];
		y_i[iy + 1] = sum[1];
		ai += incai;
		iy += incy;
	      }			/* end for */
	    } else {		/* alpha != 1 */
	      ai = 0;
	      iy = ky;
	      for (i = 0; i < leny; i++) {
		sum[0] = sum[1] = 0.0;
		sum2[0] = sum2[1] = 0.0;
		aij = ai;
		jx = kx;
		for (j = 0; j < lenx; j++) {
		  a_elem[0] = a_i[aij];
		  a_elem[1] = a_i[aij + 1];
		  a_elem[1] = -a_elem[1];
		  x_elem[0] = head_x_i[jx];
		  x_elem[1] = head_x_i[jx + 1];
		  {
		    prod[0] =
		      (double) a_elem[0] * x_elem[0] -
		      (double) a_elem[1] * x_elem[1];
		    prod[1] =
		      (double) a_elem[0] * x_elem[1] +
		      (double) a_elem[1] * x_elem[0];
		  }
		  sum[0] = sum[0] + prod[0];
		  sum[1] = sum[1] + prod[1];
		  x_elem[0] = tail_x_i[jx];
		  x_elem[1] = tail_x_i[jx + 1];
		  {
		    prod[0] =
		      (double) a_elem[0] * x_elem[0] -
		      (double) a_elem[1] * x_elem[1];
		    prod[1] =
		      (double) a_elem[0] * x_elem[1] +
		      (double) a_elem[1] * x_elem[0];
		  }
		  sum2[0] = sum2[0] + prod[0];
		  sum2[1] = sum2[1] + prod[1];
		  aij += incaij;
		  jx += incx;
		}
		{
		  tmp1[0] =
		    (double) sum[0] * alpha_i[0] -
		    (double) sum[1] * alpha_i[1];
		  tmp1[1] =
		    (double) sum[0] * alpha_i[1] +
		    (double) sum[1] * alpha_i[0];
		}
		{
		  tmp2[0] =
		    (double) sum2[0] * alpha_i[0] -
		    (double) sum2[1] * alpha_i[1];
		  tmp2[1] =
		    (double) sum2[0] * alpha_i[1] +
		    (double) sum2[1] * alpha_i[0];
		}
		tmp1[0] = tmp1[0] + tmp2[0];
		tmp1[1] = tmp1[1] + tmp2[1];
		y_i[iy] = tmp1[0];
		y_i[iy + 1] = tmp1[1];
		ai += incai;
		iy += incy;
	      }
	    }
	  } else {		/* beta != 0 */
	    if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {
	      /* save m multiplies if alpha = 1 */
	      ai = 0;
	      iy = ky;
	      for (i = 0; i < leny; i++) {
		sum[0] = sum[1] = 0.0;;
		sum2[0] = sum2[1] = 0.0;;
		aij = ai;
		jx = kx;
		for (j = 0; j < lenx; j++) {
		  a_elem[0] = a_i[aij];
		  a_elem[1] = a_i[aij + 1];
		  a_elem[1] = -a_elem[1];
		  x_elem[0] = head_x_i[jx];
		  x_elem[1] = head_x_i[jx + 1];
		  {
		    prod[0] =
		      (double) a_elem[0] * x_elem[0] -
		      (double) a_elem[1] * x_elem[1];
		    prod[1] =
		      (double) a_elem[0] * x_elem[1] +
		      (double) a_elem[1] * x_elem[0];
		  }
		  sum[0] = sum[0] + prod[0];
		  sum[1] = sum[1] + prod[1];
		  x_elem[0] = tail_x_i[jx];
		  x_elem[1] = tail_x_i[jx + 1];
		  {
		    prod[0] =
		      (double) a_elem[0] * x_elem[0] -
		      (double) a_elem[1] * x_elem[1];
		    prod[1] =
		      (double) a_elem[0] * x_elem[1] +
		      (double) a_elem[1] * x_elem[0];
		  }
		  sum2[0] = sum2[0] + prod[0];
		  sum2[1] = sum2[1] + prod[1];
		  aij += incaij;
		  jx += incx;
		}
		sum[0] = sum[0] + sum2[0];
		sum[1] = sum[1] + sum2[1];
		y_elem[0] = y_i[iy];
		y_elem[1] = y_i[iy + 1];
		{
		  tmp1[0] =
		    (double) y_elem[0] * beta_i[0] -
		    (double) y_elem[1] * beta_i[1];
		  tmp1[1] =
		    (double) y_elem[0] * beta_i[1] +
		    (double) y_elem[1] * beta_i[0];
		}
		tmp2[0] = sum[0] + tmp1[0];
		tmp2[1] = sum[1] + tmp1[1];
		y_i[iy] = tmp2[0];
		y_i[iy + 1] = tmp2[1];
		ai += incai;
		iy += incy;
	      }
	    } else {		/* alpha != 1, the most general form:
				   y = alpha*A*head_x + alpha*A*tail_x + beta*y */
	      ai = 0;
	      iy = ky;
	      for (i = 0; i < leny; i++) {
		sum[0] = sum[1] = 0.0;;
		sum2[0] = sum2[1] = 0.0;;
		aij = ai;
		jx = kx;
		for (j = 0; j < lenx; j++) {
		  a_elem[0] = a_i[aij];
		  a_elem[1] = a_i[aij + 1];
		  a_elem[1] = -a_elem[1];
		  x_elem[0] = head_x_i[jx];
		  x_elem[1] = head_x_i[jx + 1];
		  {
		    prod[0] =
		      (double) a_elem[0] * x_elem[0] -
		      (double) a_elem[1] * x_elem[1];
		    prod[1] =
		      (double) a_elem[0] * x_elem[1] +
		      (double) a_elem[1] * x_elem[0];
		  }
		  sum[0] = sum[0] + prod[0];
		  sum[1] = sum[1] + prod[1];
		  x_elem[0] = tail_x_i[jx];
		  x_elem[1] = tail_x_i[jx + 1];
		  {
		    prod[0] =
		      (double) a_elem[0] * x_elem[0] -
		      (double) a_elem[1] * x_elem[1];
		    prod[1] =
		      (double) a_elem[0] * x_elem[1] +
		      (double) a_elem[1] * x_elem[0];
		  }
		  sum2[0] = sum2[0] + prod[0];
		  sum2[1] = sum2[1] + prod[1];
		  aij += incaij;
		  jx += incx;
		}
		{
		  tmp1[0] =
		    (double) sum[0] * alpha_i[0] -
		    (double) sum[1] * alpha_i[1];
		  tmp1[1] =
		    (double) sum[0] * alpha_i[1] +
		    (double) sum[1] * alpha_i[0];
		}
		{
		  tmp2[0] =
		    (double) sum2[0] * alpha_i[0] -
		    (double) sum2[1] * alpha_i[1];
		  tmp2[1] =
		    (double) sum2[0] * alpha_i[1] +
		    (double) sum2[1] * alpha_i[0];
		}
		tmp1[0] = tmp1[0] + tmp2[0];
		tmp1[1] = tmp1[1] + tmp2[1];
		y_elem[0] = y_i[iy];
		y_elem[1] = y_i[iy + 1];
		{
		  tmp2[0] =
		    (double) y_elem[0] * beta_i[0] -
		    (double) y_elem[1] * beta_i[1];
		  tmp2[1] =
		    (double) y_elem[0] * beta_i[1] +
		    (double) y_elem[1] * beta_i[0];
		}
		tmp1[0] = tmp1[0] + tmp2[0];
		tmp1[1] = tmp1[1] + tmp2[1];
		y_i[iy] = tmp1[0];
		y_i[iy + 1] = tmp1[1];
		ai += incai;
		iy += incy;
	      }
	    }
	  }

	} else {

	  /* if beta = 0, we can save m multiplies:
	     y = alpha*A*head_x + alpha*A*tail_x  */
	  if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	    if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {
	      /* save m more multiplies if alpha = 1 */
	      ai = 0;
	      iy = ky;
	      for (i = 0; i < leny; i++) {
		sum[0] = sum[1] = 0.0;
		sum2[0] = sum2[1] = 0.0;
		aij = ai;
		jx = kx;
		for (j = 0; j < lenx; j++) {
		  a_elem[0] = a_i[aij];
		  a_elem[1] = a_i[aij + 1];

		  x_elem[0] = head_x_i[jx];
		  x_elem[1] = head_x_i[jx + 1];
		  {
		    prod[0] =
		      (double) a_elem[0] * x_elem[0] -
		      (double) a_elem[1] * x_elem[1];
		    prod[1] =
		      (double) a_elem[0] * x_elem[1] +
		      (double) a_elem[1] * x_elem[0];
		  }
		  sum[0] = sum[0] + prod[0];
		  sum[1] = sum[1] + prod[1];
		  x_elem[0] = tail_x_i[jx];
		  x_elem[1] = tail_x_i[jx + 1];
		  {
		    prod[0] =
		      (double) a_elem[0] * x_elem[0] -
		      (double) a_elem[1] * x_elem[1];
		    prod[1] =
		      (double) a_elem[0] * x_elem[1] +
		      (double) a_elem[1] * x_elem[0];
		  }
		  sum2[0] = sum2[0] + prod[0];
		  sum2[1] = sum2[1] + prod[1];
		  aij += incaij;
		  jx += incx;
		}
		sum[0] = sum[0] + sum2[0];
		sum[1] = sum[1] + sum2[1];
		y_i[iy] = sum[0];
		y_i[iy + 1] = sum[1];
		ai += incai;
		iy += incy;
	      }			/* end for */
	    } else {		/* alpha != 1 */
	      ai = 0;
	      iy = ky;
	      for (i = 0; i < leny; i++) {
		sum[0] = sum[1] = 0.0;
		sum2[0] = sum2[1] = 0.0;
		aij = ai;
		jx = kx;
		for (j = 0; j < lenx; j++) {
		  a_elem[0] = a_i[aij];
		  a_elem[1] = a_i[aij + 1];

		  x_elem[0] = head_x_i[jx];
		  x_elem[1] = head_x_i[jx + 1];
		  {
		    prod[0] =
		      (double) a_elem[0] * x_elem[0] -
		      (double) a_elem[1] * x_elem[1];
		    prod[1] =
		      (double) a_elem[0] * x_elem[1] +
		      (double) a_elem[1] * x_elem[0];
		  }
		  sum[0] = sum[0] + prod[0];
		  sum[1] = sum[1] + prod[1];
		  x_elem[0] = tail_x_i[jx];
		  x_elem[1] = tail_x_i[jx + 1];
		  {
		    prod[0] =
		      (double) a_elem[0] * x_elem[0] -
		      (double) a_elem[1] * x_elem[1];
		    prod[1] =
		      (double) a_elem[0] * x_elem[1] +
		      (double) a_elem[1] * x_elem[0];
		  }
		  sum2[0] = sum2[0] + prod[0];
		  sum2[1] = sum2[1] + prod[1];
		  aij += incaij;
		  jx += incx;
		}
		{
		  tmp1[0] =
		    (double) sum[0] * alpha_i[0] -
		    (double) sum[1] * alpha_i[1];
		  tmp1[1] =
		    (double) sum[0] * alpha_i[1] +
		    (double) sum[1] * alpha_i[0];
		}
		{
		  tmp2[0] =
		    (double) sum2[0] * alpha_i[0] -
		    (double) sum2[1] * alpha_i[1];
		  tmp2[1] =
		    (double) sum2[0] * alpha_i[1] +
		    (double) sum2[1] * alpha_i[0];
		}
		tmp1[0] = tmp1[0] + tmp2[0];
		tmp1[1] = tmp1[1] + tmp2[1];
		y_i[iy] = tmp1[0];
		y_i[iy + 1] = tmp1[1];
		ai += incai;
		iy += incy;
	      }
	    }
	  } else {		/* beta != 0 */
	    if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {
	      /* save m multiplies if alpha = 1 */
	      ai = 0;
	      iy = ky;
	      for (i = 0; i < leny; i++) {
		sum[0] = sum[1] = 0.0;;
		sum2[0] = sum2[1] = 0.0;;
		aij = ai;
		jx = kx;
		for (j = 0; j < lenx; j++) {
		  a_elem[0] = a_i[aij];
		  a_elem[1] = a_i[aij + 1];

		  x_elem[0] = head_x_i[jx];
		  x_elem[1] = head_x_i[jx + 1];
		  {
		    prod[0] =
		      (double) a_elem[0] * x_elem[0] -
		      (double) a_elem[1] * x_elem[1];
		    prod[1] =
		      (double) a_elem[0] * x_elem[1] +
		      (double) a_elem[1] * x_elem[0];
		  }
		  sum[0] = sum[0] + prod[0];
		  sum[1] = sum[1] + prod[1];
		  x_elem[0] = tail_x_i[jx];
		  x_elem[1] = tail_x_i[jx + 1];
		  {
		    prod[0] =
		      (double) a_elem[0] * x_elem[0] -
		      (double) a_elem[1] * x_elem[1];
		    prod[1] =
		      (double) a_elem[0] * x_elem[1] +
		      (double) a_elem[1] * x_elem[0];
		  }
		  sum2[0] = sum2[0] + prod[0];
		  sum2[1] = sum2[1] + prod[1];
		  aij += incaij;
		  jx += incx;
		}
		sum[0] = sum[0] + sum2[0];
		sum[1] = sum[1] + sum2[1];
		y_elem[0] = y_i[iy];
		y_elem[1] = y_i[iy + 1];
		{
		  tmp1[0] =
		    (double) y_elem[0] * beta_i[0] -
		    (double) y_elem[1] * beta_i[1];
		  tmp1[1] =
		    (double) y_elem[0] * beta_i[1] +
		    (double) y_elem[1] * beta_i[0];
		}
		tmp2[0] = sum[0] + tmp1[0];
		tmp2[1] = sum[1] + tmp1[1];
		y_i[iy] = tmp2[0];
		y_i[iy + 1] = tmp2[1];
		ai += incai;
		iy += incy;
	      }
	    } else {		/* alpha != 1, the most general form:
				   y = alpha*A*head_x + alpha*A*tail_x + beta*y */
	      ai = 0;
	      iy = ky;
	      for (i = 0; i < leny; i++) {
		sum[0] = sum[1] = 0.0;;
		sum2[0] = sum2[1] = 0.0;;
		aij = ai;
		jx = kx;
		for (j = 0; j < lenx; j++) {
		  a_elem[0] = a_i[aij];
		  a_elem[1] = a_i[aij + 1];

		  x_elem[0] = head_x_i[jx];
		  x_elem[1] = head_x_i[jx + 1];
		  {
		    prod[0] =
		      (double) a_elem[0] * x_elem[0] -
		      (double) a_elem[1] * x_elem[1];
		    prod[1] =
		      (double) a_elem[0] * x_elem[1] +
		      (double) a_elem[1] * x_elem[0];
		  }
		  sum[0] = sum[0] + prod[0];
		  sum[1] = sum[1] + prod[1];
		  x_elem[0] = tail_x_i[jx];
		  x_elem[1] = tail_x_i[jx + 1];
		  {
		    prod[0] =
		      (double) a_elem[0] * x_elem[0] -
		      (double) a_elem[1] * x_elem[1];
		    prod[1] =
		      (double) a_elem[0] * x_elem[1] +
		      (double) a_elem[1] * x_elem[0];
		  }
		  sum2[0] = sum2[0] + prod[0];
		  sum2[1] = sum2[1] + prod[1];
		  aij += incaij;
		  jx += incx;
		}
		{
		  tmp1[0] =
		    (double) sum[0] * alpha_i[0] -
		    (double) sum[1] * alpha_i[1];
		  tmp1[1] =
		    (double) sum[0] * alpha_i[1] +
		    (double) sum[1] * alpha_i[0];
		}
		{
		  tmp2[0] =
		    (double) sum2[0] * alpha_i[0] -
		    (double) sum2[1] * alpha_i[1];
		  tmp2[1] =
		    (double) sum2[0] * alpha_i[1] +
		    (double) sum2[1] * alpha_i[0];
		}
		tmp1[0] = tmp1[0] + tmp2[0];
		tmp1[1] = tmp1[1] + tmp2[1];
		y_elem[0] = y_i[iy];
		y_elem[1] = y_i[iy + 1];
		{
		  tmp2[0] =
		    (double) y_elem[0] * beta_i[0] -
		    (double) y_elem[1] * beta_i[1];
		  tmp2[1] =
		    (double) y_elem[0] * beta_i[1] +
		    (double) y_elem[1] * beta_i[0];
		}
		tmp1[0] = tmp1[0] + tmp2[0];
		tmp1[1] = tmp1[1] + tmp2[1];
		y_i[iy] = tmp1[0];
		y_i[iy + 1] = tmp1[1];
		ai += incai;
		iy += incy;
	      }
	    }
	  }

	}
      }



      break;
    }
  case blas_prec_extra:{

      int i, j;
      int iy, jx, kx, ky;
      int lenx, leny;
      int ai, aij;
      int incai, incaij;

      const double *a_i = (double *) a;
      const float *head_x_i = (float *) head_x;
      const float *tail_x_i = (float *) tail_x;
      double *y_i = (double *) y;
      double *alpha_i = (double *) alpha;
      double *beta_i = (double *) beta;
      double a_elem[2];
      float x_elem[2];
      double y_elem[2];
      double head_prod[2], tail_prod[2];
      double head_sum[2], tail_sum[2];
      double head_sum2[2], tail_sum2[2];
      double head_tmp1[2], tail_tmp1[2];
      double head_tmp2[2], tail_tmp2[2];
      FPU_FIX_DECL;

      /* all error calls */
      if (m < 0)
	BLAS_error(routine_name, -3, m, 0);
      else if (n <= 0)
	BLAS_error(routine_name, -4, n, 0);
      else if (incx == 0)
	BLAS_error(routine_name, -10, incx, 0);
      else if (incy == 0)
	BLAS_error(routine_name, -13, incy, 0);

      if ((order == blas_rowmajor) && (trans == blas_no_trans)) {
	lenx = n;
	leny = m;
	incai = lda;
	incaij = 1;
      } else if ((order == blas_rowmajor) && (trans != blas_no_trans)) {
	lenx = m;
	leny = n;
	incai = 1;
	incaij = lda;
      } else if ((order == blas_colmajor) && (trans == blas_no_trans)) {
	lenx = n;
	leny = m;
	incai = 1;
	incaij = lda;
      } else {			/* colmajor and blas_trans */
	lenx = m;
	leny = n;
	incai = lda;
	incaij = 1;
      }

      if (lda < leny)
	BLAS_error(routine_name, -7, lda, NULL);

      FPU_FIX_START;

      incx *= 2;
      incy *= 2;
      incai *= 2;
      incaij *= 2;

      if (incx > 0)
	kx = 0;
      else
	kx = (1 - lenx) * incx;
      if (incy > 0)
	ky = 0;
      else
	ky = (1 - leny) * incy;

      /* No extra-precision needed for alpha = 0 */
      if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
	if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	  iy = ky;
	  for (i = 0; i < leny; i++) {
	    y_i[iy] = 0.0;
	    y_i[iy + 1] = 0.0;
	    iy += incy;
	  }
	} else if (!(beta_i[0] == 0.0 && beta_i[1] == 0.0)) {
	  iy = ky;
	  for (i = 0; i < leny; i++) {
	    y_elem[0] = y_i[iy];
	    y_elem[1] = y_i[iy + 1];
	    {
	      /* Compute complex-extra = complex-double * complex-double. */
	      double head_t1, tail_t1;
	      double head_t2, tail_t2;
	      /* Real part */
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = y_elem[0] * split;
		a1 = con - y_elem[0];
		a1 = con - a1;
		a2 = y_elem[0] - a1;
		con = beta_i[0] * split;
		b1 = con - beta_i[0];
		b1 = con - b1;
		b2 = beta_i[0] - b1;

		head_t1 = y_elem[0] * beta_i[0];
		tail_t1 =
		  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = y_elem[1] * split;
		a1 = con - y_elem[1];
		a1 = con - a1;
		a2 = y_elem[1] - a1;
		con = beta_i[1] * split;
		b1 = con - beta_i[1];
		b1 = con - b1;
		b2 = beta_i[1] - b1;

		head_t2 = y_elem[1] * beta_i[1];
		tail_t2 =
		  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
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
	      /* Imaginary part */
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = y_elem[1] * split;
		a1 = con - y_elem[1];
		a1 = con - a1;
		a2 = y_elem[1] - a1;
		con = beta_i[0] * split;
		b1 = con - beta_i[0];
		b1 = con - b1;
		b2 = beta_i[0] - b1;

		head_t1 = y_elem[1] * beta_i[0];
		tail_t1 =
		  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = y_elem[0] * split;
		a1 = con - y_elem[0];
		a1 = con - a1;
		a2 = y_elem[0] - a1;
		con = beta_i[1] * split;
		b1 = con - beta_i[1];
		b1 = con - b1;
		b2 = beta_i[1] - b1;

		head_t2 = y_elem[0] * beta_i[1];
		tail_t2 =
		  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
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
	    y_i[iy] = head_tmp1[0];
	    y_i[iy + 1] = head_tmp1[1];
	    iy += incy;
	  }
	}
      } else {			/* alpha != 0 */
	if (trans == blas_conj_trans) {

	  /* if beta = 0, we can save m multiplies:
	     y = alpha*A*head_x + alpha*A*tail_x  */
	  if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	    if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {
	      /* save m more multiplies if alpha = 1 */
	      ai = 0;
	      iy = ky;
	      for (i = 0; i < leny; i++) {
		head_sum[0] = head_sum[1] = tail_sum[0] = tail_sum[1] = 0.0;
		head_sum2[0] = head_sum2[1] = tail_sum2[0] = tail_sum2[1] =
		  0.0;
		aij = ai;
		jx = kx;
		for (j = 0; j < lenx; j++) {
		  a_elem[0] = a_i[aij];
		  a_elem[1] = a_i[aij + 1];
		  a_elem[1] = -a_elem[1];
		  x_elem[0] = head_x_i[jx];
		  x_elem[1] = head_x_i[jx + 1];
		  {
		    double cd[2];
		    cd[0] = (double) x_elem[0];
		    cd[1] = (double) x_elem[1];
		    {
		      /* Compute complex-extra = complex-double * complex-double. */
		      double head_t1, tail_t1;
		      double head_t2, tail_t2;
		      /* Real part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[0] * split;
			a1 = con - a_elem[0];
			a1 = con - a1;
			a2 = a_elem[0] - a1;
			con = cd[0] * split;
			b1 = con - cd[0];
			b1 = con - b1;
			b2 = cd[0] - b1;

			head_t1 = a_elem[0] * cd[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[1] * split;
			a1 = con - a_elem[1];
			a1 = con - a1;
			a2 = a_elem[1] - a1;
			con = cd[1] * split;
			b1 = con - cd[1];
			b1 = con - b1;
			b2 = cd[1] - b1;

			head_t2 = a_elem[1] * cd[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_prod[0] = head_t1;
		      tail_prod[0] = tail_t1;
		      /* Imaginary part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[1] * split;
			a1 = con - a_elem[1];
			a1 = con - a1;
			a2 = a_elem[1] - a1;
			con = cd[0] * split;
			b1 = con - cd[0];
			b1 = con - b1;
			b2 = cd[0] - b1;

			head_t1 = a_elem[1] * cd[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[0] * split;
			a1 = con - a_elem[0];
			a1 = con - a1;
			a2 = a_elem[0] - a1;
			con = cd[1] * split;
			b1 = con - cd[1];
			b1 = con - b1;
			b2 = cd[1] - b1;

			head_t2 = a_elem[0] * cd[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_prod[1] = head_t1;
		      tail_prod[1] = tail_t1;
		    }
		  }
		  {
		    double head_t, tail_t;
		    double head_a, tail_a;
		    double head_b, tail_b;
		    /* Real part */
		    head_a = head_sum[0];
		    tail_a = tail_sum[0];
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
		    head_sum[0] = head_t;
		    tail_sum[0] = tail_t;
		    /* Imaginary part */
		    head_a = head_sum[1];
		    tail_a = tail_sum[1];
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
		    head_sum[1] = head_t;
		    tail_sum[1] = tail_t;
		  }
		  x_elem[0] = tail_x_i[jx];
		  x_elem[1] = tail_x_i[jx + 1];
		  {
		    double cd[2];
		    cd[0] = (double) x_elem[0];
		    cd[1] = (double) x_elem[1];
		    {
		      /* Compute complex-extra = complex-double * complex-double. */
		      double head_t1, tail_t1;
		      double head_t2, tail_t2;
		      /* Real part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[0] * split;
			a1 = con - a_elem[0];
			a1 = con - a1;
			a2 = a_elem[0] - a1;
			con = cd[0] * split;
			b1 = con - cd[0];
			b1 = con - b1;
			b2 = cd[0] - b1;

			head_t1 = a_elem[0] * cd[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[1] * split;
			a1 = con - a_elem[1];
			a1 = con - a1;
			a2 = a_elem[1] - a1;
			con = cd[1] * split;
			b1 = con - cd[1];
			b1 = con - b1;
			b2 = cd[1] - b1;

			head_t2 = a_elem[1] * cd[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_prod[0] = head_t1;
		      tail_prod[0] = tail_t1;
		      /* Imaginary part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[1] * split;
			a1 = con - a_elem[1];
			a1 = con - a1;
			a2 = a_elem[1] - a1;
			con = cd[0] * split;
			b1 = con - cd[0];
			b1 = con - b1;
			b2 = cd[0] - b1;

			head_t1 = a_elem[1] * cd[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[0] * split;
			a1 = con - a_elem[0];
			a1 = con - a1;
			a2 = a_elem[0] - a1;
			con = cd[1] * split;
			b1 = con - cd[1];
			b1 = con - b1;
			b2 = cd[1] - b1;

			head_t2 = a_elem[0] * cd[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_prod[1] = head_t1;
		      tail_prod[1] = tail_t1;
		    }
		  }
		  {
		    double head_t, tail_t;
		    double head_a, tail_a;
		    double head_b, tail_b;
		    /* Real part */
		    head_a = head_sum2[0];
		    tail_a = tail_sum2[0];
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
		    head_sum2[0] = head_t;
		    tail_sum2[0] = tail_t;
		    /* Imaginary part */
		    head_a = head_sum2[1];
		    tail_a = tail_sum2[1];
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
		    head_sum2[1] = head_t;
		    tail_sum2[1] = tail_t;
		  }
		  aij += incaij;
		  jx += incx;
		}
		{
		  double head_t, tail_t;
		  double head_a, tail_a;
		  double head_b, tail_b;
		  /* Real part */
		  head_a = head_sum[0];
		  tail_a = tail_sum[0];
		  head_b = head_sum2[0];
		  tail_b = tail_sum2[0];
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
		  head_sum[0] = head_t;
		  tail_sum[0] = tail_t;
		  /* Imaginary part */
		  head_a = head_sum[1];
		  tail_a = tail_sum[1];
		  head_b = head_sum2[1];
		  tail_b = tail_sum2[1];
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
		  head_sum[1] = head_t;
		  tail_sum[1] = tail_t;
		}
		y_i[iy] = head_sum[0];
		y_i[iy + 1] = head_sum[1];
		ai += incai;
		iy += incy;
	      }			/* end for */
	    } else {		/* alpha != 1 */
	      ai = 0;
	      iy = ky;
	      for (i = 0; i < leny; i++) {
		head_sum[0] = head_sum[1] = tail_sum[0] = tail_sum[1] = 0.0;
		head_sum2[0] = head_sum2[1] = tail_sum2[0] = tail_sum2[1] =
		  0.0;
		aij = ai;
		jx = kx;
		for (j = 0; j < lenx; j++) {
		  a_elem[0] = a_i[aij];
		  a_elem[1] = a_i[aij + 1];
		  a_elem[1] = -a_elem[1];
		  x_elem[0] = head_x_i[jx];
		  x_elem[1] = head_x_i[jx + 1];
		  {
		    double cd[2];
		    cd[0] = (double) x_elem[0];
		    cd[1] = (double) x_elem[1];
		    {
		      /* Compute complex-extra = complex-double * complex-double. */
		      double head_t1, tail_t1;
		      double head_t2, tail_t2;
		      /* Real part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[0] * split;
			a1 = con - a_elem[0];
			a1 = con - a1;
			a2 = a_elem[0] - a1;
			con = cd[0] * split;
			b1 = con - cd[0];
			b1 = con - b1;
			b2 = cd[0] - b1;

			head_t1 = a_elem[0] * cd[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[1] * split;
			a1 = con - a_elem[1];
			a1 = con - a1;
			a2 = a_elem[1] - a1;
			con = cd[1] * split;
			b1 = con - cd[1];
			b1 = con - b1;
			b2 = cd[1] - b1;

			head_t2 = a_elem[1] * cd[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_prod[0] = head_t1;
		      tail_prod[0] = tail_t1;
		      /* Imaginary part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[1] * split;
			a1 = con - a_elem[1];
			a1 = con - a1;
			a2 = a_elem[1] - a1;
			con = cd[0] * split;
			b1 = con - cd[0];
			b1 = con - b1;
			b2 = cd[0] - b1;

			head_t1 = a_elem[1] * cd[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[0] * split;
			a1 = con - a_elem[0];
			a1 = con - a1;
			a2 = a_elem[0] - a1;
			con = cd[1] * split;
			b1 = con - cd[1];
			b1 = con - b1;
			b2 = cd[1] - b1;

			head_t2 = a_elem[0] * cd[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_prod[1] = head_t1;
		      tail_prod[1] = tail_t1;
		    }
		  }
		  {
		    double head_t, tail_t;
		    double head_a, tail_a;
		    double head_b, tail_b;
		    /* Real part */
		    head_a = head_sum[0];
		    tail_a = tail_sum[0];
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
		    head_sum[0] = head_t;
		    tail_sum[0] = tail_t;
		    /* Imaginary part */
		    head_a = head_sum[1];
		    tail_a = tail_sum[1];
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
		    head_sum[1] = head_t;
		    tail_sum[1] = tail_t;
		  }
		  x_elem[0] = tail_x_i[jx];
		  x_elem[1] = tail_x_i[jx + 1];
		  {
		    double cd[2];
		    cd[0] = (double) x_elem[0];
		    cd[1] = (double) x_elem[1];
		    {
		      /* Compute complex-extra = complex-double * complex-double. */
		      double head_t1, tail_t1;
		      double head_t2, tail_t2;
		      /* Real part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[0] * split;
			a1 = con - a_elem[0];
			a1 = con - a1;
			a2 = a_elem[0] - a1;
			con = cd[0] * split;
			b1 = con - cd[0];
			b1 = con - b1;
			b2 = cd[0] - b1;

			head_t1 = a_elem[0] * cd[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[1] * split;
			a1 = con - a_elem[1];
			a1 = con - a1;
			a2 = a_elem[1] - a1;
			con = cd[1] * split;
			b1 = con - cd[1];
			b1 = con - b1;
			b2 = cd[1] - b1;

			head_t2 = a_elem[1] * cd[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_prod[0] = head_t1;
		      tail_prod[0] = tail_t1;
		      /* Imaginary part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[1] * split;
			a1 = con - a_elem[1];
			a1 = con - a1;
			a2 = a_elem[1] - a1;
			con = cd[0] * split;
			b1 = con - cd[0];
			b1 = con - b1;
			b2 = cd[0] - b1;

			head_t1 = a_elem[1] * cd[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[0] * split;
			a1 = con - a_elem[0];
			a1 = con - a1;
			a2 = a_elem[0] - a1;
			con = cd[1] * split;
			b1 = con - cd[1];
			b1 = con - b1;
			b2 = cd[1] - b1;

			head_t2 = a_elem[0] * cd[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_prod[1] = head_t1;
		      tail_prod[1] = tail_t1;
		    }
		  }
		  {
		    double head_t, tail_t;
		    double head_a, tail_a;
		    double head_b, tail_b;
		    /* Real part */
		    head_a = head_sum2[0];
		    tail_a = tail_sum2[0];
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
		    head_sum2[0] = head_t;
		    tail_sum2[0] = tail_t;
		    /* Imaginary part */
		    head_a = head_sum2[1];
		    tail_a = tail_sum2[1];
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
		    head_sum2[1] = head_t;
		    tail_sum2[1] = tail_t;
		  }
		  aij += incaij;
		  jx += incx;
		}
		{
		  /* Compute complex-extra = complex-extra * complex-double. */
		  double head_a0, tail_a0;
		  double head_a1, tail_a1;
		  double head_t1, tail_t1;
		  double head_t2, tail_t2;
		  head_a0 = head_sum[0];
		  tail_a0 = tail_sum[0];
		  head_a1 = head_sum[1];
		  tail_a1 = tail_sum[1];
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
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

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
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

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
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

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
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

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
		  head_a0 = head_sum2[0];
		  tail_a0 = tail_sum2[0];
		  head_a1 = head_sum2[1];
		  tail_a1 = tail_sum2[1];
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
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

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
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

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
		    con = alpha_i[0] * split;
		    b1 = con - alpha_i[0];
		    b1 = con - b1;
		    b2 = alpha_i[0] - b1;

		    c11 = head_a1 * alpha_i[0];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

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
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

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
		y_i[iy] = head_tmp1[0];
		y_i[iy + 1] = head_tmp1[1];
		ai += incai;
		iy += incy;
	      }
	    }
	  } else {		/* beta != 0 */
	    if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {
	      /* save m multiplies if alpha = 1 */
	      ai = 0;
	      iy = ky;
	      for (i = 0; i < leny; i++) {
		head_sum[0] = head_sum[1] = tail_sum[0] = tail_sum[1] = 0.0;;
		head_sum2[0] = head_sum2[1] = tail_sum2[0] = tail_sum2[1] =
		  0.0;;
		aij = ai;
		jx = kx;
		for (j = 0; j < lenx; j++) {
		  a_elem[0] = a_i[aij];
		  a_elem[1] = a_i[aij + 1];
		  a_elem[1] = -a_elem[1];
		  x_elem[0] = head_x_i[jx];
		  x_elem[1] = head_x_i[jx + 1];
		  {
		    double cd[2];
		    cd[0] = (double) x_elem[0];
		    cd[1] = (double) x_elem[1];
		    {
		      /* Compute complex-extra = complex-double * complex-double. */
		      double head_t1, tail_t1;
		      double head_t2, tail_t2;
		      /* Real part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[0] * split;
			a1 = con - a_elem[0];
			a1 = con - a1;
			a2 = a_elem[0] - a1;
			con = cd[0] * split;
			b1 = con - cd[0];
			b1 = con - b1;
			b2 = cd[0] - b1;

			head_t1 = a_elem[0] * cd[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[1] * split;
			a1 = con - a_elem[1];
			a1 = con - a1;
			a2 = a_elem[1] - a1;
			con = cd[1] * split;
			b1 = con - cd[1];
			b1 = con - b1;
			b2 = cd[1] - b1;

			head_t2 = a_elem[1] * cd[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_prod[0] = head_t1;
		      tail_prod[0] = tail_t1;
		      /* Imaginary part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[1] * split;
			a1 = con - a_elem[1];
			a1 = con - a1;
			a2 = a_elem[1] - a1;
			con = cd[0] * split;
			b1 = con - cd[0];
			b1 = con - b1;
			b2 = cd[0] - b1;

			head_t1 = a_elem[1] * cd[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[0] * split;
			a1 = con - a_elem[0];
			a1 = con - a1;
			a2 = a_elem[0] - a1;
			con = cd[1] * split;
			b1 = con - cd[1];
			b1 = con - b1;
			b2 = cd[1] - b1;

			head_t2 = a_elem[0] * cd[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_prod[1] = head_t1;
		      tail_prod[1] = tail_t1;
		    }
		  }
		  {
		    double head_t, tail_t;
		    double head_a, tail_a;
		    double head_b, tail_b;
		    /* Real part */
		    head_a = head_sum[0];
		    tail_a = tail_sum[0];
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
		    head_sum[0] = head_t;
		    tail_sum[0] = tail_t;
		    /* Imaginary part */
		    head_a = head_sum[1];
		    tail_a = tail_sum[1];
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
		    head_sum[1] = head_t;
		    tail_sum[1] = tail_t;
		  }
		  x_elem[0] = tail_x_i[jx];
		  x_elem[1] = tail_x_i[jx + 1];
		  {
		    double cd[2];
		    cd[0] = (double) x_elem[0];
		    cd[1] = (double) x_elem[1];
		    {
		      /* Compute complex-extra = complex-double * complex-double. */
		      double head_t1, tail_t1;
		      double head_t2, tail_t2;
		      /* Real part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[0] * split;
			a1 = con - a_elem[0];
			a1 = con - a1;
			a2 = a_elem[0] - a1;
			con = cd[0] * split;
			b1 = con - cd[0];
			b1 = con - b1;
			b2 = cd[0] - b1;

			head_t1 = a_elem[0] * cd[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[1] * split;
			a1 = con - a_elem[1];
			a1 = con - a1;
			a2 = a_elem[1] - a1;
			con = cd[1] * split;
			b1 = con - cd[1];
			b1 = con - b1;
			b2 = cd[1] - b1;

			head_t2 = a_elem[1] * cd[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_prod[0] = head_t1;
		      tail_prod[0] = tail_t1;
		      /* Imaginary part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[1] * split;
			a1 = con - a_elem[1];
			a1 = con - a1;
			a2 = a_elem[1] - a1;
			con = cd[0] * split;
			b1 = con - cd[0];
			b1 = con - b1;
			b2 = cd[0] - b1;

			head_t1 = a_elem[1] * cd[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[0] * split;
			a1 = con - a_elem[0];
			a1 = con - a1;
			a2 = a_elem[0] - a1;
			con = cd[1] * split;
			b1 = con - cd[1];
			b1 = con - b1;
			b2 = cd[1] - b1;

			head_t2 = a_elem[0] * cd[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_prod[1] = head_t1;
		      tail_prod[1] = tail_t1;
		    }
		  }
		  {
		    double head_t, tail_t;
		    double head_a, tail_a;
		    double head_b, tail_b;
		    /* Real part */
		    head_a = head_sum2[0];
		    tail_a = tail_sum2[0];
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
		    head_sum2[0] = head_t;
		    tail_sum2[0] = tail_t;
		    /* Imaginary part */
		    head_a = head_sum2[1];
		    tail_a = tail_sum2[1];
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
		    head_sum2[1] = head_t;
		    tail_sum2[1] = tail_t;
		  }
		  aij += incaij;
		  jx += incx;
		}
		{
		  double head_t, tail_t;
		  double head_a, tail_a;
		  double head_b, tail_b;
		  /* Real part */
		  head_a = head_sum[0];
		  tail_a = tail_sum[0];
		  head_b = head_sum2[0];
		  tail_b = tail_sum2[0];
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
		  head_sum[0] = head_t;
		  tail_sum[0] = tail_t;
		  /* Imaginary part */
		  head_a = head_sum[1];
		  tail_a = tail_sum[1];
		  head_b = head_sum2[1];
		  tail_b = tail_sum2[1];
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
		  head_sum[1] = head_t;
		  tail_sum[1] = tail_t;
		}
		y_elem[0] = y_i[iy];
		y_elem[1] = y_i[iy + 1];
		{
		  /* Compute complex-extra = complex-double * complex-double. */
		  double head_t1, tail_t1;
		  double head_t2, tail_t2;
		  /* Real part */
		  {
		    /* Compute double_double = double * double. */
		    double a1, a2, b1, b2, con;

		    con = y_elem[0] * split;
		    a1 = con - y_elem[0];
		    a1 = con - a1;
		    a2 = y_elem[0] - a1;
		    con = beta_i[0] * split;
		    b1 = con - beta_i[0];
		    b1 = con - b1;
		    b2 = beta_i[0] - b1;

		    head_t1 = y_elem[0] * beta_i[0];
		    tail_t1 =
		      (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
		  }
		  {
		    /* Compute double_double = double * double. */
		    double a1, a2, b1, b2, con;

		    con = y_elem[1] * split;
		    a1 = con - y_elem[1];
		    a1 = con - a1;
		    a2 = y_elem[1] - a1;
		    con = beta_i[1] * split;
		    b1 = con - beta_i[1];
		    b1 = con - b1;
		    b2 = beta_i[1] - b1;

		    head_t2 = y_elem[1] * beta_i[1];
		    tail_t2 =
		      (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
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
		  /* Imaginary part */
		  {
		    /* Compute double_double = double * double. */
		    double a1, a2, b1, b2, con;

		    con = y_elem[1] * split;
		    a1 = con - y_elem[1];
		    a1 = con - a1;
		    a2 = y_elem[1] - a1;
		    con = beta_i[0] * split;
		    b1 = con - beta_i[0];
		    b1 = con - b1;
		    b2 = beta_i[0] - b1;

		    head_t1 = y_elem[1] * beta_i[0];
		    tail_t1 =
		      (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
		  }
		  {
		    /* Compute double_double = double * double. */
		    double a1, a2, b1, b2, con;

		    con = y_elem[0] * split;
		    a1 = con - y_elem[0];
		    a1 = con - a1;
		    a2 = y_elem[0] - a1;
		    con = beta_i[1] * split;
		    b1 = con - beta_i[1];
		    b1 = con - b1;
		    b2 = beta_i[1] - b1;

		    head_t2 = y_elem[0] * beta_i[1];
		    tail_t2 =
		      (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
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
		  double head_t, tail_t;
		  double head_a, tail_a;
		  double head_b, tail_b;
		  /* Real part */
		  head_a = head_sum[0];
		  tail_a = tail_sum[0];
		  head_b = head_tmp1[0];
		  tail_b = tail_tmp1[0];
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
		  head_tmp2[0] = head_t;
		  tail_tmp2[0] = tail_t;
		  /* Imaginary part */
		  head_a = head_sum[1];
		  tail_a = tail_sum[1];
		  head_b = head_tmp1[1];
		  tail_b = tail_tmp1[1];
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
		  head_tmp2[1] = head_t;
		  tail_tmp2[1] = tail_t;
		}
		y_i[iy] = head_tmp2[0];
		y_i[iy + 1] = head_tmp2[1];
		ai += incai;
		iy += incy;
	      }
	    } else {		/* alpha != 1, the most general form:
				   y = alpha*A*head_x + alpha*A*tail_x + beta*y */
	      ai = 0;
	      iy = ky;
	      for (i = 0; i < leny; i++) {
		head_sum[0] = head_sum[1] = tail_sum[0] = tail_sum[1] = 0.0;;
		head_sum2[0] = head_sum2[1] = tail_sum2[0] = tail_sum2[1] =
		  0.0;;
		aij = ai;
		jx = kx;
		for (j = 0; j < lenx; j++) {
		  a_elem[0] = a_i[aij];
		  a_elem[1] = a_i[aij + 1];
		  a_elem[1] = -a_elem[1];
		  x_elem[0] = head_x_i[jx];
		  x_elem[1] = head_x_i[jx + 1];
		  {
		    double cd[2];
		    cd[0] = (double) x_elem[0];
		    cd[1] = (double) x_elem[1];
		    {
		      /* Compute complex-extra = complex-double * complex-double. */
		      double head_t1, tail_t1;
		      double head_t2, tail_t2;
		      /* Real part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[0] * split;
			a1 = con - a_elem[0];
			a1 = con - a1;
			a2 = a_elem[0] - a1;
			con = cd[0] * split;
			b1 = con - cd[0];
			b1 = con - b1;
			b2 = cd[0] - b1;

			head_t1 = a_elem[0] * cd[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[1] * split;
			a1 = con - a_elem[1];
			a1 = con - a1;
			a2 = a_elem[1] - a1;
			con = cd[1] * split;
			b1 = con - cd[1];
			b1 = con - b1;
			b2 = cd[1] - b1;

			head_t2 = a_elem[1] * cd[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_prod[0] = head_t1;
		      tail_prod[0] = tail_t1;
		      /* Imaginary part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[1] * split;
			a1 = con - a_elem[1];
			a1 = con - a1;
			a2 = a_elem[1] - a1;
			con = cd[0] * split;
			b1 = con - cd[0];
			b1 = con - b1;
			b2 = cd[0] - b1;

			head_t1 = a_elem[1] * cd[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[0] * split;
			a1 = con - a_elem[0];
			a1 = con - a1;
			a2 = a_elem[0] - a1;
			con = cd[1] * split;
			b1 = con - cd[1];
			b1 = con - b1;
			b2 = cd[1] - b1;

			head_t2 = a_elem[0] * cd[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_prod[1] = head_t1;
		      tail_prod[1] = tail_t1;
		    }
		  }
		  {
		    double head_t, tail_t;
		    double head_a, tail_a;
		    double head_b, tail_b;
		    /* Real part */
		    head_a = head_sum[0];
		    tail_a = tail_sum[0];
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
		    head_sum[0] = head_t;
		    tail_sum[0] = tail_t;
		    /* Imaginary part */
		    head_a = head_sum[1];
		    tail_a = tail_sum[1];
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
		    head_sum[1] = head_t;
		    tail_sum[1] = tail_t;
		  }
		  x_elem[0] = tail_x_i[jx];
		  x_elem[1] = tail_x_i[jx + 1];
		  {
		    double cd[2];
		    cd[0] = (double) x_elem[0];
		    cd[1] = (double) x_elem[1];
		    {
		      /* Compute complex-extra = complex-double * complex-double. */
		      double head_t1, tail_t1;
		      double head_t2, tail_t2;
		      /* Real part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[0] * split;
			a1 = con - a_elem[0];
			a1 = con - a1;
			a2 = a_elem[0] - a1;
			con = cd[0] * split;
			b1 = con - cd[0];
			b1 = con - b1;
			b2 = cd[0] - b1;

			head_t1 = a_elem[0] * cd[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[1] * split;
			a1 = con - a_elem[1];
			a1 = con - a1;
			a2 = a_elem[1] - a1;
			con = cd[1] * split;
			b1 = con - cd[1];
			b1 = con - b1;
			b2 = cd[1] - b1;

			head_t2 = a_elem[1] * cd[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_prod[0] = head_t1;
		      tail_prod[0] = tail_t1;
		      /* Imaginary part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[1] * split;
			a1 = con - a_elem[1];
			a1 = con - a1;
			a2 = a_elem[1] - a1;
			con = cd[0] * split;
			b1 = con - cd[0];
			b1 = con - b1;
			b2 = cd[0] - b1;

			head_t1 = a_elem[1] * cd[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[0] * split;
			a1 = con - a_elem[0];
			a1 = con - a1;
			a2 = a_elem[0] - a1;
			con = cd[1] * split;
			b1 = con - cd[1];
			b1 = con - b1;
			b2 = cd[1] - b1;

			head_t2 = a_elem[0] * cd[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_prod[1] = head_t1;
		      tail_prod[1] = tail_t1;
		    }
		  }
		  {
		    double head_t, tail_t;
		    double head_a, tail_a;
		    double head_b, tail_b;
		    /* Real part */
		    head_a = head_sum2[0];
		    tail_a = tail_sum2[0];
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
		    head_sum2[0] = head_t;
		    tail_sum2[0] = tail_t;
		    /* Imaginary part */
		    head_a = head_sum2[1];
		    tail_a = tail_sum2[1];
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
		    head_sum2[1] = head_t;
		    tail_sum2[1] = tail_t;
		  }
		  aij += incaij;
		  jx += incx;
		}
		{
		  /* Compute complex-extra = complex-extra * complex-double. */
		  double head_a0, tail_a0;
		  double head_a1, tail_a1;
		  double head_t1, tail_t1;
		  double head_t2, tail_t2;
		  head_a0 = head_sum[0];
		  tail_a0 = tail_sum[0];
		  head_a1 = head_sum[1];
		  tail_a1 = tail_sum[1];
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
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

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
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

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
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

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
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

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
		  head_a0 = head_sum2[0];
		  tail_a0 = tail_sum2[0];
		  head_a1 = head_sum2[1];
		  tail_a1 = tail_sum2[1];
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
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

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
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

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
		    con = alpha_i[0] * split;
		    b1 = con - alpha_i[0];
		    b1 = con - b1;
		    b2 = alpha_i[0] - b1;

		    c11 = head_a1 * alpha_i[0];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

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
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

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
		y_elem[0] = y_i[iy];
		y_elem[1] = y_i[iy + 1];
		{
		  /* Compute complex-extra = complex-double * complex-double. */
		  double head_t1, tail_t1;
		  double head_t2, tail_t2;
		  /* Real part */
		  {
		    /* Compute double_double = double * double. */
		    double a1, a2, b1, b2, con;

		    con = y_elem[0] * split;
		    a1 = con - y_elem[0];
		    a1 = con - a1;
		    a2 = y_elem[0] - a1;
		    con = beta_i[0] * split;
		    b1 = con - beta_i[0];
		    b1 = con - b1;
		    b2 = beta_i[0] - b1;

		    head_t1 = y_elem[0] * beta_i[0];
		    tail_t1 =
		      (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
		  }
		  {
		    /* Compute double_double = double * double. */
		    double a1, a2, b1, b2, con;

		    con = y_elem[1] * split;
		    a1 = con - y_elem[1];
		    a1 = con - a1;
		    a2 = y_elem[1] - a1;
		    con = beta_i[1] * split;
		    b1 = con - beta_i[1];
		    b1 = con - b1;
		    b2 = beta_i[1] - b1;

		    head_t2 = y_elem[1] * beta_i[1];
		    tail_t2 =
		      (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
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
		  /* Imaginary part */
		  {
		    /* Compute double_double = double * double. */
		    double a1, a2, b1, b2, con;

		    con = y_elem[1] * split;
		    a1 = con - y_elem[1];
		    a1 = con - a1;
		    a2 = y_elem[1] - a1;
		    con = beta_i[0] * split;
		    b1 = con - beta_i[0];
		    b1 = con - b1;
		    b2 = beta_i[0] - b1;

		    head_t1 = y_elem[1] * beta_i[0];
		    tail_t1 =
		      (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
		  }
		  {
		    /* Compute double_double = double * double. */
		    double a1, a2, b1, b2, con;

		    con = y_elem[0] * split;
		    a1 = con - y_elem[0];
		    a1 = con - a1;
		    a2 = y_elem[0] - a1;
		    con = beta_i[1] * split;
		    b1 = con - beta_i[1];
		    b1 = con - b1;
		    b2 = beta_i[1] - b1;

		    head_t2 = y_elem[0] * beta_i[1];
		    tail_t2 =
		      (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
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
		y_i[iy] = head_tmp1[0];
		y_i[iy + 1] = head_tmp1[1];
		ai += incai;
		iy += incy;
	      }
	    }
	  }

	} else {

	  /* if beta = 0, we can save m multiplies:
	     y = alpha*A*head_x + alpha*A*tail_x  */
	  if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	    if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {
	      /* save m more multiplies if alpha = 1 */
	      ai = 0;
	      iy = ky;
	      for (i = 0; i < leny; i++) {
		head_sum[0] = head_sum[1] = tail_sum[0] = tail_sum[1] = 0.0;
		head_sum2[0] = head_sum2[1] = tail_sum2[0] = tail_sum2[1] =
		  0.0;
		aij = ai;
		jx = kx;
		for (j = 0; j < lenx; j++) {
		  a_elem[0] = a_i[aij];
		  a_elem[1] = a_i[aij + 1];

		  x_elem[0] = head_x_i[jx];
		  x_elem[1] = head_x_i[jx + 1];
		  {
		    double cd[2];
		    cd[0] = (double) x_elem[0];
		    cd[1] = (double) x_elem[1];
		    {
		      /* Compute complex-extra = complex-double * complex-double. */
		      double head_t1, tail_t1;
		      double head_t2, tail_t2;
		      /* Real part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[0] * split;
			a1 = con - a_elem[0];
			a1 = con - a1;
			a2 = a_elem[0] - a1;
			con = cd[0] * split;
			b1 = con - cd[0];
			b1 = con - b1;
			b2 = cd[0] - b1;

			head_t1 = a_elem[0] * cd[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[1] * split;
			a1 = con - a_elem[1];
			a1 = con - a1;
			a2 = a_elem[1] - a1;
			con = cd[1] * split;
			b1 = con - cd[1];
			b1 = con - b1;
			b2 = cd[1] - b1;

			head_t2 = a_elem[1] * cd[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_prod[0] = head_t1;
		      tail_prod[0] = tail_t1;
		      /* Imaginary part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[1] * split;
			a1 = con - a_elem[1];
			a1 = con - a1;
			a2 = a_elem[1] - a1;
			con = cd[0] * split;
			b1 = con - cd[0];
			b1 = con - b1;
			b2 = cd[0] - b1;

			head_t1 = a_elem[1] * cd[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[0] * split;
			a1 = con - a_elem[0];
			a1 = con - a1;
			a2 = a_elem[0] - a1;
			con = cd[1] * split;
			b1 = con - cd[1];
			b1 = con - b1;
			b2 = cd[1] - b1;

			head_t2 = a_elem[0] * cd[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_prod[1] = head_t1;
		      tail_prod[1] = tail_t1;
		    }
		  }
		  {
		    double head_t, tail_t;
		    double head_a, tail_a;
		    double head_b, tail_b;
		    /* Real part */
		    head_a = head_sum[0];
		    tail_a = tail_sum[0];
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
		    head_sum[0] = head_t;
		    tail_sum[0] = tail_t;
		    /* Imaginary part */
		    head_a = head_sum[1];
		    tail_a = tail_sum[1];
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
		    head_sum[1] = head_t;
		    tail_sum[1] = tail_t;
		  }
		  x_elem[0] = tail_x_i[jx];
		  x_elem[1] = tail_x_i[jx + 1];
		  {
		    double cd[2];
		    cd[0] = (double) x_elem[0];
		    cd[1] = (double) x_elem[1];
		    {
		      /* Compute complex-extra = complex-double * complex-double. */
		      double head_t1, tail_t1;
		      double head_t2, tail_t2;
		      /* Real part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[0] * split;
			a1 = con - a_elem[0];
			a1 = con - a1;
			a2 = a_elem[0] - a1;
			con = cd[0] * split;
			b1 = con - cd[0];
			b1 = con - b1;
			b2 = cd[0] - b1;

			head_t1 = a_elem[0] * cd[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[1] * split;
			a1 = con - a_elem[1];
			a1 = con - a1;
			a2 = a_elem[1] - a1;
			con = cd[1] * split;
			b1 = con - cd[1];
			b1 = con - b1;
			b2 = cd[1] - b1;

			head_t2 = a_elem[1] * cd[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_prod[0] = head_t1;
		      tail_prod[0] = tail_t1;
		      /* Imaginary part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[1] * split;
			a1 = con - a_elem[1];
			a1 = con - a1;
			a2 = a_elem[1] - a1;
			con = cd[0] * split;
			b1 = con - cd[0];
			b1 = con - b1;
			b2 = cd[0] - b1;

			head_t1 = a_elem[1] * cd[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[0] * split;
			a1 = con - a_elem[0];
			a1 = con - a1;
			a2 = a_elem[0] - a1;
			con = cd[1] * split;
			b1 = con - cd[1];
			b1 = con - b1;
			b2 = cd[1] - b1;

			head_t2 = a_elem[0] * cd[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_prod[1] = head_t1;
		      tail_prod[1] = tail_t1;
		    }
		  }
		  {
		    double head_t, tail_t;
		    double head_a, tail_a;
		    double head_b, tail_b;
		    /* Real part */
		    head_a = head_sum2[0];
		    tail_a = tail_sum2[0];
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
		    head_sum2[0] = head_t;
		    tail_sum2[0] = tail_t;
		    /* Imaginary part */
		    head_a = head_sum2[1];
		    tail_a = tail_sum2[1];
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
		    head_sum2[1] = head_t;
		    tail_sum2[1] = tail_t;
		  }
		  aij += incaij;
		  jx += incx;
		}
		{
		  double head_t, tail_t;
		  double head_a, tail_a;
		  double head_b, tail_b;
		  /* Real part */
		  head_a = head_sum[0];
		  tail_a = tail_sum[0];
		  head_b = head_sum2[0];
		  tail_b = tail_sum2[0];
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
		  head_sum[0] = head_t;
		  tail_sum[0] = tail_t;
		  /* Imaginary part */
		  head_a = head_sum[1];
		  tail_a = tail_sum[1];
		  head_b = head_sum2[1];
		  tail_b = tail_sum2[1];
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
		  head_sum[1] = head_t;
		  tail_sum[1] = tail_t;
		}
		y_i[iy] = head_sum[0];
		y_i[iy + 1] = head_sum[1];
		ai += incai;
		iy += incy;
	      }			/* end for */
	    } else {		/* alpha != 1 */
	      ai = 0;
	      iy = ky;
	      for (i = 0; i < leny; i++) {
		head_sum[0] = head_sum[1] = tail_sum[0] = tail_sum[1] = 0.0;
		head_sum2[0] = head_sum2[1] = tail_sum2[0] = tail_sum2[1] =
		  0.0;
		aij = ai;
		jx = kx;
		for (j = 0; j < lenx; j++) {
		  a_elem[0] = a_i[aij];
		  a_elem[1] = a_i[aij + 1];

		  x_elem[0] = head_x_i[jx];
		  x_elem[1] = head_x_i[jx + 1];
		  {
		    double cd[2];
		    cd[0] = (double) x_elem[0];
		    cd[1] = (double) x_elem[1];
		    {
		      /* Compute complex-extra = complex-double * complex-double. */
		      double head_t1, tail_t1;
		      double head_t2, tail_t2;
		      /* Real part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[0] * split;
			a1 = con - a_elem[0];
			a1 = con - a1;
			a2 = a_elem[0] - a1;
			con = cd[0] * split;
			b1 = con - cd[0];
			b1 = con - b1;
			b2 = cd[0] - b1;

			head_t1 = a_elem[0] * cd[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[1] * split;
			a1 = con - a_elem[1];
			a1 = con - a1;
			a2 = a_elem[1] - a1;
			con = cd[1] * split;
			b1 = con - cd[1];
			b1 = con - b1;
			b2 = cd[1] - b1;

			head_t2 = a_elem[1] * cd[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_prod[0] = head_t1;
		      tail_prod[0] = tail_t1;
		      /* Imaginary part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[1] * split;
			a1 = con - a_elem[1];
			a1 = con - a1;
			a2 = a_elem[1] - a1;
			con = cd[0] * split;
			b1 = con - cd[0];
			b1 = con - b1;
			b2 = cd[0] - b1;

			head_t1 = a_elem[1] * cd[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[0] * split;
			a1 = con - a_elem[0];
			a1 = con - a1;
			a2 = a_elem[0] - a1;
			con = cd[1] * split;
			b1 = con - cd[1];
			b1 = con - b1;
			b2 = cd[1] - b1;

			head_t2 = a_elem[0] * cd[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_prod[1] = head_t1;
		      tail_prod[1] = tail_t1;
		    }
		  }
		  {
		    double head_t, tail_t;
		    double head_a, tail_a;
		    double head_b, tail_b;
		    /* Real part */
		    head_a = head_sum[0];
		    tail_a = tail_sum[0];
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
		    head_sum[0] = head_t;
		    tail_sum[0] = tail_t;
		    /* Imaginary part */
		    head_a = head_sum[1];
		    tail_a = tail_sum[1];
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
		    head_sum[1] = head_t;
		    tail_sum[1] = tail_t;
		  }
		  x_elem[0] = tail_x_i[jx];
		  x_elem[1] = tail_x_i[jx + 1];
		  {
		    double cd[2];
		    cd[0] = (double) x_elem[0];
		    cd[1] = (double) x_elem[1];
		    {
		      /* Compute complex-extra = complex-double * complex-double. */
		      double head_t1, tail_t1;
		      double head_t2, tail_t2;
		      /* Real part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[0] * split;
			a1 = con - a_elem[0];
			a1 = con - a1;
			a2 = a_elem[0] - a1;
			con = cd[0] * split;
			b1 = con - cd[0];
			b1 = con - b1;
			b2 = cd[0] - b1;

			head_t1 = a_elem[0] * cd[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[1] * split;
			a1 = con - a_elem[1];
			a1 = con - a1;
			a2 = a_elem[1] - a1;
			con = cd[1] * split;
			b1 = con - cd[1];
			b1 = con - b1;
			b2 = cd[1] - b1;

			head_t2 = a_elem[1] * cd[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_prod[0] = head_t1;
		      tail_prod[0] = tail_t1;
		      /* Imaginary part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[1] * split;
			a1 = con - a_elem[1];
			a1 = con - a1;
			a2 = a_elem[1] - a1;
			con = cd[0] * split;
			b1 = con - cd[0];
			b1 = con - b1;
			b2 = cd[0] - b1;

			head_t1 = a_elem[1] * cd[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[0] * split;
			a1 = con - a_elem[0];
			a1 = con - a1;
			a2 = a_elem[0] - a1;
			con = cd[1] * split;
			b1 = con - cd[1];
			b1 = con - b1;
			b2 = cd[1] - b1;

			head_t2 = a_elem[0] * cd[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_prod[1] = head_t1;
		      tail_prod[1] = tail_t1;
		    }
		  }
		  {
		    double head_t, tail_t;
		    double head_a, tail_a;
		    double head_b, tail_b;
		    /* Real part */
		    head_a = head_sum2[0];
		    tail_a = tail_sum2[0];
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
		    head_sum2[0] = head_t;
		    tail_sum2[0] = tail_t;
		    /* Imaginary part */
		    head_a = head_sum2[1];
		    tail_a = tail_sum2[1];
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
		    head_sum2[1] = head_t;
		    tail_sum2[1] = tail_t;
		  }
		  aij += incaij;
		  jx += incx;
		}
		{
		  /* Compute complex-extra = complex-extra * complex-double. */
		  double head_a0, tail_a0;
		  double head_a1, tail_a1;
		  double head_t1, tail_t1;
		  double head_t2, tail_t2;
		  head_a0 = head_sum[0];
		  tail_a0 = tail_sum[0];
		  head_a1 = head_sum[1];
		  tail_a1 = tail_sum[1];
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
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

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
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

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
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

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
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

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
		  head_a0 = head_sum2[0];
		  tail_a0 = tail_sum2[0];
		  head_a1 = head_sum2[1];
		  tail_a1 = tail_sum2[1];
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
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

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
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

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
		    con = alpha_i[0] * split;
		    b1 = con - alpha_i[0];
		    b1 = con - b1;
		    b2 = alpha_i[0] - b1;

		    c11 = head_a1 * alpha_i[0];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

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
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

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
		y_i[iy] = head_tmp1[0];
		y_i[iy + 1] = head_tmp1[1];
		ai += incai;
		iy += incy;
	      }
	    }
	  } else {		/* beta != 0 */
	    if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {
	      /* save m multiplies if alpha = 1 */
	      ai = 0;
	      iy = ky;
	      for (i = 0; i < leny; i++) {
		head_sum[0] = head_sum[1] = tail_sum[0] = tail_sum[1] = 0.0;;
		head_sum2[0] = head_sum2[1] = tail_sum2[0] = tail_sum2[1] =
		  0.0;;
		aij = ai;
		jx = kx;
		for (j = 0; j < lenx; j++) {
		  a_elem[0] = a_i[aij];
		  a_elem[1] = a_i[aij + 1];

		  x_elem[0] = head_x_i[jx];
		  x_elem[1] = head_x_i[jx + 1];
		  {
		    double cd[2];
		    cd[0] = (double) x_elem[0];
		    cd[1] = (double) x_elem[1];
		    {
		      /* Compute complex-extra = complex-double * complex-double. */
		      double head_t1, tail_t1;
		      double head_t2, tail_t2;
		      /* Real part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[0] * split;
			a1 = con - a_elem[0];
			a1 = con - a1;
			a2 = a_elem[0] - a1;
			con = cd[0] * split;
			b1 = con - cd[0];
			b1 = con - b1;
			b2 = cd[0] - b1;

			head_t1 = a_elem[0] * cd[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[1] * split;
			a1 = con - a_elem[1];
			a1 = con - a1;
			a2 = a_elem[1] - a1;
			con = cd[1] * split;
			b1 = con - cd[1];
			b1 = con - b1;
			b2 = cd[1] - b1;

			head_t2 = a_elem[1] * cd[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_prod[0] = head_t1;
		      tail_prod[0] = tail_t1;
		      /* Imaginary part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[1] * split;
			a1 = con - a_elem[1];
			a1 = con - a1;
			a2 = a_elem[1] - a1;
			con = cd[0] * split;
			b1 = con - cd[0];
			b1 = con - b1;
			b2 = cd[0] - b1;

			head_t1 = a_elem[1] * cd[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[0] * split;
			a1 = con - a_elem[0];
			a1 = con - a1;
			a2 = a_elem[0] - a1;
			con = cd[1] * split;
			b1 = con - cd[1];
			b1 = con - b1;
			b2 = cd[1] - b1;

			head_t2 = a_elem[0] * cd[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_prod[1] = head_t1;
		      tail_prod[1] = tail_t1;
		    }
		  }
		  {
		    double head_t, tail_t;
		    double head_a, tail_a;
		    double head_b, tail_b;
		    /* Real part */
		    head_a = head_sum[0];
		    tail_a = tail_sum[0];
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
		    head_sum[0] = head_t;
		    tail_sum[0] = tail_t;
		    /* Imaginary part */
		    head_a = head_sum[1];
		    tail_a = tail_sum[1];
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
		    head_sum[1] = head_t;
		    tail_sum[1] = tail_t;
		  }
		  x_elem[0] = tail_x_i[jx];
		  x_elem[1] = tail_x_i[jx + 1];
		  {
		    double cd[2];
		    cd[0] = (double) x_elem[0];
		    cd[1] = (double) x_elem[1];
		    {
		      /* Compute complex-extra = complex-double * complex-double. */
		      double head_t1, tail_t1;
		      double head_t2, tail_t2;
		      /* Real part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[0] * split;
			a1 = con - a_elem[0];
			a1 = con - a1;
			a2 = a_elem[0] - a1;
			con = cd[0] * split;
			b1 = con - cd[0];
			b1 = con - b1;
			b2 = cd[0] - b1;

			head_t1 = a_elem[0] * cd[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[1] * split;
			a1 = con - a_elem[1];
			a1 = con - a1;
			a2 = a_elem[1] - a1;
			con = cd[1] * split;
			b1 = con - cd[1];
			b1 = con - b1;
			b2 = cd[1] - b1;

			head_t2 = a_elem[1] * cd[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_prod[0] = head_t1;
		      tail_prod[0] = tail_t1;
		      /* Imaginary part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[1] * split;
			a1 = con - a_elem[1];
			a1 = con - a1;
			a2 = a_elem[1] - a1;
			con = cd[0] * split;
			b1 = con - cd[0];
			b1 = con - b1;
			b2 = cd[0] - b1;

			head_t1 = a_elem[1] * cd[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[0] * split;
			a1 = con - a_elem[0];
			a1 = con - a1;
			a2 = a_elem[0] - a1;
			con = cd[1] * split;
			b1 = con - cd[1];
			b1 = con - b1;
			b2 = cd[1] - b1;

			head_t2 = a_elem[0] * cd[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_prod[1] = head_t1;
		      tail_prod[1] = tail_t1;
		    }
		  }
		  {
		    double head_t, tail_t;
		    double head_a, tail_a;
		    double head_b, tail_b;
		    /* Real part */
		    head_a = head_sum2[0];
		    tail_a = tail_sum2[0];
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
		    head_sum2[0] = head_t;
		    tail_sum2[0] = tail_t;
		    /* Imaginary part */
		    head_a = head_sum2[1];
		    tail_a = tail_sum2[1];
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
		    head_sum2[1] = head_t;
		    tail_sum2[1] = tail_t;
		  }
		  aij += incaij;
		  jx += incx;
		}
		{
		  double head_t, tail_t;
		  double head_a, tail_a;
		  double head_b, tail_b;
		  /* Real part */
		  head_a = head_sum[0];
		  tail_a = tail_sum[0];
		  head_b = head_sum2[0];
		  tail_b = tail_sum2[0];
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
		  head_sum[0] = head_t;
		  tail_sum[0] = tail_t;
		  /* Imaginary part */
		  head_a = head_sum[1];
		  tail_a = tail_sum[1];
		  head_b = head_sum2[1];
		  tail_b = tail_sum2[1];
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
		  head_sum[1] = head_t;
		  tail_sum[1] = tail_t;
		}
		y_elem[0] = y_i[iy];
		y_elem[1] = y_i[iy + 1];
		{
		  /* Compute complex-extra = complex-double * complex-double. */
		  double head_t1, tail_t1;
		  double head_t2, tail_t2;
		  /* Real part */
		  {
		    /* Compute double_double = double * double. */
		    double a1, a2, b1, b2, con;

		    con = y_elem[0] * split;
		    a1 = con - y_elem[0];
		    a1 = con - a1;
		    a2 = y_elem[0] - a1;
		    con = beta_i[0] * split;
		    b1 = con - beta_i[0];
		    b1 = con - b1;
		    b2 = beta_i[0] - b1;

		    head_t1 = y_elem[0] * beta_i[0];
		    tail_t1 =
		      (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
		  }
		  {
		    /* Compute double_double = double * double. */
		    double a1, a2, b1, b2, con;

		    con = y_elem[1] * split;
		    a1 = con - y_elem[1];
		    a1 = con - a1;
		    a2 = y_elem[1] - a1;
		    con = beta_i[1] * split;
		    b1 = con - beta_i[1];
		    b1 = con - b1;
		    b2 = beta_i[1] - b1;

		    head_t2 = y_elem[1] * beta_i[1];
		    tail_t2 =
		      (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
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
		  /* Imaginary part */
		  {
		    /* Compute double_double = double * double. */
		    double a1, a2, b1, b2, con;

		    con = y_elem[1] * split;
		    a1 = con - y_elem[1];
		    a1 = con - a1;
		    a2 = y_elem[1] - a1;
		    con = beta_i[0] * split;
		    b1 = con - beta_i[0];
		    b1 = con - b1;
		    b2 = beta_i[0] - b1;

		    head_t1 = y_elem[1] * beta_i[0];
		    tail_t1 =
		      (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
		  }
		  {
		    /* Compute double_double = double * double. */
		    double a1, a2, b1, b2, con;

		    con = y_elem[0] * split;
		    a1 = con - y_elem[0];
		    a1 = con - a1;
		    a2 = y_elem[0] - a1;
		    con = beta_i[1] * split;
		    b1 = con - beta_i[1];
		    b1 = con - b1;
		    b2 = beta_i[1] - b1;

		    head_t2 = y_elem[0] * beta_i[1];
		    tail_t2 =
		      (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
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
		  double head_t, tail_t;
		  double head_a, tail_a;
		  double head_b, tail_b;
		  /* Real part */
		  head_a = head_sum[0];
		  tail_a = tail_sum[0];
		  head_b = head_tmp1[0];
		  tail_b = tail_tmp1[0];
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
		  head_tmp2[0] = head_t;
		  tail_tmp2[0] = tail_t;
		  /* Imaginary part */
		  head_a = head_sum[1];
		  tail_a = tail_sum[1];
		  head_b = head_tmp1[1];
		  tail_b = tail_tmp1[1];
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
		  head_tmp2[1] = head_t;
		  tail_tmp2[1] = tail_t;
		}
		y_i[iy] = head_tmp2[0];
		y_i[iy + 1] = head_tmp2[1];
		ai += incai;
		iy += incy;
	      }
	    } else {		/* alpha != 1, the most general form:
				   y = alpha*A*head_x + alpha*A*tail_x + beta*y */
	      ai = 0;
	      iy = ky;
	      for (i = 0; i < leny; i++) {
		head_sum[0] = head_sum[1] = tail_sum[0] = tail_sum[1] = 0.0;;
		head_sum2[0] = head_sum2[1] = tail_sum2[0] = tail_sum2[1] =
		  0.0;;
		aij = ai;
		jx = kx;
		for (j = 0; j < lenx; j++) {
		  a_elem[0] = a_i[aij];
		  a_elem[1] = a_i[aij + 1];

		  x_elem[0] = head_x_i[jx];
		  x_elem[1] = head_x_i[jx + 1];
		  {
		    double cd[2];
		    cd[0] = (double) x_elem[0];
		    cd[1] = (double) x_elem[1];
		    {
		      /* Compute complex-extra = complex-double * complex-double. */
		      double head_t1, tail_t1;
		      double head_t2, tail_t2;
		      /* Real part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[0] * split;
			a1 = con - a_elem[0];
			a1 = con - a1;
			a2 = a_elem[0] - a1;
			con = cd[0] * split;
			b1 = con - cd[0];
			b1 = con - b1;
			b2 = cd[0] - b1;

			head_t1 = a_elem[0] * cd[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[1] * split;
			a1 = con - a_elem[1];
			a1 = con - a1;
			a2 = a_elem[1] - a1;
			con = cd[1] * split;
			b1 = con - cd[1];
			b1 = con - b1;
			b2 = cd[1] - b1;

			head_t2 = a_elem[1] * cd[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_prod[0] = head_t1;
		      tail_prod[0] = tail_t1;
		      /* Imaginary part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[1] * split;
			a1 = con - a_elem[1];
			a1 = con - a1;
			a2 = a_elem[1] - a1;
			con = cd[0] * split;
			b1 = con - cd[0];
			b1 = con - b1;
			b2 = cd[0] - b1;

			head_t1 = a_elem[1] * cd[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[0] * split;
			a1 = con - a_elem[0];
			a1 = con - a1;
			a2 = a_elem[0] - a1;
			con = cd[1] * split;
			b1 = con - cd[1];
			b1 = con - b1;
			b2 = cd[1] - b1;

			head_t2 = a_elem[0] * cd[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_prod[1] = head_t1;
		      tail_prod[1] = tail_t1;
		    }
		  }
		  {
		    double head_t, tail_t;
		    double head_a, tail_a;
		    double head_b, tail_b;
		    /* Real part */
		    head_a = head_sum[0];
		    tail_a = tail_sum[0];
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
		    head_sum[0] = head_t;
		    tail_sum[0] = tail_t;
		    /* Imaginary part */
		    head_a = head_sum[1];
		    tail_a = tail_sum[1];
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
		    head_sum[1] = head_t;
		    tail_sum[1] = tail_t;
		  }
		  x_elem[0] = tail_x_i[jx];
		  x_elem[1] = tail_x_i[jx + 1];
		  {
		    double cd[2];
		    cd[0] = (double) x_elem[0];
		    cd[1] = (double) x_elem[1];
		    {
		      /* Compute complex-extra = complex-double * complex-double. */
		      double head_t1, tail_t1;
		      double head_t2, tail_t2;
		      /* Real part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[0] * split;
			a1 = con - a_elem[0];
			a1 = con - a1;
			a2 = a_elem[0] - a1;
			con = cd[0] * split;
			b1 = con - cd[0];
			b1 = con - b1;
			b2 = cd[0] - b1;

			head_t1 = a_elem[0] * cd[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[1] * split;
			a1 = con - a_elem[1];
			a1 = con - a1;
			a2 = a_elem[1] - a1;
			con = cd[1] * split;
			b1 = con - cd[1];
			b1 = con - b1;
			b2 = cd[1] - b1;

			head_t2 = a_elem[1] * cd[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_prod[0] = head_t1;
		      tail_prod[0] = tail_t1;
		      /* Imaginary part */
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[1] * split;
			a1 = con - a_elem[1];
			a1 = con - a1;
			a2 = a_elem[1] - a1;
			con = cd[0] * split;
			b1 = con - cd[0];
			b1 = con - b1;
			b2 = cd[0] - b1;

			head_t1 = a_elem[1] * cd[0];
			tail_t1 =
			  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute double_double = double * double. */
			double a1, a2, b1, b2, con;

			con = a_elem[0] * split;
			a1 = con - a_elem[0];
			a1 = con - a1;
			a2 = a_elem[0] - a1;
			con = cd[1] * split;
			b1 = con - cd[1];
			b1 = con - b1;
			b2 = cd[1] - b1;

			head_t2 = a_elem[0] * cd[1];
			tail_t2 =
			  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) +
			  a2 * b2;
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
		      head_prod[1] = head_t1;
		      tail_prod[1] = tail_t1;
		    }
		  }
		  {
		    double head_t, tail_t;
		    double head_a, tail_a;
		    double head_b, tail_b;
		    /* Real part */
		    head_a = head_sum2[0];
		    tail_a = tail_sum2[0];
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
		    head_sum2[0] = head_t;
		    tail_sum2[0] = tail_t;
		    /* Imaginary part */
		    head_a = head_sum2[1];
		    tail_a = tail_sum2[1];
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
		    head_sum2[1] = head_t;
		    tail_sum2[1] = tail_t;
		  }
		  aij += incaij;
		  jx += incx;
		}
		{
		  /* Compute complex-extra = complex-extra * complex-double. */
		  double head_a0, tail_a0;
		  double head_a1, tail_a1;
		  double head_t1, tail_t1;
		  double head_t2, tail_t2;
		  head_a0 = head_sum[0];
		  tail_a0 = tail_sum[0];
		  head_a1 = head_sum[1];
		  tail_a1 = tail_sum[1];
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
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

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
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

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
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

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
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

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
		  head_a0 = head_sum2[0];
		  tail_a0 = tail_sum2[0];
		  head_a1 = head_sum2[1];
		  tail_a1 = tail_sum2[1];
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
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

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
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

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
		    con = alpha_i[0] * split;
		    b1 = con - alpha_i[0];
		    b1 = con - b1;
		    b2 = alpha_i[0] - b1;

		    c11 = head_a1 * alpha_i[0];
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

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
		    c21 =
		      (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

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
		y_elem[0] = y_i[iy];
		y_elem[1] = y_i[iy + 1];
		{
		  /* Compute complex-extra = complex-double * complex-double. */
		  double head_t1, tail_t1;
		  double head_t2, tail_t2;
		  /* Real part */
		  {
		    /* Compute double_double = double * double. */
		    double a1, a2, b1, b2, con;

		    con = y_elem[0] * split;
		    a1 = con - y_elem[0];
		    a1 = con - a1;
		    a2 = y_elem[0] - a1;
		    con = beta_i[0] * split;
		    b1 = con - beta_i[0];
		    b1 = con - b1;
		    b2 = beta_i[0] - b1;

		    head_t1 = y_elem[0] * beta_i[0];
		    tail_t1 =
		      (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
		  }
		  {
		    /* Compute double_double = double * double. */
		    double a1, a2, b1, b2, con;

		    con = y_elem[1] * split;
		    a1 = con - y_elem[1];
		    a1 = con - a1;
		    a2 = y_elem[1] - a1;
		    con = beta_i[1] * split;
		    b1 = con - beta_i[1];
		    b1 = con - b1;
		    b2 = beta_i[1] - b1;

		    head_t2 = y_elem[1] * beta_i[1];
		    tail_t2 =
		      (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
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
		  /* Imaginary part */
		  {
		    /* Compute double_double = double * double. */
		    double a1, a2, b1, b2, con;

		    con = y_elem[1] * split;
		    a1 = con - y_elem[1];
		    a1 = con - a1;
		    a2 = y_elem[1] - a1;
		    con = beta_i[0] * split;
		    b1 = con - beta_i[0];
		    b1 = con - b1;
		    b2 = beta_i[0] - b1;

		    head_t1 = y_elem[1] * beta_i[0];
		    tail_t1 =
		      (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
		  }
		  {
		    /* Compute double_double = double * double. */
		    double a1, a2, b1, b2, con;

		    con = y_elem[0] * split;
		    a1 = con - y_elem[0];
		    a1 = con - a1;
		    a2 = y_elem[0] - a1;
		    con = beta_i[1] * split;
		    b1 = con - beta_i[1];
		    b1 = con - b1;
		    b2 = beta_i[1] - b1;

		    head_t2 = y_elem[0] * beta_i[1];
		    tail_t2 =
		      (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
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
		y_i[iy] = head_tmp1[0];
		y_i[iy + 1] = head_tmp1[1];
		ai += incai;
		iy += incy;
	      }
	    }
	  }

	}
      }

      FPU_FIX_STOP;
    }
    break;
  }
}
