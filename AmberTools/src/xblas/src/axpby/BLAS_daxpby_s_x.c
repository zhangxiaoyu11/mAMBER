#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_daxpby_s_x(int n, double alpha, const float *x, int incx,
		     double beta, double *y,
		     int incy, enum blas_prec_type prec)

/*
 * Purpose
 * =======
 *
 * This routine computes:
 *
 *      y <- alpha * x + beta * y.
 *
 * Arguments
 * =========
 * 
 * n         (input) int
 *           The length of vectors x and y.
 * 
 * alpha     (input) double
 *
 * x         (input) const float*
 *           Array of length n.
 *
 * incx      (input) int
 *           The stride used to access components x[i].
 * 
 * beta      (input) double
 *
 * y         (input) double*
 *           Array of length n.
 * 
 * incy      (input) int
 *           The stride used to access components y[i].
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
  static const char routine_name[] = "BLAS_daxpby_s_x";

  switch (prec) {
  case blas_prec_single:
  case blas_prec_double:
  case blas_prec_indigenous:
    {
      int i, ix = 0, iy = 0;
      const float *x_i = x;
      double *y_i = y;
      double alpha_i = alpha;
      double beta_i = beta;
      float x_ii;
      double y_ii;
      double tmpx;
      double tmpy;


      /* Test the input parameters. */
      if (incx == 0)
	BLAS_error(routine_name, -4, incx, NULL);
      else if (incy == 0)
	BLAS_error(routine_name, -7, incy, NULL);

      /* Immediate return */
      if (n <= 0 || (alpha_i == 0.0 && beta_i == 1.0))
	return;





      if (incx < 0)
	ix = (-n + 1) * incx;
      if (incy < 0)
	iy = (-n + 1) * incy;

      for (i = 0; i < n; ++i) {
	x_ii = x_i[ix];
	y_ii = y_i[iy];
	tmpx = alpha_i * x_ii;	/* tmpx  = alpha * x[ix] */
	tmpy = beta_i * y_ii;	/* tmpy = beta * y[iy] */
	tmpy = tmpy + tmpx;
	y_i[iy] = tmpy;
	ix += incx;
	iy += incy;
      }				/* endfor */


    }
    break;
  case blas_prec_extra:
    {
      int i, ix = 0, iy = 0;
      const float *x_i = x;
      double *y_i = y;
      double alpha_i = alpha;
      double beta_i = beta;
      float x_ii;
      double y_ii;
      double head_tmpx, tail_tmpx;
      double head_tmpy, tail_tmpy;
      FPU_FIX_DECL;

      /* Test the input parameters. */
      if (incx == 0)
	BLAS_error(routine_name, -4, incx, NULL);
      else if (incy == 0)
	BLAS_error(routine_name, -7, incy, NULL);

      /* Immediate return */
      if (n <= 0 || (alpha_i == 0.0 && beta_i == 1.0))
	return;

      FPU_FIX_START;



      if (incx < 0)
	ix = (-n + 1) * incx;
      if (incy < 0)
	iy = (-n + 1) * incy;

      for (i = 0; i < n; ++i) {
	x_ii = x_i[ix];
	y_ii = y_i[iy];
	{
	  double dt = (double) x_ii;
	  {
	    /* Compute double_double = double * double. */
	    double a1, a2, b1, b2, con;

	    con = alpha_i * split;
	    a1 = con - alpha_i;
	    a1 = con - a1;
	    a2 = alpha_i - a1;
	    con = dt * split;
	    b1 = con - dt;
	    b1 = con - b1;
	    b2 = dt - b1;

	    head_tmpx = alpha_i * dt;
	    tail_tmpx =
	      (((a1 * b1 - head_tmpx) + a1 * b2) + a2 * b1) + a2 * b2;
	  }
	}			/* tmpx  = alpha * x[ix] */
	{
	  /* Compute double_double = double * double. */
	  double a1, a2, b1, b2, con;

	  con = beta_i * split;
	  a1 = con - beta_i;
	  a1 = con - a1;
	  a2 = beta_i - a1;
	  con = y_ii * split;
	  b1 = con - y_ii;
	  b1 = con - b1;
	  b2 = y_ii - b1;

	  head_tmpy = beta_i * y_ii;
	  tail_tmpy = (((a1 * b1 - head_tmpy) + a1 * b2) + a2 * b1) + a2 * b2;
	}			/* tmpy = beta * y[iy] */
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_tmpy + head_tmpx;
	  bv = s1 - head_tmpy;
	  s2 = ((head_tmpx - bv) + (head_tmpy - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_tmpy + tail_tmpx;
	  bv = t1 - tail_tmpy;
	  t2 = ((tail_tmpx - bv) + (tail_tmpy - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_tmpy = t1 + t2;
	  tail_tmpy = t2 - (head_tmpy - t1);
	}
	y_i[iy] = head_tmpy;
	ix += incx;
	iy += incy;
      }				/* endfor */

      FPU_FIX_STOP;
    }
    break;
  }
}				/* end BLAS_daxpby_s_x */
