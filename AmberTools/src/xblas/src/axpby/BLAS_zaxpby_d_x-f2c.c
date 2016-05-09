
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zaxpby_d_x(int n, const void *alpha, const double *x, int incx,
		     const void *beta, void *y,
		     int incy, enum blas_prec_type prec);


extern void FC_FUNC_(blas_zaxpby_d_x, BLAS_ZAXPBY_D_X)
 
  (int *n, const void *alpha, const double *x, int *incx, const void *beta,
   void *y, int *incy, int *prec) {
  BLAS_zaxpby_d_x(*n, alpha, x, *incx, beta, y, *incy,
		  (enum blas_prec_type) *prec);
}
