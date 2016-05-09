
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_caxpby_x(int n, const void *alpha, const void *x, int incx,
		   const void *beta, void *y,
		   int incy, enum blas_prec_type prec);


extern void FC_FUNC_(blas_caxpby_x, BLAS_CAXPBY_X)
 
  (int *n, const void *alpha, const void *x, int *incx, const void *beta,
   void *y, int *incy, int *prec) {
  BLAS_caxpby_x(*n, alpha, x, *incx, beta, y, *incy,
		(enum blas_prec_type) *prec);
}
