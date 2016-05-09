
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_chemv_c_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, const void *alpha, const void *a, int lda,
		      const float *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);


extern void FC_FUNC_(blas_chemv_c_s_x, BLAS_CHEMV_C_S_X)
 
  (int *uplo, int *n, const void *alpha, const void *a, int *lda,
   const float *x, int *incx, const void *beta, void *y, int *incy,
   int *prec) {
  BLAS_chemv_c_s_x(blas_colmajor, (enum blas_uplo_type) *uplo, *n, alpha, a,
		   *lda, x, *incx, beta, y, *incy,
		   (enum blas_prec_type) *prec);
}
