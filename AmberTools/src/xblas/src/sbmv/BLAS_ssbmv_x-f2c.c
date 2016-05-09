
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_ssbmv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  int n, int k, float alpha, const float *a, int lda,
		  const float *x, int incx, float beta,
		  float *y, int incy, enum blas_prec_type prec);


extern void FC_FUNC_(blas_ssbmv_x, BLAS_SSBMV_X)
 
  (int *uplo, int *n, int *k, float *alpha, const float *a, int *lda,
   const float *x, int *incx, float *beta, float *y, int *incy, int *prec) {
  BLAS_ssbmv_x(blas_colmajor, (enum blas_uplo_type) *uplo, *n, *k, *alpha, a,
	       *lda, x, *incx, *beta, y, *incy, (enum blas_prec_type) *prec);
}
