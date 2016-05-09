
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zsbmv_z_c(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, int k, const void *alpha, const void *a, int lda,
		    const void *x, int incx, const void *beta,
		    void *y, int incy);


extern void FC_FUNC_(blas_zsbmv_z_c, BLAS_ZSBMV_Z_C)
 
  (int *uplo, int *n, int *k, const void *alpha, const void *a, int *lda,
   const void *x, int *incx, const void *beta, void *y, int *incy) {
  BLAS_zsbmv_z_c(blas_colmajor, (enum blas_uplo_type) *uplo, *n, *k, alpha, a,
		 *lda, x, *incx, beta, y, *incy);
}
