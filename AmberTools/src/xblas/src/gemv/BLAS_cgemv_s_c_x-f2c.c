
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_cgemv_s_c_x(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, const void *alpha, const float *a,
		      int lda, const void *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);


extern void FC_FUNC_(blas_cgemv_s_c_x, BLAS_CGEMV_S_C_X)
 
  (int *trans, int *m, int *n, const void *alpha, const float *a, int *lda,
   const void *x, int *incx, const void *beta, void *y, int *incy,
   int *prec) {
  BLAS_cgemv_s_c_x(blas_colmajor, (enum blas_trans_type) *trans, *m, *n,
		   alpha, a, *lda, x, *incx, beta, y, *incy,
		   (enum blas_prec_type) *prec);
}
