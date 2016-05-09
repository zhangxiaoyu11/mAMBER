
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zgbmv2_d_z(enum blas_order_type order, enum blas_trans_type trans,
		     int m, int n, int kl, int ku, const void *alpha,
		     const double *a, int lda, const void *head_x,
		     const void *tail_x, int incx, const void *beta,
		     void *y, int incy);


extern void FC_FUNC_(blas_zgbmv2_d_z, BLAS_ZGBMV2_D_Z)
 
  (int *trans, int *m, int *n, int *kl, int *ku, const void *alpha,
   const double *a, int *lda, const void *head_x, const void *tail_x,
   int *incx, const void *beta, void *y, int *incy) {
  BLAS_zgbmv2_d_z(blas_colmajor, (enum blas_trans_type) *trans, *m, *n, *kl,
		  *ku, alpha, a, *lda, head_x, tail_x, *incx, beta, y, *incy);
}
