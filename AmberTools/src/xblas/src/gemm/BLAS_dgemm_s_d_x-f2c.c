
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dgemm_s_d_x(enum blas_order_type order, enum blas_trans_type transa,
		      enum blas_trans_type transb, int m, int n, int k,
		      double alpha, const float *a, int lda, const double *b,
		      int ldb, double beta, double *c, int ldc,
		      enum blas_prec_type prec);


extern void FC_FUNC_(blas_dgemm_s_d_x, BLAS_DGEMM_S_D_X)
 
  (int *transa, int *transb, int *m, int *n, int *k, double *alpha,
   const float *a, int *lda, const double *b, int *ldb, double *beta,
   double *c, int *ldc, int *prec) {
  BLAS_dgemm_s_d_x(blas_colmajor, (enum blas_trans_type) *transa,
		   (enum blas_trans_type) *transb, *m, *n, *k, *alpha, a,
		   *lda, b, *ldb, *beta, c, *ldc,
		   (enum blas_prec_type) *prec);
}
