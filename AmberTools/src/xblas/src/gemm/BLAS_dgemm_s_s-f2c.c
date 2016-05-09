
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dgemm_s_s(enum blas_order_type order, enum blas_trans_type transa,
		    enum blas_trans_type transb, int m, int n, int k,
		    double alpha, const float *a, int lda, const float *b,
		    int ldb, double beta, double *c, int ldc);


extern void FC_FUNC_(blas_dgemm_s_s, BLAS_DGEMM_S_S)
 
  (int *transa, int *transb, int *m, int *n, int *k, double *alpha,
   const float *a, int *lda, const float *b, int *ldb, double *beta,
   double *c, int *ldc) {
  BLAS_dgemm_s_s(blas_colmajor, (enum blas_trans_type) *transa,
		 (enum blas_trans_type) *transb, *m, *n, *k, *alpha, a, *lda,
		 b, *ldb, *beta, c, *ldc);
}
