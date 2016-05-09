
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_chemm_x(enum blas_order_type order, enum blas_side_type side,
		  enum blas_uplo_type uplo, int m, int n,
		  const void *alpha, const void *a, int lda,
		  const void *b, int ldb, const void *beta,
		  void *c, int ldc, enum blas_prec_type prec);


extern void FC_FUNC_(blas_chemm_x, BLAS_CHEMM_X)
 
  (int *side, int *uplo, int *m, int *n, const void *alpha, const void *a,
   int *lda, const void *b, int *ldb, const void *beta, void *c, int *ldc,
   int *prec) {
  BLAS_chemm_x(blas_colmajor, (enum blas_side_type) *side,
	       (enum blas_uplo_type) *uplo, *m, *n, alpha, a, *lda, b, *ldb,
	       beta, c, *ldc, (enum blas_prec_type) *prec);
}
