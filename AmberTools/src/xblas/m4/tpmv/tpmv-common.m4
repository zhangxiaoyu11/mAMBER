dnl
dnl
define(`TPMV_NAME', 
  `ifelse(`$2', `$1', `BLAS_$1tpmv$3', `BLAS_$1tpmv_$2$3')')dnl
dnl
dnl
define(`TPMV_PARAMS', 
  `enum blas_order_type order, enum blas_uplo_type uplo, 
   enum blas_trans_type trans, enum blas_diag_type diag, 
   int n, $1_scalar alpha, const $2_array tp, 
   $1_array x, int incx`'ifelse($3, _x, `, enum blas_prec_type prec')')dnl
dnl
dnl
define(`TPMV_HEAD',
  `void TPMV_NAME($1, $2, $3)(TPMV_PARAMS($1, $2, $3))')dnl
dnl
dnl
define(`TPMV_ARGS', 
`if_blas(``s, s', `d, d', `c, c', `z, z',')dnl
`d, s', `z, c', `c, s', `z, d',
 `s, s, _x', `d, d, _x', `c, c, _x', `z, z, _x',
 `d, s, _x', `z, c, _x', `c, s, _x', `z, d, _x'')dnl
dnl
dnl
