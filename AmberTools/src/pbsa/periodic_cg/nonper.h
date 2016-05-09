!---------------- nonper.h -----------------
#define BC_NONPER 3
!    parameters for nonperiodic box specs
_REAL_ extraboxdim
parameter (extraboxdim=30.d0)
!   Center of box coords used for centering imagcrds in ew_direct.f
_REAL_ xbox0,ybox0,zbox0
common /nonper_real/xbox0,ybox0,zbox0
!----------------end of nonper.h ---------
