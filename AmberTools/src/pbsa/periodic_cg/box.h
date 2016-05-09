!     common block sizes:

integer bc_boxi,bc_boxr
parameter(bc_boxi=7)
parameter(bc_boxr=613)

! ... floats:

#include "dprec.h"
_REAL_ box,cut,scnb,scee,dielc,rad,wel,radhb,welhb, &
      cutcap,xcap,ycap,zcap,fcap,rwell
common/boxr/box(3),cut,scnb,scee,dielc, &
      cutcap,xcap,ycap,zcap,fcap,rwell, &
      rad(100),wel(100),radhb(200),welhb(200)

! ... integers:

integer ntb,ifbox,numpk,nbit,ifcap,natcap,isftrp

common/boxi/ntb,ifbox,numpk,nbit,ifcap,natcap,isftrp
