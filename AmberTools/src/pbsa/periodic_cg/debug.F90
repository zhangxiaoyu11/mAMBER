! <compile=optimized>
#include "copyright.h"
#include "dprec.h"
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine load_debug here]
subroutine load_debug(nf)
   implicit none
   integer nf
#  include "flocntrl.h"
#  include "debug.h"
#  include "md.h"

   integer ifind,j
   namelist /debugf/do_dir,do_rec,do_adj,do_self,do_bond, &
         do_angle,do_ephi,doxconst,do_cap,do_14, &
         do_debugf,neglgdel,zerochg,zerovdw,zerodip, &
         atomn,nranatm, &
         ranseed,chkvir,dumpfrc,rmsfrc,do_tgt

   ! default flow control all force routines turned on
   do_dir = 1
   do_rec = 1
   do_adj = 1
   do_self = 1
   do_bond = 1
   do_angle = 1
   do_ephi = 1
   doxconst = 1
   do_cap = 1
   do_14 = 1
   do_tgt = 1

   ! debug off by default
   do_debugf = 0
   neglgdel = 5
   zerochg = 0
   zerovdw = 0
   zerodip = 0
   nranatm = 0
   ranseed = 71277
   chkvir = 0
   dumpfrc = 0
   rmsfrc = 0
   do j = 1,natomn
      atomn(j) = 0
   end do
   call nmlsrc('debugf',nf,ifind)
   if ( ifind /= 0)then
      read(5,nml=debugf,err=190)
   end if
   if ( do_debugf == 0 )then
      do_dir = 1
      do_rec = 1
      do_adj = 1
      do_self = 1
      do_bond = 1
      do_angle = 1
      do_ephi = 1
      doxconst = 1
      do_tgt = 1
      do_cap = 1

      neglgdel = 5
      zerochg = 0
      zerovdw = 0
      zerodip = 0
      nranatm = 0
      ranseed = 71277
      chkvir = 0
      dumpfrc = 0
      rmsfrc = 0
   end if
   if ( do_debugf == 1 )then
      ! disable list calls, use old list logic
      ! this prevents strange scaling behavior when
      !nbflag = 0
      !ntnb = 0
      !nbfilter = 1.5d0*cutlist
   end if
   return
   190 write(6,*)'Error in &debugf namelist'
   call mexit(6,1)
end subroutine load_debug 
