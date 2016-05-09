#include "copyright.h"
#include "dprec.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine minrit here]
subroutine minrit(x)
   
   implicit _REAL_ (a-h,o-z)
   
   !     ----- ROUTINE TO WRITE FINAL COORDINATES -----
   
#  include "md.h"
#  include "files.h"
#  include "box.h"
   dimension x(*)
   
   nr = nrp
   nr3 = 3*nr
   tt = 0.0d0

   if (ntxo <= 0) then
      call amopen(16,restrt,owrite,'U','W')
   else
      call amopen(16,restrt,owrite,'F','W')
   end if
   
   if (ntxo /= 0) then
      write(16,40) title
      if( nr > 100000 ) then
         write(16,221) nr
      else
         write(16,220) nr
      end if
      write(16,292) (x(i),i=1,nr3)
   else
      write(16) title
      dumm = 0.0d0
      write(16) nr,dumm
      write(16) (x(i),i = 1,nr3)
   end if
   close(unit=16)
   40 format(20a4)
   220 format(i5,1x,f10.5,i5)
   221 format(i6,f10.5,i5)
   292 format(6f12.7)
   return
end subroutine minrit 
