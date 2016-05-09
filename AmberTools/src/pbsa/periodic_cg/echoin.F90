#include "copyright.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine echoin here]
subroutine echoin(in,iout)
   
   
   ! Subroutine ECHO INput
   
   ! This routine is called by the MDREAD routine for an Amber suite program.
   ! It reads and echos the information in the input file to the user.
   
   ! Special provision is made for the lines written to the file by the
   ! Amber/Interface to log the Interface lines used to generate the input
   ! file.
   
   ! Before this routine returns, the file on unit IN will be rewound.
   
   ! Author: David A. Pearlman
   ! Date: 8/92
   
   ! IN: The unit the input file is read from.
   ! IOUT: The unit to which the echoed information is to be written.
   
   character(len=80) aa
   character(len=14) bgflg,enflg
   data bgflg/'!! BEGIN Amber'/
   data enflg/'!! END   Amber'/
   
   ! First echo the Interface script, if any:
   
   inter = 0
   intlin = 0
   do 10 i = 1,999999
      read(in,500,end=20) aa
      if (aa(1:14) == enflg) inter = 0
      if (inter == 1) write(iout,500) aa(1:79)
      if (aa(1:14) == bgflg) then
         inter = 1
         if (intlin == 0) then
            write(iout,*)
            write(iout,1000)
            write(iout,*)
         end if
         intlin = 1
      end if
   10 continue
   
   ! Now echo the standard Amber input lines:
   
   20 if (intlin > 0) write(iout,501)
   rewind(in)
   write(iout,*)
   write(iout,1010)
   write(iout,*)
   inter = 0
   do 30 i = 1,999999
      read(in,500,end=40) aa
      if (aa(1:14) == bgflg) inter = 1
      if (inter == 0) write(iout,500) aa(1:79)
      if (aa(1:14) == enflg) inter = 0
   30 continue
   
   40 continue
   rewind(in)
   return
   
   ! Format statements
   
   500 format(a)
   501 format(79('-'))
   1000 format(' The Interface script used to generate the input file:')
   1010 format(' Here is the input file:')
end subroutine echoin 
