! <compile=optimized>
#include "copyright.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine get_istack here]
subroutine get_istack(ipointer,isize)
   implicit none
   integer ipointer,isize

#  include "memory.h"
#  include "istack.h"

   ipointer=itop_stk+1
   if(itop_stk+isize > lastist )then
      write(6,*) "Exceeding lastist in get_istack"
      write(6,*) "  lastist = ",lastist
      write(6,*) "  top_stk= ",itop_stk
      write(6,*) "  isize  = ",isize
      write(6,*) "  request= ",itop_stk+isize
      write(6,*) " Increase lastist in the &cntrl namelist"
      call mexit(6,1)
   end if
   itop_stk = itop_stk+isize
   inum_stkptrs=inum_stkptrs+1
   ilast_stkptr=ipointer
   if(inum_stkptrs > max_istack_ptrs)then
      write(6,*) "Exceeding MAX_ISTACK_PTRS in get_istack() "
      write(6,*) "increase max_istack_ptrs in istack.h"
      call mexit(6,1)
   end if
   istk_ptr(inum_stkptrs)=ipointer
   ihighest_stk=max(ihighest_stk,itop_stk)
!  write(6,'(a,2i10)') '  get_istack: ', itop_stk, isize
   return
end subroutine get_istack 
!------------------------------------------------------
!         FREE_STACK
!------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine free_istack here]
subroutine free_istack(ipointer)
   implicit none
   integer ipointer

#  include "istack.h"
   if(ipointer /= ilast_stkptr)then
      write(6,*)"Trying to free int stack from the middle"
      write(6,*)"highest pointer is ",ilast_stkptr
      write(6,*)"free req ptr is    ",ipointer
      write(6,*)"   Must free from top (last one made freed first"
      call mexit(6,1)
   end if
   inum_stkptrs=inum_stkptrs-1
   if(inum_stkptrs < 0)then
      write(6,*)"Problem in int stack, freed too much."
      write(6,*)" There are zero pointers left on the stack"
      write(6,*)" cannot free anymore"
      call mexit(6,1)
   end if
   if (inum_stkptrs > 0)ilast_stkptr=istk_ptr(inum_stkptrs)
   itop_stk = ipointer-1
   if(itop_stk < ibot_stk)then
      write(6,*)"Problem in integer stack, freed too much."
      write(6,*)"top of stack is below the top of heap"
      call mexit(6,1)
   end if
!  write(6,'(a,i10)') ' free_istack: ', itop_stk
   return
end subroutine free_istack 
!------------------------------------------------------
!         ISTACK_SETUP
!------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine istack_setup here]
subroutine istack_setup()
   implicit none
#  include "istack.h"
   
   itop_stk=0
   ibot_stk=0
   ihighest_stk=0
   inum_stkptrs=0
   return
end subroutine istack_setup 
