! <compile=optimized>
#include "copyright.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine get_stack here]
subroutine get_stack(ipointer,isize)
   implicit none
   integer ipointer,isize
#  include "memory.h"
#  include "rstack.h"
   ipointer=top_stk+1
   if(top_stk+isize > lastrst )then
      write(6,*) "Exceeding lastrst in get_stack"
      write(6,*) "  lastrst = ",lastrst
      write(6,*) "  top_stk= ",top_stk
      write(6,*) "  isize  = ",isize
      write(6,*) "  request= ",top_stk+isize
      write(6,*) " Increase lastrst in the &cntrl namelist"
      call mexit(6,1)
   end if
   top_stk = top_stk+isize
   num_stkptrs=num_stkptrs+1
   last_stkptr=ipointer
   if(num_stkptrs > max_stack_ptrs)then
      write(6,*) "Exceeding MAX_STACK_PTRS in get_stack() "
      write(6,*) "increase MAX_STACK_PTRS in rstack.h"
      call mexit(6,1)
   end if
   rstk_ptr(num_stkptrs)=ipointer
   highest_stk=max(highest_stk,top_stk)
   return
end subroutine get_stack 

!------------------------------------------------------
!         FREE_STACK
!------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine free_stack here]
subroutine free_stack(ipointer)
   implicit none
#  include "rstack.h"
   integer ipointer
   if(ipointer /= last_stkptr)then
      write(6,*) "Trying to free stack from the middle"
      write(6,*) "highest pointer is ",last_stkptr
      write(6,*) "free req ptr is    ",ipointer
      write(6,*) "   Must free from top (last one made freed first)"
      call mexit(6,1)
   end if
   num_stkptrs=num_stkptrs-1
   if(num_stkptrs < 0)then
      write(6,*) "Problem in stack, freed too much."
      write(6,*) " There are zero pointers left on the stack"
      write(6,*) " cannot free anymore"
      call mexit(6,1)
   end if
   if(num_stkptrs > 0)last_stkptr=rstk_ptr(num_stkptrs)
   top_stk = ipointer-1
   if(top_stk < bot_stk)then
      write(6,*) "Problem in stack, freed too much, corrupt pointer"
      write(6,*) "top of stack is below the top of heap"
      call mexit(6,1)
   end if
   return
end subroutine free_stack 

!------------------------------------------------------
!         RSTACK_SETUP
!------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine rstack_setup here]
subroutine rstack_setup()
   implicit none
   integer top_heap
#  include "rstack.h"
   top_stk=0
   bot_stk=0
   highest_stk=0
   num_stkptrs=0
   return
end subroutine rstack_setup 
