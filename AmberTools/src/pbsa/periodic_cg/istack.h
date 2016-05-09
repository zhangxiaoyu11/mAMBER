
!-------------BEGIN istack.h  ------------------------------------------------
integer max_istack_ptrs
parameter( max_istack_ptrs=100 )
integer ibot_stk,itop_stk,inum_stkptrs, &
      ilast_stkptr,ihighest_stk,istk_ptr(max_istack_ptrs)
common /int_stk/ibot_stk,itop_stk,inum_stkptrs, &
      ilast_stkptr,ihighest_stk,istk_ptr
!-------------END   istack.h  ------------------------------------------------

