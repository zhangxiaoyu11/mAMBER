
!-------------BEGIN rstack.h  ------------------------------------------------
integer max_stack_ptrs
parameter( max_stack_ptrs=100 )
integer bot_stk,top_stk,num_stkptrs, &
      last_stkptr,highest_stk,rstk_ptr(max_stack_ptrs)
common /real_stk/bot_stk,top_stk,num_stkptrs, &
      last_stkptr,highest_stk,rstk_ptr
!-------------END   rstack.h  ------------------------------------------------

