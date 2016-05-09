#include "copyright.h"
#include "dprec.h"
#ifndef MKL
!  vectorized function calls

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ vectorized exponential
subroutine vdexp( n, x, y )
   
   implicit none
   integer  i,n
   _REAL_   x(n), y(n)
   
#  ifdef MASSLIB
      call vexp ( y, x, n )
#  else
      y(1:n) = exp(x(1:n))
#  endif
   
   return
end subroutine vdexp 
!--------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ vectorized logarithm
subroutine vdln( n, x, y )
   
   implicit none
   integer  i,n
   _REAL_   x(n), y(n)
   
#  ifdef MASSLIB
      call vlog ( y, x, n )
#  else
      y(1:n) = log(x(1:n))
#  endif
   
   return
end subroutine vdln 
!--------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ vectorized inverse square root
subroutine vdinvsqrt( n, x, y )
   
   implicit none
   integer  i,n
   _REAL_   x(n), y(n)
#  include "constants.h"
   
#  ifdef MASSLIB
      call vrsqrt ( y, x, n )
#  else
      y(1:n) = one/sqrt( x(1:n) )
#  endif
   
   return
end subroutine vdinvsqrt 
!--------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ vectorized inverse
subroutine vdinv( n, x, y )
   
   implicit none
   integer  i,n
   _REAL_   x(n), y(n)
#  include "constants.h"
   
#  ifdef MASSLIB
      call vrec ( y, x, n )
#  else
      y(1:n) = one/x(1:n)
#  endif
   
   return
end subroutine vdinv 
#else

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ dummy subroutine to avoid compiling an empty file
subroutine vd_dummy_to_avoid_empty_file( )
   
   implicit none
   return
end subroutine vd_dummy_to_avoid_empty_file 
!--------------------------------------------------------------

#endif
