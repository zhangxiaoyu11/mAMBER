#include "copyright.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine mexit here]
subroutine mexit(output_unit, status)
   
   !  mexit() - machine-dependent exit() procedure, designed to return an
   !            appropriate (success/failure) value to the operating system.
   
   implicit none
   integer output_unit  ! close this unit if greater than zero
   integer status       ! exit status; error if non-zero

#ifdef MPI
#  include "mpif.h"
#  include "parallel.h"
#endif

#ifdef MPI
   
   !       ...status .gt. 0 implies an error condition, therefore
   !       kill all the nodes.  This is accomplished with mpi_abort followed
   !       by a call to MPI_FINALIZE.
   
   if (status /= 0) then
      call mpi_abort(mpi_comm_world, status, ierr)
      call exit(status)
   else
      call mpi_finalize(ierr)
   end if
#endif

#ifdef IBM_VM_CMS
   if (output_unit /= 0 .and. output_unit /= 6) then
      close(unit=output_unit)
   end if
#else
   if (output_unit > 0) then
      close(unit=output_unit)
   end if
#endif

#if XLF90 || IBM3090 || F2C
   if (status /= 0) then
      stop 1
   else
      stop 0
   end if
#else
#  ifdef UXP
   call setrcd(status)
#  else
   call exit(status)
#  endif
#endif
end subroutine mexit 
