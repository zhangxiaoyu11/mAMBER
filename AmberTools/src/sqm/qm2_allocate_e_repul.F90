! <compile=optimized>
#include "copyright.h"
#include "../include/assert.fh"
!-----------------------------------------------------------------
!Written by Ross Walker (TSRI, 2005)
!Allocates memory for 1-electron repulsion integrals
!-----------------------------------------------------------------

subroutine qm2_allocate_qmqm_e_repul(n2el)

!This routine allocates space for the QM-QM 2 electron repulsion integrals and
!optionally the one electron repulsion integrals. It should only be called once
!per sander run on the first call to QM_MM.

  use qmmm_module, only : qmmm_nml,qm2_struct, qmmm_struct, qmmm_mpi
  implicit none

!Passed in
  integer, intent(in) :: n2el

!Local
  integer :: ier=0

  allocate ( qm2_struct%qm_qm_2e_repul(n2el), stat=ier ) 
  REQUIRE ( ier == 0 )

  if (qmmm_nml%qmqm_erep_incore) then
     !only need the QM-QM electron repulsion integrals stored if we
     !are doing analytical QM-QM derivatives and qmqm_erep_incore = true.
     !qmqm_e_repul needs to be nquant_nlink*(nquant_nlink-1)/2 x 22 long
     !In parallel it only needs to be 
     !qmmm_mpi%nquant_nlink_loop_extent_end-qmmm_mpi%nquant_nlink_loop_extent_begin+1
     !long.
#ifdef MPI
     allocate ( qm2_struct%qm_qm_e_repul(22,qmmm_mpi%nquant_nlink_loop_extent_end- &
                                            qmmm_mpi%nquant_nlink_loop_extent_begin+1), stat=ier )
#else
     allocate ( qm2_struct%qm_qm_e_repul(22, &
              (qmmm_struct%nquant_nlink * (qmmm_struct%nquant_nlink-1)/2)), stat=ier )
#endif
     REQUIRE ( ier == 0 )
  end if
  return
end subroutine qm2_allocate_qmqm_e_repul

