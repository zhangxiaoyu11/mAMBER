! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"

!Subroutines for the calculation of charges on QM
!atoms in QMMM calculations.

subroutine qm2_calc_mulliken(iqm,mul_chg)

! Calculates the Mulliken atom charge for QM atom iqm regions.
! Written by Ross Walker (TSRI, 2005)

! iqm should be specified as a number from 1 to nquant_nlink.
! Mulliken charge in electron units is returned in mul_chg

! Requires a converged density matrix stored in qm2_struct%den_matrix

      use qmmm_module, only : qm2_params, qm2_struct
      implicit none

!Passed in
      integer, intent(in) :: iqm
      _REAL_, intent(out) :: mul_chg

!Local
      integer :: loop_count, orb_beg, orb_end, tri
      _REAL_ :: density_sum

      density_sum = 0.0d0
      !Find the beginning and ending orbital locations for this atom.
      orb_beg=qm2_params%orb_loc(1,iqm)
      orb_end=qm2_params%orb_loc(2,iqm)

      do loop_count=orb_beg,orb_end
        tri = qm2_params%pascal_tri2(loop_count)
        density_sum = density_sum + qm2_struct%den_matrix(tri)
      end do

      mul_chg = qm2_params%core_chg(iqm) - density_sum

      return
end subroutine qm2_calc_mulliken


