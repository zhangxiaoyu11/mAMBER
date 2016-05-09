#include "copyright.h"
#include "dprec.h"


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Emit the final minimization report 
!-----------------------------------------------------------------------
!     --- REPORT_MIN_RESULTS ---
!-----------------------------------------------------------------------
! Find the maximum component of the gradient (grdmax),
! emit this and other gradient details (printe),
! optionally do pseudocontact shift constraints,
! optionally decompose energies for mm_pbsa,
! and emit nmr related information (nmrptx).

subroutine report_min_results( nstep, gradient_rms, coordinates, &
      forces, energies, igraph, xx, ix, ih )

   use decomp, only : checkdec, printdec, printpdec
   implicit none

   integer, intent(in)    :: nstep
   _REAL_,  intent(in)    :: gradient_rms
   _REAL_,  intent(in)    :: coordinates(*)
   _REAL_,  intent(in)    :: forces(*)
   _REAL_,  intent(in)    :: energies(51)
   character(len=4), intent(in)    :: igraph(*)      ! atom name map
   _REAL_,  intent(inout) :: xx(*)          ! real dynamic memory
   integer, intent(inout) :: ix(*)          ! integer dynamic memory
   character(len=4), intent(inout) :: ih(*) ! hollerith dynamic memory

#  include "box.h"
#  include "extra.h"
#  include "files.h"
#  include "md.h"
#  include "memory.h"

   character(len=4) :: atom_name_of_gmax
   integer :: atom_number_of_gmax
   _REAL_  :: gradient_max

   if (master) then
      write(6, '(/ /20x,a,/)' ) 'FINAL RESULTS'
      call grdmax( forces, gradient_max, atom_number_of_gmax )
      atom_name_of_gmax = igraph(atom_number_of_gmax)
      if (imin /= 5) rewind(MDINFO_UNIT)
      call printe( nstep, gradient_rms, gradient_max, energies, &
             atom_number_of_gmax, atom_name_of_gmax )
      if (idecomp > 0) call checkdec()
      if (idecomp == 1 .or. idecomp == 2) call printdec(ix)
      if (idecomp == 3 .or. idecomp == 4) call printpdec(ix)
      if (imin /= 5) call amflsh(MDINFO_UNIT)
   end if

   return
end subroutine report_min_results


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Emit a minimization progress report 
!-----------------------------------------------------------------------
!     --- REPORT_MIN_PROGRESS ---
!-----------------------------------------------------------------------
! Find the maximum component of the gradient (grdmax),
! emit this and other gradient details (printe),
! and emit nmr related information (nmrptx).

subroutine report_min_progress( nstep, gradient_rms, forces, energies, igraph )

   implicit none

   integer, intent(in)    :: nstep
   _REAL_,  intent(in)    :: gradient_rms
   _REAL_,  intent(in)    :: forces(*)
   _REAL_,  intent(in)    :: energies(51)
   character(len=4), intent(in)    :: igraph(*)    ! atom name map

#  include "extra.h"
#  include "files.h"
#  include "md.h"

   character(len=4) :: atom_name_of_gmax
   integer :: atom_number_of_gmax
   _REAL_  :: gradient_max

   if (master) then
      call grdmax( forces, gradient_max, atom_number_of_gmax )
      atom_name_of_gmax = igraph(atom_number_of_gmax)
      if (imin /= 5) rewind(MDINFO_UNIT)
      call printe( nstep, gradient_rms, gradient_max, energies, &
             atom_number_of_gmax, atom_name_of_gmax )
      if (imin /= 5) call amflsh(MDINFO_UNIT)
   end if

   return
end subroutine report_min_progress


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Compute the maximum gradient component and the corresponding atom
!-----------------------------------------------------------------------
!     --- GRDMAX ---
!-----------------------------------------------------------------------

subroutine grdmax( gradient, gradient_max, atom_number_of_gmax )

   implicit none
   _REAL_,  intent(in)  :: gradient(*)
   !     This is actually two-dimensional (3,natoms), but to enable
   !     vectorization on IA32 SSE platforms they are treated as
   !     one-dimensional; this may also improve software pipelining !

   _REAL_,  intent(out) :: gradient_max
   integer, intent(out) :: atom_number_of_gmax

#  include "constants.h"
#  include "memory.h"

   integer :: i
   _REAL_  :: gi

   gradient_max = ZERO
   atom_number_of_gmax = 1
   do i = 1,3*natom
      gi = abs(gradient(i))
      if (gi > gradient_max) then
         gradient_max = gi
         atom_number_of_gmax = i
      end if
   end do
   atom_number_of_gmax = (atom_number_of_gmax - 1)/3 + 1

   return
end subroutine grdmax


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Print status of a minimization calculation.
!-----------------------------------------------------------------------
!     --- PRINTE ---
!-----------------------------------------------------------------------
! Output the step number, the gradient rms, the max gradient component,
! and the atom label and atom number of the max gradient component.

subroutine printe( nstep, gradient_rms, gradient_max, energies, &
      atom_number_of_gmax, atom_name_of_gmax )

   implicit none
   integer, intent(in)    :: nstep
   _REAL_,  intent(in)    :: gradient_rms
   _REAL_,  intent(in)    :: gradient_max
   _REAL_,  intent(in)    :: energies(51)
   integer, intent(in)    :: atom_number_of_gmax
   character(len=4), intent(in)    :: atom_name_of_gmax

#  include "md.h"

   _REAL_  epot
   _REAL_  enonb
   _REAL_  enele
   _REAL_  ehbond
   _REAL_  ebond
   _REAL_  eangle
   _REAL_  edihed
   _REAL_  enb14
   _REAL_  eel14
   _REAL_  econst
   _REAL_  epolar
   _REAL_  aveper
   _REAL_  aveind
   _REAL_  avetot
   _REAL_  esurf
   _REAL_  edisp
   _REAL_  e3bod
   _REAL_  diprms
   _REAL_  dipiter
   _REAL_  dipole_temp

   !     ----- Extract Energies. -----

   epot   = energies(1)
   enonb  = energies(2)
   enele  = energies(3)
   ehbond = energies(4)
   ebond  = energies(5)
   eangle = energies(6)
   edihed = energies(7)
   enb14  = energies(8)
   eel14  = energies(9)
   econst = energies(10)
   epolar = energies(11)
   aveper = energies(12)
   aveind = energies(13)
   avetot = energies(14)
   esurf  = energies(15)
   e3bod  = energies(16)
   edisp  = energies(18)
   diprms = energies(22)
   dipiter = energies(23)
   dipole_temp = energies(24)

   write(6,9018)
   write(6,9028) nstep, epot, gradient_rms, gradient_max, &
         atom_name_of_gmax, atom_number_of_gmax
   write(6,9038) ebond,eangle,edihed
   if( igb == 0 ) then
      write(6,9048) enonb,enele,ehbond
   elseif (igb >= 10) then
      write(6,9050) enonb,enele,ehbond
   else
      write(6,9049) enonb,enele,ehbond
   end if
   write(6,9058) enb14,eel14,econst
#ifdef QMMM
   write(6,9079) esurf
#else
   if (igb >= 10) write(6,9074) esurf,edisp
   if( gbsa > 0 ) write(6,9077) esurf
#endif
   if (abs(epolar) > 0.00001 .or. abs(e3bod) > 0.00001) &
         write(6,9068) epolar,e3bod
   if (econst /= 0.0) write(6,9078) epot-econst

   call amflsh(6)

   !     ----- SEND IT TO THE INFO FILE -----

   if (imin /= 5) then
      write(7,9018)
      write(7,9028) nstep, epot, gradient_rms, gradient_max, &
            atom_name_of_gmax, atom_number_of_gmax
      write(7,9038) ebond,eangle,edihed
      if( igb == 0 ) then
         write(7,9048) enonb,enele,ehbond
      elseif (igb >= 10) then
         write(7,9050) enonb,enele,ehbond
      else
         write(7,9049) enonb,enele,ehbond
      end if
      write(7,9058) enb14,eel14,econst
      if( gbsa > 0 ) write(7,9077) esurf
      if (igb >= 10) write(7,9074) esurf,edisp
      if (epolar /= 0.0 .or. e3bod /= 0.0) &
            write(7,9068) epolar,e3bod
      if (econst /= 0.0) write(7,9078) epot-econst
   end if

   9018 format (/ /,3x,'NSTEP',7x,'ENERGY',10x,'RMS',12x,'GMAX',9x, &
         'NAME',4x,'NUMBER')
   9028 format(1x,i6,2x,3(2x,1pe13.4),5x,a4,2x,i7,/)
   9038 format (1x,'BOND    = ',f13.4,2x,'ANGLE   = ',f13.4,2x, &
         'DIHED      = ',f13.4)
   9048 format (1x,'VDWAALS = ',f13.4,2x,'EEL     = ',f13.4,2x, &
         'HBOND      = ',f13.4)
   9049 format (1x,'VDWAALS = ',f13.4,2x,'EEL     = ',f13.4,2x, &
         'EGB        = ',f13.4)
   9050 format (1x,'VDWAALS = ',f13.4,2x,'EEL     = ',f13.4,2x, &
         'EPB        = ',f13.4)
   9058 format (1x,'1-4 VDW = ',f13.4,2x,'1-4 EEL = ',f13.4,2x, &
         'RESTRAINT  = ',f13.4)
   9068 format (1x,'EPOLAR  = ',f13.4,2x,'ETHREEB = ',f13.4)
   9074 format (1x,'ECAVITY = ',f13.4,2x,'EDISPER = ',f13.4)
   9077 format (1x,'ESURF   = ',f13.4)
   9078 format (1x,'EAMBER  = ',f13.4)
   9079 format (1x,'ESCF    = ',f13.4)
   9090 format(1x,'Dipole convergence: rms = ',e10.3,' iters = ',f6.2)
   9091 format(1x,'Dipole convergence: rms = ',e10.3, &
         ' temperature = ',f6.2)
   return
end subroutine printe
