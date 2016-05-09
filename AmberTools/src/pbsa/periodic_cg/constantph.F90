#include "copyright.h"
#include "dprec.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Select random new state for random titrating residue, set charges
subroutine cnstphbeginstep(stateinf,resstate,trescnt, &
      statene,chrgdat,dcharge, iselres, iselstat, gen)
   implicit none
#  include "dynph.h"
#  include "random.h"
#  include "md.h"   
   type(const_ph_info), intent(in) :: stateinf(0:*)
   integer, intent(in) :: resstate(0:*), trescnt
   integer, intent(out) :: iselres, iselstat
   _REAL_, intent(in) ::  statene(0:*), chrgdat(0:*)
   _REAL_, intent(out) :: dcharge(1:*)
   type (rand_gen_state), intent(inout) :: gen
   integer :: i
   _REAL_ randval

   call amrand_gen(gen,randval)         !Select residue
   iselres = int((randval*0.9999999d0)*trescnt)
   call amrand_gen(gen,randval)         !Select new state (always different from current)
   iselstat = int((randval*0.9999999d0) &
         *(stateinf(iselres)%num_states - 1))
   if (iselstat >= resstate(iselres)) then
      iselstat = iselstat + 1
   end if
   do i = 0, stateinf(iselres)%num_atoms - 1 !Set charges for new state
      dcharge(i+stateinf(iselres)%first_atom) &
            =chrgdat(stateinf(iselres)%first_charge + iselstat &
            * stateinf(iselres)%num_atoms &
            + i)
   end do
end subroutine cnstphbeginstep 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Evaluate transition energy for proposed protonation state, accept or revert
subroutine cnstphendstep (stateinf,resstate,protcnt, &
      trescnt,statene,chrgdat,dcharge,charge, dvdl, &
      iselres,iselstat, gen)
   implicit none
#  include "dynph.h"
#  include "random.h"
#  include "md.h"   
   type(const_ph_info), intent(in) :: stateinf(0:*)
   integer, intent(inout) :: resstate(0:*)
   _REAL_, intent(inout) :: dvdl
   integer, intent(in) :: protcnt(0:*)
   integer, intent(in) :: trescnt, iselres,iselstat
   _REAL_, intent(in) :: statene(0:*), chrgdat(0:*)
   _REAL_, intent(out) :: charge(1:*),dcharge(1:*)
   type (rand_gen_state), intent(inout) :: gen
   _REAL_ randval, deltae
   integer i, statebase

   
   statebase = stateinf(iselres)%first_state
   !     delta energy = E(proposed state) - E(current state)
   deltae = statene(iselstat+statebase) &
         -statene(resstate(iselres)+statebase)
   !     Correct for pH (delta protons * pH * 2.303RT)
   deltae = deltae -(protcnt(iselstat+statebase) &
         -protcnt(resstate(iselres)+statebase)) &
         *solvph*2.303d0*2.d-3*temp0
   call amrand_gen(gen,randval)

   dvdl = dvdl - deltae         !Adjust transition energy for non-elec factors

   if ((dvdl < 0) .or. (randval <= exp(-dvdl/(2.d-3*temp0)))) then
      !     Make new charges permanent
      do i = 0, stateinf(iselres)%num_atoms-1
         charge(i+stateinf(iselres)%first_atom) &
               =chrgdat(stateinf(iselres)%first_charge &
               + iselstat * stateinf(iselres)%num_atoms &
               + i)
      end do
      resstate(iselres) = iselstat
   else
      !     Revert dcharge to old charges
      do i = 0, stateinf(iselres)%num_atoms-1
         dcharge(i+stateinf(iselres)%first_atom) &
               =chrgdat(stateinf(iselres)%first_charge &
               + resstate(iselres) * stateinf(iselres)%num_atoms &
               + i)
      end do
   end if
end subroutine cnstphendstep 

subroutine cnstphwriterestart(stateinf, resstate, protcnt, trescnt, statene, inchrgdat)
#  include "dynph.h"
#  include "files.h"
#  include "constants.h"  
   
   type (const_ph_info), intent (in) :: stateinf(0:TITR_RES_C-1)
   integer, intent(in) :: resstate(0:TITR_RES_C-1), protcnt(0:TITR_STATES_C-1)
   integer, intent(in) :: trescnt
   _REAL_, intent(in) ::  statene(0:TITR_STATES_C-1), inchrgdat(0:ATOM_CHRG_C-1)
   _REAL_ :: chrgdat(0:ATOM_CHRG_C-1)
   integer :: iatom
   logical, save :: first = .true.
   character(len=40) :: resname(0:TITR_RES_C)
   character(len=7) :: stat
   
   common /cnstphresname/ resname
  
   namelist /cnstph/ stateinf, resstate, protcnt, chrgdat, statene,trescnt, resname

  do iatom=0, ATOM_CHRG_C-1
      chrgdat(iatom) = inchrgdat(iatom) / AMBER_ELECTROSTATIC
   end do
   
   if (first) then              ! Need DELIM, so can't use amopen
      if (owrite == 'N') then
         stat = 'NEW'
      else if (owrite == 'O') then
         stat = 'OLD'
      else if (owrite == 'R') then
         stat = 'REPLACE'
      else if (owrite == 'U') then
         stat = 'UNKNOWN'
      end if
      open(unit=CNSTPH_UNIT, file=cprestrt, status=stat, form='FORMATTED', DELIM='APOSTROPHE')
      first = .false.
   else
      open(unit=CNSTPH_UNIT, file=cprestrt, status='OLD', form='FORMATTED', DELIM='APOSTROPHE')
   end if
   write(CNSTPH_UNIT, nml=cnstph)
   close(CNSTPH_UNIT)

   call amflsh(CPOUT_UNIT)       ! Make sure all cpout data up to this restart point is on disk
   
end subroutine cnstphwriterestart

subroutine cnstphwrite(resstate,iselres,trescnt)
#  include "dynph.h"
#  include "md.h"
#  include "extra.h"
#  include "files.h"
  
   integer, intent(in) :: resstate(0:TITR_RES_C-1)
   integer, intent(in) :: trescnt, iselres
   
   logical :: full
   logical, save :: first = .true.
   integer :: i
   
   if (.not. master) then 
      return
   end if

   full = (first .or. irespa == nstlim .or. mod(irespa,ntwr)==0)
   
   if (full) then
       write(CPOUT_UNIT, '(a,f8.5)') 'Solvent pH: ',solvph
      write(CPOUT_UNIT, '(a,i8)') 'Monte Carlo step size: ',ntcnstph
      write(CPOUT_UNIT, '(a,i8)') 'Time step: ',irespa
      write(CPOUT_UNIT, '(a,f10.3)') 'Time: ',t
      do i=0,trescnt-1
         write(CPOUT_UNIT, '(a,i4,a,i2)') 'Residue ',i,' State: ',resstate(i)
      end do
   else
      write(CPOUT_UNIT, '(a,i4,a,i2)') 'Residue ',iselres,' State: ',resstate(iselres)
   end if
   write(CPOUT_UNIT, '()')
   first = .false.
end subroutine cnstphwrite

  
