! <compile=optimized>
#include "copyright.h"
#include "dprec.h"
!------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fill_timer_array here]
subroutine fill_timer_array()
   implicit none
#  include "def_time.h"

   ! after adding the defined integer to "def_time.h"
   ! add a call to add_timer
   ! the arguments are the index of the timer, the
   ! index of its parent timer (must have been added before
   ! unless its -1) and an output string <= 20 characters

   call add_timer(TIME_TOTAL,-1,'Total time')
   call add_timer(TIME_RDPARM, TIME_TOTAL,'Read topology time')
   call add_timer(TIME_RDCRD,  TIME_TOTAL,'Read coords time')
   call add_timer(TIME_FASTWT, TIME_TOTAL,'Fast Water setup')
   call add_timer(TIME_RUNMD,  TIME_TOTAL,'Runmd Time')
#ifdef PSANDER
       call add_timer(TIME_DCMPSETUP,  TIME_RUNMD,'Decomp setup time')
         call add_timer(TIME_DIRECTFILL, TIME_DCMPSETUP, &
              'fill_spatial_direct')
         call add_timer(TIME_RECIPFILL, TIME_DCMPSETUP, &
              'fill_spatial_recip')
         call add_timer(TIME_SHAKEUPDATE, TIME_DCMPSETUP, &
              'shake_spatial_update')
         call add_timer(TIME_NEEDFRCFILL, TIME_DCMPSETUP, &
              'dcmp_fill_needfrc')
#endif
   call add_timer(TIME_FORCE,  TIME_RUNMD,'Force time')
   call add_timer(TIME_NONBON,  TIME_FORCE,'Nonbond force')
   !-- list
   call add_timer(TIME_LIST,   TIME_NONBON,'List time')
   call add_timer(TIME_BLDLST, TIME_LIST,'Build the list')
   !--ewald
   call add_timer(TIME_EWALD,  TIME_NONBON,'Ewald time')
   call add_timer(TIME_DIR,    TIME_EWALD,'Direct Ewald time')
   call add_timer(TIME_SHORT_ENE, TIME_DIR,'Short_ene time')
   call add_timer(TIME_ADJ,   TIME_EWALD,'Adjust Ewald time')
   call add_timer(TIME_SELF,  TIME_EWALD,'Self Ewald time')
   call add_timer(TIME_REC,   TIME_EWALD,'Recip Ewald time')
   call add_timer(TIME_EWFRCADJ, TIME_EWALD,"Force Adjust")
   call add_timer(TIME_EWVIRIAL, TIME_EWALD,"Virial junk")
   call add_timer(TIME_EWFSTRT, TIME_EWALD,"Start sycnronization")
   call add_timer(TIME_BSPL,TIME_REC,'Fill Bspline coeffs')
   call add_timer(TIME_FILLG,    TIME_REC,'Fill charge grid')
   call add_timer(TIME_SCSUM,    TIME_REC,'Scalar sum')
   call add_timer(TIME_GRADS,    TIME_REC,'Grad sum')
   call add_timer(TIME_FFT,      TIME_REC,'FFT time')
   call add_timer(TIME_FFTCOMM, &
         TIME_FFT,'FFT communication time')
   call add_timer(TIME_EWFACTORS,TIME_REC,'Calc trig tables')
   call add_timer(TIME_EWRECSUM, &
         TIME_REC,'Regular Ewald sum')
   call add_timer(TIME_DISTDIP, &
         TIME_EWALD,'dipole distribute time')
   call add_timer(TIME_COLLFIELD, &
         TIME_EWALD,'Field Collect time')
   call add_timer(TIME_LESADJ,TIME_EWALD,'LES adjust time')
   !-- egb
   call add_timer(TIME_EGB,TIME_NONBON,'Gen Born time')
   call add_timer(TIME_GBRAD1,TIME_EGB,'Calc gb radii')
   call add_timer(TIME_GBRADDIST, &
         TIME_EGB,'Communicate gb radii')
   call add_timer(TIME_GBRAD2,TIME_EGB,'Calc gb diag')
   call add_timer(TIME_GBFRC,TIME_EGB,'Calc gb off-diag')
   call add_timer(TIME_GBSA,TIME_EGB,'Surface area energy')
   !-- pb
   call add_timer(TIME_PBFORCE,TIME_NONBON,'PB Nonbond')
   call add_timer(TIME_PBLIST,TIME_PBFORCE,'PB NB list')
   call add_timer(TIME_PBSETUP,TIME_PBFORCE,'PB FD grid')
   call add_timer(TIME_PBSAS,TIME_PBFORCE,'PB Sasa')
   call add_timer(TIME_PBFDFRC,TIME_PBFORCE,'PB FD Force')
   call add_timer(TIME_PBEPS,TIME_PBFDFRC,'PB Epsmap')
   call add_timer(TIME_PBSOLV,TIME_PBFDFRC,'PB Solver')
   call add_timer(TIME_PBDBE,TIME_PBFDFRC,'PB DB Force')
   call add_timer(TIME_PBDIRECT,TIME_PBFORCE,'PB Direct')
   call add_timer(TIME_PBMP,TIME_PBFORCE,'PB Multiple')
   !-- np
   call add_timer(TIME_NPFORCE,TIME_NONBON,'NP Nonbond')
   call add_timer(TIME_NPSAS,TIME_NPFORCE,'NP Sasa')
   call add_timer(TIME_NPCAV,TIME_NPFORCE,'NP Cavity')
   call add_timer(TIME_NPDIS,TIME_NPFORCE,'NP Dispersion')
   !-- misc nonbond
   call add_timer(TIME_TRIPDIP,TIME_NONBON,'Triple dipole time')
   call add_timer(TIME_YAMMPNB,TIME_NONBON,'Yammp nonbond time')
   !-- bonded
   call add_timer(TIME_QMMM,   TIME_FORCE,'QMMM energy')
   call add_timer(TIME_BOND,   TIME_FORCE,'Bond/Angle/Dihedral')
   !-- force communications
   call add_timer(TIME_COLLFRC,TIME_FORCE,'FRC Collect time')
#ifdef PSANDER
   call add_timer(TIME_FPCKSND,TIME_COLLFRC,'Pack n send')
   call add_timer(TIME_FSNDRCV,TIME_COLLFRC,'send recv')
   call add_timer(TIME_FWTRCV,TIME_COLLFRC,'wait recv')
#endif
   !-- more runmd stuff
   call add_timer(TIME_SHAKE,   TIME_RUNMD,'Shake time')
#ifdef PSANDER
   call add_timer(TIME_DCMPFSTWAT,   TIME_RUNMD,'Decomp fstwat time')
#endif
   call add_timer(TIME_VERLET,  TIME_RUNMD,'Verlet update time')
   call add_timer(TIME_DIPUP,   TIME_RUNMD,'Dipole update time')
   call add_timer(TIME_EKCMR,   TIME_RUNMD,'Ekcmr time')
   !-- Coord communications
   call add_timer(TIME_DISTCRD, TIME_RUNMD,'CRD distribute time')
#ifdef PSANDER
   call add_timer(TIME_SPATCRD, TIME_DISTCRD,'Spatial CRD dist')
      call add_timer(TIME_CPCKSND, TIME_SPATCRD,'CRD packsend')
      call add_timer(TIME_CWTRCV,  TIME_SPATCRD,'CRD wt/rcv/unpack')
   call add_timer(TIME_PSSKIN,  TIME_DISTCRD,'Psander check skin')
   call add_timer(TIME_LISTCRD, TIME_DISTCRD,'List CRD+V dist')
      call add_timer(TIME_LCRDSETUP,   TIME_LISTCRD,'Find moved atm')
      call add_timer(TIME_LCRDSENDMV,  TIME_LISTCRD,'Send mvd crd')
      call add_timer(TIME_LCRDUPDTCRD, TIME_LISTCRD,'Update moved crd')
      call add_timer(TIME_LCRDPOST,    TIME_LISTCRD,'Redo spatial setup')
         call add_timer(TIME_FILLRECIP,   TIME_LCRDPOST,'send atom lists')
         call add_timer(TIME_SHKUPDATE,   TIME_LCRDPOST,'shake_spatial_update')
         call add_timer(TIME_FILLNEEDFRC, TIME_LCRDPOST,'dcmp_fill_needfrc')
#endif

   !-- NOE
   call add_timer(TIME_NOE,TIME_FORCE,'Noe calc time ')
   call add_timer(TIME_CALRATE,TIME_NOE,'Calrate time')
   call add_timer(TIME_DSPEV,TIME_NOE,'Dspev time')
   call add_timer(TIME_SHFDER,TIME_NOE,'shift der time')
   call add_timer(TIME_KMAT,TIME_NOE,'Kmat time')
   call add_timer(TIME_NOECALC1,TIME_NOE,'Noecalc1 time')
   call add_timer(TIME_NOECALC2,TIME_NOE,'Noecalc2 time')
   call add_timer(TIME_CALDIS,TIME_NOE,'Calc dis time')
   call add_timer(TIME_REMARC,TIME_NOE,'Remarc time')
   call add_timer(TIME_DINTEN,TIME_NOE,'Dinten time')
   call add_timer(TIME_RINGCURR,TIME_NOE,'Ringcurr time')
   call add_timer(TIME_ELECNOE,TIME_NOE,'Electro. noe time')
   call add_timer(TIME_ANISO,TIME_NOE,'Anisotr. noe time')
   call add_timer(TIME_DRATES,TIME_NOE,'Anisotr. noe time')

!cnt N.TAKADA
#ifdef MDM_MD
   call add_timer(TIME_INIT_MDM,TIME_NONBON,'Init_mdm ')
   call add_timer(TIME_E_ELE_MDM,TIME_NONBON,'Eele_mdm ')
   call add_timer(TIME_F_ELE_MDM,TIME_NONBON,'Fele_mdm ')
   call add_timer(TIME_E_VDW_MDM,TIME_NONBON,'Evdw_mdm ')
   call add_timer(TIME_F_VDW_MDM,TIME_NONBON,'Fvdw_mdm ')
   call add_timer(TIME_EMU_MDM,TIME_NONBON,'MDMemu_mdm ')
#endif
!cnt N.TAKADA

   return
end subroutine fill_timer_array
!------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine init_timers here]
subroutine init_timers()
   implicit none
#  include "new_time.h"
#  include "def_time.h"
   integer i,success
   do i = 1,maxtime
      tpar_p(i) = 0
      tchild_p(i) = 0
      tsib_p(i) = 0
      t_level(i) = 0
      t_state(i) = 0
      t_added(i) = 0
      tnext_p(i) = 0
      tprev_p(i) = 0
   end do
   do i = 1,maxtime
      t_accum(i) = 0.d0
      tch_acc(i) = 0.d0
      t_curr(i) = 0.d0
   end do
   do i = 1,maxtime
      t_string(i) = ' '
   end do
   call fill_timer_array()

   ! establish what level each timer is at

   do i = 1,10000
      call get_level(success)
      if ( success == 1 )goto 100
   end do
   100 continue

   ! build the lists(s) for output
   ! should imitate a depth-first tree walk

   call  build_lists()
   return
end subroutine init_timers
!----------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine add_timer here]
subroutine add_timer(index,par_index,string)
   implicit none
#  include "new_time.h"
   integer index,par_index
   character(len=*) string
   integer i

   if ( index > maxtime )then
      write(6,*)'index ',index,' bigger than MAXTIME ',maxtime
      write(6,*)'attempt to add timer ',string
      call mexit(6,1)
   end if
   if ( t_added(index) /= 0 )then
      write(6,8)index,string
      8 format(1x,'add_timer: timer (',i2,1x,a,') already added! ')
      call mexit(6,1)
   end if
   if (par_index > 0 )then
      if(t_added(par_index) == 0 )then
         write(6,*)'add_timer: parent timer not yet added! '
         write(6,9)index,par_index,string
         9 format(1x,'index,par_index,string = ',i2,1x,i2,1x,a)
         call mexit(6,1)
      end if
   end if
   t_added(index) = 1
   t_string(index) = string
   tpar_p(index) = par_index
   if ( par_index > 0 )then

      ! fortran linked list of sibs

      if ( tchild_p(par_index) == 0 )then
         tchild_p(par_index) = index
      else
         i = tchild_p(par_index)
         10 continue
         if ( tsib_p(i) > 0 )then
            i = tsib_p(i)
            goto 10
         else
            tsib_p(i) = index
         end if
      end if
   end if
   return
end subroutine add_timer
!---------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Synchronize processors to ensure timer accuracy.
subroutine timer_barrier( communicator )

   implicit none
   integer communicator
   logical barrier_active
   integer ierr

#ifdef MPI
#  include "mpif.h"
#endif

   barrier_active = .true.
   if ( barrier_active ) then
#ifdef MPI
      call mpi_barrier( communicator, ierr )
#endif
   end if
   return
end subroutine timer_barrier
!---------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine timer_start here]
subroutine timer_start(istart)
   implicit none
#  include "new_time.h"
   _REAL_ tim
   integer istart
   call wallclock(tim)
   if (istart > 0 .and. istart < maxtime)then
      if ( t_state(istart) /= 0 )then
         write(6,*)'Invalid start to time: ', &
               istart,' ',t_string(istart)
         call mexit(6,1)
      end if
      t_curr(istart) = tim
      t_state(istart) = 1
   end if
   return
end subroutine timer_start
!---------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine timer_stop here]
subroutine timer_stop(istop)
   implicit none
#  include "new_time.h"
   _REAL_ tim
   integer istop
   call wallclock(tim)
   if (istop > 0 .and. istop < maxtime)then
      if ( t_state(istop) == 0 )then
         write(6,*)'Invalid stop to time: ', &
               istop,' ',t_string(istop)
         call mexit(6,1)
      end if
      t_accum(istop) = t_accum(istop) + tim - t_curr(istop)
      ! reset, so can't stop again without start
      t_state(istop) = 0
   end if
   return
end subroutine timer_stop
!---------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine timer_stop_start here]
subroutine timer_stop_start(istop,istart)
   implicit none
#  include "new_time.h"
   _REAL_ tim
   integer istop,istart

   ! TEC suggestion... combine above 2 to save calls

   call wallclock(tim)
   if (istop > 0 .and. istop < maxtime)then
      if ( t_state(istop) == 0 )then
         write(6,*)'Invalid stop to time: ',t_string(istop)
         call mexit(6,1)
      end if
      t_accum(istop) = t_accum(istop) + tim - t_curr(istop)

      ! reset, so can't stop again without start

      t_state(istop) = 0
   end if
   if (istart > 0 .and. istart < maxtime)then
      if ( t_state(istart) /= 0 )then
         write(6,*)'Invalid start to time: ',t_string(istart)
         call mexit(6,1)
      end if
      t_curr(istart) = tim
      t_state(istart) = 1
   end if
   return
end subroutine timer_stop_start
!---------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine get_level here]
subroutine get_level(success)
   implicit none
   integer success
#  include "new_time.h"
   integer i
   do i = 1,maxtime

      ! only look at undetermined ones

      if ( t_added(i) == 1)then
         if ( t_level(i) == 0 )then
            if ( tpar_p(i) == -1 )then
               t_level(i) = 1
            else if ( t_level(tpar_p(i)) /= 0 )then
               t_level(i) = t_level(tpar_p(i)) + 1
            end if
         end if
      end if
   end do

   ! check if all defined

   success = 1
   do i = 1,maxtime
      if ( t_added(i) == 1 .and. t_level(i) == 0 )success = 0
   end do
   return
end subroutine get_level
!---------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine build_lists here]
subroutine build_lists()
   implicit none
#  include "new_time.h"
   integer i,j,numfirst,numlast
   do i = 1,maxtime
      if ( t_added(i) == 1 )then
         if ( tsib_p(i) > 0 )then
            if ( tchild_p(tsib_p(i)) > 0 )then

               ! follow the children to get Tnext_P(i)

               j = tchild_p(tsib_p(i))
               5 continue
               if ( tchild_p(j) /= 0 )then
                  j = tchild_p(j)
                  goto 5
               end if
               tnext_p(i) = j
            else
               tnext_p(i) = tsib_p(i)
            end if
         else if ( tpar_p(i) > 0 )then
            tnext_p(i) = tpar_p(i)
         end if
         if ( tnext_p(i) /= 0 )then
            tprev_p(tnext_p(i)) = i
         end if
      end if
   end do
   t_numtree = 0
   do i = 1,maxtime
      if ( t_added(i) == 1 )then
         if ( tnext_p(i) == 0 )then
            t_numtree = t_numtree + 1
            if ( t_numtree > t_maxtree )then
               write(6,*)'Too many timer trees. Check T_MAXTREE'
               call mexit(6,1)
            end if
            t_last(t_numtree) = i
         end if
      end if
   end do

   ! follow ends back to get start of lists

   do j = 1,t_numtree
      i = t_last(j)
      10 continue
      if ( tprev_p(i) /= 0 )then
         i = tprev_p(i)
         goto 10
      end if
      t_first(j) = i
   end do
   return
end subroutine build_lists
!---------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine finish_timers here]
subroutine finish_timers(all)
   implicit none
   _REAL_ all
#  include "new_time.h"
   _REAL_ t1,p1,p2,p3,minp
   integer i,j,k,max

   ! get the max level

   max = 0
   do i = 1,maxtime
      if ( t_level(i) > max )max = t_level(i)
   end do

   ! accumulate child times
   ! do from maxlevel up to level 2 (level 1 has no parent)
   ! first zero them

   do i = 1,maxtime
      tch_acc(i) = 0.d0
   end do
   do j = 1, max-1
      k = max - j + 1
      do i = 1,maxtime
         if ( t_level(i) == k )then

            ! make sure accumulated time is at least child times (roundoff)

            if ( tch_acc(i) > t_accum(i) )t_accum(i) = tch_acc(i)
            if ( tpar_p(i) > 0 )then
               tch_acc(tpar_p(i)) = tch_acc(tpar_p(i)) + t_accum(i)
            end if
         end if
      end do
   end do

   ! get other times

   do i = 1,maxtime
      t_other(i) = t_accum(i) - tch_acc(i)
   end do

   ! minimal percent of parent to be printed

   minp = 0.005d0

   ! T_print = 0 means don't print
   ! T_print = 1 means print
   ! T_print = 2 means there is time unaccounted for by child times
   !             that needs to be printed as well
   ! find who gets printed. Need to go down by levels
   ! you don't get printed unless your parent is

   do i = 1,maxtime
      t_print(i) = 0
   end do
   do j = 1, max
      do i = 1,maxtime
         if ( t_level(i) == j )then
            if ( tpar_p(i) > 0 )then
               t1 = t_accum(tpar_p(i))
            else
               t1 = all
            end if
            if( t1 > 1.d-4 ) then
               p1 = 100.d0*t_accum(i)/t1
            else
               p1 = 0.d0
            end if
            if( t_accum(i) > 1.d-4 ) then
               p2 = 100.d0*t_other(i)/t_accum(i)
               p3 = 100.d0*tch_acc(i)/t_accum(i)
            else
               p2 = 0.d0
               p3 = 0.d0
            end if
            if ( tpar_p(i) > 0 )then

               ! is your parent printable and are you a big enough percent of it??

               if (t_print(tpar_p(i)) > 0 .and. p1 > minp)then

                  ! do you have children and is there some time not in child times??

                  if ( tchild_p(i) > 0 .and. p2 > minp &
                        .and. p3 > minp )then
                     t_print(i) = 2
                  else
                     t_print(i) = 1
                  end if
               end if
            else
               if ( p1 > minp )then
                  if ( tchild_p(i) > 0 .and. p2 > minp )then
                     t_print(i) = 2
                  else
                     t_print(i) = 1
                  end if
               end if
            end if
         end if  ! ( t_level(i) == j )
      end do  !  i = 1,maxtime
   end do  !  j = 1, max

   return
end subroutine finish_timers
!---------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine write_timers here]
subroutine write_timers(all,nf)
   implicit none
   _REAL_ all
   integer nf
#  include "new_time.h"
#  include "def_time.h"
   integer i,j,k
   integer length(10)
   _REAL_ t,p
   _REAL_ total
   character(len=35) space
   character(len=5) short
   character(len=20) other
   other = 'Other'
   space = '                         '
   length(1) = 1
   do j = 2,10
      length(j) = length(j-1) + 3
   end do
   do k = 1,t_numtree
      i = t_first(k)
      10 continue
      if ( t_print(i) > 0 )then
         if ( t_print(i) == 2 )then

            ! need to first print leftover time from sum of children

            short = t_string(i)(1:5)
            t = t_other(i)
            p = 100.d0*t_other(i)/t_accum(i)
            j = t_level(i) + 1
            if ( p > 99.9 )then
               write(nf,29)space(1:length(j)),other,t,short
            else
               write(nf,30)space(1:length(j)),other,t,p,short
            end if
         end if

         ! now print out entry for this time

         if ( tpar_p(i) > 0 )then
            short = t_string(tpar_p(i))(1:5)
            j = tpar_p(i)
            total = t_accum(j)
         else
            short = 'ALL '
            total = all
         end if
         t = t_accum(i)
         p = 100.d0*t/total
         j = t_level(i)
         if ( j > 10 )then
            write(nf,*)'level too big:',t_string(i)
         end if
         if ( p > 99.9 )then
            write(nf,29)space(1:length(j)),t_string(i),t,short
         else
            write(nf,30)space(1:length(j)),t_string(i),t,p,short
         end if
      end if  ! ( t_print(i) > 0 )
      if ( tnext_p(i) /= 0 )then
         i = tnext_p(i)
         goto 10
      end if
   end do  !  k = 1,t_numtree
   29 format('|',a,a,1x,f10.2,' (100.0% of ',a,')')
   30 format('|',a,a,1x,f10.2,' (',f5.2,'% of ',a,')')
   return
end subroutine write_timers
!------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine write_timerstats here]
subroutine write_timerstats(tim,tim2,mint,maxt, &
      to,to2,minto,maxto,nf)
   implicit none
   _REAL_ all
   _REAL_ tim(*),tim2(*),mint(*),maxt(*), &
         to(*),to2(*),minto(*),maxto(*)
   integer nf
#  include "new_time.h"
#  include "def_time.h"
   integer i,j,k
   integer length(10)
   _REAL_ total
   character(len=35) space
   character(len=20) other
   other = 'Other'
   space = '                         '
   length(1) = 1
   do j = 2,10
      length(j) = length(j-1) + 3
   end do
   do k = 1,t_numtree
      i = t_first(k)
      10 continue
      if ( t_print(i) > 0 )then
         if ( t_print(i) == 2 )then
            ! need to first print leftover time from sum of children
            j = t_level(i) + 1
            write(nf,30)space(1:length(j)),other,to(i), &
                  minto(i),maxto(i),to2(i)
         end if
         ! now print out entry for this time
         j = t_level(i)
         if ( j > 10 )then
            write(nf,*)'level too big:',t_string(i)
         end if
         write(nf,30)space(1:length(j)),t_string(i), &
               tim(i),mint(i),maxt(i),tim2(i)
      end if
      if ( tnext_p(i) /= 0 )then
         i = tnext_p(i)
         goto 10
      end if
   end do
   30 format('|',a,a,1x,f10.2,' (',3(f10.2),')')
   return
end subroutine write_timerstats
!------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine profile_time here]
subroutine profile_time(all)
   implicit none
   _REAL_ all
#  include "new_time.h"
#  include "extra.h"
#  include "istack.h"
#  include "rstack.h"

#ifdef MPI
#  include "parallel.h"
#ifdef MPI_DOUBLE_PRECISION
#undef MPI_DOUBLE_PRECISION
#endif
#  include "mpif.h"
#ifdef CRAY_PVP
#define MPI_DOUBLE_PRECISION MPI_REAL8
#endif
   integer ist(mpi_status_size),ier
#endif
   _REAL_ buf(maxtime+1),tim(maxtime), &
         to(maxtime),to2(maxtime),minto(maxtime), &
         maxto(maxtime), &
         tim2(maxtime),mint(maxtime),maxt(maxtime), &
         allave
   integer i,j

   if( master ) write(6,'(/80(1H-)/,''   5.  TIMINGS'',/80(1H-)/)')
#ifndef MPI
   ! single processor code
   call finish_timers(all)
   call write_timers(all,6)
   write(6,'(/,a,i10)') '| Highest rstack allocated: ',highest_stk
   write(6,'(a,i10)') '| Highest istack allocated: ',ihighest_stk
   return
#else
   if ( numtasks == 1 .or. mpi_orig )then
      call finish_timers(all)
      call write_timers(all,6)
      write(6,'(/,a,i10)') '| Highest rstack allocated: ',highest_stk
      write(6,'(a,i10)') '| Highest istack allocated: ',ihighest_stk
      return
   end if

   ! numtasks > 1...get processor by processor info, do statistics

   if ( .not. master)then
      do i = 1,maxtime
         buf(i) = t_accum(i)
      end do
      buf(maxtime+1) = all

      ! send to master

      call mpi_send(buf,maxtime+1,MPI_DOUBLE_PRECISION, &
            0,mytaskid,commsander,ier)
      return
   else

      ! go through processors. first yourself

      open(unit=8,file='profile_mpi')
      j = 0
      write(8,16)j
      call finish_timers(all)
      call write_timers(all,8)

      ! update for average,stats

      do i = 1,maxtime
         tim(i) = t_accum(i)
         tim2(i) = t_accum(i)**2
         mint(i) = t_accum(i)
         maxt(i) = t_accum(i)
         to(i) = t_other(i)
         to2(i) = t_other(i)**2
         minto(i) = t_other(i)
         maxto(i) = t_other(i)
      end do
      allave = all
      do j = 1,numtasks-1
         call mpi_recv(buf,maxtime+1,MPI_DOUBLE_PRECISION, &
               j,j,commsander,ist,ierr)
         do i = 1,maxtime
            t_accum(i) = buf(i)
         end do
         all = buf(maxtime+1)
         write(8,16)j
         call finish_timers(all)
         call write_timers(all,8)

         ! update for average,stats

         do i = 1,maxtime
            tim(i) = tim(i) + t_accum(i)
            tim2(i) = tim2(i) + t_accum(i)**2
            if ( t_accum(i) < mint(i) )mint(i)=t_accum(i)
            if ( t_accum(i) > maxt(i) )maxt(i)=t_accum(i)
            to(i) = to(i) + t_other(i)
            to2(i) = to2(i) + t_other(i)**2
            if ( t_other(i) < minto(i) )minto(i) = t_other(i)
            if ( t_other(i) > maxto(i) )maxto(i) = t_other(i)
         end do
         allave = allave + all
      end do
      do i = 1,maxtime
         t_accum(i) = tim(i)/numtasks
      end do
      all = allave/numtasks

      ! the averages are run through finish_timers.. this sets print options

      write(6,17)
      call finish_timers(all)
      call write_timers(all,6)

      ! do the statistics

      do i = 1,maxtime
         tim(i) = tim(i)/numtasks
         tim2(i) = sqrt(tim2(i)/numtasks - tim(i)**2)
         to(i) = to(i)/numtasks
         to2(i) = sqrt(to2(i)/numtasks - to(i)**2)
      end do
      write(8,18)
      write(8,19)
      call write_timerstats(tim,tim2,mint,maxt,to,to2,minto,maxto,8)
      close(unit=8)
      write(6,'(/,a,i10)') '| Highest rstack allocated: ',highest_stk
      write(6,'(a,i10)') '| Highest istack allocated: ',ihighest_stk
   end if  ! ( .not. master)
   return
#endif
   16 format('|>>>>>>>>PROFILE of TIMES  for process ',i4)
   17 format('|>>>>>>>>PROFILE of Average TIMES>>>>>>>>> ')
   18 format('|>>>>>>>>Statistics of TIMES>>>>>>>>> ')
   19 format('|>>>>>>>>Printed as average time ', &
         '(min,max,sd) >>>>>>>>> ')
end subroutine profile_time
