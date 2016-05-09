! <compile=optimized>
#include "copyright.h"
#include "dprec.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ main driver routine to compute energies and forces
subroutine force(xx,ix,ih,ipairs,x,f,ener,vir,r_stack,i_stack, &
      fs, rborn, reff ,do_list_update )

   use genborn
   use poisson_boltzmann, only : pb_force
   use dispersion_cavity, only : npopt, np_force

   implicit none
   logical do_list_update
   integer   ipairs(*)
   _REAL_ xx(*)
   integer   ix(*)
   character(len=4) ih(*)
   _REAL_ fs(*),rborn(*),reff(*),dvdl
#  include "constants.h"
#  include "def_time.h"
#ifdef MPI
#  include "parallel.h"
#ifdef MPI_DOUBLE_PRECISION
#undef MPI_DOUBLE_PRECISION
#endif
#  include "mpif.h"
#endif
   logical belly,nocrst
#  include "md.h"
#  include "box.h"
#  include "memory.h"
#  include "parms.h"
#  include "files.h"
#  include "extra.h"
#  include "flocntrl.h"
#  include "pb_md.h"

   _REAL_  r_stack(*)
   integer i_stack(*)

   _REAL_  enmr(3),devdis(4),devang(4),devtor(4),entr

   _REAL_  x(*),f(*),ene(30),vir(*)
   _REAL_  ener(*) ! offsets in this ener array = offsets in runmd ener - 22
   save ene

   integer i,m,nttyp,npair,nhb
   _REAL_  virvsene,evdw,eelt,e3bod,epol,esurf,edisp
   _REAL_  epolar,aveper,aveind,avetot,emtot,dipiter,dipole_temp
   integer l_r2x,l_rjx,l_tmp1,l_tmp2,l_tmp3,l_tmp4,l_tmp5
   integer l_tmp6,l_tmp7,l_tmp8,l_jj,l_skipv, l_kvls,l_jvls,l_psi
   integer l_da,l_sumd
   integer newbalance
   save newbalance

   call timer_start(TIME_FORCE)
   ene(:) = ZERO

   belly = ibelly > 0
   nocrst = .false.
   nttyp = ntypes*(ntypes+1)/2
#ifdef MPI
   if (mpi_orig) then
      !     =========================== AMBER/MPI ===========================

      !     Check to see if we are done yet in mpi_orig case (tec3).
      !     This is done by monitoring the status of an integer notdone.
      !     If notdone .eq. 1 then we keep going.  notdone is set to zero
      !     when we no longer want to call force().  This perhaps is not the
      !     most efficient means to implement this check...

      call mpi_bcast(notdone,1,mpi_integer,0,commsander,ierr)
      if (notdone /= 1) return

      !       Send copies of xyz coords, setbox common block, vir array
      !       and NTNB value to all nodes from master with a broadcast.

      if (numtasks > 1) then

         call mpi_bcast(box,bc_boxr,MPI_DOUBLE_PRECISION,0, &
               commsander,ierr)
         call mpi_bcast(ntb,bc_boxi,mpi_integer,0,commsander,ierr)
         call mpi_bcast(vir,3,MPI_DOUBLE_PRECISION,0, &
               commsander,ierr)
         call mpi_bcast(xx(lcrd),3*natom,MPI_DOUBLE_PRECISION, &
               0,commsander,ierr)
         call mpi_bcast(ntnb,1,mpi_integer,0,commsander,ierr)
         if (iabs(ntb) >= 2) then
            call mpi_bcast(xx(l45),3*natom,MPI_DOUBLE_PRECISION, &
                  0,commsander,ierr)
         end if
      end if
   end if

   !     ========================= END AMBER/MPI =========================
#endif /* MPI */

   !     ----- ZERO OUT THE ENERGIES AND FORCES -----

   aveper=0.d0
   aveind=0.d0
   avetot=0.d0
   dipiter=0.d0
   dvdl=0.d0
   dipole_temp=0.d0
   do i=1,3
      enmr(i) = 0.d0
   end do
   do i=1,4
      vir(i) = 0.d0
   end do
   virvsene = 0.d0
   do i=1,3*natom
      f(i) = 0.d0
   end do

   epolar = 0.d0
   e3bod = 0.d0

   call timer_start(TIME_BOND)
   if( ntf < 2 ) then

      !*****************************
      !          Bonds with H
      !*****************************

      call bond(nbonh,ix(iibh),ix(ijbh),ix(iicbh),x,xx,ix,f,ene(6),nocrst)
   end if

   if( ntf < 3 ) then

      !*****************************
      !          Bonds without H
      !*****************************

      call bond(nbona+nbper,ix(iiba),ix(ijba),ix(iicba),x,xx,ix,f,ene(7),nocrst)

   end if

   if( ntf < 4 ) then

      !***************************
      !          ANGLES with H
      !***************************

      call angl(ntheth,ix(i24),ix(i26),ix(i28),ix(i30),x,xx,ix,f,ene(8))
   end if

   if( ntf < 5 ) then

      !*****************************
      !          ANGLES without H
      !*****************************

      call angl(ntheta+ngper,ix(i32),ix(i34),ix(i36),ix(i38),x,xx,ix,f,ene(9))
   end if

   if( ntf < 6 ) then

      !******************************
      !          DIHEDRALS with H
      !******************************

      call ephi(nphih,ix(i40),ix(i42),ix(i44),ix(i46),ix(i48), &
         xx(l15),ix(i04),x,xx,ix,f,dvdl,ene(10),ene(11),ene(12),xx(l190))
   end if

   if( ntf < 7 ) then

      !*********************************
      !          DIHEDRALS without H
      !*********************************

      call ephi(nphia+ndper,ix(i50),ix(i52),ix(i54),ix(i56),ix(i58), &
         xx(l15),ix(i04),x,xx,ix,f,dvdl,ene(13),ene(14),ene(15),xx(l190))
   end if
   call timer_stop(TIME_BOND)

   !     ----- CALCULATE THE POSITION CONSTRAINT ENERGY -----

   if(natc > 0) then

      ! ---could be for ntr=1 (positional restraints) or targeted MD (itgtmd=1)

      if (ntr == 1) then
         call xconst(natc,entr,ix(icnstrgp),x,f,xx(lcrdr),xx(l60))
         ene(20) = entr
      end if

   end if

   if(ifcap > 0) call capwat(natom,x,f)

   !     ---- get the noesy volume penalty energy: ------

   ene(22) = 0.0

   !     -- when igb!=0 and igb!=10, all nonbonds are done in routine egb:

   esurf = 0.d0
   if( igb /= 0 .and. igb < 10) then

      call timer_start(TIME_EGB)
      call egb( x,xx,ix,f,rborn,fs,reff,xx(l15),ix(i04),ix(i06), &
            ix(i08),ix(i10),xx(l190), &
            cut,ntypes,natom,natbel,epol,eelt,evdw, &
            esurf,dvdl,xx(l165),ix(i82),xx(l170),xx(l175), &
            xx(l180),xx(l185), xx(l186),xx(l187),xx(l188),xx(l189))

      ene(2) = evdw
      ene(3) = eelt
      ene(4) = epol
      ene(23) = esurf
      ene(21) = dvdl
      call timer_stop(TIME_EGB)

   end if  ! ( igb /= 0 .and. igb < 10)
   

   !     -- when igb==10, all nonbonds are done in routine pb_force:
   !                      all nonpolar interactions are done in np_force:

#ifdef MPI
   if(mytaskid == 0)then
#endif
      if( igb >= 10 ) then

         call timer_start(TIME_PBFORCE)
         call pb_force(natom,nres,ntypes,ix(i02),ix(i04),ix(i06),ix(i10), &
                       cn1,cn2,xx(l15),x,f,evdw,eelt,epol)
         if ( pbgrid ) pbgrid = .false.
         if ( pbinit ) pbinit = .false.
         ene(2) = evdw
         ene(3) = eelt
         ene(4) = epol
         call timer_stop(TIME_PBFORCE)

         call timer_start(TIME_NPFORCE)
         esurf = 0.0d0; edisp = 0.0d0
         if ( ifcap == 0 .and. npopt /= 0 ) &
         call np_force(natom,nres,ntypes,ix(i02),ix(i04),ix(i06),&
                       cn1,cn2,x,f,esurf,edisp)
         if ( pbprint ) pbprint = .false.
         ene(23) = esurf
         ene(24) = edisp
         call timer_stop(TIME_NPFORCE)

      end if  ! ( igb >= 10 )
#ifdef MPI
   end if
#endif

#ifdef MPI
   call timer_barrier( commsander )
   call timer_start(TIME_COLLFRC)

   !     add force, ene, vir, npair, nhb copies from all nodes
   !            also add up newbalance for nonperiodic.

   call fdist(f,xx(lfrctmp),ene,vir,npair,nhb,r_stack,newbalance)
   call timer_stop(TIME_COLLFRC)
#endif

   ! ---- at this point, the parallel part of the force calculation is
   !      finished, and the forces have been distributed to their needed
   !      locations.  All forces below here are computed redundantly on
   !      all processors, and added into the force vector.  Hence, below
   !      is the place to put any component of the force calculation that
   !      has not (yet) been parallelized.

   !     ----- CALCULATE TOTAL ENERGY AND GROUP THE COMPONENTS -----


   do m = 2,15
      ene(1) = ene(1) + ene(m)
   end do
   ene(1) = ene(1) + epolar + e3bod + ene(23) + ene(24)

   ene(5) = ene(6)+ene(7)
   ene(6) = ene(8)+ene(9)
   ene(7) = ene(10)+ene(13)
   ene(8) = ene(11)+ene(14)
   ene(9) = ene(12)+ene(15)
   ene(10) = ene(17)+ene(20)+enmr(1)+enmr(2)+enmr(3)
   ene(1) = ene(1)+ene(10)

   !    Here is a summary of how the ene array is used.  For parallel runs,
   !    these values get summed then rebroadcast to all nodes (via
   !    mpi_allreduce).

   !    ene(1):    total energy
   !    ene(2):    van der Waals
   !    ene(3):    electrostatic energy
   !    ene(4):    10-12 (hb) energy, or GB/PB energy when igb.gt.0
   !    ene(5):    bond energy
   !    ene(6):    angle energy
   !    ene(7):    torsion angle energy
   !    ene(8):    1-4 nonbonds
   !    ene(9):    1-4 electrostatics
   !    ene(10):   constraint energy
   !    ene(11-19):  used a scratch, but not needed further below
   !    ene(20):   position constraint energy
   !    ene(21):   charging free energy result
   !    ene(22):   noe volume penalty
   !    ene(23):   surface-area dependent solvation energy or cavity energy
   !    ene(24):   surface-area dependent dispersion energy

   !     ----- TRANSFER THE ENERGIES TO THE ARRAY ENER, USED IN PRNTMD -----

   ener(1:10) = ene(1:10)
   ener(11) = epolar
   ener(12) = aveper
   ener(13) = aveind
   ener(14) = avetot
   ener(15) = ene(23)
   ener(16) = e3bod
   ener(17) = ene(21)
   ener(18) = ene(24)

   !     ----- IF BELLY IS ON THEN SET THE BELLY ATOM FORCES TO ZERO -----

   if (belly) call bellyf(natom,ix(ibellygp),f)
   call timer_stop(TIME_FORCE)
   return
end subroutine force
