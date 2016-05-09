#include "copyright.h"
#ifdef MPI
#include "dprec.h"

!     The AMBER/MPI implementation and support routines were
!     originally and independently implemented and contributed
!     by James Vincent (JV) 7/94.  Modified by tec3, dac and JV.


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ broadcast data from master to all other nodes, at beginning
subroutine startup(xx,ix,ih)
   !************************************************************
   !     Send data needed by all nodes once at startup from master
   !     after master has read in all data
   !************************************************************

   implicit none
#  include "parallel.h"
#ifdef MPI_DOUBLE_PRECISION
#undef MPI_DOUBLE_PRECISION
#endif
#  include "mpif.h"
#  include "extra.h"
#  include "md.h"
#  include "memory.h"
#  include "parms.h"
#  include "box.h"
#  include "files.h"
#  include "nonper.h"
#  include "flocntrl.h"
#  include "debug.h"
#  include "rstack.h"
#  include "istack.h"
#  include "new_time.h"
   
   _REAL_ xx(*)
   integer ix(*), ier
   character(len=4) ih(*)
   
   !     Send and receive common blocks from the master node:
   
   !  files.h:
   
   call mpi_bcast(ntpr,BC_HULP,MPI_INTEGER,0,commsander,ierr)
   call mpi_bcast(iredir,8,MPI_INTEGER,0,commsander,ierr)
   
   !  md.h:
   
   call mpi_bcast(nrp,BC_MDI,MPI_INTEGER,0,commsander,ierr)
   call mpi_bcast(t,BC_MDR,MPI_DOUBLE_PRECISION,0,commsander,ierr)
   
   !  extra.h:
   
   call mpi_bcast(ilbopt,6,MPI_INTEGER,0,commsander,ierr)
   call mpi_bcast(lbwght,1,MPI_DOUBLE_PRECISION,0,commsander,ierr)
   
   !  box.h:
   
   call mpi_bcast(ntb,BC_BOXI,MPI_INTEGER,0,commsander,ierr)
   call mpi_bcast(box,BC_BOXR,MPI_DOUBLE_PRECISION,0,commsander,ierr)
   
   !  parms.h:
   
   call mpi_bcast(rk,BC_PARMR,MPI_DOUBLE_PRECISION,0,commsander,ierr)
   call mpi_bcast(ipn,BC_PARMI,MPI_INTEGER,0,commsander,ierr)
   call mpi_bcast(numbnd,5,MPI_INTEGER,0,commsander,ierr)
   
   !  nonper.h
   
   call mpi_bcast(xbox0,BC_NONPER,MPI_DOUBLE_PRECISION,0,commsander,ier)
   
   !  local_nb.h
   !     common/dirpars/
   !     common/nb_bound/
   !     common/nb_float/
   !     common/nb_integ/
   
   call mpi_bcast(numnptrs,BC_DIRPARS,MPI_INTEGER,0,commsander,ier)
   call mpi_bcast(maxnblst,BC_NB_BOUND,MPI_INTEGER,0,commsander,ier)
   call mpi_bcast(limgcrds,BC_NBFLOAT,MPI_INTEGER,0,commsander,ier)
   call mpi_bcast(inumatg,BC_NBINTEG,MPI_INTEGER,0,commsander,ier)
   
   ! flocntrl.h
   
   call mpi_bcast(do_dir,BC_FLOCNTRL,MPI_INTEGER,0,commsander,ier)
   
   ! debug.h
   
   call mpi_bcast(do_debugf,BC_DEBUG,MPI_INTEGER,0,commsander,ier)
   call mpi_bcast(lscg,BC_DEB_HEAP,MPI_INTEGER,0,commsander,ier)
   
   ! timer info new_time.h
   
   call mpi_bcast(tpar_p,BC_TIME_PAR,MPI_INTEGER,0,commsander,ier)
   
   ! extra points
   
   call mpi_bcast(ifrtyp,BC_EXTRA_PT,MPI_INTEGER,0,commsander,ier)
   
   !  rstack.h,istack.h
   
   call mpi_bcast(bot_stk,5,MPI_INTEGER,0,commsander,ier)
   call mpi_bcast(ibot_stk,5,MPI_INTEGER,0,commsander,ier)
   
   !     IX,XX,IH
   
   call mpi_bcast(ix(1),lasti,MPI_INTEGER,0,commsander,ierr)
   call mpi_bcast(ih(1),lasth,MPI_INTEGER,0,commsander,ierr)
   call mpi_bcast(xx(1),lastr,MPI_DOUBLE_PRECISION,0,commsander,ierr)
   
   call mpi_barrier(commsander,ierr)
   
   !   ---- divide atoms up among the processors, always splitting on
   !        residue boundaries:
   
   call setpar(ix(i02),xx(lmass))
   return
end subroutine startup 
!----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fdist here]
subroutine fdist(f,forcetmp,ene,vir,npair,nhb,r_stack,newbalance)
   !************************************************************
   !   for the "original version", (when mpi_orig is set):
   
   !     James Vincent 7/94
   !     Gather all copies of the force array with ene and vir
   !     arrays tacked onto the end, then extract ene and vir
   !     arrays back out.
   !     Input array f is current local PE copy, forcetmp is
   !     scratch space, but results are put back into f, thus
   !     overwriting what was there.
   !     f: final force array - result of reduce operation
   !     forcetmp: scratch space
   !     ene: final ene array with sum of all ene values
   !     vir: final vir array
   !     npair: final npair value
   !     nhb: final nhb values
   
   !   when mpi_orig is false, does a distributed sum of the forces,
   !     so that each node ends up with only a piece of the total
   !     force array; does a global sum on the energy and virial
   
   !************************************************************

   implicit none
#  include "memory.h"
#  include "parallel.h"
#ifdef MPI_DOUBLE_PRECISION
#undef MPI_DOUBLE_PRECISION
#endif
#  include "mpif.h"
#  include "md.h"
#  include "extra.h"
   
   !     Parameters:
   
   _REAL_ f(*),forcetmp(*),ene(*),vir(*)
   _REAL_ r_stack(*)
   integer npair,nhb,newbalance
   
   !     Local:
   
   integer i,j
   
   if (numtasks == 1) return
   
   !     Tack ene, vir, and npair, nhb onto end of f:
   
   j = 3*natom+iscale+1
   f(j)   = vir(1)
   f(j+1)   = vir(2)
   f(j+2)   = vir(3)
   do i = 2,30
      f(j+1+i) = ene(i)
   end do
   f(j+32) = npair
   f(j+33) = nhb
   f(j+34) = newbalance
   
   if (mpi_orig) then
      
      !       ---Add all copies of force array together from all nodes:
      
      call mpi_reduce(f,forcetmp,(3*natom+iscale+35), &
            MPI_DOUBLE_PRECISION,mpi_sum,0,commsander,ierr)
      vir(1) = forcetmp(j)
      vir(2) = forcetmp(j+1)
      vir(3) = forcetmp(j+2)
      do i = 2,30
         ene(i) = forcetmp(j+i+1)
      end do
      do i=1,3*natom
         f(i) = forcetmp(i)
      end do
   else
      
      !       ---Add all copies of virial and energy and put result back on ALL nodes:
      
      call mpi_allreduce(f(j),forcetmp(j),35, &
            MPI_DOUBLE_PRECISION,mpi_sum,commsander,ierr)
      vir(1) = forcetmp(j)
      vir(2) = forcetmp(j+1)
      vir(3) = forcetmp(j+2)
      do i = 2,30
         ene(i) = forcetmp(j+i+1)
      end do
      newbalance=0
      if(forcetmp(j+34) > 0.d0)newbalance=1
      
      !       --- now do a distributed sum of the force array:
      
      !       ...due to lack of parallelization in the initial parts
      !       of runmd in the init=3 case, the more efficient
      !       reduce_scatter needs to be replaced with an mpi_allreduce call
      
      if (init /= 3) then
         
         call fsum(f,forcetmp)
         
         !   use following if init=3 really needs to be used:
         
      else
         call mpi_allreduce(f, forcetmp, 3*natom, &
               MPI_DOUBLE_PRECISION,mpi_sum,commsander,ierr)
         do i=1, 3*natom
            f(i) = forcetmp(i)
         end do
      end if
   end if  ! (mpi_orig)
   
   478 format(t2,'NB-update: NPAIRS =',i8,'  HBPAIR =',i8)
   return
end subroutine fdist 
!=======================================================================

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fsum here]
subroutine fsum(f,tmp)
   
   !     equivalent to an MPI_REDUCE_SCATTER on f:  all processors contribute
   !       to f, and the appropriate part of the result winds up on each
   !       processor
   

   implicit none
   _REAL_ f(*),tmp(*)
   integer i

   
#  include "parallel.h"
#  include "extra.h"
#  include "memory.h"
#ifdef MPI_DOUBLE_PRECISION
#  undef MPI_DOUBLE_PRECISION
#endif
#  include "mpif.h"
#ifdef CRAY_PVP
#  define MPI_DOUBLE_PRECISION MPI_REAL8
#endif
   
#ifndef noBTREE
   integer other,ncyclesm1,k,bit,cs,cr,ns,nr
   integer ist(mpi_status_size)
#endif
   
#ifndef noBTREE
   if (numtasks <= 1) return
   ncyclesm1 = logtwo(numtasks) - 1
   bit=numtasks/2
   
   do k = 0,ncyclesm1
      
      other=ieor(mytaskid,bit)
      
      !        send chunk:
      
      cs = ishft(other,-((ncyclesm1)-k))*bit
      ns = iparpt3(cs+bit)-iparpt3(cs)
      
      !        recv chunk:
      
      cr = ishft(mytaskid,-((ncyclesm1)-k))*bit
      nr = iparpt3(cr+bit)-iparpt3(cr)
      
      
      call mpi_sendrecv( &
            f(iparpt3(cs)+1),ns,MPI_DOUBLE_PRECISION,other,5, &
            tmp(iparpt3(cr)+1),nr,MPI_DOUBLE_PRECISION,other,5, &
            commsander, ist, ierr )
      do i=iparpt3(cr)+1,iparpt3(cr+bit)
         f(i) = f(i) + tmp(i)
      end do
      
      bit = ishft(bit,-1)
      
   end do  !  k = 0,ncyclesm1
#else
   
   call mpi_reduce_scatter(f, tmp(iparpt3(mytaskid)+1), &
         rcvcnt3, MPI_DOUBLE_PRECISION, mpi_sum, &
         commsander, ierr)
   do i=1, 3*natom
      f(i) = tmp(i)
   end do
#endif
   return
end subroutine fsum 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Distribute the coordinates to all processors.
subroutine xdist(x)

   implicit none
   
   _REAL_ x(*)
   
#  include "parallel.h"
#    ifdef MPI_DOUBLE_PRECISION
#      undef MPI_DOUBLE_PRECISION
#    endif
#  include "mpif.h"
#    ifdef CRAY_PVP
#      define MPI_DOUBLE_PRECISION MPI_REAL8
#    endif
   
#  ifndef noBTREE
   integer other,ncyclesm1,k,bit,cs,cr,ns,nr
   integer ist(mpi_status_size),ireq
!  logical sendfirst
#endif
   
#  ifndef noBTREE
   
   if (numtasks <= 1) return
   ncyclesm1 = logtwo(numtasks) - 1
   bit=1
   do k = 0,ncyclesm1
      other=ieor(mytaskid,bit)
      cs = ishft(mytaskid,-k)*bit
      cr = ishft(other,-k)*bit
      ns = iparpt3(cs+bit)-iparpt3(cs)
      nr = iparpt3(cr+bit)-iparpt3(cr)
!     sendfirst = iand(mytaskid,bit).gt.0
      
#if 1
      call mpi_sendrecv( &
            x(iparpt3(cs)+1),ns,MPI_DOUBLE_PRECISION,other,5, &
            x(iparpt3(cr)+1),nr,MPI_DOUBLE_PRECISION,other,5, &
            commsander, ist, ierr )
#else
!        if (sendfirst) then 
!           call MPI_SEND( x(iparpt3(cs)+1),ns,MPI_DOUBLE_PRECISION,  &
!               other,25,commsander,ierr)
!           call MPI_RECV( x(iparpt3(cr)+1),nr,  &
!             MPI_DOUBLE_PRECISION,other,25,commsander,ist,ierr )
!        else
            call MPI_iRECV( x(iparpt3(cr)+1),nr,  &
              MPI_DOUBLE_PRECISION,other,25,commsander,ireq,ierr )
            call MPI_SEND( x(iparpt3(cs)+1),ns,MPI_DOUBLE_PRECISION,  &
                other,25,commsander,ierr)
            call mpi_wait( ireq, ist, ierr )
!        endif
#endif
      
      bit = ishft(bit,1)
   end do
#  else
   
   !       --- Assume an "in-place" allgatherv works: this seems(?) to
   !           be true everywhere....
   
   call mpi_allgatherv( &
         x(iparpt3(mytaskid)+1),rcvcnt3(mytaskid), &
         MPI_DOUBLE_PRECISION,x,rcvcnt3,iparpt3, &
         MPI_DOUBLE_PRECISION,commsander, ierr)
#  endif
   return
end subroutine xdist 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fgblsum here]
subroutine fgblsum(x,tmp)
   implicit none
   _REAL_ x(*),tmp(*)
   
   call fsum(x,tmp)
   call xdist(x)
   
   return
end subroutine fgblsum 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ puts two integers into a real array, two ints fit into one real
subroutine setb(b,nd,nz)
   implicit none
   integer b(2),nd,nz
   b(1)=nd
   b(2)=nz
   return
end subroutine setb 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ gets two integers from a real array, two ints fit into one real
subroutine getb(b,nd,nz)
   implicit none
   integer b(2),nd,nz
   nd=b(1)
   nz=b(2)
   return
end subroutine getb 
#else
subroutine dummy_parallel()
end subroutine dummy_parallel
#endif  /* MPI  */
