! <compile=optimized>
#include "copyright.h"
#include "dprec.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine shake here]
subroutine shake(nrr,nbonh,nbon,nbper,ib,jb,igrp,winv,conp, &
      skip,x,xp,niter,belly,ifstwt &
#ifdef PSANDER
      ,shkh,qspatial &
#endif
      )
   
   
   !     Mods for version 4.1:
   !     - Add IFSTWT(I), so that waters shaken by fast analytic 3-point
   !       shake will not be shake-n here -- dap
   !     Mods for Rev A by gls:
   !     - this is now library routine, calling args changed in gibbs.
   !     - cpp-selectable single/double prec.
   !     - removed unreferenced statement labels
   !------------------------------------------------------------------------
   implicit none
   logical skip(*),ready,first
   logical belly
   
   !     ----- ALL THE BONDS INVOLVING HYDROGEN ATOMS ARE LOADED FIRST
   !           IN ARRAYS IB AND JB FOLLOWED BY THOSE INVOLVING NON-
   !           HYDROGEN ATOMS AND THE PERTURBED ATOMS -----
   
#ifdef MPI
   !     =========================== AMBER/MPI ===========================
   
   !     ...the only MPI code specific to this routine is wrappers for
   !     all the I/O to allow only the master to write...
   
#  include "parallel.h"
   
   !     ========================= END AMBER/MPI =========================
#endif
#  include "md.h"
#  include "extra.h"
   
   integer nrr,nbonh,nbon,nbper
   _REAL_ winv(*)
   integer ib(*),jb(*),igrp(*)
   _REAL_ conp(*),x(*),xp(*)
   _REAL_ xij(3),xpij(3)
   integer ifstwt(*)
   logical qspatial


   _REAL_ winvi,winvj,diff,rrpr,acor,xh,toler,rpij2, &
         tol2,zero
   integer m,i,j,i3,j3,istart,iend,niter,nbt,nit,k,ll
   integer nbt_new
   data first /.true./
   data zero /0.0d0/
   save first,zero,nbt_new
#ifdef PSANDER
   integer shkh(nbonh),llind
#else
   qspatial=.false.
#endif
   
   if(qspatial)then
      if(ntc == 1)then
         niter=1
         return
      end if
   else
      if(ntc == 1)then
         niter=1
         return
      else if(ntc == 2)then
         nbt = nbonh
      else
         nbt = nbonh+nbon+nbper
      end if
      
      !     first time through, figure out last bond that is not a 3-point water
      !     bond:
      
      if (first) then
         nbt_new = 0
         do i=1,nbt
            if(ifstwt(i) == 0 ) nbt_new = i
         end do
         first= .false.
      end if
      if( nbt_new == 0 ) then
         niter=1
         return
      end if
   end if

   niter = 0
   tol2 = tol
   
   nit = 0
   do k = 1,nrr
      skip(k) = .true.
      skip(nrr+k) = .false.
   end do
   ready = .false.
   !         print *,mytaskid," SHAKE bonds ",nbonh
   
   !===========================================================================
   !        BIG ITERATION LOOP
   
   do while(.not.ready)
      if(nit > 3000) then
         write(6,311)
         write(6,322)
         if (imin == 1) write(6,*) &
               ' *** Especially for minimization,', &
               ' try ntc=1 (no shake)'
         call mexit(6,1)
      end if

      ready = .true.
#ifdef MPI
      
      !     ---- divide atoms up among the processors, always splitting on
      !     residue boundaries:
      
      if (ntc /= 2) then
         write(6,*) 'this parallel version only works for ntc < 3'
         write(6,*) '      only fast water and bonds with H.'
         call mexit(6,1)
      end if
      istart = iparpt(mytaskid) + 1
      iend = iparpt(mytaskid+1)
#endif
      
      !     ----- LOOP OVER ALL THE BONDS that are not 3-point waters: -----
      
#ifdef PSANDER
      do 100 llind = 1,nbonh
         ll=shkh(llind)
#else
      do 100 ll = 1,nbt_new
         
         !     Skip bonds constrained in fast 3-point routine:
         
         if (ifstwt(ll) == 1) goto 100
         
#endif
         i3 = ib(ll)
         i  = i3/3+1
         j3 = jb(ll)
         j  = j3/3+1
         if(skip(nrr+i).and.skip(nrr+j)) goto 100
#ifdef MPI
# ifndef PSANDER
         
         !     --- skip constrained bonds not destined for this processor:
         
         if( .not.mpi_orig ) then
            if (i < istart .or. i > iend) goto 100
            if (j < istart .or. j > iend) then
               write(6,*) 'partition error in shake'
               call mexit(6,1)
            end if
         end if
# endif
#endif
         
         toler = conp(ll)
         
         !     ----- IF BELLY OPTION IS ON THEN THE RESETTING OF THE FROZEN ATOM
         !              IS TO BE PREVENTED -----
         
         winvi = winv(i)
         winvj = winv(j)
         if(belly) then
            if(igrp(i) <= 0) winvi = zero
            if(igrp(j) <= 0) winvj = zero
         end if
         
         !     --- calc nominal distance squared
         
         diff = toler
         rpij2 = zero
         do m = 1,3
            xpij(m) = xp(i3+m)-xp(j3+m)
            rpij2 = rpij2+xpij(m)**2
         end do
         
         !     ----- APPLY THE CORRECTION -----
         
         diff = diff-rpij2
         if( abs(diff) < toler*tol2) goto 100
         do m = 1,3
            xij(m) = x(i3+m)-x(j3+m)
         end do
         
         !     ----- SHAKE RESETTING OF COORDINATE IS DONE HERE -----
         
         rrpr = zero
         do m = 1,3
            rrpr = rrpr+xij(m)*xpij(m)
         end do
         if(rrpr < toler*1.0d-06) then
            write(6,321) niter,nit,ll,i,j
            write(6,322)
            if (imin == 1) write(6,*) &
                  ' *** Especially for minimization,', &
                  ' try ntc=1 (no shake)'
            call mexit(6,1)
         end if
         acor = diff/(rrpr*(winv(i)+winv(j)+winv(i)+winv(j)))
         do m = 1,3
            xh = xij(m)*acor
            xp(i3+m) = xp(i3+m)+xh*winvi
            xp(j3+m) = xp(j3+m)-xh*winvj
         end do
         skip(i) = .false.
         skip(j) = .false.
         ready = .false.
      100 continue
      
      nit = nit+1
      do k = 1,nrr
         skip(nrr+k) = skip(k)
         skip(k) = .true.
      end do
      
   end do  !  100 llind = 1,nbonh
   !          END OF ITERATION LOOP
   !==================================================================
   niter = niter+nit
   return

   311 format(/5x,'Coordinate resetting (SHAKE) was not accomplished', &
         /5x,'within 3000 iterations')
   321 format(/5x,'Coordinate resetting (SHAKE) cannot be accomplished,', &
         /5x,'deviation is too large', &
         /5x,'NITER, NIT, LL, I and J are :',7i5)
   322 format(/5x,'Note: This is usually a symptom of some deeper', &
         /5x,'problem with the energetics of the system.' )
end subroutine shake 
