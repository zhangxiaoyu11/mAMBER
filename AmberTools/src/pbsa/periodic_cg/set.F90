#include "copyright.h"
#include "dprec.h"


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine setatm here]
subroutine setatm(nat,nres,natb,ipres,igrp,igres)
   implicit _REAL_ (a-h,o-z)
   logical active
   dimension ipres(*),igrp(*),igres(*)
   
   idum = 0
   do 100 i = 1,nat
      if(igrp(i) <= 0) goto 100
      idum = idum+1
   100 continue
   natb = idum
   
   do 120 i = 1,nres
      i1 = ipres(i)
      i2 = ipres(i+1)-1
      igres(i) = 0
      active = .false.
      do 140 j = i1,i2
         if(igrp(j) <= 0) goto 140
         active = .true.
      140 continue
      if(active) igres(i) = 1
   120 continue
   return
end subroutine setatm 
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine setbon here]
subroutine setbon(nb,mb,ib,jb,icb,igrp)
   
   implicit _REAL_ (a-h,o-z)
   dimension ib(*),jb(*),icb(*),igrp(*)
   
   nba = 0
   nca = 0
   do 100 i = 1,nb
      iat = ib(i)/3+1
      jat = jb(i)/3+1
      if(igrp(iat) <= 0 .and. igrp(jat) <= 0) goto 100
      if(i > mb) nca = nca+1
      nba = nba+1
      ib(nba) = ib(i)
      jb(nba) = jb(i)
      icb(nba) = icb(i)
   100 continue
   nb = nba
   mb = nba-nca
   return
end subroutine setbon 
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine setang here]
subroutine setang(nt,mt,it,jt,kt,ict,igrp)
   
   implicit _REAL_ (a-h,o-z)
   dimension it(*),jt(*),kt(*),ict(*),igrp(*)
   
   nta = 0
   mta = 0
   do 100 i = 1,nt
      iat = it(i)/3+1
      jat = jt(i)/3+1
      kat = kt(i)/3+1
      if(igrp(iat) <= 0 .and. igrp(jat) <= 0 .and. igrp(kat) <= 0) &
            goto 100
      if(i > mt) mta = mta+1
      nta = nta+1
      it(nta) = it(i)
      jt(nta) = jt(i)
      kt(nta) = kt(i)
      ict(nta) = ict(i)
   100 continue
   nt = nta
   mt = nta-mta
   return
end subroutine setang 
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine setdih here]
subroutine setdih(np,mp,ip,jp,kp,lp,icp,igrp)
   
   implicit _REAL_ (a-h,o-z)
   dimension ip(*),jp(*),kp(*),lp(*),icp(*),igrp(*)
   
   npa = 0
   mpa = 0
   do 100 i = 1,np
      iat = ip(i)/3+1
      jat = jp(i)/3+1
      kat = iabs(kp(i))/3+1
      lat = iabs(lp(i))/3+1
      if(igrp(iat) <= 0 .and. igrp(jat) <= 0 .and. igrp(kat) <= 0 &
            .and. igrp(lat) <= 0) goto 100
      if(i > mp) mpa = mpa+1
      npa = npa+1
      ip(npa) = ip(i)
      jp(npa) = jp(i)
      kp(npa) = kp(i)
      lp(npa) = lp(i)
      icp(npa) = icp(i)
   100 continue
   np = npa
   mp = npa-mpa
   return
end subroutine setdih 
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine rdrest here]
subroutine rdrest(nr,ntx,ntrx,xc)
   
   implicit _REAL_ (a-h,o-z)
   
   !     ----- ROUTINE TO READ THE REFERENCE POSITIONS FOR RESTRAINING ----
   
#  include "files.h"
   character(len=80) line
   
   dimension xc(*)
   
   nr3 = 3*nr
   write(6,191)
   
   if(ntrx /= 0) goto 185
   
   read(10) title1
   read(10) natom
   if(natom /= nr) goto 300
   read(10) (xc(i3),i3 = 1,nr3)
   goto 190
   
   185 read(10,40) title1
   read(10,'(a)') line
   if( line(6:6) == ' ' ) then   ! this is an old, i5 file
      read(line,'(i5)') natom
   else                          ! assume a new, i6 file
      read(line,'(i6)') natom
   end if
   if(natom < nr) goto 300
   read(10,92) (xc(i),i = 1,nr3)
   190 continue
   write(6,31) title1
   return
   300 continue
   write(6,1000)
   call mexit(6, 1)
   191 format(/,'   5.  REFERENCE ATOM COORDINATES',/)
   31 format(2x,20a4)
   40 format(20a4)
   92 format(6f12.7)
   94 format(10i5)
   1000 format(/2x,'FATAL: NATOM mismatch in constraint coord', &
         ' and topology files')
end subroutine rdrest 
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine dihdup here]
subroutine dihdup(nphia,ip,jp,kp,lp,icp,pn)
   
   !     duplicates pointers to multi-term dihedrals for vector ephi.
   !     H-atom diheds are duplicated at the end of the h-atom lists,
   !     but heavy atom diheds must be inserted between the end of the
   !     heavy-atom diheds, but before the constraint diheds if they are
   !     present.  In order to use this code, extra space MUST be allocated
   !     in LOCMEM
   
   !     (Note: this is only for ancient prmtop files: anything created by
   !     LEaP should report no duplications)
   
   !     Author:  George Seibel
   
   implicit _REAL_ (a-h,o-z)
   
   !     COMMON:
   
#  include "md.h"
   
   !     INPUT:
   
   integer nphia
   !        ... number of dihedrals
   integer ip(*), jp(*), kp(*), lp(*)
   !        ... pointers to atoms of dihedrals
   integer icp(*)
   !        ... pointers to dihedral parameters
   dimension pn(*)
   !        ... periodicity; negative if multi-term, read until + encountered
   
   !     INTERNAL:
   
   integer ndup
   !        ... number of diheds duplicated
   logical const
   !        ... const flag: true if constraints are present
   integer ic
   !        ... working parameter pointer
   integer i
   !        ... do loop index
   
   ndup = 0
   ierr = 0
   
   do i = 1, nphia
      ic = icp(i)
      100 continue
      if (pn(ic) < 0) then
         ndup = ndup + 1
         if (ndup > maxdup) then
            ierr = 1
         else
            
            !            --- duplicate pointers of multi-term dihedrals ---
            
            ip(nphia+ndup)  = ip(i)
            jp(nphia+ndup)  = jp(i)
            kp(nphia+ndup)  = kp(i)
            lp(nphia+ndup)  = lp(i)
            
            !            --- but use the NEXT parameter pointer ---
            
            icp(nphia+ndup) = ic + 1
         end if
         
         !            --- check for a third or higher term ---
         
         ic = ic + 1
         goto 100
      end if
   end do
   
   if (ierr == 1) then
      write(6,'(/,5x,a,i5,a)') 'MAXDUP =',maxdup,' exceeded'
      write(6,'(/,5x,a,i5,a)') 'set MAXDUP =',ndup,' in locmem.f'
      call mexit(6, 1)
   end if
   
   nphia = nphia + ndup
   write(6,'(a,i5,a)') '| Duplicated',ndup,' dihedrals'
   return
end subroutine dihdup 
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine dihpar here]
subroutine dihpar(numphi,pk,pn,phase,gamc,gams,ipn,fmn)
   
   implicit _REAL_ (a-h,o-z)
   
   !     ----- ROUTINE TO GET ADDITIONAL PARAMETERS FOR THE
   !           VECTORISED EPHI -----
   
   dimension pk(*),pn(*),phase(*),gamc(*),gams(*),ipn(*)
   dimension fmn(*)
   
   data zero,one,tenm3,tenm10/0.0d+00,1.0d+00,1.0d-03,1.0d-06/
   data four,pi/4.0d+00,3.141592653589793d+00/
   
   pim = four*atan(one)
   do 100 i = 1,numphi
      dum = phase(i)
      if(abs(dum-pi) <= tenm3) dum = sign(pim,dum)
      dumc = cos(dum)
      dums = sin(dum)
      if(abs(dumc) <= tenm10) dumc = zero
      if(abs(dums) <= tenm10) dums = zero
      gamc(i) = dumc*pk(i)
      gams(i) = dums*pk(i)
      fmn(i) = one
      if(pn(i) <= zero) fmn(i) = zero
      pn(i) = abs(pn(i))
      ipn(i) = int(pn(i)+tenm3)
   100 continue
   return
end subroutine dihpar 
!-----------------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine setgms here]
subroutine setgms(numphi)
   
   ! Subroutine SET GaMS.
   
   ! This routine sets the values of the GAMC() and GAMS() arrays, which
   ! are used in vectorized torsional energy routines. This routine is only
   ! called when the torsional force constants are changed in routine MODWT.
   ! Otherwise, these arrays are set only once, by a call to DIHPAR from RDPARM,
   ! at the start of the program. This is an abbreviated version of DIHPAR which
   ! passes most arguments by common, not call-list.
   
   ! Author: David A. Pearlman
   ! Date: 5/92
   
   implicit _REAL_ (a-h,o-z)
   
#  include "parms.h"
   
   data zero,one,tenm3,tenm10/0.0d+00,1.0d+00,1.0d-03,1.0d-06/
   data four,pi/4.0d+00,3.141592653589793d+00/
   
   pim = four*atan(one)
   do 100 i = 1,numphi
      dum = phase(i)
      if(abs(dum-pi) <= tenm3) dum = sign(pim,dum)
      dumc = cos(dum)
      dums = sin(dum)
      if(abs(dumc) <= tenm10) dumc = zero
      if(abs(dums) <= tenm10) dums = zero
      gamc(i) = dumc*pk(i)
      gams(i) = dums*pk(i)
   100 continue
   return
end subroutine setgms 
#ifdef MPI

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ DISTRIBUTE ATOMS TO PROCESSORS USING RESIDUE BOUNDARIES
subroutine setpar(ipres,amass)
   
   dimension ipres(*)
   _REAL_ amass(*)

#  include "memory.h"
#  include "parallel.h"
#ifdef MPI_DOUBLE_PRECISION
#undef MPI_DOUBLE_PRECISION
#endif
#  include "mpif.h"
#  include "extra.h"
   integer target,fraction,i,iat
   integer thismol, thismolsfirstnode, thismolslastnode
   integer thismolslastatom, group_world, group_temp

   if (mpi_orig) then
      
      !   when only master will be running the runmd,runmin routines:
      
      iparpt(0) = 0
      iparpt(1) = natom

   else
      
      !   all nodes will divide the work in runmd,runmin:
      
      if (nres < numtasks) then
         write(6,*) 'Must have more residues than processors!'
         call mexit(6,1)
      end if
      fraction = natom/numtasks
      target = 0
      iparpt(0) = 0
      iparpt(numtasks) = natom
      do 20 node=1,numtasks-1
         target = target + fraction
         do 10 ires=1,nres

            !   --- do not stop before or after a single-atom residue... This is a
            !       check to avoid the potential SHAKE failure by the parallel
            !       residue-based distribution: if there is a residue with a single
            !       atom at the parallel distribution boundary, and if this is a hydrogen
            !       which may be SHAKEn, this is an error since the hydrogen may be on
            !       a different processor than the residue to which it is connected.
            !       Granted, this is a fairly rare case since most of the available
            !       force fields do not use this peculiar and particularly hacked 
            !       topology construction to terminate a chain however the check for
            !       single atom residues goes way back to amber 4.1.  To allow for 
            !       simulations of LJ fluids, we expand this check to not only look
            !       for single atom residues at the boundary but also check the
            !       mass of the atom.  This is not fully rigorous (i.e. hydrogen
            !       like SHAKE-able atoms with a mass of > 3.0 will sneak by this
            !       check) however this *should* work in standard usage.
            !
            !       A more descriptive error message is generated than previously.
            
            if ((ipres(ires+1) - ipres(ires)) < 2) then
               iat = 0
               do i=1,ires
                  iat = iat + ipres(i+1) - ipres(i)
               enddo
               if (amass(iat) < 3.0) goto 10
            endif

            if ( ires > 1 )then
               if( (ipres(ires) - ipres(ires-1)) < 2 ) then
                  iat = 0
                  do i=1,ires-1
                     iat = iat + ipres(i+1) - ipres(i)
                  enddo
                  if (amass(iat) < 3.0) goto 10
               endif
            endif
 
            !   --- if this residue starts past the target atom, end the atom list
            !       at the end of the previous residue:
            
            if (ipres(ires) > target) then
               iparpt(node) = ipres(ires) - 1
               goto 20
            end if
         10 continue

         write(6,*) 'ERROR IN SETPAR() upon atom distribution for parallel runs.'
         write(6,*) 'Note that the atom distribution for the parallel implementation'
         write(6,*) 'is done on a residue basis.  A residue containing a single atom has'
         write(6,*) 'been detected at the boundary between processors.  As the mass of'
         write(6,*) 'this atom marks it as a SHAKE-able hydrogen, this is flagged as an'
         write(6,*) 'error since SHAKE across processor boundaries is not supported at'
         write(6,*) 'processor boundaries.  To get around this error, either run a '
         write(6,*) 'sequential run or modify the code in set.f, routine setpar().'
         write(6,*) 'Aborting!!!'
         call mexit(6,1)

      20 continue
      
      !       --- following used when forces, coords. communicated:
      
      do node=0,numtasks
         iparpt3(node) = 3*iparpt(node)
      end do
      
      !  --- include "extra" dynamical variables for final processor:
      
      iparpt3(numtasks) = iparpt3(numtasks) + iscale
      
      do i=0,numtasks-1
         rcvcnt(i) = iparpt(i+1) - iparpt(i)
         rcvcnt3(i) = iparpt3(i+1) - iparpt3(i)
      end do
      if (master) then
         write(6,'(a)') '|  Atom division among processors:'
         write(6,'(3H|  8i8)') (iparpt(j),j=0,numtasks)
      end if
   end if  !  10 ires=1,nres
   return
end subroutine setpar 
#endif
