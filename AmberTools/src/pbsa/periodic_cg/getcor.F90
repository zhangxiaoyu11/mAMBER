#include "dprec.h"
#include "copyright.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine getcor here]
subroutine getcor(nr,x,v,f,ntx,box,irest,tt)
   
   !     --- reads initial coords, vel, and box for MD.
   
   !     Rev A mods: amopen calls, EOF protection, XC() no longer
   !                 read from restart file.  improved error message.
   !                 output new box lengths.
   
   !     EWALD: read a modified form of the coordinate/restart file
   
   !      IMPLICIT _REAL_ (A-H,O-Z)
   
   implicit none
#  include "files.h"
   integer lun,nr,ntx,irest
   _REAL_ x(*),v(*),f(*),box(3),tt
   
   logical form
   _REAL_ fvar(7),a,b,c,alpha,beta,gamma
   integer i,ivar,ifld(7),ihol(1)
   integer natom,nr3,ier
   character(len=256) line_test
   
   
   nr3 = 3*nr
   lun = 9 !FIXME: use lun=-1 for new amopen()
   
   write(6,9108)
   
   !     ----- OPEN THE COORDINATE FILE -----
   
   form = (ntx == 1 .or. ntx == 5 .or. ntx == 7 )
   
   
   if(form)then
      !        subr amopen(lun,fname,fstat,fform,facc)
      call amopen(lun,inpcrd,'O','F','R')
      read(lun,9008) title1
      
      read(lun,'(a80)') line_test
      if( line_test(6:6) == ' ' ) then ! this is an old, i5 file
         read(line_test,9010) natom,tt
      else                   ! assume a new, i6 file
         read(line_test,9011) natom,tt
      end if
      9010 format(i5,e15.7)
      9011 format(i6,e15.7)
      
      if(natom /= nr) then
         write(6,9118)
         call mexit(6, 1)
      end if
      read(lun,9028) (x(i),i=1,nr3)
      
      if(ntx == 1) then
         do i = 1,nr3
            v(i) = 0.d0
         end do
         write(6,9008) title1
         write(6,9009) tt
         close(lun, iostat=ier)
         return
      end if
      
      !     ----- LOAD THE VELOCITY -----
      
      read(lun,9028,end=1010) (v(i),i=1,nr3)
      write(6,9008) title1
      write(6,9009) tt
      close(lun, iostat=ier)
      return
      
      !     ----- BINARY READING -----
      
   else
      call amopen(lun,inpcrd,'O','U','R')
      if(ntx == 2)then
         read(lun) title1
         read(lun) natom
         if(natom /= nr) then
            write(6,9118)
            call mexit(6, 1)
         end if
         read(lun,end=1000,err=1000) (x(i),i = 1,nr3)
         do i = 1,nr3
            v(i) = 0.d0
         end do
         write(6,9008) title1
         write(6,9009) tt
         close(lun, iostat=ier)
         return
      end if
      
      if(ntx == 3)then
         read(lun) title1
         read(lun) natom
         if(natom /= nr) then
            write(6,9118)
            call mexit(6, 1)
         end if
         read(lun,end=1000,err=1000) (x(i),i = 1,nr3)
         read(lun) (f(i),i = 1,nr3)
         write(6,9008) title1
         write(6,9009) tt
         close(lun, iostat=ier)
         return
      end if
      
      read(lun) title1
      read(lun) natom
      if(natom /= nr) then
         write(6,9118)
         call mexit(6, 1)
      end if
      read(lun,end=1000,err=1000) (x(i),i = 1,nr3)
      read(lun,end=1010,err=1010) (v(i),i = 1,nr3)
      if(ntx < 6) then
         write(6,9008) title1
         write(6,9009) tt
         close(lun, iostat=ier)
         return
      end if
      
      !     READ(lun,end=1020,err=1020) a,b,c,alpha,beta,gamma
      do 102 i=1,6
         ifld(i)=3
         fvar(i)=0.0d0
      102 continue
      ifld(7)=0
      ihol(1)=0
      call rfree(ifld,ihol,ivar,fvar,lun,6)
      if((fvar(4) /= 0).or.(fvar(5) /= 0).or.(fvar(6) /= 0)) then
         alpha=fvar(4)
         beta=fvar(5)
         gamma=fvar(6)
      else
         !     -- defaults:
         alpha=90.0d0
         beta=90.0d0
         gamma=90.0d0
      end if
      a=fvar(1)
      b=fvar(2)
      c=fvar(3)
      box(1) = a
      box(2) = b
      box(3) = c
      if ( alpha < 1.d0 ) then
         write(6,'(/,a)') 'EWALD: BAD BOX PARAMETERS in inpcrd!'
         call mexit(6, 1)
      end if
      !call fill_ucell(a,b,c,alpha,beta,gamma)
      write(6,9129) a,b,c,alpha,beta,gamma
   end if
   
   write(6,9008) title1
   write(6,9009) tt
   close(lun, iostat=ier)
   return
   1000 continue
   write(6,'(a,a)') 'FATAL: Could not read coords from ',inpcrd
   call mexit(6, 1)
   1010 continue
   write(6,'(a,a)') 'FATAL: Could not read velocities from ',inpcrd
   call mexit(6, 1)
   1020 continue
   write(6,'(a,a)') 'FATAL: Could not read BOX from ',inpcrd
   call mexit(6, 1)
   
   9008 format(a80)
   9009 format(t2,'begin time read from input coords =', &
         f10.3,' ps'/)
   9028 format(6f12.7)
   
   9108 format &
         (/80(1h-)/,'   3.  ATOMIC COORDINATES AND VELOCITIES',/80(1h-)/)
   9118 format(/2x,'FATAL: NATOM mismatch in coord and ', &
         'topology files')
   9128 format(t2,'NEW BOX DIMENSIONS from inpcrd file:', &
         /5x,'X =',f10.5,'  Y =',f10.5,'  Z =',f10.5,/)
   9129 format(t2,'NEW EWALD BOX PARAMETERS from inpcrd file:', &
         /5x,'A     =',f10.5,'  B    =',f10.5,'  C     =',f10.5,/, &
         /5x,'ALPHA =',f10.5,'  BETA =',f10.5,'  GAMMA =',f10.5,/)
end subroutine getcor 
