! <compile=optimized>
#include "copyright.h"

!--------------------------------------------------------------
! CHAMBER: SORT2 module
!
! Description: Module containing routines for sorting
!              arrays.
!
! Authors: Mark Williamson (SDSC, 2009)
!          Michael Crowley (NREL, 2009)
!          Ross Walker (SDSC, 2009)
!
!--------------------------------------------------------------

module sort2

contains
!============================================================================
!
!                  SORTER
!
!----------------------------------------------------------------------------
subroutine sorter(n,array0,array1)
   use psfprm, only: allint
   implicit none
   integer, intent(in) :: n
   integer, dimension(n),intent(inout) :: array0,array1
   integer, dimension(:),pointer ::   tmparr
   integer  :: nd,nd1,l1,l2,n0,n1,top, index
   logical :: flag
   real(kind=8) randnum

   if(n==0)return
   call allint(int(1+dlog(real(n,8))*4.0d0),tmparr)
   top=0
   n0=1
   nd=n

 3333 do while(nd.eq.n0) 
         if(top.eq.0) then
           deallocate(tmparr)
           return
         endif
         nd=tmparr(top)
         n0=tmparr(top-1)
         top=top-2
      enddo
      if(nd  ==  n0+1) then
         if (.not. correct_order(n0,nd,array0,array1)) then
            call exchange(n0,nd,array0,array1)
         end if
      endif

! process-partition-sort2
      call random_number(randnum)
      index=int((nd-n0)*randnum)+n0
      n1=n0-1
      nd1=nd+1

      do while(n1.lt.nd1)
         flag=.true.
         do while(flag)
            n1=n1+1
            flag=correct_order(n1,index,array0,array1) .and. (n1.ne.nd)
         enddo
         
         flag=.true.
         do while(flag)
            nd1=nd1-1
            flag=correct_order(index,nd1,array0,array1) .and. (nd1 /= n0) 
         enddo
         
         if (n1.lt.nd1) call exchange(n1,nd1,array0,array1)
      enddo
      
      if (n1.lt.index) then
         call exchange(n1,index,array0,array1)
      else
         if (index.lt.nd1) then
            call exchange(nd1,index,array0,array1)
         end if
      end if
!
! fin process-partition-sort2
!
      l1=nd1-n0
      l2=nd-n1
      if (l1.le.l2) then
       !----- l2>=l1 do l1 first
         top=top+2
         tmparr(top-1)=n1
         tmparr(top)=nd
         nd=nd1
      else
         top=top+2
         tmparr(top-1)=n0
         tmparr(top)=nd1
         n0=n1
      end if
      goto 3333
   deallocate(tmparr)
   return
end subroutine sorter
!---------------------------------------------------------------------
!            correct_order
!---------------------------------------------------------------------
logical function correct_order(k,l,x1,x2)
   implicit none
   integer,intent(in) ::   k, l
   integer, dimension(:), intent(in) ::x1(*), x2(*)

   correct_order=(x1(k).le.x1(l))
   if (x1(k).eq.x1(l)) then
      correct_order=(x2(k).le.x2(l))
      if (x2(k).eq.x2(l)) correct_order=(k.le.l)
   end if

   return
end function correct_order

!---------------------------------------------------------------------
!            exchange
!---------------------------------------------------------------------
subroutine exchange(i1,i2,a1,a2)
   implicit none 
   integer,intent(in) ::   i1, i2
   integer,intent(inout), dimension(*) :: a1, a2
   integer itmp

   itmp=a1(i1)
   a1(i1)=a1(i2)
   a1(i2)=itmp

   itmp=a2(i1)
   a2(i1)=a2(i2)
   a2(i2)=itmp

   return
end subroutine exchange

end module sort2
