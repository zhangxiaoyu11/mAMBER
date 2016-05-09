! <compile=optimized>
#include "copyright.h"

!--------------------------------------------------------------
! CHAMBER: PSF_STRINGS Module
!
! Description: Module containing all of the routines
!              for processing the strings used in psf files.
!
! Authors: Mark Williamson (SDSC, 2009)
!          Michael Crowley (NREL, 2009)
!          Ross Walker (SDSC, 2009)
!
!--------------------------------------------------------------

module psf_strings
private 

public :: getwords
public :: next_word
public :: next_int
public :: strUpCase
public :: strLowCase

character( * ), private, parameter :: lower_case = 'abcdefghijklmnopqrstuvwxyz'
character( * ), private, parameter :: upper_case = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' 

contains


!------------------- STRING FUNCTIONS -----------------------------------
integer function getwords(words,maxw,num,unit) result(retval)
   implicit none
   character(len=80) :: l,blank
   data blank/"                                                                                "/
   character(len=80), intent(out), dimension(1:maxw) :: words
   integer, intent(in) :: maxw,unit
   integer, intent(out) :: num
   
   character(len=80)  :: word
   integer ll  ,n , lw,comment

   ll=80
   retval=0
   l = " "
   read(unit,'(a80)',end=999)l
   !   write(outu,*)" ----getwords line ---- ",l(1:len_trim(l))
   comment = index(l,"!")
   if(comment .gt. 0 ) l(comment:80)=" "
   words(1:maxw)=blank
   nxtwrd: do n=1,maxw
      call next_word(l,lw,ll,word)
!      write(6,*)"     --- getwords word --- ",n,word(1:lw),lw
!     write(6,*)"     --- getwords line --- ",n,l(1:ll),ll
      words(n)(1:lw)=word(1:lw)
      if(lw == 0) exit nxtwrd
      if(words(n)(1:1) == "!") exit nxtwrd
   enddo nxtwrd
   num=n-1
!   write(6,*)"    ----- getwords returning ",num," words"
   if(lw == 0 ) retval=1
   return
 999  retval=2
      return
end function getwords

subroutine next_int(l,ll,num)
  implicit none
  character(len=80), intent(inout) :: l
  integer,intent(out) :: num
  integer,intent(inout) :: ll
  character(len=80) :: word
  character(len=80) :: fmt
  integer lw

  call next_word(l,lw,ll,word)
  write(fmt,'("(i",i1,")")')lw
  read(word,fmt)num
  return
end subroutine next_int

subroutine next_word(l,lw,ll,word)
    
   implicit none

   character(len=80), intent(inout) :: l
   character(len=80), intent(out) :: word
   integer, intent(inout) :: ll
   integer, intent(out) :: lw
   character(len=80) :: tmp_l
   character(len=11),parameter :: routine="-next_word "
   integer i,k

   word(1:80) = ' '
   ll = len_trim(l)

   !Note by Ross Walker - The original next_word code is buggy and does not
   !have correct protection to deal with situations where there is just 1 character
   !on a line. The following is a hack to deal with this.
   !--- Begin hack ---
   if (ll == 1) then
      if (l(1:1) /= " ") then
         lw = 1
         word(1:1)=l(1:1)
         l = ""
         return
      end if
   end if
   !--- End hack ---

   i=1
   k=0
   word=" "
   !--- remove leading blanks ----------------
   do while(i < ll .and. l(i:i) == " ")
      i=i+1
   enddo
   tmp_l = l
   l(1:ll-i+1)=tmp_l(i:ll)
   ll=ll-i+1
 !--- fill in next word ---------------------
   i=1
   do while(i < ll .and. l(i:i) /= " ")
      i=i+1
   enddo
   if(i > 1 .and. i < ll ) then
      word(1:i-1)=l(1:i-1)
      word(i:i)=" "
      lw=i-1
  !elseif(i > 1 .and. i == ll ) then !this will fail in the situation
                                     !when there is only one, one character
                                     !word on a line
   elseif(i .ge. 1 .and. i == ll ) then
      word(1:i)=l(1:i)
      lw=i
   else
      l = " "
      ll = 0
      lw = 0
      word=""
     return
   endif
   if(i< ll)then
      tmp_l = l
      l(1:ll-i+1)=tmp_l(i:ll)
      ll=ll-i+1
      if(ll > 0 ) ll=len_trim(l(1:ll))
   else
      l = " "
      ll = 0
   endif
 !--- remove leading blanks ----------------
   i=1
   do while(i < ll .and. l(i:i) == " ")
      i=i+1
   enddo
   tmp_l = l
   l(1:ll-i+1)=tmp_l(i:ll)
   ll=ll-i+1
   l(ll+1:80)=" "
   return
end subroutine next_word

function next_float(l,ll) result(value)
   implicit none
   real(kind=8) :: value
   character(len=ll), intent(inout) :: l
   character(len=ll) :: tmp_l
   character(len=80) :: fmt
   integer, intent(inout) :: ll
   integer i,k

   i=1
   k=0
   do while(i < ll .and. l(i:i) == " ")
      i=i+1
   enddo
   tmp_l = l
   l(1:ll-i+1)=tmp_l(i:ll)
   ll=ll-i+1
   i=0
   do while(i < ll .and. l(i:i) /= " ")
      i=i+1
   enddo
   if(i > 1) then
      if(i<10)then
         write(fmt,'("(G",i1,".",i1,")")')i,i
      else
         write(fmt,'("(G",i2,".5)")')i
      endif
      read(l(1:i),fmt)value
   else
      l = " "
      ll = 0
      value=0.d0 
     return
   endif
   if(i< ll)then
      tmp_l = l
      l(1:ll-i+1)=tmp_l(i:ll)
      ll=ll-i+1
      if(ll > 0 ) ll=len_trim(l(1:ll))
   else
      l = " "
      ll = 0
   endif
   return
end function next_float

subroutine skip_field(l,ll)
   implicit none
   character(len=ll), intent(inout) :: l
   character(len=ll) :: tmp_l
   integer, intent(inout) :: ll
   integer i,k

   i=0
   k=0
   do while(i < ll .and. l(i:i) /= " ")
      i=i+1
   enddo
   do while(i < ll .and. l(i:i) == " ")
      i=i+1
   enddo
   if(i == ll ) then
      l = " "
      ll = 0
      return
   endif
   tmp_l = l
   l(1:ll-i+1)=tmp_l(i:ll)
   ll=ll-i+1
   return
end subroutine skip_field


!These two functions are taken from:
!Figure 3.5B, pg 80, "Upgrading to Fortran 90", by Cooper Redwine,
!   1995 Springer-Verlag, New York. 

FUNCTION strUpCase ( Input_String ) RESULT ( Output_String )
  ! -- Argument and result
  CHARACTER( * ), INTENT( IN )     :: Input_String
  CHARACTER( LEN( Input_String ) ) :: Output_String
  ! -- Local variables
  INTEGER :: i, n

  ! -- Copy input string
  Output_String = Input_String
  ! -- Loop over string elements
  DO i = 1, LEN( Output_String )
    ! -- Find location of letter in lower case constant string
    n = INDEX( LOWER_CASE, Output_String( i:i ) )
    ! -- If current substring is a lower case letter, make it upper case
    IF ( n /= 0 ) Output_String( i:i ) = UPPER_CASE( n:n )
  END DO
END FUNCTION StrUpCase 

FUNCTION strLowCase ( Input_String ) RESULT ( Output_String )
  ! -- Argument and result
  CHARACTER( * ), INTENT( IN )     :: Input_String
  CHARACTER( LEN( Input_String ) ) :: Output_String
  ! -- Local variables
  INTEGER :: i, n

  ! -- Copy input string
  Output_String = Input_String
  ! -- Loop over string elements
  DO i = 1, LEN( Output_String )
    ! -- Find location of letter in upper case constant string
    n = INDEX( UPPER_CASE, Output_String( i:i ) )
    ! -- If current substring is an upper case letter, make it lower case
    IF ( n /= 0 ) Output_String( i:i ) = LOWER_CASE( n:n )
  END DO
END FUNCTION StrLowCase


end module psf_strings
