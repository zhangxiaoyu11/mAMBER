#include "copyright.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine mdfil here]
subroutine mdfil
   
   !     Author: George Seibel; many subsequent modifications by many others.
   implicit none
   
   !     Modified for multisander to allow reading command line from a string
   !     rather than the command line.  Use iargc_wrap and getarg_wrap instead
   !     of the intrinsics.
   
   !     OUTPUT: (to common)
   
#  include "files.h"
#ifdef MPI
#  include "parallel.h"
#endif
   
   !     INTERNAL:
   
   character(len=80) arg
   !         temp for each of the whitespace delimited command line arguments
   integer iarg
   !         index of the current argument
   integer iargc_wrap
   !         wrapper to intrinsic that returns the index of the last argument
   !         from either the command line or a string
   integer last_arg_index
   !         index of the last argument
   
   !     --- default file names ---
   
   mdin   = 'mdin'
   mdout  = 'mdout'
   inpcrd = 'inpcrd'
   parm   = 'prmtop'
   restrt = 'restrt'
   refc   = 'refc'
   mdvel  = 'mdvel'
   mden   = 'mden'
   mdcrd  = 'mdcrd'
   mdinfo = 'mdinfo'
   vecs   = 'vecs'
   freqe   = 'dummy'
   rstdip = 'rstdip'
   inpdip = 'inpdip'
   mddip = 'mddip'
   radii = 'radii'
   cpin = 'cpin'
   cpout = 'cpout'
   cprestrt = 'cprestrt'
#ifdef MMTSB
   mmtsb_setup_file = 'mmtsb_setup.job'
#endif
   if (numgroup == 1) groups = ' '

   
   !     --- default status of output: New
   
   owrite = 'N'

   !     --- get command line arguments ---
   
   iarg = 0
   last_arg_index = iargc_wrap()
   do while (iarg < last_arg_index)
      iarg = iarg + 1

      call getarg_wrap(iarg,arg)

      if (arg == '-O') then
#ifdef ABSOFT_WINDOWS
         ! although Replace is standard Fortran 90 apparently Absoft f90.exe 
         ! cannot handle it
         owrite = 'U' !      status of output: unknown
#else
         owrite = 'R' !      status of output: Replace
#endif
      else if (arg == '-i') then
         iarg = iarg + 1
         call getarg_wrap(iarg,mdin)
      else if (arg == '-o') then
         iarg = iarg + 1
         call getarg_wrap(iarg,mdout)
      else if (arg == '-p') then
         iarg = iarg + 1
         call getarg_wrap(iarg,parm)
      else if (arg == '-c') then
         iarg = iarg + 1
         call getarg_wrap(iarg,inpcrd)
      else if (arg == '-vecs') then
         iarg = iarg + 1
         call getarg_wrap(iarg,vecs)
      else if (arg == '-radii') then
         iarg = iarg + 1
         call getarg_wrap(iarg,radii)
      else if (arg == '-f') then
         iarg = iarg + 1
         call getarg_wrap(iarg,freqe)
      else if (arg == '-r') then
         iarg = iarg + 1
         call getarg_wrap(iarg,restrt)
      else if (arg == '-ref' .or. arg == '-z') then
         iarg = iarg + 1
         call getarg_wrap(iarg,refc)
      else if (arg == '-e') then
         iarg = iarg + 1
         call getarg_wrap(iarg,mden)
      else if (arg == '-v') then
         iarg = iarg + 1
         call getarg_wrap(iarg,mdvel)
      else if (arg == '-x'.or.arg == '-t') then
         iarg = iarg + 1
         call getarg_wrap(iarg,mdcrd)
      else if (arg == '-inf') then
         iarg = iarg + 1
         call getarg_wrap(iarg,mdinfo)
      else if (arg == '-idip') then
         iarg = iarg + 1
         call getarg_wrap(iarg,inpdip)
      else if (arg == '-rdip') then
         iarg = iarg + 1
         call getarg_wrap(iarg,rstdip)
      else if (arg == '-mdip') then
         iarg = iarg + 1
         call getarg_wrap(iarg,mddip)
      else if (arg == '-cpin') then
         iarg = iarg + 1
         call getarg_wrap(iarg,cpin)
      else if (arg == '-cpout') then
         iarg = iarg + 1
         call getarg_wrap(iarg,cpout)
      else if (arg == '-cprestrt') then
         iarg = iarg + 1
         call getarg_wrap(iarg,cprestrt)
#ifdef MMTSB
      else if (arg == '-mmtsb') then
         iarg = iarg + 1
         call getarg_wrap(iarg, mmtsb_setup_file)
#endif
#ifdef MPI
      else if (arg(1:3) == '-p4') then
         iarg = iarg+1
      else if (arg == '-np') then
         iarg = iarg+1
      else if (arg == '-mpedbg') then
         continue
      else if (arg == '-dbx') then
         continue
      else if (arg == '-gdb') then
         continue
#endif
#ifdef MPI
         !     Parse input options for multisander
         
      else if (arg == '-ng') then
         iarg = iarg + 1
         call getarg_wrap(iarg,arg)
         read(arg,'(i5)',err=91) numgroup

      else if (arg == '-ng-nonsequential') then
         ng_sequential = .false.

      else if (arg == '-groupfile') then
         iarg = iarg + 1
         call getarg_wrap(iarg,groups)

      else if (arg == '-gpes') then
         iarg = iarg + 1
         call getarg_wrap(iarg,gpes)
#endif
      else if (arg == ' ') then
         continue
      else
         write(6,'(/,5x,a,a)') 'mdfil: Error unknown flag: ',arg
         write(6,9000)
         call mexit(6, 1)
      end if  ! (arg == '-O')
   end do  !  while (iarg < last_arg_index)
   
   return
   
#ifdef MPI
   91 write(6,*)'mdfil: Error "-nrecip" and "-ng" expect integer arguments'
   write(6,*)'                   '
   call mexit(6, 1)
#endif
   9000 format(/,5x, &
         'usage: sander  [-O] -i mdin -o mdout -p prmtop -c inpcrd ', &
         '-r restrt',/19x,'[-ref refc -x mdcrd -v mdvel -e mden ', &
         '-idip inpdip -rdip rstdip -mdip mddip ', &
         '-inf mdinfo -radii radii]' &
         , /, 'Consult the manual for additional options.')
end subroutine mdfil 

!        '-O                Overwrite existing files.',
!        '-i MDIN           Namelist control input file',
!        '-o MDOUT          Output file.',
!        '-p PARM           ParmTop file.',
!        '-c INPCRD         Coordinate file.',
!        '-vecs VECS        ???',
!        '-radii RADII      ???',
!        '-f FREQE          ???',
!        '-r RESTRT         ???',
!        '-ref REFC         ???',
!        '-z REFC           alias for -ref.',
!        '-e MDEN           ???',
!        '-v MDVEL          ???',
!        '-x MDCRD          ???',
!        '-t MDCRD          alias for -x MDCRD',
!        '-inf MDINFO       ???',
!        '-idip INPDIP      ???',
!        '-rdip RSTDIP      ???',
!        '-mdip MDDIP       ???',
!        '-cpin CPDAT       Constant pH state information ',
!        '-cprstrt CPDAT    Constant pH state restart information',
!        '-cpout CPOUT      Constant pH protonation output
!#ifdef MMTSB
!        '-mmtsb MMTSB      MMTSB Setup file; contents server generated'
!#endif


!     WRAPPER FOR IARGC() TO SUPPORT READING COMMAND LINE INPUT FROM A STRING



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of integer function iargc_wrap here]
integer function iargc_wrap()

   implicit none
   integer iargc
#ifdef MPI
#  include "files.h"
#  include "parallel.h"

   integer istart, iend
   integer ia, is, ie

   if (numgroup > 1 .and. groups /= ' ') then
      ia = 0
      istart = 0
      iend = len(groupbuffer)

      do while (istart < iend)
         if ( groupbuffer(istart:istart) == ' ' ) then
            istart = istart + 1
         else
            is = istart
            ie = istart
            do while ( ie < iend .and. &
                  groupbuffer(ie:ie) /= ' ' )
               ie = ie + 1
            end do
            ia = ia + 1
            istart = ie
         end if
      end do
      !         write (6,*) 'DEBUG: IARGC_WRAP, ARGUMENT COUNT IS ', ia
      iargc_wrap = ia
   else
#endif

      iargc_wrap = iargc()

#ifdef MPI
   end if
#endif
end function iargc_wrap 



!     WRAPPER FOR GETARG TO SUPPORT GRABBING COMMAND LINE ARGUMENTS FROM A STRING


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine getarg_wrap here]
subroutine getarg_wrap(iarg, arg)

   implicit none
   integer iarg
   character(len=*) arg

#ifdef MPI
#  include "files.h"
   integer istart, iend
   integer ia, is, ie

   if (groups /= ' ') then

      ia = 0
      istart = 1
      iend = len(groupbuffer)

      do while (istart < iend)

         if ( groupbuffer(istart:istart) == ' ' ) then
            istart = istart + 1
         else

            is = istart
            ie = istart
            do while ( ie < iend .and. &
                  groupbuffer(ie:ie) /= ' ' )
               ie = ie + 1
            end do
            ia = ia + 1

            if (iarg == ia) then
               arg = groupbuffer(is:ie)
               return
            end if
            istart = ie
         end if
      end do

   else
#endif

      call getarg(iarg, arg)

#ifdef MPI
   end if
#endif
end subroutine getarg_wrap 

