! <compile=optimized>
#include "copyright.h"

!--------------------------------------------------------------
! CHAMBER
!
! Description: CHAMBER is a program for converting CHARMM
!              PSF files, along with associated topology
!              and coordinate files into prmtop and inpcrd
!              files for use in Sander and PMEMD. With 
!              suitably patched versions of Sander and PMEMD
!              this allows simulations to be run using the
!              CHARMM force field.
!
! Authors: Mark Williamson (SDSC, 2009)
!          Michael Crowley (NREL, 2009)
!          Ross Walker (SDSC, 2009)
!
!--------------------------------------------------------------

program main
!==============================================================================
!            MAIN PROGRAM
!=============================================================================
   use psfprm
   use prmwrt, only:write_prmtop_header,write_prmtop_atres,read_coordinates, &
                    print_prmtop_summary
   use molnt
   implicit none
   !Used for execution timing
   real :: start_time, stop_time
   
   psf_filename = "psf.psf"
   param_filename = "par_all27_prot_na.prm"
   topology_filename="top_all27_prot_na.rtf"
   prmtop_filename="prmtop"
   outu_name= "stdout"
   chmcrd_filename = "chmpdb.pdb"
   inpcrd_filename = "inpcrd"
   chmrst_filename = ""
   vmd_prmtop_filename="vmd_prmtop"

   call print_header() !Generic header explaining the program

   call cmd_args
   if(verbose_debug)write(outu,*) &
        "param ",param_filename(1:len_trim(param_filename))
   if(verbose_debug)write(outu,*)"psf ",psf_filename(1:len_trim(psf_filename))

   call cpu_time(start_time) !Start the clock

   call open_units !Open all files that will be used
   call read_psf_flags !Read in all the flags on the 1st line of the PSF
                       !Also checks if this appears to be a psf file.
   call set_formats !Set what fmt0{0-6} mean (this is a function of qextfmt)
   call read_title
   call read_psf_atom_data        !read psf, topolog, get nb parameters
   call read_psf_bond_angle_dihe  !get bonded lists from psf

   if( has_psf_flag("CHEQ") )  then
     ! Skip over the Charge Equilibration data which gives
     ! the total number of molecules and indirectly how many 
     ! atoms in each molecule
     call skip_over_psf_molnt
   else
     ! Skip over, information is provided via
     ! derive_molecule_information() later on
   endif

   call read_psf_lp

   if( has_psf_flag("CMAP")  )  then
     call read_psf_cmap
   elseif (CMAP_enabled) then
     !User specified cmap on command line
     !but no CMAP flag was found in the psf.
     write(outu,'(a)') "FATAL ERROR: -cmap specified on command line but no CMAP flag found in PSF header."
     call mexit(outu,1)
   endif

   ! Derive molecular topology from the bonded list.
   ! Essentially reproduces the information
   ! present in charmm's MOLNT and IMOLNT
   ! Used when writing box information and by get_gbradii()
   call derive_molecule_information 

   call print_psf_summary         !Summary of what was parsed from the PSF
 
 !--- psf file is read, now collect important information

   call atom_types           !organize atom types found in psf file 
   call get_atom_parameters  !get nonbond parameters for atom types
   call bond_types           !organize bond types found in psf file
   call angle_types          !         angle
   call dihe_types           !         dihe
   call improper_types       !         impropers
   call excluded_atoms
   call get_parameters

   call get_bonded_params

   if( has_psf_flag("CMAP") .and. CMAP_enabled )then
     !CMAP data is read from the par file even if CMAP wasn't enabled
     !on the command line. It just won't get written to the prmtop if
     !-nocmap was specified.
     call assign_cmap_types
   endif

   call print_assignment_summary !Summary of what was assigned

   call print_prmtop_summary !Summary of extra prmtop related terms

   !This is read before write_prmtop* because if there is any box information
   !it will change the output of the prmtop file.
   call read_coordinates(chmcrd_unit)

   call write_prmtop_header(prmtop_unit, .FALSE.)
   call write_prmtop_atres(prmtop_unit, .FALSE.)

   if ( want_vmd_prmtop ) then
     call write_prmtop_header(vmd_prmtop_unit, .TRUE.)
     call write_prmtop_atres(vmd_prmtop_unit, .TRUE.)
   endif

   call deallocate_everything
   call deallocate_molnt

   !stop the clock
   call cpu_time(stop_time)
   write(outu,'(a28,f8.2,a8)') "| Conversion carried out in ",stop_time - start_time," seconds"
   
   call mexit(outu,0)

end program main

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine cmd_args
  use psfprm
  implicit none
   character(len=80) arg
   !         temp for each of the whitespace delimited command line arguments
   integer iarg
   !         index of the current argument
   !integer iargc_wrap
   !         wrapper to intrinsic that returns the index of the last argument
   !         from either the command line or a string
   integer last_arg_index
   !         index of the last argument
   integer iargc

   logical :: cmap_option_specified = .false.

   iarg = 0
   last_arg_index = iargc()
   do while (iarg < last_arg_index)
      iarg = iarg + 1

      call getarg_wrap(iarg,arg)

      select case ( arg )

      case ('-box')
         ! We are making the assumption that the box is cubic at this moment:
         ifbox = 1
         box%pbcType = 1 !Cubic

         box%alpha   = 90.0d0
         box%beta    = 90.0d0
         box%gamma   = 90.0d0


         if(verbose_debug)write(outu,*)"Found command line arg -box"

         iarg = iarg + 1
         call getarg_wrap(iarg,arg)
         !Use an Internal file to convert from CHAR to real
         read( arg, '(f12.7)') box%a

         iarg = iarg + 1
         call getarg_wrap(iarg,arg)
         read( arg, '(f12.7)') box%b

         iarg = iarg + 1
         call getarg_wrap(iarg,arg)
         read( arg, '(f12.7)') box%c

         !Sanity checks
         if ( box%a <= 0 .or. box%b <= 0 .or. box%c <= 0) then
           write(outu,'(a)') "FATAL ERROR: Box length cannot be zero or negative"
           call mexit(outu,1)
         endif

         if ( box%a < 4.0 .or. box%b < 4.0 .or. box%c < 4.0) then
           write(outu,'(a)') "FATAL ERROR: Box length must be greater than 4A"
           call mexit(outu,1)
         endif

         if ( box%a > 9999.0 .or. box%b > 9999.0 .or. box%c > 9999.0) then
           write(outu,'(a)') "FATAL ERROR: Box length must be less than 9999A"
           call mexit(outu,1)
         endif

    
         if(verbose_debug) then
            write(outu,*) &
              "     using CUBIC box with the following a,b,c dimensions :"
            write(outu,'(7x,3f12.6)') box%a,box%b,box%c
         endif


      case ('-p')
         if(verbose_debug)write(outu,*)"Found command line arg -p"
         iarg = iarg + 1
         call getarg_wrap(iarg,prmtop_filename)
         if(verbose_debug)write(outu,*) &
              "     using ",prmtop_filename(1:len_trim(prmtop_filename))
         
      case ('-inpcrd') 
         if(verbose_debug)write(outu,*)"Found command line arg -inpcrd"
         iarg = iarg + 1
         call getarg_wrap(iarg,inpcrd_filename)
         if(verbose_debug)write(outu,*) &
              "     using ",inpcrd_filename(1:len_trim(inpcrd_filename))
         
      case ( '-param') 
         if(verbose_debug)write(outu,*)"Found command line arg -param"
         iarg = iarg + 1
         call getarg_wrap(iarg,param_filename)
         if(verbose_debug)write(outu,*) &
              "     using ",param_filename(1:len_trim(param_filename))
         
      case ( '-psf') 
         if(verbose_debug)write(outu,*)"Found command line arg -psf"
         iarg = iarg + 1
         call getarg_wrap(iarg,psf_filename)
         if(verbose_debug)write(outu,*) &
              "     using ",psf_filename(1:len_trim(psf_filename))
         l_xplor_psf = .false.
         
      case ( '-xpsf') 
         if(verbose_debug)write(outu,*)"Found command line arg -xpsf"
         iarg = iarg + 1
         call getarg_wrap(iarg,psf_filename)
         if(verbose_debug)write(outu,*) &
              "     using ",psf_filename(1:len_trim(psf_filename))
         l_xplor_psf = .true.
         
      case ( '-crd') 
         if(verbose_debug)write(outu,*)"Found command line arg -crd"
         iarg = iarg + 1
         call getarg_wrap(iarg,chmcrd_filename)
         if(verbose_debug)write(outu,*) &
              "     using ",chmcrd_filename(1:len_trim(chmcrd_filename))

      case ( '-top') 
         if(verbose_debug)write(outu,*)"Found command line arg -top"
         iarg = iarg + 1
         call getarg_wrap(iarg,topology_filename)
         if(verbose_debug)write(outu,*) &
              "   using ",topology_filename(1:len_trim(topology_filename))
      case ( '-tip3_flex') 
         if(verbose_debug)write(outu,'("Found command line arg -tip3_flex")')
         tip3_flex=.true.
         if(verbose_debug)write(outu, &
              '("   Will keep angle information for tip3 water")')
      case ( '-verbose') 
         verbose_debug=.true.
         verbose_debug_vdw=.true.
      case ( '-nocmap')
         if (cmap_option_specified) then
           !Problem we have already found a -cmap option - they
           !are mutually excludive.
           write(outu,'(/,5x,a,a)') ' Error -cmap and -nocmap are mutually exclusive.'
           call usage_message(1)
         else
           cmap_option_specified = .true.
           CMAP_enabled = .false.
         end if
      case ( '-cmap')
         if (cmap_option_specified) then
           !Problem we have already found a -cmap option - they
           !are mutually excludive.
           write(outu,'(/,5x,a,a)') ' Error -cmap and -nocmap are mutually exclusive.'
           call usage_message(1)
         else
           cmap_option_specified = .true.
           CMAP_enabled = .true.
         end if
      case ( '-radius_set')
         if(verbose_debug)write(outu,*)"Found command line arg -radius_set"
         iarg = iarg + 1
         call getarg_wrap(iarg,radius_set_read)
         read(radius_set_read,'(i1)') radius_set
      case ( '-vmd')
         if(verbose_debug)write(outu,*)"Found command line arg -vmd"
         want_vmd_prmtop = .true.
      case ( ' ') 
         continue
      case ( '-help', '--help', '-h', '/help')
         call usage_message(0)
      case default
         write(outu,'(/,5x,a,a)') ' Error unknown flag: ',arg
         call usage_message(1)
      end select
   end do  !  while (iarg < last_arg_index)

   !Either CMAP or NOCMAP must be specified there is no default.
   if (.not. cmap_option_specified) then
      write(outu,'(/,5x,a,a)') ' Error CMAP option was not specified.'
      write(outu,'(/,5x,a,a)') '       either -cmap or -nocmap must be'
      write(outu,'(/,5x,a,a)') '       specified. There is no default.'
      call usage_message(1)
   end if

   return
end subroutine cmd_args
 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Wrapper for GETARG to support grabbing command line arguments from a string
subroutine getarg_wrap(iarg, arg)

   implicit none
   integer iarg
   character(len=*) arg

   integer*4 which_argument

      ! Intrinsic getarg requires a 4 byte integer argument; 
      ! this guards the argument for builds with default 8 byte integers.
      which_argument = iarg
      call getarg(which_argument, arg)

 end subroutine getarg_wrap
 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Wrapper for IARGC to support reading command line input from a string
integer function iargc_wrap()

   implicit none
   integer iargc
      iargc_wrap = iargc()

end function iargc_wrap 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+
subroutine mexit(output_unit, status)
   
   !  mexit() - machine-dependent exit() procedure, designed to return an
   !            appropriate (success/failure) value to the operating system.
   
   implicit none
   integer output_unit  ! close this unit if greater than zero
   integer status       ! exit status; error if non-zero

   if (output_unit > 0) then
      close(unit=output_unit)
   end if

   call exit(status)

end subroutine mexit 

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine usage_message(status)
   use psfprm, only : outu

   implicit none
   integer status

   write(outu,'("")')
   write(outu,'(" Usage: chamber [args] ")')
   write(outu,'(" args for input are       <default> ")')
   write(outu,'("                 -top     <top_all27_prot_na.rtf>  ")')
   write(outu,'("                 -param   <par_all27_prot_na.prm>     ")')
   write(outu,'("                 -psf     <psf.psf> ")')
   write(outu,'("                 -crd     <chmpdb.pdb> ")')
   write(outu,'("")')
   write(outu,'("Note: -crd can specify a pdb, a CHARMM crd or CHARMM rst file.")')
   write(outu,'("           The filetype is auto detected.")')
   
   write(outu,'(/" args for output are      <default> ")')
   write(outu,'("                 -p       <prmtop> ")')
   write(outu,'("                 -inpcrd  <inpcrd>    ")')
   
   write(outu,'(/" args for options are: ")')
   write(outu,'("                 -cmap / -nocmap (Required option. Specifies")')
   write(outu,'("                                  whether CMAP terms should be")')
   write(outu,'("                                  included or excluded.")')
   write(outu,'("                 -tip3_flex  (allow angle in water) ")')
   write(outu,'("                 -box  a b c ")')
   write(outu,'("                     Set the Orthorhombic lattice parameters  a b c ")')
   write(outu,'("                     for the generated inpcrd file.")')
   write(outu,'("                 -verbose    (lots of progress messages) ")')
   write(outu,'("                 -vmd ")')
   write(outu,'("                     Write a VMD compatible form of the prmtop file ")')

   write(outu,'(/" -radius_set (GB radius set) options are: <default>     ")')
   write(outu,'("     0          bondi radii                 (bondi)  ")')
   write(outu,'("     1          amber6 modified Bondi radii (amber6) ")')
   write(outu,'("    <2>         modified Bondi radii        (mbondi) ")')
   write(outu,'("     6          H(N)-modified Bondi radii   (mbondi2)")')

   write(outu,'(/" arg for help (this message)")')
   write(outu,'("                 -h     ")')
   
   call mexit(outu,status)

end subroutine usage_message

subroutine print_header()

  use psfprm, only : outu

  implicit none
  
  write(outu,'(a)') " "
  write(outu,'(a)') "|                *****************************************************"
  write(outu,'(a)') "|                *   CHAMBER: Charmm psf to AMBER prmtop convertor   *"
  write(outu,'(a)') "|                *                                                   *"
  write(outu,'(a)') "|                *                       v1.0                        *"
  write(outu,'(a)') "|                *                                                   *"
  write(outu,'(a)') "|                *                    Written by:                    *"
  write(outu,'(a)') "|                *                 Mark J. Williamson                *"
  write(outu,'(a)') "|                *                 Michael F. Crowley                *"
  write(outu,'(a)') "|                *                   Ross C. Walker                  *"
  write(outu,'(a)') "|                *                                                   *"
  write(outu,'(a)') "|                *****************************************************"
  write(outu,'(a)') " "

  return  

end subroutine print_header
