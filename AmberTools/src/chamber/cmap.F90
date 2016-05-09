! <compile=optimized>
#include "copyright.h"

!--------------------------------------------------------------
! CHAMBER: CMAP MODULE
!
! Description: Module file to define types for processing
!              CHARMM cmap terms in CHARMM force field files.
!
! Authors: Mark Williamson (SDSC, 2009)
!          Michael Crowley (NREL, 2009)
!          Ross Walker (SDSC, 2009)
!
!--------------------------------------------------------------

module cmap

   !Energy correction map 

   !http://dx.doi.org/10.1002/jcc.20065
   !http://dx.doi.org/10.1021/ja036959e

   integer :: cmap_term_count !Number of Cross (CMAP) Terms from PSF
   integer :: cmap_type_count !Number of unique cmap terms found.


   integer, pointer, dimension(:,:) :: cmap_index !Contains the atom numbers of the
                                                  !atoms making up the two dihedrals
                                                  !that form the cmap term. 6th index
                                                  !ultimately gets filled with the location
                                                  !of the cmap parameters - essentially a
                                                  !re-basing of the types.
                                                  !(6,cmap_term_count)

   integer,parameter :: cmap_entries_per_line = 5


   !Experimental cmap type
   type cmapParameter

     !Atom types for a given dihedral type pair:
     !C    NH1  CT1  C    NH1  CT1  C    NH1   24
     character(len=4),dimension(8) :: label !Charmm atom type label

     ! 5    9    6    5    9    6    5    9
     ! Charmm's internal numerical atom type label
     integer,dimension(8) :: charmm_atom_type_numerical_lbl

     integer :: resolution ! The number of cmap data points on each axis
                           ! Since it is a grid, this is the same for both
                           ! axis.

     integer :: gridStepSize  ! set once number_of_grid_steps is known
                              !  360/resolution

     real(kind=8), pointer, dimension(:,:) :: grid
                           !The number of grid points for that CMAP parameter
                           !The total should be resolution**2


     !The CMAP grid is used as a basis for a spline fit to obtain a CMAP
     !energy correction for any value of phi,psi

     ! Example CMAP grid with a resolution of 5:
     !
     !   +    +    +    +    +
     !
     ! P +    +    +    +    +
     ! s
     ! i +    +    +    +    +
     !
     !   +    +    +    +    +          Hence the degree step size is
     !                                  360/5 == 72.
     !   +    +    +    +    +
     !          Phi                     The origin is always -180 degrees
     !
     !   ^    ^    ^    ^    ^
     !   |    |    |    |    |
     ! -180 -108  -36  36   108 (angle in degrees)
     integer  :: gridOrigin=-180 !Where the 2D grid starts in degrees
   end type cmapParameter

   !and now, an array of them
   type(cmapParameter),pointer,dimension(:) :: needed_cmap_types

   contains


   !Display a summary of what was assigned from the CHARMM param file
   subroutine print_cmap_summary(u)
     implicit none
     integer, intent(in) :: u

     !Local
     integer             :: i


     write(u,'(a60)') &
     "====================================================================="
     write(u,'(a50)') "         CMAP parameters assignment index"
     write(u,'(a60)')&
    "====================================================================="
     write(u,'(a)') ""
     write(u,'(a,i4,a)') "Dumping entire contents of needed_cmap_types(",&
                                                 cmap_type_count,")"

     write(u,'(a)') ""

     do i=1,cmap_type_count
       write(6,'(a,i4,a,i4,a)') "needed_cmap_type(",i,"/",cmap_type_count,")"
       write(u,'(8(a4,1x),i4,i4)') needed_cmap_types(i)%label(1:8),&
                                     needed_cmap_types(i)%resolution,&
                                     needed_cmap_types(i)%gridStepSize
       write(u,'(5(f9.5))') &
                                 needed_cmap_types(i)%grid(&
                                 1:needed_cmap_types(i)%resolution,&
                                 1:needed_cmap_types(i)%resolution)
     enddo
   end subroutine print_cmap_summary


   subroutine deallocate_cmap
     implicit none
     integer :: outu = 6
     integer :: i,ierr

     if(associated(cmap_index)) then
       deallocate(cmap_index,stat=ierr)
       if (ierr /= 0) then
         write(outu,'(a)') "ERROR Deallocating cmap_index."
         call mexit(outu,1)
       end if
     end if

     do i=1, cmap_type_count
       if(associated(needed_cmap_types(i)%grid)) then
         deallocate(needed_cmap_types(i)%grid,stat=ierr)
         if (ierr /= 0) then
           write(outu,'(a)') "ERROR Deallocating needed_cmap_types(i)%grid"
           call mexit(outu,1)
         end if
       end if
     enddo

     if(associated(needed_cmap_types)) then
       deallocate(needed_cmap_types,stat=ierr)
       if (ierr /= 0) then
         write(outu,'(a)') "ERROR Deallocating needed_cmap_types."
         call mexit(outu,1)
       end if
     end if

   end subroutine deallocate_cmap

end module cmap
