! <compile=optimized>
#include "copyright.h"

!--------------------------------------------------------------
! CHAMBER: MOLNT MODULE
!
! Description: Module file for the generation of molecular information
!              solely from the bonded list.
!
! Authors: Mark Williamson (SDSC, 2009)
!
!--------------------------------------------------------------

module molnt

  ! Derive, just from the bonded information within the system, how many
  ! molecules there are and then which atoms belong to which molecule.
  ! This is all done using the fact that if an atom, i, is bonded to atom,
  ! j, then j is in the same molecule as i.
  !
  ! CHARMM's psf *sometimes* contains this information if the CHEQ keyword
  ! has been used in its compile; the MOLNT variable and iMOLNT arrays contain
  ! this. However, this cannot be assumed to always be present, hence we need
  ! to calculate this from what we DO know to always be present: the bonding 
  ! list.
  !
  ! Ultimately, this information is used to generate the SOLVENT_POINTERS and
  ! ATOMS_PER_MOLECULE sections in the outputted prmtop file if -box is defined.

  use psfprm, only: allint, natom, nbond, ib, jb
  use psfprm, only: verbose_debug, outu !Used for logging


  !What the outside world can see....

  !Subroutines:
  public                               :: derive_molecule_information
  public                               :: deallocate_molnt

  !Variables:
  public                               :: total_number_of_molecules
                                        ! Total number of molecules 
                                        ! in the entire system.
  public                               :: atom_to_molecule_map
                                        ! Given an atom index, 
                                        ! this is what molecule it is in
  public                               :: number_of_atoms_per_molecule
                                        ! Given an molecule index, this is
                                        ! how many atoms are in it.

  public    :: bonded_list


  !The rest is private
  private

  integer, save                        :: total_number_of_molecules = 0! MOLNT
  integer, pointer, dimension(:), save :: atom_to_molecule_map !imolnt
                                            ! map between atom index and
                                            ! molecule that it is in

  integer, parameter                   :: max_no_of_bonds_per_atom = 10
  integer, pointer, dimension(:,:)     :: bonded_list
  ! A 2D array containing a list of atoms that a given atom is bonded to.
  ! Given that max_no_of_bonds_per_atom = 5 in these examples:
  !
  ! e.g.        bonded_list(1:5, 4) is 1,2,3,0,0
  ! indicates that atom 4 is bonded to atoms 1, 2, and 3.
  !
  ! e.g.        bonded_list(1:5, 3) is 1,2,3,4,0
  ! indicates that atom 3 is bonded to atoms 1, 2, 3 and 4.
  !
  ! e.g.        bonded_list(1:5, 9) is 0,0,0,0,0
  ! indicates that atom 9 is not bonded to anything, i.e. 
  ! it is probably a single atom ion


  integer, pointer, dimension(:), save :: number_of_atoms_per_molecule

  contains


! The main entry point really and should populate the following:
!
!  total_number_of_molecules
!  atom_to_molecule_map(natom)
!  number_of_atoms_per_molecule(nmol)

subroutine derive_molecule_information
  implicit none

  !local
  integer :: i
  integer :: allocate_status


  !Allocate bonded list
  allocate( bonded_list(max_no_of_bonds_per_atom, natom), stat=allocate_status )
  if (allocate_status .ne. 0) then
    write(outu,'(a)') "bonded_list allocation error"
    call mexit(outu,1)
  endif

  !Populate the bonded list
  call generate_bonded_list(ib, jb, bonded_list)

  !DEBUG
  !write(outu, '(a1)')     ""
  !write(outu, '(a30)')    "Atomic bonding topology"
  !do i=1, natom
  !  write(*,*), i, bonded_list(1:max_no_of_bonds_per_atom, i)
  !enddo
  !write(outu, '(a1)')     ""



  !Now allocate for the molecular list
  call allint(natom, atom_to_molecule_map)


  call generate_molecule_list(bonded_list,               &
                                  total_number_of_molecules, &
                                  atom_to_molecule_map)

  !DEBUG
  if(verbose_debug) then
    write(outu, '(a1)')     ""
    write(outu, '(a30)')     "Molecule information"
    write(outu, '(a30,i8)')  "Total number of molecules is: ", total_number_of_molecules
    do i=1, natom
      write(outu, '(a6,i8,a30,i8)') "Atom :",i," is part of molecule number: ", atom_to_molecule_map(i)
    enddo
    write(outu, '(a30)') "End Molecule information"
    write(outu, '(a1)') ""
  endif

  call allint(total_number_of_molecules, number_of_atoms_per_molecule)

  ! Zero this array since to avoid any
  ! conditional jump or move depending on uninitialized value in these arrays
  number_of_atoms_per_molecule = 0

  call cal_num_of_atoms_per_molecule(atom_to_molecule_map, number_of_atoms_per_molecule)

end subroutine derive_molecule_information







subroutine generate_bonded_list(ib, jb, bonded_list)
! Populates the 2D bonded_list array using the list of bonds contained
! within arrays ib() and jb()

  implicit none
  integer, intent(in)           :: ib(nbond)
  integer, intent(in)           :: jb(nbond) !List of i and j'th atoms in each bond
  integer, intent(out)          :: bonded_list(max_no_of_bonds_per_atom, natom)

  ! local
  integer                       :: i

  ! Prepare the bonded_list for population
  ! a value of zero indicates a free space
  do i=1, natom
    bonded_list(1:max_no_of_bonds_per_atom, i) = 0
  enddo

  do i=1, nbond ! Walk the list of bonds
    call add_to_bonded_list(ib(i), jb(i), bonded_list)
    ! Do the other way  
    call add_to_bonded_list(jb(i), ib(i), bonded_list)
    ! TODO think this through... this 2nd call to add_to_bonded_list is 
    ! is only needed for all the atoms bar the first in the LAST molecule
  enddo

end subroutine generate_bonded_list




subroutine add_to_bonded_list(i, j, bonded_list)
  ! Used by generate_bonded_list to
  ! add atom j to atom i's list of bonded atoms

  implicit none
  integer, intent(in)             :: i, j
  integer, intent(inout)          :: bonded_list(max_no_of_bonds_per_atom, natom)

  ! Local
  integer                         :: k

  do k=1, max_no_of_bonds_per_atom ! Walk available placeholders
    if ( bonded_list( k, i ) .eq. 0  ) then ! atom i has got a free slot
      bonded_list( k, i ) = j !now add it
      exit ! We've added the atom, hence we're done, now quit 
    elseif (k == max_no_of_bonds_per_atom) then ! We've run out of placeholders
      write(outu,'(a)') "max_no_of_bonds_per_atom exceeded"
      call mexit( outu, 1 )
      stop
    endif
  enddo

end subroutine add_to_bonded_list





subroutine generate_molecule_list(bonded_list,               &
                                  total_number_of_molecules, &
                                  atom_to_molecule_map)
! Using the bonded_list information, and the fact that if an atom in a molecule
! is bonded to another atom, then that other atom is in the same molecule, derive:
!
!   i) Which molecule an atom is in for all atoms in the system.
!  ii) Work out the total number of molecules in the system.

  implicit none
  integer, intent(in)   :: bonded_list(max_no_of_bonds_per_atom, natom)
  integer, intent(out)  :: total_number_of_molecules
  integer, intent(out)  :: atom_to_molecule_map( natom )

  ! Local
  integer               :: i
  integer               :: current_molecule_index

  ! Set the current molecule index to zero
  current_molecule_index = 0

  ! Zero the molecular list
  atom_to_molecule_map = 0

  do i=1, natom

    if ( atom_to_molecule_map(i) .eq. 0 ) then !We have an atom in need of tagging
      current_molecule_index = current_molecule_index + 1
      ! write(outu,'(a)') ""
      ! write(outu,'(a,i8)') "current_atom_index is :",i
      ! write(outu,'(a,i8)') "current_molecule_index is :",current_molecule_index

      call tag_atoms_in_same_molecule(i,                      &
                                      current_molecule_index, &
                                      bonded_list,            &
                                      atom_to_molecule_map)

    ! Since tag_atoms_in_same_molecule() is recursive, when reaching this
    ! point we will have tagged every atom in the current_molecule_index, hence
    ! now increment the current_molecule_index on the next loop.
    endif

  enddo

  total_number_of_molecules = current_molecule_index

end subroutine generate_molecule_list





recursive subroutine tag_atoms_in_same_molecule(current_atom, &
                                                current_molecule_index, &
                                                bonded_list, &
                                                atom_to_molecule_map)

  ! Follow a given atom's bonded list, tagging every atom contained in that
  ! list with the current_molecule_index via the atom_to_molecule_map.
  ! NOTE this is recursive; it will follow every bonded atom until all
  ! are tagged with a molecule_index

  implicit none
  integer, intent(in)          :: current_atom
  integer, intent(in)          :: current_molecule_index
  integer, intent(in)          :: bonded_list(max_no_of_bonds_per_atom, natom)
  integer, intent(inout)       :: atom_to_molecule_map(natom)

  !Local
  integer                      :: i, current_bonded_atom

  do i=1, max_no_of_bonds_per_atom
    current_bonded_atom = bonded_list(i, current_atom)
   ! e.g.                 bonded_list(1:5, 4) is 1,2,3,0,0
   ! indicates that atom 4 is bonded to atoms 1, 2, and 3.

    if ( (i .eq. 1 )  .and.  current_bonded_atom .eq. 0  ) then
      ! write(outu,'(a)') "We have a one atom molecule"
      ! TAG the *current_atom* itself and not the bonded one!
      atom_to_molecule_map(current_atom) = current_molecule_index
      ! write(outu,'(a,i8,a,i8)') "tag_atoms_in_same_molecule has tagged atom ",&
      !    current_atom," as being part of molecule",current_molecule_index
      return
    endif


    ! Ensure that a) there actually is a bonded atom here
    !             b) has not already been tagged as being part of a molecule
    if( current_bonded_atom .ne. 0 ) then
      if( atom_to_molecule_map(current_bonded_atom ) .eq. 0 ) then

        ! Tag the bonded atom as being in the same molecule as current_atom
        atom_to_molecule_map(current_bonded_atom) = current_molecule_index

        ! Now follow that atom RECURSIVELY
        call tag_atoms_in_same_molecule(bonded_list(i, current_atom), &
                                        current_molecule_index,       &
                                        bonded_list,                  &
                                        atom_to_molecule_map)
      endif
    endif

  enddo
end subroutine tag_atoms_in_same_molecule



subroutine cal_num_of_atoms_per_molecule(atom_to_molecule_map, number_of_atoms_per_molecule)
! Given a complete atom_to_molecule_map list, this will work out how many
! atoms there are in each molecule and add this to the number_of_atoms_per_molecule array

  implicit none
  integer, pointer, dimension(:)   :: atom_to_molecule_map
  integer, pointer, dimension(:)   :: number_of_atoms_per_molecule

  !Local
  integer                          :: i, j
  integer                          :: sum_of_atoms_in_molecules


  ! Check is number_of_atoms_per_molecule is actually allocated
  if ( associated(number_of_atoms_per_molecule) .eqv. .FALSE. ) then
    write(outu,'(a)') "In cal_num_of_atoms_per_molecule ()"
    write(outu,'(a)') "number_of_atoms_per_molecule has not been allocated"
    call mexit( outu, 1 )
  endif

  ! This requires total_number_of_molecules * natom loops 
  ! whereas the old method just required one loop over natom
  do i=1, total_number_of_molecules
    do j=1, natom
       if ( atom_to_molecule_map(j) == i) then
         number_of_atoms_per_molecule(i) = number_of_atoms_per_molecule(i) + 1
       endif
    enddo
  enddo


  ! Sanity checks:

  ! 1) Sum of all atoms in all molecules must equal total number of atoms
  sum_of_atoms_in_molecules = sum(number_of_atoms_per_molecule)

  if ( natom .ne. sum_of_atoms_in_molecules ) then
    write(outu,'(a50)') "FATAL ERROR:"
    write(outu,'(a50, i8)') &
     "Sum of all atoms in all molecules:", sum_of_atoms_in_molecules
    write(outu,'(a50, i8)') &
     "does NOT equal to the total number of atoms:", natom
    call mexit(outu,1)
  endif


 ! debug
 if(verbose_debug) then
   do i=1, total_number_of_molecules
     write(outu,'(a9,i8,a8,i8,a12)') "molecule(",i,") has ", &
                  number_of_atoms_per_molecule(i), "atoms in it"
   enddo
 endif

end subroutine cal_num_of_atoms_per_molecule



subroutine deallocate_molnt

  integer            :: ierr

  if ( associated ( atom_to_molecule_map ) ) then
    deallocate( atom_to_molecule_map, stat=ierr)
    if (ierr /= 0) then
      write(outu,'(a)') "ERROR Deallocating atom_to_molecule_map."
      call mexit(outu,1)
     endif
  endif

  if ( associated ( bonded_list ) ) then
    deallocate( bonded_list, stat=ierr)
    if (ierr /= 0) then
      write(outu,'(a)') "ERROR Deallocating bonded_list."
      call mexit(outu,1)
     endif
  endif

  if ( associated ( number_of_atoms_per_molecule ) ) then
    deallocate( number_of_atoms_per_molecule, stat=ierr)
    if (ierr /= 0) then
      write(outu,'(a)') "ERROR Deallocating number_of_atoms_per_molecule."
      call mexit(outu,1)
     endif
  endif

end subroutine deallocate_molnt




end module molnt

