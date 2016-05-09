! <compile=optimized>
#include "copyright.h"

!--------------------------------------------------------------
! CHAMBER: PRMWRT module
!
! Description: Routines for writing the prmtop file.
!              And for reading coordinates.
!
! Authors: Mark Williamson (SDSC, 2009)
!          Michael Crowley (NREL, 2009)
!          Ross Walker (SDSC, 2009)
!
!--------------------------------------------------------------

module prmwrt

character(len=12) :: flag_format='("%FLAG ",a)'
character(len=14) :: form_format='("%FORMAT",a)'
integer :: maxexcl


   integer, pointer, dimension(:) :: nonbond_excl_14,nonbond_excl_pointers_14,ico
   logical, parameter :: verbose= .false.

contains
subroutine write_prmtop_header(unit, is_vmd_compatible)
   use psfprm
   implicit none
   character(len=21) :: routine="<write_prmtop_header>"
   integer, intent(in) :: unit
   logical, intent(in) :: is_vmd_compatible ! Write a prmtop file that can be read by VMD
   character(len=80) :: fmt
   integer :: ntheta_am,nbona_am,nphia_am
   logical :: done

   maxexcl=natom*8

   write(unit,'(2a)')"%VERSION  VERSION_STAMP = V0001.000  ", &
                               "DATE = 06/11/03  12:02:08"
   if (is_vmd_compatible) then
     write(unit,flag_format)"TITLE"
     ! Bug in VMD - does not currently read the format string for the title
     ! it always expects 20a4. Additionally it then expects the title to be blank.
     write(unit,form_format) "(20a4)"
     write(unit,'(20a4)') "" 
   else
     ! Note CTITLE used here instead of TITLE to deliberately prevent old versions
     ! of PMEMD or sander that do not support the CHARMM force field from being used.
     write(unit,flag_format)"CTITLE"
     write(unit,form_format) "(a80)"
     write(unit,'(a80)')title
   endif
   write(unit,flag_format)"POINTERS"
   fmt="(10I8)"
   write(unit,form_format)fmt

   ntheta_am=mtheta
   nbona_am = mbona
   nphia_am = mphia
   done=.false.
   do while(.not. done)
      allocate(nonbond_excl_14(maxexcl),nonbond_excl_pointers_14(natom))
      call make_exclusion_list(nonbond_excl_14,nonbond_excl_pointers_14,maxexcl,done) 
      if (.not. done)then
         deallocate(nonbond_excl_14,nonbond_excl_pointers_14)
         maxexcl=maxexcl+maxexcl
         write(outu,'(2a,i12)')routine," make_exclusion_list reallocating ",maxexcl
      endif
   enddo
   natyp = ntypes
   nmxrs = max_res_at
   write(outu,'(//2a,i6)')routine," NPHB ",nphb
   write(unit,fmt)natom,  ntypes, nbonh,  mbona,  ntheth, mtheta, &
                  nphih,  mphia,  nhparm, nparm,  &
                  nnb,    nres, &
                  nbona_am,  ntheta_am, nphia_am,  numbnd, numang, nptra, &
                  natyp,  nphb,   &
                  ipert,  nbper,  ngper,  ndper, &
                  mbper,  mgper,  mdper,  ifbox,  nmxrs,  ifcap, &
                  nextra

   write(unit,flag_format)"FORCE_FIELD_TYPE"
   write(unit,form_format) "(i2,a78)"
   write(unit,'(" 1 CHARMM",i4,2x,a65)')chmff_verno,chmff_type
   

   return
end subroutine write_prmtop_header

subroutine write_prmtop_atres(unit, is_vmd_compatible)
   use psfprm 
   use cmap
   use molnt !SOLVENT_POINTERS + ATOMS_PER_MOLECULE information
   implicit none
   integer, intent(in) :: unit
   logical, intent(in) :: is_vmd_compatible ! Write a prmtop file that can be read by VMD
   integer :: i,j,index
   character(len=20) :: fmt
   character(len=2) :: word
   integer, dimension(ntypes,ntypes) :: tmpico


   write(unit,flag_format)"ATOM_NAME"
   fmt="(20a4)"
   write(unit,form_format)fmt
   write(unit,fmt)(atom_label(i),i=1,natom)

   if (is_vmd_compatible) then
     fmt="(5e16.8)"
   else
     fmt="(3e24.16)"
   endif
   write(unit,flag_format)"CHARGE"
   if (is_vmd_compatible) then
     ! VMD cannot deal with the valid %COMMENT section
   else
     write(unit,'(a)')"%COMMENT  Atomic charge multiplied by sqrt(332.0716D0) (CCELEC) "
   endif
   write(unit,form_format)fmt
   write(unit,fmt)(cg(i),i=1,natom)
   
   write(unit,flag_format)"MASS"
   fmt="(5e16.8)"
   write(unit,form_format)fmt
   write(unit,fmt)(amass(i),i=1,natom)
   
   write(unit,flag_format)"ATOM_TYPE_INDEX"
   fmt="(10i8)"
   write(unit,form_format)fmt
   write(unit,fmt)(locattype(i),i=1,natom)
   
   write(unit,flag_format)"NUMBER_EXCLUDED_ATOMS"
   fmt="(10i8)"
   write(unit,form_format)fmt
   write(unit,fmt)(nonbond_excl_pointers_14(i),i=1,natom)

   write(unit,flag_format)"EXCLUDED_ATOMS_LIST"
   write(unit,form_format)fmt
   write(unit,fmt)(nonbond_excl_14(i),i=1,nnb)

   deallocate(nonbond_excl_14,nonbond_excl_pointers_14)

   
   call allint(ntypes*ntypes,ico)
   index=0
   do i=1,ntypes
      do j=1,i
         index=index+1
         tmpico(i,j)=index
         tmpico(j,i)=index
      enddo
   enddo
   index=0
   do i=1,ntypes
      do j=1,ntypes
         index=index+1
         ico(index)=tmpico(i,j)
      enddo
   enddo
   write(unit,flag_format)"NONBONDED_PARM_INDEX"
   fmt="(10i8)"
   write(unit,form_format)fmt
   write(unit,fmt)(ico(i),i=1,ntypes*ntypes)
   call deallocate_int(ico)

   write(unit,flag_format)"RESIDUE_LABEL"
   fmt="(20a4)"
   write(unit,form_format)fmt
   write(unit,fmt)(lres(i),i=1,nres)

   write(unit,flag_format)"RESIDUE_POINTER"
   fmt="(10i8)"
   write(unit,form_format)fmt
   write(unit,fmt)(ipres(i),i=1,nres)

   write(unit,flag_format)"BOND_FORCE_CONSTANT"
   fmt="(5e16.8)"
   write(unit,form_format)fmt
   write(unit,fmt)(rk(i),i=1,numbnd)

   write(unit,flag_format)"BOND_EQUIL_VALUE"
   write(unit,form_format)fmt
   write(unit,fmt)(req(i),i=1,numbnd)

   write(unit,flag_format)"ANGLE_FORCE_CONSTANT"
   write(unit,form_format)fmt
   write(unit,fmt)(tk(i),i=1,numang)

   fmt="(3e25.17)"
   write(unit,flag_format)"ANGLE_EQUIL_VALUE"
   write(unit,form_format)fmt
   write(unit,fmt)(teq(i),i=1,numang)

   write(unit,flag_format)"CHARMM_UREY_BRADLEY_COUNT"
   write(unit,'(a)')"%COMMENT  V(ub) = K_ub(r_ik - R_ub)**2"
   write(unit,'(a)')"%COMMENT  Number of Urey Bradley terms and types"
   fmt="(2i8)"
   write(unit,form_format)fmt
   write(unit,fmt)nub, nubtypes

   write(unit,flag_format)"CHARMM_UREY_BRADLEY"
   write(unit,'(a)')"%COMMENT  List of the two atoms and its parameter index"
   write(unit,'(a)')"%COMMENT  in each UB term: i,k,index"
   fmt="(10i8)"
   write(unit,form_format)fmt
   write(unit,fmt)(ub_atm_i(i), ub_atm_k(i), ub_idx(i),i=1,nub)

   write(unit,flag_format)"CHARMM_UREY_BRADLEY_FORCE_CONSTANT"
   write(unit,'(a)')"%COMMENT  K_ub: kcal/mole/A**2"
   fmt="(5e16.8)"
   write(unit,form_format)fmt
   write(unit,fmt)(kub(i),i=1,nubtypes)

   write(unit,flag_format)"CHARMM_UREY_BRADLEY_EQUIL_VALUE"
   write(unit,'(a)')"%COMMENT  r_ub: A "
   write(unit,form_format)fmt
   write(unit,fmt)(rub(i),i=1,nubtypes)

   write(unit,flag_format)"DIHEDRAL_FORCE_CONSTANT"
   write(unit,form_format)fmt
   write(unit,fmt)(pk(i),i=1,ndihtypes)

   write(unit,flag_format)"DIHEDRAL_PERIODICITY"
   write(unit,form_format)fmt
   write(unit,fmt)(pn(i),i=1,ndihtypes)

   write(unit,flag_format)"DIHEDRAL_PHASE"
   write(unit,form_format)fmt
   write(unit,fmt)(phase(i),i=1,ndihtypes)

   write(unit,flag_format)"SCEE_SCALE_FACTOR"
   write(unit,form_format)fmt
   write(unit,fmt)(SCEE_SCALE_FACTOR(i),i=1,ndihtypes)

   write(unit,flag_format)"SCNB_SCALE_FACTOR"
   write(unit,form_format)fmt
   write(unit,fmt)(SCNB_SCALE_FACTOR(i),i=1,ndihtypes)


   write(unit,flag_format)"CHARMM_NUM_IMPROPERS"
   write(unit,'(a)')"%COMMENT  Number of terms contributing to the"
   write(unit,'(a)')"%COMMENT  quadratic four atom improper energy term:"
   write(unit,'(a)')"%COMMENT  V(improper) = K_psi(psi - psi_0)**2"
   fmt="(10i8)"
   write(unit,form_format)fmt
   write(unit,fmt)nimphi

   write(unit,flag_format)"CHARMM_IMPROPERS"
   write(unit,'(a)')"%COMMENT  List of the four atoms in each improper term"
   write(unit,'(a)')"%COMMENT  i,j,k,l,index  i,j,k,l,index"
   write(unit,'(a)')"%COMMENT  where index is into the following two lists:"
   write(unit,'(a)')"%COMMENT  CHARMM_IMPROPER_{FORCE_CONSTANT,IMPROPER_PHASE}"
   write(unit,form_format)fmt
   write(unit,fmt)(im(i),jm(i),km(i),lm(i),imp(i),i=1,nimphi)

   write(unit,flag_format)"CHARMM_NUM_IMPR_TYPES"
   write(unit,'(a)')"%COMMENT  Number of unique parameters contributing to the"
   write(unit,'(a)')"%COMMENT  quadratic four atom improper energy term"
   fmt="(i8)"
   write(unit,form_format)fmt
   write(unit,fmt)nimprtypes

   write(unit,flag_format)"CHARMM_IMPROPER_FORCE_CONSTANT"
   write(unit,'(a)')"%COMMENT  K_psi: kcal/mole/rad**2 "
   fmt="(5e16.8)"
   write(unit,form_format)fmt
   write(unit,fmt)(pk_impr(i),i=1,nimprtypes)

   write(unit,flag_format)"CHARMM_IMPROPER_PHASE"
   write(unit,'(a)')"%COMMENT  psi: degrees"
   write(unit,form_format)fmt
   write(unit,fmt)(phase_impr(i),i=1,nimprtypes)

   write(unit,flag_format)"SOLTY"
   write(unit,form_format)fmt
   write(unit,fmt)(0.d0,i=1,natyp)

   if (is_vmd_compatible) then
     fmt="(5e16.8)"
   else
     fmt="(3e24.16)"
   endif
   write(unit,flag_format)"LENNARD_JONES_ACOEF"
   write(unit,form_format)fmt
   write(unit,fmt)(cn1(i),i=1,ntypes*(ntypes+1)/2)

   write(unit,flag_format)"LENNARD_JONES_BCOEF"
   write(unit,form_format)fmt
   write(unit,fmt)(cn2(i),i=1,ntypes*(ntypes+1)/2)

   write(unit,flag_format)"LENNARD_JONES_14_ACOEF"
   write(unit,form_format)fmt
   write(unit,fmt)(cn114(i),i=1,ntypes*(ntypes+1)/2)

   write(unit,flag_format)"LENNARD_JONES_14_BCOEF"
   write(unit,form_format)fmt
   write(unit,fmt)(cn214(i),i=1,ntypes*(ntypes+1)/2)

   write(unit,flag_format)"BONDS_INC_HYDROGEN"
   fmt="(10i8)"
   write(unit,form_format)fmt
   write(unit,fmt)(ibh(i),jbh(i),icbh(i),i=1,nbonh)

   write(unit,flag_format)"BONDS_WITHOUT_HYDROGEN"
   write(unit,form_format)fmt
   write(unit,fmt)(iba(i),jba(i),icba(i),i=1,mbona)

   write(unit,flag_format)"ANGLES_INC_HYDROGEN"
   write(unit,form_format)fmt
   write(unit,fmt)(ith(i),jth(i),kth(i),icth(i),i=1,ntheth)

   write(unit,flag_format)"ANGLES_WITHOUT_HYDROGEN"
   write(unit,form_format)fmt
   write(unit,fmt)(ita(i),jta(i),kta(i),icta(i),i=1,mtheta)

   write(unit,flag_format)"DIHEDRALS_INC_HYDROGEN"
   write(unit,form_format)fmt
   write(unit,fmt)(iph2(i),jph2(i),kph2(i),lph2(i),icph2(i),i=1,nphih)

   write(unit,flag_format)"DIHEDRALS_WITHOUT_HYDROGEN"
   write(unit,form_format)fmt
   write(unit,fmt)(ipa2(i),jpa2(i),kpa2(i),lpa2(i),icpa2(i),i=1,mphia)

   write(unit,flag_format)"HBOND_ACOEF"
   fmt="(5e16.8)"
   write(unit,form_format)fmt
   ! Nothing to write, this is not used unless there are 10-12 potentials
   write(unit,fmt)

   write(unit,flag_format)"HBOND_BCOEF"
   write(unit,form_format)fmt
   ! Nothing to write, this is not used unless there are 10-12 potentials
   write(unit,fmt)

   write(unit,flag_format)"HBCUT"
   write(unit,form_format)fmt
   ! Nothing to write, this is not used unless there are 10-12 potentials
   write(unit,fmt)

   write(unit,flag_format)"AMBER_ATOM_TYPE"
   fmt="(20a4)"
   write(unit,form_format)fmt
   write(unit,fmt)(attype_name(locattype(i)),i=1,natom)

   write(unit,flag_format)"TREE_CHAIN_CLASSIFICATION"
   write(unit,form_format)fmt
   write(unit,fmt)("BLA ",i=1,natom)

   write(unit,flag_format)"JOIN_ARRAY"
   fmt="(10i8)"
   write(unit,form_format)fmt
   write(unit,fmt)(0,i=1,natom)

   write(unit,flag_format)"IROTAT"
   write(unit,form_format)fmt
   write(unit,fmt)(0,i=1,natom)

   write(unit,flag_format)"RADIUS_SET"
   fmt='(1a80)'
   write(unit,form_format)fmt
   write(unit,fmt) radius_set_name

   write(unit,flag_format)"RADII"
   fmt='(5E16.8)'
   write(unit,form_format)fmt
   write(unit,fmt)(radii(i),i=1,natom)

   write(unit,flag_format)"SCREEN"
   write(unit,form_format)fmt
   write(unit,fmt)(screen(i),i=1,natom)

   if( has_psf_flag("CMAP") .and. CMAP_enabled )then
     !Only write CMAP related data if enabled
     write(unit,flag_format)"CHARMM_CMAP_COUNT"
     write(unit,'(a)')"%COMMENT  Number of CMAP terms, number of unique CMAP parameters"
     fmt='(2I8)'
     write(unit,form_format)fmt
     write(unit,fmt)cmap_term_count,cmap_type_count

     write(unit,flag_format)"CHARMM_CMAP_RESOLUTION"
     write(unit,'(a)')"%COMMENT  Number of steps along each phi/psi CMAP axis"
     write(unit,'(a)')"%COMMENT  for each CMAP_PARAMETER grid"
     fmt='(20I4)'
     write(unit,form_format)fmt
     write(unit,fmt)( needed_cmap_types(i)%resolution,i=1,cmap_type_count)

     !Loop over the unique CMAP parameters
     do i=1,cmap_type_count
       write(word,'(i2.2)')i
       write(unit,flag_format)"CHARMM_CMAP_PARAMETER_" // word
       write(unit,'(a,8(a4,1x))')"%COMMENT       ",needed_cmap_types(i)%label

       fmt='(8(F9.5))'
       write(unit,form_format)fmt
       write(unit,fmt) needed_cmap_types(i)%grid(&
                         1:needed_cmap_types(i)%resolution, &
                         1:needed_cmap_types(i)%resolution)

     enddo !i=1,cmap_type_count


     write(unit,flag_format)"CHARMM_CMAP_INDEX"
     write(unit,'(a)')"%COMMENT  Atom index i,j,k,l,m of the cross term"
     write(unit,'(a)')"%COMMENT  and then pointer to CHARMM_CMAP_PARAMETER_n"
     fmt='(6I8)'
     write(unit,form_format)fmt
     !Only write out i,j,k,l,m since these are the only ones needed
     !to calculate the dihedral angle between two adjacent dihedrals.
     ! hence (note the duplication of terms 2:4 and 5:7 
     !
     !  11      13      15      21      13      15      21      23     1
     ! becomes
     !
     !  11      13      15      21      23       1
     write(unit,fmt)(cmap_index(1:4,i),cmap_index(8:9,i),i=1,cmap_term_count)

     !End CMAP related data
   endif

   if ( ifbox > 0) then
     write(unit,flag_format)"SOLVENT_POINTERS"
     !write(unit,'(a)')"%COMMENT  TODO"
     fmt="(3i8)"
     write(unit,form_format)fmt
        !write(unit,fmt) IPTRES, NSPM, NSPSOL
        !FORMAT(12I6)  IPTRES, NSPM, NSPSOL
        !  IPTRES : final residue that is considered part of the solute,
        !           reset in sander and gibbs
        !  NSPM   : total number of molecules
        !  NSPSOL : the first solvent "molecule"

     !Note, PMEMD, uses NSPM (probably for next array (NSP) FLAG assignment),
     !but throws away IPTRES and NSPSOL. Hence we will set these to zero.
     write(unit,fmt) 0, total_number_of_molecules, 0 

     write(unit,flag_format)"ATOMS_PER_MOLECULE"
     fmt="(12i6)"
     write(unit,form_format)fmt
       !write(unit,fmt)(NSP(i), i=1,NSPM)
       !  NSP    : the total number of atoms in each molecule,
       !           necessary to correctly perform the pressure
       !           scaling.
     write(unit, fmt)(number_of_atoms_per_molecule(i), &
                                    i=1, total_number_of_molecules)

   endif



   return

end subroutine write_prmtop_atres


!--------------- READ COORDINATES -----------------------------
!
subroutine read_coordinates(u)

! This routine can process 3 types of files.
!
! 1) PDB File.
! 2) Charmm coordinate file.
! 3) Charmm restart file.
!
! In the case of PDB and coordinate files no velocity
! information is written. Additionally there is no
! support for determining box sizes and so if the
! user wants a periodic system they need to manually
! specify this.
!
! In the case of restart files both coordinate and
! velocity data is read and written to the restart
! file. Periodicity is not currently supported.
!
! The file type is auto detected.
!

  use psfprm, only: atom_label,outu,natom,title,inpcrd_unit,chmcrd_filename
  use psfprm, only: fmt00,fmtCor,box,ifbox
  use psf_strings, only:getwords
   implicit none
   integer, intent(in) :: u
   character(len=18) :: routine="<read_coordinates>"
   character(len=80),dimension(20) :: words
   character(len=128) :: line
   integer :: stat,i,atnum,numw
   real(kind=8), dimension(3,natom) :: crd
   real(kind=8), dimension(3,natom) :: vel
   real(kind=8), dimension(43) :: pbc
   real(kind=8) :: current_time
   logical :: crd_file = .false.
   logical :: pdb_file = .false.
   logical :: rst_file = .false.
   logical :: restart  = .false. !is this a restart file - write velocities?
   !Temporary variables used to pad the read format of the coordinate file
   integer      :: tmpI1,tmpI2
   character    :: tmpA1,tmpA2,tmpA3,tmpA4
   real(kind=8) :: tmpF1


! Determine the file type.

   write(outu,'(2a)') " Determining filetype of coordinate file: ", &
                       chmcrd_filename(1:len_trim(chmcrd_filename))

! 1) Is this a charmm coordinate file?
! 
!    If it is then it should begin with:
!
!    * TITLE1
!    * TITLE2
!    * ...
!    *
!    ATOMCOUNT
!
!  Check that the first line is a *
   rewind(u)
   stat=getwords(words,20,numw,u)
   if ( stat /= 2) then
     !We have not hit end of file.
     !Check if the first word is a *
     if ( index(words(1),"*") /= 0 ) then
       ! We found a * so this could be a crd file.
       ! Now read until we no longer get *'s
       do
         stat=getwords(words,20,numw,u)
         if ( stat == 2 ) exit !End of file
         if ( index(words(1),"*") == 0 ) then
            !This line does not begin with a *.
            !Break the loop and then we will 
            !check if it contains the correct atom count.
            exit
         end if
       end do
       !Check the atom count matches.
       read(words(1),'(I8)') atnum !Not right, should be a function of EXT?
       !read(words(1),fmt00) atnum !This should work, but it does not
       if ( atnum == natom ) then
         crd_file = .true. 
         restart = .false.
         write(outu,'(a)') " Identified CHARMM CRD File."
       else
         ! This is not a charmm crd file.
         crd_file = .false.
       end if
     else
       crd_file = .false.
     end if
   end if

! 2) Is this a charmm rst file?
!
!    If it is then it should begin with for example
!
!REST    35     1
!
!       3 !NTITLE followed by title
!* POLY PROLINE IN GAS PHASE
!* BY ROSS WALKER & MARK WILLIAMSON (SDSC)
!*  DATE:     1/23/ 9     15:49: 8      CREATED BY USER: rcw

!NATOM,NPRIV,NSTEP,NSAVC,NSAVV,JHSTRT,NDEGF,SEED,NSAVL
!          65        1000        1000          10           0           0         189 0.422033020000000D+09           0
!
!  Check that the first word of line 1 is "REST"

   if (.not. crd_file) then
     rewind(u)
     stat=getwords(words,20,numw,u)
     if ( stat /= 2) then
       !We have not hit end of file.
       !Check if the first word is REST
       if ( index(words(1),"REST") /= 0 ) then
         ! We found REST so this could be a RST file.
         ! Now look forward until we find the marker for NATOM
         !  !NATOM,NPRIV,NSTEP,NSAVC,NSAVV,JHSTRT,NDEGF,SEED,NSAVL
         ! This is actually a comment so we can't use the getwords code for this.

         if ( index(words(4),"CUBI") /= 0 ) then
           ! Update the box information with what we know from this keyword
           ! We are making the assumption that the box is cubic at this moment:
           ifbox = 1
           box%pbcType = 1 !Cubic

           box%alpha   = 90.0d0
           box%beta    = 90.0d0
           box%gamma   = 90.0d0
         endif
         !CHECK is word is something else and then quit with an error

         do 
           read(u,'(a128)',iostat=stat) line
           if ( stat < 0 ) exit !End of file
           if (line(3:7) == "NATOM") exit !We found the natom
         end do

         !Next line should be the atom count
         stat=getwords(words,20,numw,u)
         if ( stat == 2 ) then  !End of file
            rst_file = .false.
         else
            !Check the atom count matches.
            read(words(1),'(I8)') atnum
            if ( atnum == natom ) then
              rst_file = .true.
              restart = .true.
              write(outu,'(a)') " Identified CHARMM RST File."

              !The current release does not support CHARMM restart files.
              write(outu,'(a)')  "ERROR: CHARMM RST Files are not currently supported."
              call mexit(outu,1)
            else
              ! This is not a charmm rst file.
              rst_file = .false.
            end if
         end if
       else
         rst_file = .false.
       end if
     end if
   end if !if not crd_file

! 3) Is this a pdb file?
!
!    For the moment if it is not a rst or crd file assume it is a pdb file.
!
   if (.not. rst_file .and. .not. crd_file) then
     write(outu,'(a)') " Assuming PDB File."
     pdb_file = .true.
     restart = .false.
   end if

   if (crd_file) then
     do atnum=1,natom
       !Read in each coordinate line using a fixed format
       !This is used here because getwords cannot deal with lines
       !greater than 80 char at the moment and hence breaks when
       !reading in an EXT coordinate file 
       read(u,fmtCor) tmpI1,tmpI2,tmpA1,tmpA2,&
                      crd(1,atnum),crd(2,atnum),crd(3,atnum),&
                      tmpA3,tmpA4,tmpF1
     enddo

   else if (rst_file) then

     !Look for the current MD Time
     current_time = 0.0d0

     if (ifbox > 0 ) then
       !Search through the file and find " !CRYSTAL PARAMETERS"
       !Note getwords() ignores comments so we can't use it.
       rewind(u)
       do
          read(u,'(a128)',iostat=stat) line
          if ( stat < 0 ) then
            !End of file
            write(outu,'(2a)') &
              "FATAL ERROR: Reached end of file on CHARMM RST: ", &
            chmcrd_filename(1:len_trim(chmcrd_filename))
            write(outu,'(a)') &  
              "             While searching for  !CRYSTAL PARAMETERS"
            call mexit(outu,1)
          end if
          if (line(3:20) == "CRYSTAL PARAMETERS") exit !We found the FLAG
       end do

       !We are now at the position where the PBC dimensions should be
       !Loop over the number of atoms - it should be 3 coordinates per line
       do i = 1, 15, 3
          ! !CRYSTAL PARAMETERS
          ! 0.309000000000000D+02 0.000000000000000D+00 0.309000000000000D+02
          ! 0.000000000000000D+00 0.000000000000000D+00 0.309000000000000D+02
          ! 0.000000000000000D+00 0.000000000000000D+00 0.000000000000000D+00
          ! etc...
          read(u,'(3D22.15)',iostat=stat) pbc(i), pbc(i+1), pbc(i+2)
          if ( stat < 0 ) then
            !End of file
            write(outu,'(2a)') &
               "FATAL ERROR: Reached end of file on CHARMM RST: ", &
            chmcrd_filename(1:len_trim(chmcrd_filename))
            write(outu,'(a)')&
               "             While searching for PBC parameters "
            call mexit(outu,1)
          end if
       end do

       box%a = pbc(1)
       box%b = pbc(3)
       box%c = pbc(6) 

     endif !ifbox



     !Search through the file and find " !XOLD, YOLD, ZOLD"
     !Note getwords ignores comments so we can't use it.
     do
        read(u,'(a128)',iostat=stat) line
        if ( stat < 0 ) then
          !End of file
          write(outu,'(2a)') "FATAL ERROR: Reached end of file on CHARMM RST: ", &
          chmcrd_filename(1:len_trim(chmcrd_filename))
          write(outu,'(a)')  "             While searching for !XOLD, YOLD, ZOLD"
          call mexit(outu,1)
        end if
        if (line(3:18) == "XOLD, YOLD, ZOLD") exit !We found the FLAG
     end do

     !We are now at the position where the coordinates are.
     !Loop over the number of atoms - it should be 3 coordinates per line
     do i = 1, natom
        read(u,'(3D22.15)',iostat=stat) crd(1,i),crd(2,i),crd(3,i)
        if ( stat < 0 ) then
          !End of file
          write(outu,'(2a)') "FATAL ERROR: Reached end of file on CHARMM RST: ", &
          chmcrd_filename(1:len_trim(chmcrd_filename))
          write(outu,'(a,i6)')  "             While searching for XOLD, YOLD, ZOLD of atom: ",i
          call mexit(outu,1)
        end if
     end do

     !Next we need to read the velocities
     !Search through the file and find " !VX, VY, VZ"
     !Note getwords ignores comments so we can't use it.
     do
        read(u,'(a128)',iostat=stat) line
        if ( stat < 0 ) then
          !End of file
          write(outu,'(2a)') "FATAL ERROR: Reached end of file on CHARMM RST: ", &
          chmcrd_filename(1:len_trim(chmcrd_filename))
          write(outu,'(a)')  "             While searching for !VX, VY, VZ"
          call mexit(outu,1)
        end if
        if (line(3:12) == "VX, VY, VZ") exit !We found the FLAG
     end do

     !We are now at the position where the velocities are.
     !Loop over the number of atoms - it should be 3 velocities per line
     do i = 1, natom
        read(u,'(3D22.15)',iostat=stat) vel(1,i),&
                                        vel(2,i),&
                                        vel(3,i)
        if ( stat < 0 ) then
          !End of file
          write(outu,'(2a)') "FATAL ERROR: Reached end of file on CHARMM RST: ", &
          chmcrd_filename(1:len_trim(chmcrd_filename))
          write(outu,'(a,i6)')  "             While searching for VX, VY, VZ of atom: ",i
          call mexit(outu,1)
        end if
     end do

     !Rescale the velocities
     !THIS IS WIP; I cannot find the correct scaling; 
     !this value here was reverse engineered by comparing totKE 
     !between CHARMM outputs and AMBER restarts
     !I DO NOT BELIEVE THAT IT IS A SIMPLE SCALING FACTOR. THAT MAKES
     !NO SENSE AT ALL.
     do i = 1,natom
       !1773.6276
       vel(1,i) = vel(1,i) * 0.933136d0
       vel(2,i) = vel(2,i) * 0.933136d0
       vel(3,i) = vel(3,i) * 0.933136d0
     enddo
 
   else if (pdb_file) then
     rewind(u)
     atnum=0
     do 
        stat=getwords(words,20,numw,u)
        if ( stat == 2 ) exit !End of file
        if( (index(words(1),"ATOM") /= 0) .or. (index(words(1),"HETA") /= 0) ) then
           atnum=atnum+1

           !Check we don't exceed the number of atoms
           if (atnum > natom) then
              write(outu,'(a)')    'FATAL ERROR: Atom count in pdb file exceeds natom in psf file.'
              write(outu,'(a,i8)') '             natom = ',natom
              call mexit(outu,1)
           end if

           read(words(6),'(f12.5)')crd(1,atnum)
           read(words(7),'(f12.5)')crd(2,atnum)
           read(words(8),'(f12.5)')crd(3,atnum)
           if(verbose) write(outu,'(a,i9,3f12.5)') &
                routine,atnum,crd(1,atnum),crd(2,atnum),crd(3,atnum)
        endif
     enddo
     !Check we found the correct number of atoms
     if (atnum /= natom) then
        write(outu,'(a)')    'FATAL ERROR: Atom count in pdb does not match psf file.'
        write(outu,'(a,i8)') '             pdb atom count = ',atnum
        write(outu,'(a,i8)') '             psf atom count = ',natom
        call mexit(outu,1)
     end if
   end if

   if (restart) then
     !Read the PBC details
     !readAndSetPbcDimensions(box)
     write(inpcrd_unit,'(a37,a43)')"restrt generated from psfprm utility:",title(1:43)
     if (natom>99999) then
       write(inpcrd_unit,'(i6,e15.7)')natom,current_time
     else
       write(inpcrd_unit,'(i5,e15.7)')natom,current_time
     end if
     write(inpcrd_unit,'(6f12.7)')(crd(1,i),crd(2,i),crd(3,i),i=1,natom)
     write(inpcrd_unit,'(6f12.7)')(vel(1,i),vel(2,i),vel(3,i),i=1,natom)
     if (ifbox > 0 ) then
       write(inpcrd_unit,'(6f12.7)') box%a,    box%b,   box%c, &
                                     box%alpha,box%beta,box%gamma !box info
     endif
   else
     write(inpcrd_unit,'(a37,a43)') "inpcrd generated from psfprm utility:",title(1:43)
     if (natom>99999) then
       write(inpcrd_unit,'(i6)')natom
     else
       write(inpcrd_unit,'(i5)')natom
     end if
     write(inpcrd_unit,'(6f12.7)')(crd(1,i),crd(2,i),crd(3,i),i=1,natom)
     if (ifbox > 0 ) then
       write(inpcrd_unit,'(6f12.7)') box%a,    box%b,   box%c, &
                                     box%alpha,box%beta,box%gamma !box info
     endif
   end if

   return

end subroutine read_coordinates


subroutine print_prmtop_summary
   use psfprm, only:outu,nphih,mphia,mbona,nbonh,mtheta,ntheth
   implicit none
   character(len=80) :: print_fmt = '(A50,1X,I8)'

   write(outu,'(a)') ""
   write(outu,'(a)') &
   "==========================================================="
   write(outu,'(a)') "        Created prmtop summary"
   write(outu,'(a)') &
   "==========================================================="
   write(outu,'(a)') ""
   write(outu,print_fmt)  "Number of bonds with hydrogen:", nbonh
   write(outu,print_fmt)  "Number of bonds without hydrogen:", mbona
   write(outu,'(a)') ""
   write(outu,print_fmt)  "Number of angles with hydrogen:", ntheth
   write(outu,print_fmt)  "Number of angles without hydrogen:", mtheta
   write(outu,'(a)') ""
   write(outu,print_fmt)  "Number of dihedrals with hydrogen:", nphih
   write(outu,print_fmt)  "Number of dihedrals without hydrogen:", mphia
   write(outu,'(a)') &
   "==========================================================="
   write(outu,'(a)') ""


   !write(outu,print_fmt), "1-4 NB pairs",
   !write(outu,print_fmt), "1-4 Pairs",
   return

end subroutine print_prmtop_summary

!TO BE IMPLEMENTED FOR THE RESTART
subroutine readAndSetPbcDimensions(box)
 !CHARMM conveys information about the periodic boundary
 !conditions in a
 use psfprm, only: natom, boxType
 implicit none
 type(boxType), intent(out) :: box


 !Search for  !CRYSTAL PARAMETERS section in the restart

  box%pbcType = 1 !Cubic

  box%alpha   = 90.0
  box%beta    = 90.0
  box%gamma   = 90.0

  box%a   = 90.0
  box%b   = 90.0
  box%c   = 90.0

end subroutine readAndSetPbcDimensions



end module prmwrt
