! <compile=optimized>
#include "copyright.h"

!--------------------------------------------------------------
! CHAMBER: PSFPRM Module
!
! Description: Main PSF passing routines.
!
! Authors: Mark Williamson (SDSC, 2009)
!          Michael Crowley (NREL, 2009)
!          Ross Walker (SDSC, 2009)
!
!--------------------------------------------------------------

module psfprm
   implicit none

   integer,save :: chmff_verno
   character(len=60),save :: chmff_type

   logical,save :: CMAP_enabled !True if cmap terms should be included.
   logical,save :: want_vmd_prmtop = .false. ! Write an additional prmtop
                                ! that VMD can read
   logical,save :: verbose_debug=.false.,verbose_debug_vdw=.false.

   !Flags read in from the first line of the PSF: PSF, CMAP, CHEQ etc
   character(len=8), dimension(:), allocatable :: psf_flags 

   character(len=256) :: outu_name,psf_filename,prmtop_filename, &
                        param_filename, topology_filename, &
                        chmcrd_filename, chmrst_filename, inpcrd_filename, &
                        vmd_prmtop_filename
   integer,parameter :: maxwrd=40
   character(len=80),dimension(maxwrd) :: words


   integer, parameter :: psf_unit=20, &
                         outu=6, &
                         prmtop_unit=21, &
                         prm_unit=22, &
                         top_unit=23, &
                         chmcrd_unit=24, &
                         chmrst_unit=25, &
                         inpcrd_unit=26, &
                         vmd_prmtop_unit=27
              
   integer max_res_at
   integer :: idLen !Length of the CHARMM atom type name
   logical :: tip3_flex = .false.
   logical :: l_xplor_psf = .false.
   character(len=80) :: fmt00,fmt01,fmt03,fmt04,fmt05,fmt06,fmtCor

   character(len=80) :: title,title2

   integer, parameter :: maxattypes=400

   integer :: nbond,  ndon,   nacc,ngrp,nst2

   integer :: natom,  ntypes, nbonh,  mbona,  ntheth, mtheta
   integer :: nphih,  mphia,  nhparm, nparm,  nnb,    nres
   integer :: nbona,  ntheta, nphia,  numbnd, numang, nptra
   integer :: natyp,  nphb,   ipert,  nbper,  ngper,  ndper
   integer :: mbper,  mgper,  mdper,  ifbox,  nmxrs,  ifcap


   integer :: iptres,nspm,nspsol,natcap,nextra
   real(kind=8) :: cutcap,xcap,ycap,zcap

   integer :: nattypes,& !number of atom types in this system
              tip3_htype,tip3_otype
   character(len=8), dimension(maxattypes) :: att_txt
   character(len=8), dimension(:), allocatable :: &
    atom_label, &!atom label name from the PSF
    atom_type, &!atom type name from the PSF
    lresat,    &!residue name from the PSF
    lres,      &!
    attype_name !Character value of the atom type obtained from the CHARMM
                !topology file. e.g. :
                !
                !MASS    22 CT1   12.01100 C ! aliphatic sp3 C for CH
                !           ^
                !Runs over ntypes

   !Map between locattype and the atom's Z value
   !Populated by chemicalSymbolToZ() and this array is ultimately
   !need by the gb_radii assignment routines.
   integer, dimension(:), allocatable :: element

   integer, pointer, dimension(:) :: iresid,ipres,&
    iac, &      !numerical value of the atom type obtained from the PSF file
                !this is a function of atom index (natom)
    imove, &
    im,jm,km,lm, &
    idon,ihd1,iacc,iac1, &
    igpbs,igptyp,imoveg, &
    join,irotat,nsp, &
    attypes,&   !Numerical value of the atom type obtained from the CHARMM
                !topology file. e.g. :
                !
                !MASS    22 CT1   12.01100 C ! aliphatic sp3 C for CH
                !        ^
                !attypes will point the local atom type 
                !back to the charmm atom type
                !Runs over ntypes
    locattype,& !Map between natom value ---> attypes()
                !                            \--> attype_name()
                !                             \-> element()
                !Runs over natom
    attypeh     !flags whether it is a heavy atom type (2) or H type (1)

! bonds
   integer :: nbndtypes
   integer, pointer, dimension(:) :: ib,jb, &
                                         bndtypes,icb,ibtyp,jbtyp, &
                                         ibh,jbh,iba,jba,icbh,icba
   real(kind=8),pointer,dimension(:) :: rk,req

! angles
   integer :: nangtypes !numang is the amber total number of angle types
   integer, pointer, dimension(:) :: &
      it,jt,kt,   &     !Index to atoms i,j,k in angle term
      ict,        &     !Index to that angle term's parameter
                        !All four of these arrays extend as far as ntheta
      angtypes,   &
      ith,jth,kth,&     !Angles with hydrogen
      icth,       &     !Index to angles with hydrogen
                        !associated parameters
      ita,jta,kta,&     !Angles *without* hydrogen
      icta,       &     !Index to angles *without* hydrogen's 
                        !associated parameters
                        !The above four lines are all allocated as far
                        !as ntheta
      ittyp,jttyp,kttyp !Contain the integer local types of atoms in angle
                        !Is mapped to atom types via the attype_name() array
   real(kind=8),pointer,dimension(:) :: &
      tk,teq            !Parameters (r_k and r_eq) for the angle type.
                        !Extends as far as nangtypes

! angles: Urey Bradley
   integer, save                     :: &
      nub, nubtypes     !Number of Urey Bradley terms, number of UB types
   integer, pointer, dimension(:)    :: &
      ub_atm_i, ub_atm_k,& !Index of atom i,k in the UB term
      ub_idx               !Index into rub and kub for
                           !that UB term's parameters  
                           !All 3 should extend to nub
   real(kind=8),pointer,dimension(:) :: rub_tmp, kub_tmp
   real(kind=8),pointer,dimension(:) :: &
      rub, kub          !Parameters (r_k and r_eq) for the UB type.
                        !Extends as far as nubtypes


! dihedrals
   integer :: nphi, ndihtypes
   integer :: nextra_dih
   integer, pointer, dimension(:) :: ip,jp,kp,lp,icp, &
                                         iph,jph,kph,lph,icph, &
                                         iph2,jph2,kph2,lph2,icph2, &
                                         ipa2,jpa2,kpa2,lpa2,icpa2, &
                                         ipa,jpa,kpa,lpa,icpa,dihtypes, &
                                         idtyp,jdtyp,kdtyp,ldtyp

   integer :: Num_Expl_Dihed_Types_Assigned = 0
   !Counter for the logging output of assigned dihedral types

!impropers
   integer,save :: nimphi,nimprtypes
   real(kind=8),pointer,dimension(:),save :: pk_impr,phase_impr
   integer, pointer, dimension(:),save :: impropertypes, &
                                         ipm,jpm,kpm,lpm,icpm, &
                                         imtyp,jmtyp,kmtyp,lmtyp, &
                                         imp
                                         

   integer, parameter :: maxextradih=1000
   real(kind=8), dimension(maxextradih) :: pkx,phasex
   integer     , dimension(maxextradih) :: ipnx,dupx
   
! excluded atoms
   integer, pointer, dimension(:) :: inb,iblo,natex

   character(len=4), allocatable, dimension(:) ::labres,isymbl,itree

   real(kind=8),pointer, dimension(:) :: cg,amass,pk,pn, &
                                            phase,solty, &
                                            cn1,cn2,cn114,cn214, &
                                            asol,bsol, &
                                            hbcut

   real(kind=8),pointer, dimension(:) :: scnb_scale_factor,&
                                        scee_scale_factor

! GB radius set, radii and screen
   integer                           :: radius_set = 2 !(defaults to mbondi)
   character(len=1)                  :: radius_set_read
   character(len=80)                 :: radius_set_name
   real(kind=8),pointer, dimension(:) :: radii,screen


   
   real(kind=8),parameter :: pi=3.14159265358979323846d0, &
      two=2.d0

!PBC Box settings
   type boxType
     integer      :: pbcType !1 = octahedral 
                             !2 = truncated octahedral
     real(kind=8) :: a,b,c !box dimensions
     real(kind=8) :: alpha,beta,gamma !box dimensions
   end type boxType

   type(boxType)  :: box

!Other
   integer :: numlp,numlph

 contains
!-----------------------------------------------------------------
!------ GET ATOM_PARAMETERS---------------------------------------
!  1)  Read topology file to identify FF version and description
!  2)  Atom types are in the psf file as integers, but the parameters
!           for bonds, angles, nonbonds, etc, are identified in the parameter
!           file from the atom type strings. The map of the atom type string 
!           to the atom type integer is in the top of the rtf file that is a
!           partner to the parameter file. That is all we need from the 
!           topology file.
!       Read atom types out of topology file to map atom type integers to
!           atom type strings so atom types match
!           the correct parameters in the parameter file.
!  3)  Fill GB radii array
!  4)  Read the nonbond parameters for the atom types from the parameter file
!  5)  Create combined nonbond parameters
!
!  TODO:
!     1) Check version of topology file since some may not print the
!        element at the end of the mass array. Needed for GB radii.
!
!           22    1
!           MASS     4 HT     1.00800  H ! TIP3P water hydrogen
!
!           19    1                   ?????
!           MASS     4 HT     1.00800 ?????
!
!     2) Assign default radii for non supported elements and print
!        warning about GB not being supported for this psf due to
!        missing GB radii parameters. Change element identifying code
!        to warn instead of fail.
!
!        Can we assign a fake radii that will crash sander / PMEMD?
!
!-----------------------------------------------------------------
   subroutine get_atom_parameters
     use psf_strings
     implicit none
     real(kind=8),pointer,dimension(:) :: kb,b0
     character(len=80) :: line,word_type,word_elem,word
     character(len=8) :: tip3_htype_word,tip3_otype_word
     integer :: itype,n,nfound,linlen,wlen,i,j
     integer :: numwrd,ierr
     real(kind=8) :: e,r
     real(kind=8),pointer,dimension(:) :: ee,rr,ee14,rr14
     integer alen,wrdlen
     logical chartest,wrdmatch

     line="    "
     tip3_htype_word="HT      "
     tip3_otype_word="OT      "
     rewind(top_unit)
     rewind(prm_unit)
     allocate(attype_name(ntypes))
     allocate(element(ntypes))
     call allreal(ntypes*(ntypes+1)/2,cn1,cn2,cn114,cn214)
     call allreal(ntypes,kb,b0)

     !-------Read the CHARMM topology file to get FF version number  and name
     read(top_unit,'(a60)',end=990)chmff_type
     read(top_unit,'(a80)',end=990)line
     do while(line(1:1) == "*" ) 
        read(top_unit,'(a80)',end=990)line
     enddo
     read(line,'(i4)')chmff_verno

     !-------Read the CHARMM topology file to find the atom types it knows-----
     do while(line(1:4) /= "MASS" ) 
        read(top_unit,'(a80)',end=990)line
     enddo
     nfound = 0
     !------- look for the MASS lines and stop when reach the DECL lines
     do while ( (line(1:4) /= "DECL") .and.  (line(1:4) /= "RESI") )
        call next_word(line,wlen,linlen,word)
        if(word(1:wlen) == "MASS")then
           !MASS     1 H      1.00800 H ! polar H
           !OR
           !MASS     1 H      1.00800 H
           call next_int(line,linlen,itype)

           ! Check here to see if itype > maxattypes
           if ( itype > maxattypes) then
              write(outu,'(a,i5,a,i5)'),"chmtype", itype, &
                 " is greater than maxattypes hardlimit of ",maxattypes
              write(outu,'(a)'),"Recompile with a larger value if maxattypes"
              call err_exit("Cannot continue")
           endif

           call next_word(line,wlen,linlen,word_type)
           if(word_type(1:idlen) == tip3_htype_word(1:idlen))tip3_htype=itype
           if(word_type(1:idlen) == tip3_otype_word(1:idlen))tip3_otype=itype
           if(l_xplor_psf) then
              do n=1,ntypes
                 if(word_type(1:idlen) == att_txt(n)(1:idlen)) then
                    nfound=nfound+1
                    attypes(n)=itype 
                    attype_name(n)(1:idlen) = att_txt(n)(1:idlen)
                    !Get the element of the atom type
                    call next_word(line,wlen,linlen,word_elem)
                    !Hack to skip over mass value since next_float seems not to work
                    call next_word(line,wlen,linlen,word_elem)
                    element(n)=chemicalSymbolToZ(word_elem(1:2))

                    if(verbose_debug)write(outu,'(a,a,i3,a,i3,1x,a,i2)') &
                         "<get_atom_parameters> ", &
                         "  Found type in topology ",n, &
                         "  chmtype",attypes(n),attype_name(n)(1:idlen),element(n)
                 endif
              enddo
              do j=1,natom
                 do n=1,ntypes
                    if(atom_type(j)(1:idlen) == attype_name(n)(1:idlen)) iac(j)=attypes(n)
                 enddo
              enddo
           else
              do n=1,ntypes
                 if(itype == attypes(n))then
                    nfound=nfound+1
                    attype_name(n)(1:idlen)=word_type(1:idLen)
                    !Get the element of the atom type
                    call next_word(line,wlen,linlen,word)
                    !Hack to skip over mass value since next_float seems not to work
                    call next_word(line,wlen,linlen,word)
                    element(n)=chemicalSymbolToZ(word(1:2))

                    if(verbose_debug)write(outu,'(a,a,i3,a,i3,1x,a,i2)') &
                         "<get_atom_parameters> ", &
                         "  Found type in topology ",n, &
                         "  chmtype",attypes(n),attype_name(n)(1:idlen),element(n)
                 endif
              enddo
           endif
        endif
        read(top_unit,'(a80)',end=991)line
     enddo
     if(l_xplor_psf)then
        !--- This work was done in subroutine atom_types for non-xplor psf's
        do i=1,natom
           do j=1,ntypes
              if(iac(i) == attypes(j) ) then
                 locattype(i)=j
                 if(amass(i) < 1.5 ) then
                    if(attypeh(j).eq.2)then
                       write(outu,*)"<get_atom_parameters>", &
                            "ERROR atom with mass 1 has type already", &
                            "designated as non-hydrogen type"
                       write(outu,*)"<get_atom_parameters>"," atom ",i,amass(i),atom_label(i)
                       call err_exit("Cannot continue")
                    endif
                    attypeh(j)=1
                 else
                    if(attypeh(j).eq.1)then
                       write(outu,*)"<get_atom_parameters>", &
                            "ERROR atom with mass >1 has type already", &
                            "designated as non-hydrogen type"
                       write(outu,*)"<get_atom_parameters>"," atom ",i,amass(i),atom_label(i)
                       call err_exit("Cannot continue")
                    endif
                    attypeh(j)=2
                 endif
                 if(verbose_debug) &
                      write(outu,'(a,i10,a,i3,a,i3)') &
                      "<get_atom_parameters> atom ",i, &
                      "  chmtype ",attypes(j), &
                      "    local type ",locattype(i)
              endif
           enddo
        enddo
     endif
     if(nfound /= ntypes) goto 991
     if(verbose_debug)write(outu,'("<get_atom_parameters> found all types")')
     if(verbose_debug)write(outu,'("<get_atom_parameters> found tip3 H type",i6)') &
          tip3_htype
     if(verbose_debug)write(outu,'("<get_atom_parameters> found tip3 O type",i6)') &
          tip3_otype
     write(outu,'("   --- Found all ATOM TYPES in topology file"/)')

     !-------  Go get some GB radii, just default amber radii for now -----
     if(verbose_debug)write(outu,'("calling get_gbradii")')
     call get_gbradii(radius_set)

     !------------------------------------------------------------------
     !   Look for the nonbond parameters in the parameter file
     !
     !   read(prm_unit,'(a80)',end=992)line
     ierr = getwords(words,maxwrd,numwrd,prm_unit)
     if(ierr == 2 ) then
        write(outu,*) "ERROR end of file for reading prm_unit, retval=",ierr
        call err_exit("Cannot continue")
        stop
     endif

     !  look for line starting with NONBONDED, if not there or not in col 1
     !  we have a problem.
     do while(index(words(1),"NONBONDED") /= 1 )
        if(getwords(words,maxwrd,numwrd,prm_unit) == 2 ) &
             call err_exit("ERROR end of file for reading prm_unit, No NONBONDED")
     enddo
     if(getwords(words,maxwrd,numwrd,prm_unit) == 2 ) &
          call err_exit( &
          "ERROR end of file for reading prm_unit, No NONBONDED")
     if(verbose_debug) &
          write(outu,'(/"<get_atom_parameters> reading nonbonded parameters")')
     call allreal(ntypes,ee,rr,ee14,rr14)
     nfound=0
     ee=0.d0
     rr=0.d0
     ee14=0.d0
     rr14=0.d0
     do while(nfound < ntypes)
        if(numwrd >= 4)then
           wrdlen=len_trim(words(1))
           do n=1,ntypes
              alen=len_trim(attype_name(n))
              !DEBUG
              !  write(outu,*) words(1)(1:6)," " ,attype_name(n)(1:6), &
              !                 index(words(1),attype_name(n)(1:6))
              if(attype_name(n)(1:idLen) == words(1)(1:idLen))then
                 if(ee(n) == 0.d0)then
                    nfound=nfound+1
                    if(verbose_debug)write(outu, &
                         '("        Found specific ",a8)')attype_name(n)(1:idLen)
                 else
                    if(verbose_debug)write(outu, &
                         '("        Found specific replaces generic  ",a8)') &
                         attype_name(n)(1:idLen)
                 endif

                 read(words(3),"(f10.6)")ee(n)
                 read(words(4),"(f10.6)")rr(n)
                 if(numwrd >= 7) then
                    read(words(6),"(f10.6)")ee14(n)
                    read(words(7),"(f10.6)")rr14(n)
                 else
                    ee14(n)=ee(n)
                    rr14(n)=rr(n)
                 endif
              endif
           enddo
           !--- Check for wildcarding ----------------------
           if( (words(1)(wrdlen:wrdlen) == "%") .or. &
                (words(1)(wrdlen:wrdlen) == "*" )) then
              do n=1,ntypes
                 if(ee(n) == 0.d0)then
                    !write(outu,*)"Testing ",attype_name(n)(1:idLen)," against ",words(1)(1:idLen)
                    alen=len_trim(attype_name(n))
                    wrdmatch=.true.
                    i=0
                    do while(i<wrdlen)
                       i=i+1
                       chartest=words(1)(i:i) == attype_name(n)(i:i)
                       chartest=chartest .or. (words(1)(i:i) == "*")
                       chartest=chartest .or. (words(1)(i:i) == "%")
                       wrdmatch = wrdmatch .and. chartest
                       !write(outu,*)attype_name(n)(i:i)," ",words(1)(i:i),chartest,wrdmatch
                    enddo
                    !--- if last char of parameter spec is % then lengths must match
                    if(words(1)(i:i) == "%") then
                       wrdmatch = wrdmatch .and. (alen == wrdlen)
                    else
                       wrdmatch = wrdmatch .and. (wrdlen <= alen)
                    endif
                    if(wrdmatch)then
                       nfound=nfound+1
                       if(verbose_debug)write(6,*) &
                            "Found generic ",attype_name(n)(1:idLen),words(1)(1:idLen)
                       read(words(3),"(f10.6)")ee(n)
                       read(words(4),"(f10.6)")rr(n)
                       if(numwrd >= 7) then
                          read(words(6),"(f10.6)")ee14(n)
                          read(words(7),"(f10.6)")rr14(n)
                       else
                          ee14(n)=ee(n)
                          rr14(n)=rr(n)
                       endif
                    endif
                 endif
              enddo
           endif
        endif
        if(getwords(words,maxwrd,numwrd,prm_unit) == 2 ) then
           if(verbose_debug)write(outu,*) "Found ",nfound," types  out of ",ntypes
           do i = 1,ntypes
              if(ee(i) /= 0.d0 .and. verbose_debug)  &
                   write(outu,*)attype_name(i),ee(i),rr(i)
           enddo
           call err_exit( &
                "ERROR did not find all nonbond parameters")
        endif
     enddo
     write(outu,'("   --- Found all nonbond parameters")')

     !---- Create the combined parameters ------------------------------------
     if(verbose_debug) &
          write(outu,'(/"--- Creating combined nonbond parameters")')
     n=0
     do i=1,ntypes
        do j=1,i
           n=n+1
           e=sqrt(ee(i)*ee(j))
           r=rr(i)+rr(j)
           cn1(n)=e*r**12
           cn2(n)=two*e*r**6
           e=sqrt(ee14(i)*ee14(j))
           r=rr14(i)+rr14(j)
           cn114(n)=e*r**12
           cn214(n)=two*e*r**6
           if(verbose_debug_vdw) &
                write(outu,'("<get_atom_parameters> ",2a8," cn1,cn2  ",8e10.3)') &
                attype_name(i)(1:idLen),attype_name(j)(1:idLen), &
                cn1(n),cn2(n),cn114(n),cn214(n)
           if(verbose_debug_vdw) &
                write(outu,'("                      ",2a8," ee       ",8e10.3)') &
                attype_name(i)(1:idLen),attype_name(j)(1:idLen), &
                ee(i),ee(j),ee14(i),ee14(j)
           if(verbose_debug_vdw) &
                write(outu,'("                      ",2a8," rr       ",8e10.3)') &
                attype_name(i)(1:idLen),attype_name(j)(1:idLen), &
                rr(i),rr(j),rr14(i),rr14(j)
           !         if(verbose_debug_vdw) &
           !               write(outu,'("                      ",2a8," e,r      ",8e10.3)') &
           !               attype_name(i)(1:idLen),attype_name(j)(1:idLen), &
           !               e,r
        enddo
     enddo
     call deallocate_real(kb,b0,ee,rr,ee14,rr14)
     if(verbose_debug) &
          write(outu,'("   --- Assigned NONBONDED parameters ")')
     return
990  write(outu,*)"<get_atom_parameters>","ERROR, came to end of topology file without finding MASS lines"
     return
991  write(outu,'(a,a,/a,i3,a,i3)')"<get_atom_parameters> ", &
          "ERROR, topology file ends without finding all atom types,", &
          "number found:",nfound, &
          " needed  ntypes:",ntypes
     return
     !992  write(outu,*)"<get_atom_parameters>","ERROR, came to end of parameter file without finding NONBONDED"
     !     return
     !993  write(outu,*)"<get_atom_parameters>","ERROR, came to end of parameter file without finding all atom type vdW parameters"
     !     return

   end subroutine get_atom_parameters

!-----------------------------------------------------------------
!------ GET BONDED_PARAMS---------------------------------------
!   1) Get all the needed parameters from the CHARMM parameter file
!      Bonds, angles, dihedrals
!-----------------------------------------------------------------
   subroutine get_bonded_params()
     use psf_strings, only: getwords

     implicit none

     character(len=19) :: routine="<get_bonded_params>"
     character(len=80),dimension(10) :: words
     logical :: stillreading=.true., &
          lookingforsection=.true., &
          have_dih
     real(kind=8),dimension(ntypes,ntypes) :: b00,rk0
     real(kind=8),dimension(ntypes,ntypes,ntypes) :: tk0,teq0,rub0,kub0
     real(kind=8),dimension(ndihtypes,0:10) :: pk0,phase0
     integer     ,dimension(ndihtypes,0:10) :: ipn0
     integer itt,jtt,ktt,ltt,len,n,i,nfound,linenum
     integer, allocatable, dimension(:,:) :: newtyp
     real(kind=8),parameter :: PI = 3.1415926535897932384626433832795d0 

     real(kind=8) :: pk00,phase00
     integer :: ipn00,totdihtyp,totimprtyp,ntype,oldtyp,nphih_new,mphia_new
     integer :: numwrd

     logical :: wildCard=.false.
     ! Flag to indicate the the dihedral parameter being examined
     ! has come from a wildcard entry as opposed to an explicit one

     logical, allocatable, dimension(:) :: alreadyAssigned
     ! Flag to indicate if the dihedral type has been assigned
     ! a parameter

     logical, allocatable, dimension(:) :: alreadyAssignedExplicit
     ! Flag on a known dihedral type to indicate it has ALREADY
     ! been matched against an explicit dihedral parameter:
     !            A-B-C-D
     !
     ! as opposed to a wildcarded dihedral parameter:
     !
     !            X-B-C-X
     !
     ! This should be of length ndihtypes and its purpose is to stop
     ! wildcarded dihedral parameters displacing explicit
     ! dihedral parameters that have already been read in
     ! The rules governing this are covered in detail in CHARMM's
     ! parmfile.doc in section:
     ! "Rules for the use of multiple dihedrals in CHARMM24"
     ! (specifically point no.6)


     nfound=0
     linenum=0
     rewind(prm_unit)
     rk=-1.
     req=-1.
     if(verbose_debug)write(outu,*) 
     !============ BONDS ============================
     do while(lookingforsection)
        if(getwords(words,4,numwrd,prm_unit) == 2) call err_exit( &
             "<get_bonded_params> ERROR end of file param file A")
        if(words(1)(1:3) == "BON") then
           if(verbose_debug)write(outu,'(a,a,a,a)') routine, &
                "  found ",words(1)(1:len_trim(words(1))), " Section"
           lookingforsection=.false.
        endif
     enddo

     do while(stillreading)
        if(getwords(words,4,numwrd,prm_unit) == 2)call err_exit( &
             "<get_bonded_params> ERROR end of file in param file B ")
        linenum=linenum+1
        do n=1,nbndtypes
           itt=ibtyp(n)
           jtt=jbtyp(n)
           if(itt > jtt) call swap(itt,jtt,ktt,ltt)
           if(words(1)(1:idLen) == attype_name(itt)(1:idLen)) then
              if( words(2)(1:idLen) == attype_name(jtt)(1:idLen)) then
                 len=80
                 read(words(3),'(g20.12)')rk0(itt,jtt)
                 read(words(4),'(g20.12)')b00(itt,jtt)
                 if(verbose_debug)write(outu,'(a,i4,a,a8,1x,a8,1x,2f10.3,3i4)') &
                      routine,linenum," found bond ", &
                      attype_name(itt),attype_name(jtt), &
                      b00(itt,jtt),rk0(itt,jtt),itt,jtt,n
                 nfound=nfound+1
              endif
           elseif(words(1)(1:idLen) ==  attype_name(jtt)(1:idLen)) then
              if(words(2)(1:idLen) == attype_name(itt)(1:idLen)) then
                 len=80
                 read(words(3),'(g20.12)')rk0(itt,jtt)
                 read(words(4),'(g20.12)')b00(itt,jtt)
                 if(verbose_debug)write(outu,'(a,i4,a,a8,1x,a8,1x,2f10.3,3i4)') &
                      routine,linenum," found bond ", &
                      attype_name(itt),attype_name(jtt), &
                      b00(itt,jtt),rk0(itt,jtt), &
                      itt,jtt,n
                 nfound=nfound+1
              endif
           endif
        enddo
        stillreading = nfound < nbndtypes
     enddo
     do n=1,nbndtypes
        if( (b00(ibtyp(n),jbtyp(n)) < 0.d0 ) .or. &
             (rk0(ibtyp(n),jbtyp(n)) < 0.d0 ))then 
           if(verbose_debug)write(outu,'(2a,a8,1x,a8)')routine, &
                "Cannot fine bond type:", &
                attype_name(ibtyp(n)),attype_name(jbtyp(n))
        else
           req(n)=b00(ibtyp(n),jbtyp(n))
           rk(n)=rk0(ibtyp(n),jbtyp(n))
        endif
     enddo

     if(verbose_debug)write(outu,'(2a,i6,a,i6/)')routine," BONDS FOUND: ", &
          nfound," out of ", nbndtypes

     !============ ANGLES ============================
     if(nangtypes > 0 ) then
        rewind( prm_unit)   
        lookingforsection=.true.
        do while(lookingforsection)
           if(getwords(words,4,numwrd,prm_unit) == 2) call err_exit( &
                "<get_bonded_params> ERROR end of file param file A")
           if(words(1)(1:3) == "ANG") then
              if(verbose_debug)write(outu,'(a,a,a,a)') routine, &
                   " found ",words(1)(1:len_trim(words(1))), " Section"
              lookingforsection=.false.
           endif
        enddo

        stillreading=.true.

        nfound=0
        linenum=0
        tk0  = -1.d0
        teq0 = -1.d0
        rub0 =  0.d0
        kub0 =  0.d0
        do while(stillreading)
           if(getwords(words,7,numwrd,prm_unit) == 2) &
                call err_exit( &
                "<get_bonded_params> ERROR end param file reading ANGLES")
           linenum=linenum+1
           if(words(1)(1:3) == "DIH" .or. (words(1)(1:3) == "PHI")) then
              stillreading=.false.
           else
              do n=1,nangtypes 

                 itt=ittyp(n)
                 jtt=jttyp(n)
                 ktt=kttyp(n)
                 if(angle_match(words(1:3)(1:idLen),attype_name(itt), &
                      attype_name(jtt),attype_name(ktt))) then
                    len=80
                    read(words(4),'(g20.12)')tk0(itt,jtt,ktt)
                    read(words(5),'(g20.12)')teq0(itt,jtt,ktt)

                    !This means we've found a UB term within the ANGLE term:
                    if(numwrd > 5 ) read(words(6),'(g20.12)')kub0(itt,jtt,ktt) !Kub
                    if(numwrd > 6 ) read(words(7),'(g20.12)')rub0(itt,jtt,ktt) !S0
                    !This is denoted by an extra two terms (Kub, S0)
                    !on the angle parameter line
                    ! V(angle) = Ktheta(Theta - Theta0)**2
                    !
                    ! V(Urey-Bradley) = Kub(S - S0)**2
                    !
                    !!  atom types     Ktheta    Theta0   Kub     S0
                    !!CT1  CT1  CT1   53.350    111.00    8.00   2.56100 ! ALLOW ALI

                    if(verbose_debug)write(outu,'(2a,i5,a,3i5,1x,a8,1x,a8,1x,a8,4f10.3)')&
                         routine," ANGLE param line:", &
                         linenum,"  found angle ",&
                         itt,jtt,ktt,&
                         attype_name(itt), &
                         attype_name(jtt), &
                         attype_name(ktt), &
                         tk0(itt,jtt,ktt), &
                         teq0(itt,jtt,ktt), &
                         rub0(itt,jtt,ktt), &
                         kub0(itt,jtt,ktt)
                    nfound=nfound+1
                 endif !if(angle_match(words(1:3)(1:idLen),attype_name(itt)
              enddo
           endif
           stillreading =  stillreading .and. ( nfound < nangtypes )
        enddo
        if(verbose_debug)write(outu,'(2a,i6,a,i6)')routine," ANGLES FOUND: ", &
             nfound," out of ",nangtypes
        call allreal(nangtypes,tk,teq)

        call allreal(ntheta,rub_tmp,kub_tmp)
        rub_tmp = 0.d0
        kub_tmp = 0.d0

        do n=1,numang
           if(tk0(ittyp(n),jttyp(n),kttyp(n)) < 0.d0)then
              write(outu,*)"<get_bonded_params>","Cannot fine angle type:", &
                   attype_name(ittyp(n)),attype_name(jttyp(n)),attype_name(kttyp(n))
           else
              tk(n)=tk0(  ittyp(n),jttyp(n),kttyp(n))
              teq(n)=teq0(ittyp(n),jttyp(n),kttyp(n))*PI/180.d0
              rub_tmp(n)=rub0(ittyp(n),jttyp(n),kttyp(n))
              kub_tmp(n)=kub0(ittyp(n),jttyp(n),kttyp(n))
           endif

        enddo
        if ( nfound < nangtypes) call err_exit("Cannot continue, angle types missing from parameter file, check for lower case.")
        nfound=0
        linenum=0
     endif
     if(verbose_debug)write(outu,*)


     !Overallocate ub type related arrays, since we don't know nubtypes yet
     call allreal(numang,rub,kub)

     !Work out how many nubtypes there are
     nubtypes = 0
     do n=1,numang
        if( kub0(ittyp(n),jttyp(n),kttyp(n))  /= 0.d0 ) then
           nubtypes = nubtypes + 1
           rub(nubtypes) = rub0(ittyp(n),jttyp(n),kttyp(n))
           kub(nubtypes) = kub0(ittyp(n),jttyp(n),kttyp(n))
        endif
     enddo

     !Overallocate ub term related arrays, since we don't know nub yet
     call allint(ntheta,ub_atm_i,ub_atm_k,ub_idx)

     nub = 0
     do n=1, ntheta !Look up each angle term's kub value to
        !determine if it has an associated UB term
        if ( (kub_tmp( ict(n) ) /= 0.d0) .and. (ict(n) .ne. 0) ) then
           !the ict(n) .ne. 0 check is to guard against a situation
           !in the water box example where ict(n) always points to zero!
           if(verbose_debug)write(outu,*) "Adding nub term"
           if(verbose_debug)write(outu,*) "n, kub_tmp( ict(n) )",n,kub_tmp( ict(n) )
           nub = nub + 1
           ub_atm_i(nub) = it(n)
           ub_atm_k(nub) = kt(n)

           !link up mapping
           do i=1, nubtypes
              if ( (kub_tmp( ict(n) ) ==  kub(i) ) .and. &
                   (rub_tmp( ict(n) ) ==  rub(i) ) ) then
                 ub_idx(nub) = i
              endif
           enddo

           !write(outu,*) it(n),jt(n),kt(n),kub_tmp( ict(n) ),rub_tmp( ict(n) )
        endif
     enddo


     !These have served their purpose; deallocate
     if(associated(rub_tmp))deallocate(rub_tmp)
     if(associated(kub_tmp))deallocate(kub_tmp)


     !============ DIHEDRALS ============================
     if(ndihtypes > 0 ) then
        rewind( prm_unit)   
        lookingforsection=.true.
        do while(lookingforsection)
           if(getwords(words,4,numwrd,prm_unit) == 2) call err_exit( &
                "<get_bonded_params> ERROR end of file param file A")
           if(words(1)(1:3) == "DIH") then
              if(verbose_debug)write(outu,'(4a)') routine, &
                   " found ",words(1)(1:len_trim(words(1))), " Section"
              lookingforsection=.false.
           endif
        enddo !do while(lookingforsection)
        stillreading=.true.

        ! Allocate an array of True or False to indicate whether
        ! a corresponding known dihedral type has been already assigned 
        ! an explicit type. It should be of length ndihtypes
        allocate(alreadyAssignedExplicit(ndihtypes)) 
        allocate(alreadyAssigned(nphi)) 

        ! Assume all dihedrals types are wildcards since we want
        ! to be able to match them after we have exhausted matching 
        ! all explicit types
        alreadyAssignedExplicit=.false.
        alreadyAssigned=.false.

        nfound=0
        linenum=0
        pk0=-1.
        ipn0=-1
        phase0=-1.


        if(verbose_debug)then
           !Show a preview of what types we're going to be on 
           !the look out dihedral parameter-wise for as we step 
           !through the CHARMM prn file
           write(outu, '(32x,a9,2x,a20,a17)') "ndihtypes",&
                "dihedral atom types",&
                "charmm atom type"
           do n=1,ndihtypes      
              itt=idtyp(n)
              jtt=jdtyp(n)
              ktt=kdtyp(n)
              ltt=ldtyp(n)
              write(outu,'(2a,i6,4(1x,a8),2x,4(1x,i2)  )')routine,&
                   " Looking for type ",n,attype_name(itt), &
                   attype_name(jtt),attype_name(ktt), &
                   attype_name(ltt),&
                   attypes(itt),attypes(jtt),attypes(ktt),attypes(ltt)

           enddo !do n=1,ndihtypes
        endif !verbose_debug

        do while(stillreading)
           if(getwords(words,7,numwrd,prm_unit) == 2)call err_exit( &
                "<get_bonded_params> ERROR end param file reading DIHEDRALS")
           linenum=linenum+1

           if(words(1)(1:4) == "NBON" .or. words(1)(1:4) == "NONB" &
                .or. words(1)(1:4) == "IMPR" ) then
              stillreading=.false.
           else
              do n=1,ndihtypes
                 itt=idtyp(n)
                 jtt=jdtyp(n)
                 ktt=kdtyp(n)
                 ltt=ldtyp(n)
                 if(jtt > ktt)call swap(jtt,ktt,itt,ltt)
                 if(jtt == ktt .and. itt > ltt)call swap(jtt,ktt,itt,ltt)

                 !Are we dealing with a read-in parameter that is a wildcard?
                 ! X  -  CT1 - CT2 - X
                 wildCard = (words(1) == "X   " .and. words(4) == "X   ")

                 if(dihe_match(words(1:4)(1:idLen),attype_name(itt), &
                      attype_name(jtt),attype_name(ktt), &
                      attype_name(ltt))) then
                    read(words(5),'(g20.12)')pk00    !barrier height
                    read(words(6),'(i20)')ipn00      !perodicity
                    read(words(7),'(g20.12)')phase00 !phase angle

                    if( (ipn00 > 10) .or. (ipn00 < 1) ) then
                       write(outu,*)"<get_bonded_parameters>", &
                            "DIHEDRAL has multiplicity > 10", &
                            " cannot be handled in this translator"
                       write(outu,*)"<get_bonded_parameters>", &
                            words(1:4)(1:idLen), pk00,ipn00,phase00
                       call err_exit("Dihedral reading bombs")
                    endif

                    ! We do not want to overwrite an explicit dihedral
                    ! with a wildcarded one
                    if ( wildCard .and. alreadyAssignedExplicit(n) ) then
                       if (verbose_debug) then
                          write(outu,'(a20,i4,4x,4(a8,1x),5i3)') &
                               "Not overwriting",&
                               n,&
                               attype_name(itt), &
                               attype_name(jtt), &
                               attype_name(ktt), &
                               attype_name(ltt), &
                               itt,jtt,ktt,ltt
                          write(outu,*) "alreadyAssignedExplicit(n) ",alreadyAssignedExplicit(n)

                          write(outu,'(a20,8x,4(a8,1x),f7.3,i4,f8.3)') &
                               "with parameter  ",words(1:4),&
                               pk0(n,ipn00), &
                               ipn0(n,ipn00), &
                               phase0(n,ipn00)
                          write(outu,*) "wildCard ",wildCard
                          write(outu,'(a)')""
                       endif !verbose_debug

                       !not sure about this counter? Mike?
                       !nfound=nfound-1
                    else
                       if(verbose_debug) then
                          write(outu,'(a)') "----"
                          write(outu,'(a,i4)')&
                               "Processing dihedral type number: ",n

                          write(outu,'(a33,4(a8,1x))')&
                               "type: ",&
                               attype_name(itt), &
                               attype_name(jtt), &
                               attype_name(ktt), &
                               attype_name(ltt)

                          write(outu,'(a,2x,L1)') " Has it already been assigned? ",&
                               alreadyAssigned(n)

                          if (alreadyAssigned(n)) then
                             write(outu,'(a,2x,L1)') &
                                  " Has it been assigned with an explicit parameter?",&
                                  alreadyAssignedExplicit(n)

                             write(outu,'(a)') &
                                  " Currently assigned parameters:"
                             do i=1,10
                                if (ipn0(n,i) .ne. -1) then
                                   write(outu,'(a32,5x,i4,f7.3,f8.3)') &
                                        "periodicity,height,phase:",&
                                        ipn0(n,i),pk0(n,i),phase0(n,i)
                                endif
                             enddo

                             write(outu,'(a)')""
                          endif

                          write(outu,'(a)')""

                          write(outu,'(a)') &
                               " Selecting parameter from param file:"

                          write(outu,'(a33,4(a8,1x))') &
                               "type: ",words(1:4)

                          write(outu,'(a32,5x,i4,f7.3,f8.3)') &
                               "periodicity,height,phase:",&
                               ipn00,pk00,phase00



                          write(outu,'(a,2x,L1)') &
                               " Is this selected parameter a wildcard? ",wildCard
                          write(outu,'(a)')""

                          !Get rid of any previous wildcard assignments
                          if (              alreadyAssigned(n) &
                               .and. (.not. alreadyAssignedExplicit(n)) &
                               .and. (.not. wildCard) ) then 
                             write(outu,'(a)') &
                                  " Removing any existing stale parameters for this dihedral"
                             do i=1,10
                                ipn0(n,i)   = -1
                                pk0(n,i)    = -1
                                phase0(n,i) = -1
                             enddo

                          endif

                          !show what is being overwritten
                          if (alreadyassigned(n)) then
                             if (ipn0(n,ipn00) .eq. -1) then
                                write(outu,'(a)') &
                                     " **Appending this as a sub parameter"
                             else
                                write(outu,'(a)') &
                                     " **Overwriting a sub parameter with this"
                             endif
                          else
                             write(outu,'(a)') " **Assigning this parameter"
                          endif

                       endif !if(verbose_debug)

                       !Actually add the parameter found, into our array
                       pk0(n,ipn00)=pk00
                       ipn0(n,ipn00)=ipn00
                       phase0(n,ipn00)=phase00

                       if(verbose_debug) then

                          write(outu,'(a)') &
                               " Currently assigned parameters:"
                          do i=1,10
                             if (ipn0(n,i) .ne. -1) then
                                write(outu,'(a32,5x,i4,f7.3,f8.3)') &
                                     "periodicity,height,phase:",&
                                     ipn0(n,i),pk0(n,i),phase0(n,i)
                             endif
                          enddo
                       endif !if(verbose_debug)

                       alreadyAssigned(n) = .true.
                       ! and tag it as being explicit if it's not come from a wildcard
                       alreadyAssignedExplicit(n) = .not. wildCard

                       if ( alreadyAssignedExplicit(n) ) then
                          Num_Expl_Dihed_Types_Assigned = &
                               Num_Expl_Dihed_Types_Assigned + 1
                       endif
                       ! not only is this simpler but it avoids a bug when
                       ! alreadyAssignedExplicit(n) was not initialized to false

                       nfound=nfound+1

                       if(verbose_debug) then
                          !we've done on this cycle
                          write(outu,'(a,i4)') &
                               "Finished processing dihedral type number: ",n
                          write(outu,'(a)') "----"
                          write(outu,'(a)') ""
                          write(outu,'(a)') ""
                       endif !if(verbose_debug)

                    endif !if (wildCard==.true. .and. alreadyAssignedExplicit(n)==.true.)
                 endif !if(dihe_match(words(1:4)(1:idLen),attype_name(itt), &
              enddo !do n=1,ndihtypes
           endif !if(words(1)(1:4) == "NBON" .or. words(1)(1:4) == "NONB" &
        enddo !do while(stillreading)
        if(verbose_debug) then
           write(outu,'(2a,i6,a,i6/)')routine," DIHEDRALS FOUND: ", &
                nfound," out of ",ndihtypes
        endif !verbose_debug
        totdihtyp=0
        do n=1,ndihtypes
           have_dih = .false.
           do i=1,10
              if(ipn0(n,i) > 0 )totdihtyp=totdihtyp+1
              have_dih = have_dih .or. (ipn0(n,i)>0)
           enddo
           if(.not. have_dih)then
              write(outu,'(i6,a,i5,4a8)') &
                   linenum," MISSING dihedral: ", n, &
                   attype_name(locattype(idtyp(n))), &
                   attype_name(jdtyp(n)), &
                   attype_name(kdtyp(n)), &
                   attype_name(ldtyp(n))
              call err_exit("Cannot continue")
           endif
        enddo
     endif !if(ndihtypes > 0 ) then
     !--- all dihedrals are accounted for now make a new list -------------------
     allocate(newtyp(ndihtypes,10) )
     newtyp=0
     ntype=0
     call allreal(nfound,pk,pn,phase)
     call allreal(nfound,scnb_scale_factor,scee_scale_factor)

     ! Since improper terms are in a separate list, we can assign
     ! a CHARMM default value of 1.0 for the scnb and scee scale factor
     ! to all dihedral types
     scnb_scale_factor = 1.0
     scee_scale_factor = 1.0

     do n=1,ndihtypes
        do i=1,10
           if(ipn0(n,i) > 0 ) then
              ntype=ntype+1
              newtyp(n,i)=ntype
              pk(ntype)=pk0(n,i)
              pn(ntype)=real(ipn0(n,i))
              phase(ntype)=pi*phase0(n,i)/180.d0
           endif
        enddo
     enddo

     nphih_new=0
     do n=1,nphih
        oldtyp=icph(n)
        do i=1,10
           if(newtyp(oldtyp,i)>0)then
              nphih_new=nphih_new+1
           endif
        enddo
     enddo

     call allint(nphih_new,iph2,jph2,kph2,lph2,icph2)

     iph2=0; jph2=0; kph2=0; lph2=0; icph2=0
     nphih_new=0
     do n=1,nphih
        oldtyp=icph(n)
        do i=1,10
           if(newtyp(oldtyp,i)>0)then
              nphih_new=nphih_new+1
              iph2(nphih_new)=iph(n)
              jph2(nphih_new)=jph(n)
              kph2(nphih_new)=kph(n)
              lph2(nphih_new)=lph(n)
              icph2(nphih_new)=newtyp(icph(n),i)
           endif
        enddo
     enddo

     if(verbose_debug)write(outu,'(/2a,2i6)')routine,"NPHIH ",nphih,nphih_new

     mphia_new=0
     do n=1,mphia
        oldtyp=icpa(n)
        do i=1,10
           if(newtyp(oldtyp,i)>0)then
              mphia_new=mphia_new+1
           endif
        enddo
     enddo

     call allint(mphia_new,ipa2,jpa2,kpa2,lpa2,icpa2)
     ipa2=0; jpa2=0; kpa2=0; lpa2=0; icpa2=0
     mphia_new=0
     do n=1,mphia
        oldtyp=icpa(n)
        do i=1,10
           if(newtyp(oldtyp,i)>0)then
              mphia_new=mphia_new+1
              ipa2(mphia_new)=ipa(n)
              jpa2(mphia_new)=jpa(n)
              kpa2(mphia_new)=kpa(n)
              lpa2(mphia_new)=lpa(n)
              icpa2(mphia_new)=newtyp(icpa(n),i)
           endif
        enddo
     enddo

     if(verbose_debug)write(outu,'(2a,2i6)')routine,"MPHIA ",mphia,mphia_new
     deallocate(newtyp)
     mphia=mphia_new
     nphih=nphih_new
     if(verbose_debug)write(outu,'(2a,2i6//)')routine,"Dih TYPES ",nptra,ntype
     nptra=ntype
     ndihtypes=ntype


     !Reverses the order of the four integers (i,j,k,l) that make up a 
     !torsion index to preempt tag_14_duplicates() writing a negative zero on k
     !With hydrogens
     call reorder_dihedrals_to_remove_0_k(nphih,iph2,jph2,kph2,lph2)
     !Without hydrogens
     call reorder_dihedrals_to_remove_0_k(mphia,ipa2,jpa2,kpa2,lpa2)

     !Identification of duplicate 1-4's due to multiple periodicities and or being part
     !of a ring system. This requires duplicates to have the 3rd entry of their atom id
     !set negative. They are then skipped in the 1-4 calculation.
     !With hydrogens
     call tag_14_duplicates(nphih,iph2,kph2,lph2,nbonh,ibh,jbh,ntheth,ith,kth)
     !Without hydrogens
     call tag_14_duplicates(mphia,ipa2,kpa2,lpa2,mbona,iba,jba,mtheta,ita,kta)

     !Clean up
     if(allocated(alreadyAssignedExplicit)) deallocate(alreadyAssignedExplicit)
     if(allocated(alreadyAssigned)) deallocate(alreadyAssigned)

     !============ IMPROPERS ============================
     if(nimprtypes > 0 ) then
        rewind( prm_unit)   
        lookingforsection=.true.
        do while(lookingforsection)
           if(getwords(words,4,numwrd,prm_unit) == 2) call err_exit( &
                "<get_bonded_params> ERROR end of file param file A")
           if(words(1)(1:4) == "IMPR") then
              if(verbose_debug)write(outu,'(4a)') routine, &
                   " found ",words(1)(1:len_trim(words(1))), " Section"
              lookingforsection=.false.
           endif
        enddo
        stillreading=.true.

        nfound=0
        linenum=0
        pk0=-1.
        if(nimprtypes > ndihtypes) &
             call err_exit("<get_bonded_params> impr types > ndihtypes, need to rewrite")
        ipn0=-1
        phase0=-1.


        if(verbose_debug)then
           do n=1,nimprtypes      
              itt=imtyp(n)
              jtt=jmtyp(n)
              ktt=kmtyp(n)
              ltt=lmtyp(n)
              write(outu,'(2a,i6,2x,4a)')routine," Looking for type ",n,attype_name(itt), &
                   attype_name(jtt),attype_name(ktt),attype_name(ltt)
           enddo
        endif

        do while(stillreading)
           if(getwords(words,7,numwrd,prm_unit) == 2)call err_exit( &
                "<get_bonded_params> ERROR end param file reading DIHEDRALS")
           linenum=linenum+1
           if(words(1)(1:4) == "NBON" .or. words(1)(1:4) == "NONB" &
                .or. words(1)(1:4) == "CMAP" ) then
              stillreading=.false.
           else
              do n=1,nimprtypes
                 itt=imtyp(n)
                 jtt=jmtyp(n)
                 ktt=kmtyp(n)
                 ltt=lmtyp(n)
                 if(jtt > ktt)call swap(jtt,ktt,itt,ltt)
                 if(jtt == ktt .and. itt > ltt)call swap(jtt,ktt,itt,ltt)
                 if(impr_match(words(1:4)(1:idLen),attype_name(itt), &
                      attype_name(jtt),attype_name(ktt), &
                      attype_name(ltt))) then

                    if(verbose_debug)write(outu,'(2a,i6,2x, 9a)')routine,"found ",n, &
                         words(1:4)(1:idLen),"  matches  ",attype_name(itt), &
                         attype_name(jtt),attype_name(ktt), &
                         attype_name(ltt)

                    read(words(5),'(g20.12)')pk00
                    read(words(6),'(i20)')ipn00
                    read(words(7),'(g20.12)')phase00
                    if( ipn00 /= 0 ) then
                       write(outu,'(2a,i3,a,a)')routine, &
                            "IMPROPER has multiplicity",ipn00," (should be 0)", &
                            " Error reading improper parameters"
                       write(outu,*) "line ",linenum," of param file improper section"
                       write(outu,*)"      ",words(1:4)(1:idLen), pk00,ipn00,phase00
                       call err_exit("Bomb during reading of Improper dihedral parameters")
                    endif

                    if(ipn0(n,ipn00) >= 0)then
                       write(outu,'(a,a,4a,f8.3,i4,f8.3,a,f8.3,i4,f8.3)')routine, &
                            " found duplicate dihedral:", &
                            words(1:4)(1:4), pk00,ipn00,phase00, " had", &
                            pk0(n,ipn00), ipn00,phase0(n,ipn00)
                       write(outu,'(2a)')routine,"    Replacing with latest read values"
                       nfound=nfound-1
                    endif

                    pk0(n,ipn00)=pk00
                    ipn0(n,ipn00)=ipn00
                    phase0(n,ipn00)=phase00
                    if(verbose_debug)write(outu,'(2a,4(a4,1x),f7.3,i4,f8.3,5i3)') &
                         routine," found improper ", &
                         attype_name(itt), &
                         attype_name(jtt), &
                         attype_name(ktt), &
                         attype_name(ltt), &
                         pk0(n,ipn00), &
                         ipn0(n,ipn00), &
                         phase0(n,ipn00), &
                         itt,jtt,ktt,ltt !,n
                    nfound=nfound+1
                 endif
              enddo
           endif
        enddo
        if(verbose_debug)write(outu,'(/2a,i6,a,i6/)')routine," IMPROPERS FOUND: ", &
             nfound," out of ",nimprtypes
        if(nfound /= nimprtypes) &
             call err_exit("Bomb reading improper parameters nfound /= nimprtypes")
        totimprtyp=0
        do n=1,nimprtypes
           have_dih = .false.
           if(ipn0(n,0) == 0 )totimprtyp=totimprtyp+1
           have_dih = have_dih .or. (ipn0(n,0) == 0)

           if(.not. have_dih)then
              write(outu,'(a,i6,a,i5,2x,4a4)') routine, &
                   linenum," MISSING dihedral: ", n, &
                   attype_name(imtyp(n)), &
                   attype_name(jmtyp(n)), &
                   attype_name(kmtyp(n)), &
                   attype_name(lmtyp(n))
              call err_exit("Cannot continue")
           endif
        enddo
     endif

     !--- all impropers are accounted for now make a new list -------------------
     call allreal(nfound,pk_impr,phase_impr)
     do n=1,nimprtypes
        pk_impr(n)=pk0(n,0)
        phase_impr(n)=pi*phase0(n,0)/180.d0
     enddo

     if(verbose_debug)write(outu,'(2a,2i6)')routine," Impropers, TYPES ",nimphi,nimprtypes

     return
   end subroutine get_bonded_params



!-----------------------------------------------------------------
!---------------MATCH FUNCTIONS----------------------------------------------
!-----------------------------------------------------------------
logical function angle_match(words,ati,atj,atk) result(match)
   implicit none
   character(len=idLen),dimension(3),intent(in) :: words
   character(len=idLen),intent(in) :: ati,atj,atk
   
   match=.false.
   match = match .or. ((words(1) == ati) .and. &
                      ((words(2) == atj) .and. &
                       (words(3) == atk)))
   match = match .or. ((words(1) == atk) .and. &
                      ((words(2) == atj) .and. &
                       (words(3) == ati)))
   return
end function angle_match

logical function dihe_match(words,ati,atj,atk,atl) result(match)
   implicit none
   character(len=idLen),dimension(4),intent(in) :: words
   character(len=idLen),intent(in) :: ati,atj,atk,atl
   
   match=.false.
   match = match .or. ((words(2) == atj) .and. &
                      ((words(3) == atk) .and. &
                      ((words(1) == "X   ") .and. &
                       (words(4) == "X   "))))
   match = match .or. ((words(2) == atk) .and. &
                      ((words(3) == atj) .and. &
                      ((words(1) == "X   ") .and. &
                       (words(4) == "X   "))))
   match = match .or. ((words(2) == atj) .and. &
                      ((words(3) == atk) .and. &
                      ((words(1) == ati) .and. &
                       (words(4) == atl))))
!   match = match .or. ((words(2) == atj) .and. &
!                      ((words(3) == atk) .and. &
!                      ((words(1) == atl) .and. &
!                       (words(4) == ati))))
!   match = match .or. ((words(2) == atk) .and. &
!                      ((words(3) == atj) .and. &
!                      ((words(1) == ati) .and. &
!                       (words(4) == atl))))
   match = match .or. ((words(2) == atk) .and. &
                      ((words(3) == atj) .and. &
                      ((words(1) == atl) .and. &
                       (words(4) == ati))))
   return
end function dihe_match


logical function dih_nmatch(i,j,k,l,ii,jj,kk,ll) result(match)
   implicit none
   integer, intent(in) :: i,j,k,l,ii,jj,kk,ll
   
   match=.false.
   match = match .or. ((j == jj) .and. &
                      ((k == kk) .and. &
                      ((i == ii) .and. &
                       (l == ll))))
   match = match .or. ((j == jj) .and. &
                      ((k == kk) .and. &
                      ((l == ii) .and. &
                       (i == ll))))
   match = match .or. ((k == jj) .and. &
                      ((j == kk) .and. &
                      ((i == ii) .and. &
                       (l == ll))))
   match = match .or. ((k == jj) .and. &
                      ((j == kk) .and. &
                      ((l == ii) .and. &
                       (i == ll))))
   return
end function dih_nmatch

logical function impr_match(words,ati,atj,atk,atl) result(match)
   implicit none
   character(len=idLen),dimension(4),intent(in) :: words
   character(len=idLen),intent(in) :: ati,atj,atk,atl
   
   match=.false.
   match = match .or. ((words(1) == ati) .and. &
                      ((words(4) == atl) .and. &
                      ((words(2) == "X   ") .and. &
                       (words(3) == "X   "))))
   match = match .or. ((words(1) == atl) .and. &
                      ((words(4) == ati) .and. &
                      ((words(2) == "X   ") .and. &
                       (words(3) == "X   "))))
   match = match .or. ((words(1) == ati) .and. &
                      ((words(4) == atl) .and. &
                      ((words(2) == atj) .and. &
                       (words(3) == atk))))
!   match = match .or. ((words(2) == atj) .and. &
!                      ((words(3) == atk) .and. &
!                      ((words(1) == atl) .and. &
!                       (words(4) == ati))))
!   match = match .or. ((words(2) == atk) .and. &
!                      ((words(3) == atj) .and. &
!                      ((words(1) == ati) .and. &
!                       (words(4) == atl))))
   match = match .or. ((words(1) == atl) .and. &
                      ((words(4) == ati) .and. &
                      ((words(2) == atk) .and. &
                       (words(3) == atj))))
   return
end function impr_match





!-----------------------------------------------------------------
!------ GET PARAMETERS---------------------------------------
!-----------------------------------------------------------------
subroutine get_parameters
   implicit none
   rewind(prm_unit)
   call allreal(nbndtypes,rk,req)
   

end subroutine get_parameters


!-----------------------------------------------------------------
!--------- EXCLUDED_ATOMS----------------------------------------
!-----------------------------------------------------------------
subroutine excluded_atoms

   return
end subroutine excluded_atoms

!-----------------------------------------------------------------
!------ IMPROPER_TYPES-----------------------------------------------------
!-----------------------------------------------------------------
!  go through atoms and collect dihedral types
subroutine improper_types
   implicit none
   logical, allocatable, dimension(:,:,:,:) :: limpt
   integer :: n,i,j,k,l,ityp,jtyp,ktyp,ltyp,itl,jtl,ktl,ltl,dt=0
   integer :: it
   logical :: l1,l2,have_impr
   
   allocate(limpt(nattypes,nattypes,nattypes,nattypes))
   if(verbose_debug)write(outu,'(/"<impr_types>",i6)')nimphi
   limpt=.false.
   nimprtypes=0

  !------------------------------------------------------------------------------
  !------- Determine how many dihedral types there are and set mask of the types 
  !   These are set in this loop:
  !      dt = number of dihe types in this system
  !      limpt(i,l,j,k) logical true if this type dihedral exists in this system
   do n=1,nimphi
      
      i=im(n)
      j=jm(n)
      k=km(n)
      l=lm(n)
      ityp=iac(i)
      jtyp=iac(j)
      ktyp=iac(k)
      ltyp=iac(l)
      itl=locattype(i)
      jtl=locattype(j)
      ktl=locattype(k)
      ltl=locattype(l)
      
      if(itl > ltl)then
         call swap(jtl,ktl,jtyp,ktyp)
         call swap(itl,ltl,ityp,ltyp)
      endif
      if( (itl == ltl) .and. (jtl > ktl))then
         call swap(jtl,ktl,jtyp,ktyp)
      endif
      if(.not.limpt(itl,ltl,jtl,ktl) )then
         dt=dt+1
         if(verbose_debug)write(outu,'(a,i3,":",10i4,4(1x,a4))') &
              "imp search",dt,itl,ltl,jtl,ktl,i,l,j,k
      endif
      limpt(itl,ltl,jtl,ktl)=.true.
   enddo

   call allint(dt,impropertypes,imtyp,jmtyp,kmtyp,lmtyp)
   call allint(nimphi,imp)
   impropertypes=0; imtyp=0; jmtyp=0; kmtyp=0; lmtyp=0; imp=0   

  !---- Assign improper type to each improper

   do n=1,nimphi
      i=im(n)
      j=jm(n)
      k=km(n)
      l=lm(n)
      ityp=iac(i)
      jtyp=iac(j)
      ktyp=iac(k)
      ltyp=iac(l)
      itl=locattype(i)
      jtl=locattype(j)
      ktl=locattype(k)
      ltl=locattype(l)
      
      if(itl > ltl)then
         call swap(jtl,ktl,jtyp,ktyp)
         call swap(itl,ltl,ityp,ltyp)
      endif
      if( (itl == ltl) .and. (jtl > ktl))then
         call swap(jtl,ktl,jtyp,ktyp)
      endif
      have_impr=.false.
      check_known: do it = 1,nimprtypes
         l1 = (itl == imtyp(it)) .and.  (ltl == lmtyp(it))
         l2 = (jtl == jmtyp(it)) .and.  (ktl == kmtyp(it))
         if( l1 .and. l2)then
            imp(n)=it
            have_impr= .true.
            exit check_known
         endif
      enddo check_known
      if(.not. have_impr)then
         nimprtypes= nimprtypes+1
         imtyp(nimprtypes)=itl
         jmtyp(nimprtypes)=jtl
         kmtyp(nimprtypes)=ktl
         lmtyp(nimprtypes)=ltl
         imp(n)=nimprtypes
      endif
   enddo
   if (verbose_debug) then
    write(outu,'(1a8,i8,1a16)') "Found ",nimprtypes," improper types"
   endif !verbose_debug

  return
end subroutine improper_types

!--------------------------------------------------------------------------+
!------ DIHEDRAL_TYPES-----------------------------------------------------|
!--------------------------------------------------------------------------+
!  go through atoms and collect dihedral types
subroutine dihe_types
   implicit none
   logical, dimension(nattypes,nattypes,nattypes,nattypes) :: diht
   integer, dimension(nattypes,nattypes,nattypes,nattypes) :: idiht
   integer :: n,i,j,k,l,ityp,jtyp,ktyp,ltyp,itl,jtl,ktl,ltl,dt=0
   integer :: it
   logical :: dih_h,have_dihe,l1,l2

   if(verbose_debug)write(outu,'(/a)')"<dihe_types>"
   diht=.false.
   idiht=0
   ndihtypes=0
  !------------------------------------------------------------------------------
  !------- Determine how many dihedral types there are and set mask of the types 
  !   These are set in this loop:
  !      dt = number of dihe types in this system
  !      diht(i,j,k,l) logical true if this type dihedral exists in this system
   do n=1,nphi
      i=ip(n)
      j=jp(n)
      k=kp(n)
      l=lp(n)
      ityp=iac(i)
      jtyp=iac(j)
      ktyp=iac(k)
      ltyp=iac(l)
      itl=locattype(i)
      jtl=locattype(j)
      ktl=locattype(k)
      ltl=locattype(l)
      if(jtl > ktl)then
         call swap(jtl,ktl,jtyp,ktyp)
         call swap(itl,ltl,ityp,ltyp)
      endif
      if( (jtl == ktl) .and. (itl > ltl))then
         call swap(itl,ltl,ityp,ltyp)
      endif
      if(.not.diht(jtl,ktl,itl,ltl) )then
         dt=dt+1
         !write(outu,'("phi search ",i3,":",4i4,4(1x,a4))') &
         !     dt,itl,jtl,ktl,ltl  , &
         !           attype_name(itl), &
         !           attype_name(jtl),attype_name(ktl), &
         !           attype_name(ltl)
      endif
      diht(jtl,ktl,itl,ltl)=.true.
   enddo

  !------------------------------------------------------------------------------
  !------   Assign dihedral type to each dihedral
   call allint(dt,dihtypes,idtyp,jdtyp,kdtyp,ldtyp)
   call allint(nphi,icp)
   dihtypes=0; idtyp=0; jdtyp=0; kdtyp=0; ldtyp=0; icp=0
   do n=1,nphi
      i=ip(n)
      j=jp(n)
      k=kp(n)
      l=lp(n)
      ityp=iac(i)
      jtyp=iac(j)
      ktyp=iac(k)
      ltyp=iac(l)
      itl=locattype(i)
      jtl=locattype(j)
      ktl=locattype(k)
      ltl=locattype(l)
      
      if(jtl > ktl)then
         call swap(jtl,ktl,jtyp,ktyp)
         call swap(itl,ltl,ityp,ltyp)
      endif
      if( (jtl == ktl) .and. (itl > ltl))then
         call swap(itl,ltl,ityp,ltyp)
      endif
      have_dihe=.false.
      !Look through our current list of known dihedral types
      check_known: do it = 1,ndihtypes 
         l1 = (itl == idtyp(it)) .and.  (ltl == ldtyp(it))
         l2 = (jtl == jdtyp(it)) .and.  (ktl == kdtyp(it))
         if( l1 .and. l2)then
            !Map the dihedral to its type
            icp(n)=it
            have_dihe= .true.
            exit check_known
         endif
      enddo check_known
      !We don't have this dihedral type, hence add it
      if(.not. have_dihe)then
         ndihtypes= ndihtypes+1
         idtyp(ndihtypes)=itl
         jdtyp(ndihtypes)=jtl
         kdtyp(ndihtypes)=ktl
         ldtyp(ndihtypes)=ltl
         icp(n)=ndihtypes
      endif
   enddo


   if(verbose_debug) then
    !Dump the list of the known dihedral types atom type indexes and associated types
    write(6,'(a)')""
    write(outu,'(1a30,i8,1a30)')&
    "List of the         ",ndihtypes," known dihedral types          "
    write(outu,'(a)')&
    "===================================================================="
    write(outu, '(1a16,4a8,1a16)')&
     "Dihedral type     ",&
     "idtyp()",&
     "jdtyp()",&
     "kdtyp()",&
     "ldtyp()",&
     "   Atom types     "

    write(outu, '(1a16,4a8)')&
     "    index         "

    write(outu,'(a)')&
    "===================================================================="

    write(outu, '(i8,4x,4i8,4x, a8,1x,a8,1x,a8,1x,a8 )' )&
    (i,idtyp(i),jdtyp(i),kdtyp(i),ldtyp(i),&
    attype_name(idtyp(i)),attype_name(jdtyp(i)),&
    attype_name(kdtyp(i)),attype_name(ldtyp(i)),&
    i=1,ndihtypes) 
    write(6,'(a)')""

    write(outu,'(1a30,i8,1a30)')&
    "icp(): Mapping between the ",nphi,"dihedrals in system"
    write(outu,'(1a30,i8,1a30)')&
    "and the ",ndihtypes," known dihedral types"
    write(outu,'(a)')&
    "=================================="
    write(outu, '(1a16,1a16)')&
    !1234567890123456
     "  Dihedral      ",&
     "Dihedral index  "
    write(outu,'(a)')&
    "=================================="

    write(outu,'(i8,4x,1i8)') (i,icp(i),i=1,nphi)

    write(outu,'(a)')&
    "================================================================"
    write(6,'(a)')""
    endif !if (verbose_debug)

   


  !----------------------------------------------------------------------------------
  !----- make separate lists of dihedrals with H and without H ----------------------
   call allint(nphi,iph,jph,kph,lph,ipa,jpa,kpa,lpa)
   call allint(nphi,icph,icpa)
   iph=0; jph=0; kph=0; lph=0; ipa=0; jpa=0; kpa=0; lpa=0; icph=0; icpa=0
   nphih=0
   mphia=0
   do n=1,nphi
      itl=locattype(ip(n))
      jtl=locattype(jp(n))
      ktl=locattype(kp(n))
      ltl=locattype(lp(n))
      dih_h=.false.
      if( ( attypeh(jtl) == 1 ) .or. (attypeh(ktl) == 1)) then
         write(outu,*)"dihe_types>","ERROR dihedral has H in middle ",n
         call err_exit("Cannot continue")
      endif
      dih_h= dih_h .or. (( attypeh(itl) == 1 ) &
                         .or. (attypeh(ltl) == 1)) 
      if(dih_h)then
         nphih=nphih+1
         iph(nphih)=ip(n)*3-3
         jph(nphih)=jp(n)*3-3
         kph(nphih)=kp(n)*3-3
         lph(nphih)=lp(n)*3-3
         icph(nphih)=icp(n)
      else
         mphia=mphia+1
         ipa(mphia)=ip(n)*3-3
         jpa(mphia)=jp(n)*3-3
         kpa(mphia)=kp(n)*3-3
         lpa(mphia)=lp(n)*3-3
         icpa(mphia)=icp(n)
      endif
   enddo

   if(verbose_debug)then
      write(outu,'(a,a,3i8)') &
        "<dihe_types>","nphi,nphih,mphia : ",nphi,nphih,mphia
   endif

   return
   
end subroutine dihe_types



!-----------------------------------------------------------------
!------ ANGLE_TYPES-----------------------------------------------------
!-----------------------------------------------------------------
!  go through atoms and collect atom types
subroutine angle_types
   implicit none
   logical, dimension(nattypes,nattypes,nattypes) :: langt
   integer, dimension(nattypes,nattypes,nattypes) :: &
     iangt !contains the local angle type of the local atom type combination
   integer :: n,m,i,j,k,ityp,jtyp,ktyp,itl,jtl,ktl,at=0, ioff
   integer :: center,ntheta_nowat,ictn
   logical :: ang_h
   integer,parameter :: zero=1.d0

   langt=.false.
   iangt=0
   if(verbose_debug)write(outu,*)"<angle_types>"
   !--- Count the angle types, avoid overcounting by keeping itl < ktl
   !--- flag which types are present in array langt
   !--- outer loop over center atom type
   do center=1,nattypes
      do n=1,ntheta
         j=jt(n)
         jtyp=iac(j)
         jtl=locattype(j)
         if(jtl == center)then
            i=it(n)
            k=kt(n)
            ityp=iac(i)
            ktyp=iac(k)
            itl=locattype(i)
            ktl=locattype(k)
            if(itl > ktl)call swap(i,k,itl,ktl)
            if(.not.langt(itl,ktl,center) )then
               at=at+1
            endif
            langt(itl,ktl,center)=.true.
         endif  !----- jtl == center 
      enddo
   enddo

   call allint(at,angtypes,ittyp,jttyp,kttyp)
   call allint(ntheta,ict)
   angtypes=0; ittyp=0; jttyp=0; kttyp=0; ict=0

   !--- iangt will contain the local angle type of the local atom type combination
   !--- ittyp,jttyp,kttyp contain the local types of atoms in angle
   !--- numang is the amber total number of angle types
   nangtypes=0
   do center=1,nattypes
      do itl=1,nattypes
         do ktl=itl,nattypes
            if(langt(itl,ktl,center))then
               nangtypes=nangtypes+1
               iangt(itl,ktl,center)= wat_check(nangtypes)
               ittyp(nangtypes)=itl
               jttyp(nangtypes)=center
               kttyp(nangtypes)=ktl
               if(verbose_debug) write(outu,'(2a,2i4,a,3i4)')" <angle_types>", &
                     "angle type ",nangtypes,iangt(itl,center,ktl), &
                     " angle ",itl,center,ktl
                     
            endif
         enddo
      enddo
   enddo
   numang=nangtypes

   !--- Assign angles to their local angle types based on their local atom types
   ntheta_nowat=0
   do n=1,ntheta
      jtl=locattype(jt(n))
      itl=locattype(it(n))
      ktl=locattype(kt(n))
      if(itl > ktl)call swap(i,k,itl,ktl)

      ictn=iangt(itl,ktl,jtl)

      if(ictn == 0 ) then
         write(outu,'("error in finding bond type",10i4)') &
                          n,itl,jtl,ktl,ioff,iangt(itl,ktl,center)
         call err_exit("Cannot continue")
      elseif(ictn == -1) then
         !this is a WAT angle, get rid of it
      else
         ntheta_nowat=ntheta_nowat+1
         ict(ntheta_nowat)=ictn
         jt(ntheta_nowat) = jt(n)
         it(ntheta_nowat) = it(n)
         kt(ntheta_nowat) = kt(n)
      endif
   enddo
!
!  !----- list angles with and without H ----------------------
   call allint(ntheta,ith,jth,kth,ita,jta,kta)
   call allint(ntheta,icth,icta)
   ith=0; jth=0; kth=0 ; ita=0; jta=0; kta=0
   icth=zero; icta=zero
   ntheth=0
   mtheta=0
   do n=1,ntheta
      itl=locattype(it(n))
      jtl=locattype(jt(n))
      ktl=locattype(kt(n))
      ang_h=.false.
      !--- Check if any of the atoms in bond are H
      do m=1,nattypes
         ang_h= ang_h .or. (( attypeh(itl) == 1 ) &
                         .or. (attypeh(jtl) == 1) &
                         .or. (attypeh(ktl) == 1)) 
      enddo

      if(ang_h)then
         ntheth=ntheth+1
         ith(ntheth)=it(n)*3-3
         jth(ntheth)=jt(n)*3-3
         kth(ntheth)=kt(n)*3-3
         icth(ntheth)=ict(n)
      else
         mtheta=mtheta+1
         ita(mtheta)=it(n)*3-3
         jta(mtheta)=jt(n)*3-3
         kta(mtheta)=kt(n)*3-3
         icta(mtheta)=ict(n)
      endif
   enddo
   if(verbose_debug)write(outu,*)"<angle_types>","ntheta,ntheth,mtheta : ",ntheta,ntheth,mtheta
   return

   contains
     !----------------WAT_CHECK -------------------------------------
     integer function wat_check(n) result(t)
       integer::n
       logical::iswat
       t = n
       if(tip3_flex)return

       if(verbose_debug) &
            write(outu,'("checking for water",3i6)')n,iac(it(n)),iac(jt(n))
       iswat = iac(it(n)) == tip3_htype
       iswat = iswat .and. (iac(kt(n)) == tip3_htype)
       iswat = iswat .and. (iac(jt(n)) == tip3_otype)
       if(verbose_debug) write(outu,'("     Result: ",l3)')iswat
       if(iswat)t = -1
       return
     end function wat_check
       
       
end subroutine angle_types

!-----------------------------------------------------------------
!------ Bond_TYPES-----------------------------------------------------
!-----------------------------------------------------------------
!  go through bonds and collect bond types
subroutine bond_types
   implicit none
   logical, dimension(nattypes,nattypes) :: bndt
   integer, dimension(nattypes,nattypes) :: ibndt
   integer :: n,m,i,j,it,jt,itl,jtl,bt=0, ioff
   logical :: bond2h

   !---------------------------------------------------
   ! Counts bond types : nbndtypes
   ! Fills   ibndtyp(it,jt)  = bond typd for types it,jt
   !         ibtyp(n)        = atom type for atom i in bond type n
   !         jbtyp(n)        = atom type for atom j in bond type n
   !         
   bndt=.false.
   ibndt=0
   do n=1,nbond
      i=ib(n)
      j=jb(n)
      it=iac(i)
      jt=iac(j)
      itl=locattype(i)
      jtl=locattype(j)
      if(itl > jtl)call swap(i,j,itl,jtl)
      if(.not.bndt(itl,jtl) )then
         bt=bt+1
      endif
      bndt(itl,jtl)=.true.
   enddo

   nbndtypes=0

   call allint(bt,bndtypes,ibtyp,jbtyp)
   call allint(nbond,icb)
   bndtypes=0; ibtyp=0; jbtyp=0; icb=0

   do it=1,nattypes
      do jt=it,nattypes
         if(bndt(it,jt))then
            nbndtypes=nbndtypes+1 
            ibtyp(nbndtypes)=it
            jbtyp(nbndtypes)=jt
            ibndt(it,jt)=nbndtypes
         endif
      enddo
   enddo
   numbnd=nbndtypes

   do n=1,nbond
      i=ib(n)
      j=jb(n)
      it=locattype(i)
      jt=locattype(j)
      if(it > jt)call swap(i,j,it,jt)
      icb(n)=ibndt(it,jt)
      if(icb(n) == 0 ) then
         write(outu,*)"<bond_types>","error in finding bond type",n,it,jt,ioff
         call err_exit("Cannot continue")
      endif
   enddo

  !----- list bonds to H ----------------------
   call allint(nbond,ibh,jbh,iba,jba,icbh,icba)
   ibh=0; jbh=0; iba=0; jba=0; icbh=0; icba=0

   nbonh=0
   mbona=0
   do n=1,nbond
      i=ib(n)
      j=jb(n)
      itl=locattype(i)
      jtl=locattype(j)
      bond2h=.false.
      do m=1,nattypes
         bond2h= bond2h .or. (( attypeh(itl) == 1 ) .or. (attypeh(jtl) == 1)) 
      enddo
      if(bond2h)then
         nbonh=nbonh+1
         ibh(nbonh)=(i-1)*3
         jbh(nbonh)=(j-1)*3
         icbh(nbonh)=icb(n)
         if(verbose_debug)write(outu,*)"<bond_types>","with H ",nbonh,i,j,icb(n)
      else
         mbona=mbona+1
         iba(mbona)=i*3-3
         jba(mbona)=j*3-3
         icba(mbona)=icb(n)
         if(verbose_debug)write(outu,*)"<bond_types>","without H (A) ",mbona,i,j,icb(n)
      endif
   enddo
   deallocate(icb,bndtypes)
   if(verbose_debug)write(outu,*)"<bond_types>","nbond,nbonh,mbona : ",nbond,nbonh,mbona
   return

end subroutine bond_types

subroutine swap(i,j,it,jt)
   implicit none
   integer :: i,j,it,jt
   integer :: ii,iit

   ii=j
   iit=jt
   j=i
   jt=it
   i=ii
   it=iit
   return
end subroutine swap


!-----------------------------------------------------------------
!------ ATOM_TYPES-----------------------------------------------------
!-----------------------------------------------------------------
!  go through atoms and collect atom types
!  as they were found in the psf file
!
subroutine atom_types
  implicit none
  logical, dimension(maxattypes) :: att
  integer :: i,j,nt=0
  logical :: counted

  att_txt(1:maxattypes)(1:8) = '        '
  if(verbose_debug)write(outu,'(a)')"<atom_types>"

  !-------- XPLOR format----------------------------
  if(l_xplor_psf) then
     nt=0
     !-- count unique atom type character strings
     do i=1,natom
        counted=.false.
        do j=1,nt
           counted = counted .or. (atom_type(i)(1:idlen) == att_txt(j)(1:idlen))
        enddo
        if(.not. counted) then
           nt=nt+1
           att_txt(nt)(1:idlen) = atom_type(i)(1:idlen)
        endif
     enddo
  !------- CHARMM format-------------------------------------
  else
     !--- count the atom types in this system  nt ----------------
     att=.false.
     do i=1,natom
        if(.not. att(iac(i)))nt=nt+1
        att(iac(i))=.true.
     enddo
  endif

  call allint(nt,attypes,attypeh)
  attypes=0
  attypeh=0
  call allint(natom,locattype)
  locattype=0

  !------- For XPLOR psf types, we do not have the type number yet so we have to do the
  !        rest of the work after reading the topology file to match up the atom type strings
  !        with the atom type numbers. The rest of the work in this routine will happen in 
  !        get_atom_parameters
  if(l_xplor_psf) then
     nattypes = nt
     ntypes=nt
     return
  endif


  !--- fill attypes array with type number for parameter file -----------
  !--- locattype array will have the "local atom type" for each atom
  !--- attypes will point the local atom type back to the charmm atom type
  !--- nattypes has the number of atom types in this system
  !--- attypeh array flags whether it is a heavy atom type (2) or H type (1)
  nattypes=0
  do i=1,maxattypes
     if(att(i))then
        nattypes=nattypes+1
        attypes(nattypes)=i
     endif
  enddo
  if(nattypes /= nt) print *,"Trouble in river city, nt /= nattypes"

  do i=1,natom
     do j=1,nattypes
        if(iac(i) == attypes(j) ) then
           locattype(i)=j
           if(amass(i) < 1.5 ) then
              if(attypeh(j).eq.2)then
                 write(outu,*)"<atom_types>", &
                      "ERROR atom with mass 1 has type already", &
                      "designated as non-hydrogen type"
                 write(outu,*)"<atom_types>"," atom ",i,amass(i),atom_label(i)
                 call err_exit("Cannot continue")
              endif
              attypeh(j)=1
           else
              if(attypeh(j).eq.1)then
                 write(outu,*)"<atom_types>", &
                      "ERROR atom with mass >1 has type already", &
                      "designated as non-hydrogen type"
                 write(outu,*)"<atom_types>"," atom ",i,amass(i),atom_label(i)
                 call err_exit("Cannot continue")
              endif
              attypeh(j)=2
           endif
           if(verbose_debug) &
                write(outu,'(a,i10,a,i3,a,i3)') &
                "<atom_types> atom ",i, &
                "  chmtype ",attypes(j), &
                "    local type ",locattype(i)
        endif
     enddo
  enddo
  !--- ntypes is the amber number of atom types --------------
  ntypes=nattypes
  return
end subroutine atom_types

!-----------------------------------------------------------------
!----------- READ_BOND_ANGLE_DIHE--------------------------------------
!-----------------------------------------------------------------
subroutine read_psf_bond_angle_dihe()
   implicit none
   integer :: i,u

   u=psf_unit
   if(verbose_debug)write(outu,'(a)')"<read_bond_angle_dihe> bonds"
          read(u,fmt00) nbond
          if(nbond > 0 )then
             call allint(nbond,ib,jb)
             read(u,fmt03) (ib(i),jb(i),i=1,nbond)
          else
             read(u,fmt03)
          endif

   if(verbose_debug)write(outu,'(a)')"<read_bond_angle_dihe> angles"
          read(u,fmt00) ntheta
          if(ntheta > 0 )then
             call allint(ntheta,it,jt,kt)
             read(u,fmt04) (it(i),jt(i),kt(i),i=1,ntheta)
             !Remember, there is no forth index term here;
             !that is assigned later.
          else
             read(u,fmt04)
          endif

   if(verbose_debug)write(outu,'(a)')"<read_bond_angle_dihe> dihedrals"
          read(u,fmt00) nphi
          if(nphi > 0 )then
             call allint(nphi,ip,jp,kp,lp)
             read(u,fmt03) (ip(i),jp(i),kp(i),lp(i),i=1,nphi)
             if(verbose_debug) then
              write(*,*) "I should find ",nphi," dihedrals :"
              write(*,fmt03) (ip(i),jp(i),kp(i),lp(i),i=1,nphi)
             endif
          else
             read(u,fmt03)
          endif

   if(verbose_debug)write(outu,'(a)')"<read_bond_angle_dihe> impropers"
          read(u,fmt00) nimphi
          if(nimphi > 0 )then
             call allint(nimphi,im,jm,km,lm)
             read(u,fmt03)(im(i),jm(i),km(i),lm(i),i=1,nimphi)
          else
             read(u,fmt03)
          endif

   if(verbose_debug)write(outu,'(a)')"<read_bond_angle_dihe> donors"
          read(u,fmt00) ndon
          if(ndon > 0 )then
             call allint(ndon,idon,ihd1)
             read(u,fmt03) (idon(i),ihd1(i),i=1,ndon)
          else
             read(u,fmt03)
          endif

   if(verbose_debug)write(outu,'(a)')"<read_bond_angle_dihe> acceptors"
          read(u,fmt00) nacc
          if(nacc > 0 )then
             call allint(nacc,iacc,iac1)
             read(u,fmt03) (iacc(i),iac1(i),i=1,nacc)
          else
             read(u,fmt03)
          endif

   if(verbose_debug)write(outu,'(a)')"<read_bond_angle_dihe> nonbond"
          read(u,fmt00) nnb
          if(nnb.gt.0)then
             call allint(nnb,inb)
             read(u,fmt03) (inb(i),i=1,nnb)
          else
             read(u,fmt03)
          endif
          call allint(natom,iblo)
          read(u,fmt03) (iblo(i),i=1,natom)

   if(verbose_debug)write(outu,'(a)')"<read_bond_angle_dihe> group"
          read(u,fmt05) ngrp,nst2
          call allint(ngrp,igpbs,igptyp,imoveg)
          read(u,fmt04) (igpbs(i),igptyp(i),imoveg(i),i=1,ngrp)
          
   if(verbose_debug)write(outu,'(a)')"<read_bond_angle_dihe> Done"
   return
end subroutine read_psf_bond_angle_dihe

subroutine skip_over_psf_molnt()
! Skips over the CHEQ related MOLNT and IMOLNT information
! This is routine is kind of pointless, but needed since 
! there is not an effective mechanism to navigate the psf
! internally.

   implicit none
   integer                        :: dummy_molnt
   integer, pointer, dimension(:) :: dummy_imolnt
   integer                        :: u,i

   u=psf_unit
!   !Snag the following:
!   !       1 !MOLNT
   read(u,fmt00) dummy_molnt

!   !imolnt should extend as far as natom

   call allint(natom, dummy_imolnt)
   read(u,fmt03) (dummy_imolnt(i), i=1, natom)
!     216 !MOLNT
!       1       1       1       2       2       2       3       3
!       3       4       4       4       5       5       5       6
  deallocate( dummy_imolnt )

end subroutine skip_over_psf_molnt

subroutine read_psf_lp
!Reads the number of lone pairs (NUMLP)
!and the number of lone pairs in the host table (NUMLPH)

   implicit none
   integer :: u

   u=psf_unit
   !Snag the following
   !       0       0 !NUMLP NUMLPH
   read(u,fmt05) numlp,numlph

   return

end subroutine read_psf_lp

subroutine read_psf_cmap
!Reads the number of cmap terms and the atom numbers of the associated 
!dihedral pairs from the psf. 

!       ____________ Total number of cross terms (i.e. cmap_term_count)
!      |
!      |
!     157 !NCRTERM: cross-terms
!      18      20      22      37      20      22      37      39  <- 1st term
!      37      39      41      48      39      41      48      50  <- 2nd term

!      |-------------------------|    |--------------------------|
!            Dihedral 1                       Dihedral 2


!It will fill an array cmap_index() with the atom index of the two 
!dihedrals making up this cmap term. The cmap_index() will have
!cmap_term_count members.

!The 9th term in this array will eventually be a pointer into a cmap type array. 
!A value of 0 indicated an unassigned cmap type for that cmap term. During this
!initial read in phase, this value is set to 0 for all terms; only later in the
!subroutine cmap_types is this 9th value actually set.
   use cmap
   implicit none
   integer :: i,ierr

   !     157 !NCRTERM: cross-terms
   read(psf_unit,fmt00) cmap_term_count
   
   if (verbose_debug) write(outu,'(i8,a)') cmap_term_count," CMAP cross terms found"

   !Allocate space for the CMAP indexes bases on already read in cmap_term_count
   allocate(cmap_index(9,cmap_term_count),stat=ierr)
   if ( ierr /= 0 ) then
     write(outu,'(a,i8,a)') "FATAL ERROR: Allocation of cmap_index(9,",cmap_term_count,") Failed."
     call mexit(outu,1)
   end if

   ! Zero out the cmap_index array - we do this so we can check the array for
   ! zeros as a sanity check to find any unassigned terms.
   cmap_index(1:9,1:cmap_term_count) = 0

   ! Read the cmap terms from the psf file.
   ! The procedure is:
   !  1) Read the 8 x cmap_term_count atom numbers.
   ! fmt03='(8I10)'
   !      18      20      22      37      20      22      37      39 
   !      37      39      41      48      39      41      48      50

   read(psf_unit,fmt03) (cmap_index(1:8,i),i=1,cmap_term_count)

   !Debug
   if (verbose_debug) write(outu,fmt03) (cmap_index(1:8,i),i=1,cmap_term_count)


   return

end subroutine read_psf_cmap

subroutine assign_cmap_types
!Initially, this works out how how many unique cmap types there are
!in the cmap_index() array.

!Once it has these unique CMAP types it will assigns numbers to them.

!Next it will look in the CHARMM parameter file and find the associated
!parameter.

   use psf_strings, only: getwords
   use cmap
   implicit none

   integer :: current_cmap_term_type(8)
   logical :: new_cmap_type, match_found, lookingforsection, stillreading

   character(len=80),dimension(10) :: words !needed by getwords
   character(len=idLen),dimension(8) :: current_cmap_param !temporary

   integer :: ierr,i,j,k,numwrd
   integer :: k_n_plus_1, remaining_k_left_to_read
   integer :: statall

   integer, pointer, dimension(:,:) :: cmap_type_index !Contains the atom types for the atom
                                                       !atoms that make up each unique cmap term.
                                                       !The 6th entry of cmap_index points into
                                                       !this array.


   ! cmap_type_index - we do not know ahead of time how many unique types there will be. But we
   !                   know it will be <= cmap_term_count so for simplicity we will just allocate
   !                   it to be cmap_term_count long here.
   !                   MJW- currently in par_all22_prot.inp there are ONLY 6

   allocate(cmap_type_index(8,cmap_term_count),stat=ierr)
   if ( ierr /= 0 ) then
     write(outu,'(a,i8,a)') "FATAL ERROR: Allocation of cmap_type_index(8,",cmap_term_count,") Failed."
     call mexit(outu,1)
   end if


   ! The procedure is:
   !     Translate these to atom types using locattype() and work out how many unique
   !     cmap types there are and assign them type index numbers that
   !     are then added to column 9 of the cmap_index.


   !Next we need to loop through each cmap term convert the atom numbers to types
   !and then check if the type is a duplicate and if it isn't add it to a type array.
   cmap_type_count = 0
 
   do i = 1, cmap_term_count
     !Step 1 - obtain the atom types for this cmap term.
     current_cmap_term_type(1:8) = locattype(cmap_index(1:8,i))
     
     !Step 2 - now we have the types - in integer representation we
     !         loop through our current list of types looking for
     !         duplicates.
     new_cmap_type=.true.
     do j = 1,cmap_type_count

       !Check if all 8 terms in this dihedral pair match the current type.
       match_found = .true.
       do k = 1,8
         if (current_cmap_term_type(k) /= cmap_type_index(k,j)) then
           match_found = .false.
           exit
         end if
       end do

       if (match_found) then
         !We found a matching cmap type already defined for this cmap term.
         new_cmap_type = .false.
         cmap_index(9,i) = j
         !exit the search loop.
         exit
       end if

     end do !j = 1,cmap_type_count

     if (new_cmap_type) then
       !We did not find a duplicate already so this is a new type.
       cmap_type_count = cmap_type_count + 1
       !Store the list of 8 atom type numbers making up this type.
       cmap_type_index(1:8,cmap_type_count) = current_cmap_term_type(1:8)

       !Store this type id as the type to be used for the current cmap
       !term we are looking at.
       cmap_index(9,i) = cmap_type_count
     end if

   end do !i = 1, cmap_term_count


   !Once we get here we should have an array (cmap_type_index) that is cmap_type_count long and
   !contains the atom types of each unique cmap term. Additionally we should have the 9th element
   !of each cmap_index list set to the type that it maps to.


   !Sanity check - make sure no value in cmap_index is zero. This would imply it was not set.
   !Not really needed since it should never be able to have a situation where something in cmap_index
   !was not set but we will test anyway for debugging purposed.
   do i = 1, cmap_term_count
     do j = 1,9
       if ( cmap_index(j,i) == 0 ) then
          write(outu,'(a,i6,a,i6,a)') "FATAL ERROR: Element: ',j,',',i,' of cmap_index is 0."
          write(outu,'(a)') "             This should never occur - logic failure in"
          write(outu,'(a)') "             CMAP identification code."
          call mexit(outu,1)
       end if
     end do
   end do !i = 1, cmap_term_count

   !If verbose debugging print some info about the cmap terms.
   if (verbose_debug) then
     write(outu,'(a,i6,a)') "Found ", cmap_type_count," unique cmap types:"
     do i = 1,cmap_type_count
       write(outu,'(a,i4,a,8(I4,1X))') "Type: ",i," - ",cmap_type_index(1:8,i)
       write(outu,'(13X,8(A4,1X))') attype_name(cmap_type_index(1:8,i))
       write(outu,'(a)') ""
     end do
   end if

   !We now know how large to make the needed_cmap_types array
   !from the unique number of cmap types found ( cmap_type_count )
   !hence allocate it now
   allocate(needed_cmap_types(cmap_type_count),stat=statall)
   if ( statall /= 0 ) then
     call allocate_error("needed_cmap_types",&
                   "cmapParam",cmap_type_count)
   endif


   !Partially populate each cmapParameter type within the known_cmap_types array
   !with the known label information. Later, this will be used to populate the rest
   !of the cmapParameter type
   do i = 1,cmap_type_count
     needed_cmap_types(i)%label(1:8) = &
                                   attype_name(cmap_type_index(1:8,i))(1:4)

     needed_cmap_types(i)%charmm_atom_type_numerical_lbl(1:8) = &
                                                     cmap_type_index(1:8,i)

      !write(outu,'(8(A4X))'),needed_cmap_types(i)%label(1:8)
   enddo

   !We have finished with cmap_type_index; deallocate
   if(associated(cmap_type_index)) then
     deallocate(cmap_type_index,stat=ierr)
     if (ierr /= 0) then
       write(outu,'(a)') "ERROR Deallocating cmap_type_index."
       call mexit(outu,1)
     end if
   end if

   !Now start the navigation of the CHARMM parameter file
   rewind(prm_unit)

   !Get to the CMAP section
   lookingforsection=.true. !Still navigating to CMAP section
   do while(lookingforsection)
     if(getwords(words,1,numwrd,prm_unit) == 2 ) &  !This 2nd parameter *must* 
                                                    !not be greater than 35
      call err_exit("ERROR end of file for reading prm_unit, No CMAP")
      if(words(1)(1:4) == "CMAP") then
        lookingforsection=.false.
      endif
   enddo 
    
   stillreading=.true. !We're currently in the CMAP section

   do while(stillreading)
     if(getwords(words,9,numwrd,prm_unit) == 2 ) &
           call err_exit("ERROR end of file for reading CMAP")

     !Attempt to find a CMAP parameter by looking for
     !a line in the CMAP section with nine words:
     !C    NH1  CT1  C    NH1  CT1  C    NH1   24
     if(numwrd == 9)then
       read(words(1:8),'(a4,1x)') current_cmap_param(1:8)

       do i = 1,cmap_type_count

         !Check current_cmap_param against needed_cmap_types array
         !This next expression forms a logical array and then uses the all
         !intrinsic to to see if they are all True
         if( all (needed_cmap_types(i)%label ==  current_cmap_param) ) then
           if(verbose_debug) then
             write(outu,'(a,i4)') "Match found in param file for index: ",i
             write(outu,'(8(a4,1x))') current_cmap_param(1:8)
           endif !if(verbose_debug)

           !Now read in the resolution associated with this parameter:
           !C    NH1  CT1  C    NH1  CT1  C    NH1   24
           !                                         ^
           read(words(9),'(i4)') needed_cmap_types(i)%resolution
           !and set the step size 
           needed_cmap_types(i)%gridStepSize = &
                             360/needed_cmap_types(i)%resolution

           !TODO check if the resolution is sane
           ! 360/grid_resolution should *not* have a remainder

           ! From the needed_cmap_types(i)%resolution we now know how big the
           ! temporary array for storing the cmap data points should be
           ! hence allocate it now

           allocate(needed_cmap_types(i)%grid(&
                                    needed_cmap_types(i)%resolution, &
                                    needed_cmap_types(i)%resolution),&
                    stat=statall)

           if ( statall /= 0 ) then
             call allocate_error("needed_cmap_types(i)%grid",&
                              "allreal",needed_cmap_types(i)%resolution)
           endif

           ! Step through block of grid points corresponding
           ! to 360/needed_cmap_types(i)%resolution
           do j=1, needed_cmap_types(i)%resolution

             !The point of the next do-while loop is to step over the following
             !marked lines in the charmm parameter file:

             ! C    NH1  CT1  C    NH1  CT1  C    NH1   24
             !                            <-----------
             ! !-180                      <-----------
             ! 0.126790 0.768700 0.971260 1.250970 2.121010
             do while( len_trim(words(1)) < 5 ) !Blank line, hence skip
               if(getwords(words,5,numwrd,prm_unit) == 2 ) &
               call err_exit("ERROR end of file for reading CMAP")
             enddo


             !Now we start to parse the following section:
             !
             !0.126790 0.768700 0.971260 1.250970 2.121010
             !2.695430 2.064440 1.764790 0.755870 -0.713470
             !0.976130 -2.475520 -5.455650 -5.096450 -5.305850
             !-3.975630 -3.088580 -2.784200 -2.677120 -2.646060
             !-2.335350 -2.010440 -1.608040 -0.482250
             !
             ! into the 
             !needed_cmap_types(i)%grid(i->needed_cmap_types(i)%resolution,j)

             ! We assume there are 5 data points per line. This is set by the
             ! cmap_entries_per_line parameter
             do k=1, needed_cmap_types(i)%resolution, cmap_entries_per_line
               ! Generally, needed_cmap_types(i)%resolution = 24
               ! hence k will run, 1,6,11,16,21 
               k_n_plus_1 = k+cmap_entries_per_line
               if (k_n_plus_1 <= needed_cmap_types(i)%resolution) then
                 read(words(1:cmap_entries_per_line),'(f9.6)') &
                  needed_cmap_types(i)%grid(k:k+(cmap_entries_per_line - 1),j)
               else
                 remaining_k_left_to_read = needed_cmap_types(i)%resolution - k
                 read(words(1:(remaining_k_left_to_read+1) ),'(f9.6)') &
                     needed_cmap_types(i)%grid(k:(k+remaining_k_left_to_read),j)
               endif

               if(getwords(words,5,numwrd,prm_unit) == 2 ) &
               call err_exit("ERROR end of file for reading CMAP")

             enddo

           enddo ! j=1, needed_cmap_types(i)%resolution

           !The next loop over will get the 2nd block of CMAP data
           !associated with the current_cmap_param :
           !
           !!-165
           !-0.802290 1.377090 1.577020 1.872290 2.398990
           !2.461630 2.333840 1.904070 1.061460 0.518400
           !-0.116320 -3.575440 -5.284480 -5.160310 -4.196010
           !etc


           !DEBUG; show what we have from this loop
           if (verbose_debug) then
             write(outu,'(a,i4,a)') "Harvesting for cmap_type  ",i,&
                               " has found the following:"
             write(outu,'(8(a4,1x),i4)') needed_cmap_types(i)%label(1:8),&
                                 needed_cmap_types(i)%resolution
             write(outu,'(5(f9.6))') &
                needed_cmap_types(i)%grid(&
                                  1:needed_cmap_types(i)%resolution,&
                                  1:needed_cmap_types(i)%resolution)
           endif !(verbose_debug) then

         endif !if( all (needed_cmap_types(i)%label ==  current_cmap_param) ) then

         !do nothing here 
      
       enddo !i = 1,cmap_type_count
     endif !if(numwrd == 9)then

     !Break out clause
     if(words(1)(1:4) == "NONB") then
       stillreading=.false. !We've left the CMAP section
     endif

   enddo !do while(stillreading)





   !Debug
   if (verbose_debug) then
     call print_cmap_summary(outu)
   endif !(verbose_debug) then


end subroutine assign_cmap_types



subroutine print_psf_summary
!Give a summary to the end user of what was parsed in the PSF
   use cmap, only : cmap_term_count, cmap_type_count
   implicit none
   character(len=80) :: print_fmt = '(A50,1X,I8)'

   write(outu,'(a)') ""  
   write(outu,'(a)') &
   "==========================================================="
   write(outu,'(a)') "          PSF input parsing summary"
   write(outu,'(a)') &
   "==========================================================="
   write(outu,'(a)') ""
   write(outu,print_fmt)  "Number of PSF flags found:", size(psf_flags)
   write(outu,'(a)') ""
   write(outu,print_fmt)  "Number of atoms found:", natom
   write(outu,print_fmt)  "Number of residues found:", nres
   write(outu,'(a)') ""
   write(outu,print_fmt)  "Number of bonds found:", nbond
   write(outu,print_fmt)  "Number of angles found:", ntheta
   write(outu,print_fmt)  "Number of dihedrals found:", nphi
   write(outu,print_fmt)  "Number of impropers found:", nimphi
   write(outu,'(a)') ""
   write(outu,print_fmt)  "Number of donors found:", ndon
   write(outu,print_fmt)  "Number of acceptors found:", nacc
   write(outu,print_fmt)  "Number of explicit nonbonded exclusions found:", nnb 
   write(outu,'(a)') ""
   write(outu,print_fmt)  "Number of groups found:", ngrp
   write(outu,print_fmt)  "Number of ST2 waters found:", nst2
   write(outu,'(a)') ""
   if( has_psf_flag("CMAP") .and. CMAP_enabled )then
     write(outu,print_fmt)  "Number of cross terms found:", cmap_term_count
   endif

   write(outu,'(a)') &
   "==========================================================="
   write(outu,'(a)') ""
   
   return

end subroutine print_psf_summary

subroutine print_assignment_summary
!Give a summary to the end user of what assignments were made using
!CHARMM'S top and parm files
   use cmap, only : cmap_type_count
   implicit none
   character(len=80) :: print_fmt = '(A50,1X,I8)'

   write(outu,'(a)') ""
   write(outu,'(a)') &
   "==========================================================="
   write(outu,'(a)') "          Parameter assignment summary             "
   write(outu,'(a)') &
   "==========================================================="

   write(outu,'(a)') ""
   write(outu,'(A18,1X,I4)')  "CHARMM ff version:", chmff_verno
   write(outu,'(A18,1X,A40)')  "CHARMM ff type:", chmff_type
   write(outu,'(a)') ""
   write(outu,print_fmt)  "Number of atom types assigned:", ntypes
   write(outu,'(a)') ""
   write(outu,print_fmt)  "Number of bond types assigned:", nbndtypes
   write(outu,print_fmt)  "Number of angle types assigned:", nangtypes
   write(outu,print_fmt)  "Number of UB angle terms found:", nub
   write(outu,print_fmt)  "Number of UB angle types assigned:", nubtypes
   write(outu,print_fmt)  "Number of dihedrals types assigned:", ndihtypes
   write(outu,'(A40,3X,I8)')  "Explicit:", Num_Expl_Dihed_Types_Assigned
   write(outu,'(A40,3X,I8)')  "WildCard:", ndihtypes - Num_Expl_Dihed_Types_Assigned
   write(outu,print_fmt)  "Number of improper types assigned:", nimprtypes
   write(outu,'(a)') ""
   if( has_psf_flag("CMAP") .and. CMAP_enabled )then
     write(outu,print_fmt)  "Number of unique cross types found:", cmap_type_count
     !write(outu,print_fmt), "Number of unique cross types assigned:", cmap_type_count
   endif

   write(outu,'(a)') &
   "==========================================================="
   write(outu,'(a)') ""

   return

end subroutine print_assignment_summary

!-----------------------------------------------------------------
!--------- READ_ATOM_DATA ----------------------------------
!    Read first part of psf after the title:
!           natom 
!           atom number, resname, resnum, restype
!           atom name, atom_type, charge, mass, CHEQ stuff
!    Allocate and fill arrays, convert charges to internal Amber chgs
!-----------------------------------------------------------------
subroutine read_psf_atom_data
  implicit none
  integer :: statall
  integer ii,i,this_resid,this_resnum
  character(len=8) :: lsegid
  integer :: n_res_at

  read(psf_unit,fmt00)natom

  allocate(atom_label(natom),lresat(natom),stat=statall)
  if ( statall /= 0 ) call allocate_error("atom_label","read_atom_data",natom)

  allocate(atom_type(natom),stat=statall)
  if ( statall /= 0 ) call allocate_error("atom_type","read_atom_data",natom)

  allocate(iac(natom),imove(natom),stat=statall)
  if ( statall /= 0 ) call allocate_error("iac, imove ","read_atom_data", &
       natom)

  allocate(cg(natom),amass(natom),stat=statall)
  if ( statall /= 0 ) call allocate_error("cg, amass ","read_atom_data",natom)

  call allint(natom,iresid)
  call allreal(natom,radii,screen)   

  if(verbose_debug) write(outu,'(a,i8)')"<read_atom_data> natom",natom

  nres=0
  this_resid=0
  n_res_at=1
  max_res_at = 1

  do i=1,natom
     if( l_xplor_psf)then
        !        1 GLUC 1    BGLC C1   CBS    0.200000       12.0110           0   0.00000     -0.301140E-02
        read(psf_unit,fmt01) ii,           & ! 1
             lsegid,        & ! GLUC
             iresid(i),     & ! 1
             lresat(i),     & ! BGLC
             atom_label(i), & ! C1
             atom_type(i),  & ! CBS
             cg(i),         & ! 0.200000
             amass(i),      & ! 12.0110
             imove(i)        ! 0
     else
        !       1 AAL  1    ASN  N      54  -0.470000       14.0070           0
        read(psf_unit,fmt01) ii,           & ! 1
             lsegid,       & ! AAL
             iresid(i),    & ! 1
             lresat(i),    & ! ASN
             atom_label(i), & ! N
             iac(i),       & ! 54
             cg(i),        & ! -0.47000
             amass(i),     & ! 14.0070
             imove(i)        ! 0

        if(atom_label(i) == "OH2") atom_label(i) = "O" !TIP3+TP3M
        !Rename atomtype to ensure SANDER's fastwat routine picks these up
        !as waters:
        !
        ! Waters are defined as residues that meet the following criteria:
        !   1) residue name = watnam
        !   2) atom names are owatnm, hwatnm(1) and hwatnm(2)
     endif

     if(iresid(i) /= this_resid)then
        nres=nres+1
        this_resid=iresid(i)
        n_res_at=1
     else
        n_res_at=n_res_at+1
     endif
     max_res_at = max(max_res_at,n_res_at)
  enddo

  allocate(lres(nres),stat=statall)
  if ( statall /= 0 ) call allocate_error("lres ","read_atom_data",natom)

  this_resid=0
  ii=0
  do i=1,natom
     if(iresid(i) /= this_resid)then
        ii=ii+1
        this_resid=iresid(i)
        lres(ii)=lresat(i)
        if (lresat(i) == "TIP3") lres(ii) = "WAT" !TIP3_TP3M mod: See above
     endif
  enddo


  !--- convert to internal amber charges (whole cg array)-------
  cg=cg * sqrt(332.0716D0) !CCELEC

  !--- Fill ipres -------
  this_resid=0
  this_resnum=0
  call allint(nres,ipres)
  ipres=0
  do i=1,natom
     if(this_resid /= iresid(i))then
        this_resid=iresid(i)
        this_resnum=this_resnum+1
        ipres(this_resnum)=i
        if(verbose_debug) &
             write(outu,'(a,i6,i6)')"<read_atom_data> ipres ",this_resnum,i
     endif
  enddo

  return
end subroutine read_psf_atom_data


!-----------------------------------------------------------------
!--------- OPEN_UNITS ---------------------------------------
!-----------------------------------------------------------------
subroutine open_units()
   implicit none
   logical :: file_exists
   
   inquire(file=psf_filename,exist=file_exists)
   if(file_exists)then
      open(unit=psf_unit,file=psf_filename)
   else
      write(outu,*)"ERROR in opening psf file named ",psf_filename(1:len_trim(psf_filename))
      call usage_message(1)
   endif

   inquire(file=param_filename,exist=file_exists)
    if( file_exists)then
      open(unit=prm_unit,file=param_filename)
   else
      write(outu,*)"ERROR in opening param file named: ",param_filename
      call usage_message(1)
   endif

   inquire(file=topology_filename,exist=file_exists)
    if( file_exists)then
      open(unit=top_unit,file=topology_filename)
   else
      write(outu,*)"ERROR in opening topology file named: ",topology_filename
      call usage_message(1)
   endif

   inquire(file=chmcrd_filename,exist=file_exists)
   if( file_exists)then
      open(unit=chmcrd_unit,file=chmcrd_filename)
   else
      write(outu,*)"ERROR in opening file specified by -crd named: ",chmcrd_filename
      call usage_message(1)
   endif

   if ( want_vmd_prmtop ) then
      open(unit=vmd_prmtop_unit,file=vmd_prmtop_filename, status="replace")
   endif

   open(unit=prmtop_unit,file=prmtop_filename, status="replace")

   open(unit=inpcrd_unit,file=inpcrd_filename, status="REPLACE")

   if(index(outu_name,"stdout") .ne. 1) then
      inquire(file=outu_name,exist=file_exists)
      if(.not. file_exists)then
         open(outu,file=outu_name, status="REPLACE")
      else
         write(outu,*)"ERROR in opening ",outu_name
         call usage_message(1)
      endif
   endif
   return
end subroutine open_units

!-----------------------------------------------------------------
!--------- READ_PSF_FLAGS ---------------------------------------
!-----------------------------------------------------------------
! Reads the first line of a PSF:
! 
! PSF CMAP CHEQ
!
! and parses these flags into the psf_flags array

subroutine read_psf_flags
   use psf_strings, only: getwords
   implicit none
   integer :: numwrd,ierr,i
   
   ierr = getwords(words,maxwrd,numwrd,psf_unit)
   if(ierr == 2 ) then
      write(outu,*) "ERROR end of file for reading psf_unit, retval=",ierr
         call err_exit("Cannot continue")
      stop
   endif

   if(verbose_debug) then
    write(outu,'(i4,a)') numwrd," FLAGs found in the PSF file header"
    do i=1,numwrd
     write(outu,'(a)') words(i)
    enddo
   endif
 
   !We now know how many flags we have, hence we can allocate the space
   allocate(psf_flags(numwrd))

   do i=1,numwrd
    !Add to psfprm flag array
    psf_flags(i) = words(i)(1:8)
   enddo

   !The first flag MUST be PSF
   if(index(psf_flags(1),"PSF").ne.1)then
      write(outu,'(a,a)')"PSF file is not a psf, label is: ",psf_flags(1)
      write(outu,'(a)') "Check the psf file, exiting."
      call mexit(outu,1)
   else
      if(verbose_debug)write(outu,'(a,a)') "PSF file has label ",psf_flags(1)
   endif

end subroutine read_psf_flags


function has_psf_flag(flag)
!Check to see if a supplied flag is in the psf_flags array:
!This is used to dictate the program's execution direction during
!the parsing of the PSF file

   logical :: has_psf_flag !Returns true if flag exists.
   character(len=*), intent(in) :: flag
   integer :: i

   has_psf_flag = .false.

   if( .not. allocated(psf_flags) ) then
     write(outu,'(a)') "FATAL ERROR in has_psf_flag: psf_flags is not been allocated."
     call mexit(outu,1)
   endif

   if (verbose_debug) write(outu,'(a,a)')"Checking for flag: ",flag
   do i=1,size(psf_flags)
     if (verbose_debug) write(outu,'(a,a)')"Examining potential match: ", psf_flags(i)
     if(index(psf_flags(i),flag) == 1) then
       has_psf_flag = .true.
       if (verbose_debug) write(outu,'(a,a)')flag," Found"
       return
     endif
   enddo 

   return
end function has_psf_flag

!-----------------------------------------------------------------
!------------- SET_FORMATS --------------------------------
!-----------------------------------------------------------------
subroutine set_formats
   implicit none

   if (  has_psf_flag("EXT") ) then
      !http://www.charmm.org/package/changelogs/c31log.html
      write(6,'(a)') "EXTented file format is being used"
      idLen=8 !Length of the CHARMM atom type name
      fmt00='(/I10)'
      fmt01='(I10,1X,A8,1X,I8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8)'
      fmt03='(8I10)'
      fmt04='(9I10)'
      fmt05='(/2I10)'
      fmt06='(2I10,3X,L1,3G14.6)'
      fmtCor='(I10,I10,2X,A8,2X,A8,3F20.10,2X,A8,2X,A8,F20.10)'
   else
      idLen=4 !Length of the CHARMM atom type name
      fmt00='(/I8)'
      fmt01='(I8,1X,A4,1X,I4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8)'
      fmt03='(8I8)'
      fmt04='(9I8)'
      fmt05='(/2I8)'
      fmt06='(2I8,3X,L1,3G14.6)'
      fmtCor='(I5,I5,1X,A4,1X,A4,3F10.5,1X,A4,1X,A4,F10.5)'
   endif
   if(l_xplor_psf) then
      if (  has_psf_flag("EXT") ) then
         fmt01='(I10,1X,A8,1X,I8,1X,A8,1X,A8,1X,a6,1X,2G14.6,I8)'
      else
         fmt01='(I8,1X,A4,1X,I4,1X,A4,1X,A4,1X,a4,1X,2G14.6,I8)'
      endif
   endif
   return
end subroutine set_formats

!-----------------------------------------------------------------
!---------- READ_TITLE ------------------------------------------
!-----------------------------------------------------------------
subroutine read_title
   implicit none
   integer :: ntitle,i

   read(psf_unit,fmt00)ntitle
   read(psf_unit,'(a)')title
   if(verbose_debug)write(outu,'(2a)')"PSF Title: ",title
   do i=2,ntitle
      read(psf_unit,'(a)')title2
      if(verbose_debug)write(outu,'(2a)')"PSF Title: ",title2
   enddo

   return
end subroutine read_title

!-----------------------------------------------------------------
!--------------- ALLINT -------------------------------------------
!-----------------------------------------------------------------
subroutine allint(size,a0,a1,a2,a3,a4,a5,a6,a7,a8)
   implicit none
   integer, intent(in) :: size
   integer, optional, dimension(:), pointer :: &
                     a0,a1,a2,a3,a4,a5,a6,a7,a8
   integer statall

   if(.not. present(a0)) return
   allocate(a0(size),stat=statall)
   if ( statall /= 0 ) call allocate_error("a0","allint",size)

   if(.not. present(a1)) return
   allocate(a1(size),stat=statall)
   if ( statall /= 0 ) call allocate_error("a1","allint",size)

   if(.not. present(a2)) return
   allocate(a2(size),stat=statall)
   if ( statall /= 0 ) call allocate_error("a2","allint",size)

   if(.not. present(a3)) return
   allocate(a3(size),stat=statall)
   if ( statall /= 0 ) call allocate_error("a3","allint",size)

   if(.not. present(a4)) return
   allocate(a4(size),stat=statall)
   if ( statall /= 0 ) call allocate_error("a4","allint",size)

   if(.not. present(a5)) return
   allocate(a5(size),stat=statall)
   if ( statall /= 0 ) call allocate_error("a5","allint",size)

   if(.not. present(a6)) return
   allocate(a6(size),stat=statall)
   if ( statall /= 0 ) call allocate_error("a6","allint",size)

   if(.not. present(a7)) return
   allocate(a7(size),stat=statall)
   if ( statall /= 0 ) call allocate_error("a7","allint",size)

   if(.not. present(a8)) return
   allocate(a8(size),stat=statall)
   if ( statall /= 0 ) call allocate_error("a8","allint",size)

   return
 end subroutine allint

!----------------------------------------------------------------- 
!--------------- ALLREAL -------------------------------------------
!-----------------------------------------------------------------

subroutine allreal(size,a0,a1,a2,a3,a4,a5,a6,a7,a8)
   implicit none
   integer, intent(in) :: size
   real(kind=8), optional, dimension(:), pointer :: &
                     a0,a1,a2,a3,a4,a5,a6,a7,a8
   integer statall

   if(.not. present(a0)) return
   allocate(a0(size),stat=statall)
   if ( statall /= 0 ) call allocate_error("a0","allreal",size)

   if(.not. present(a1)) return
   allocate(a1(size),stat=statall)
   if ( statall /= 0 ) call allocate_error("a1","allreal",size)

   if(.not. present(a2)) return
   allocate(a2(size),stat=statall)
   if ( statall /= 0 ) call allocate_error("a2","allreal",size)

   if(.not. present(a3)) return
   allocate(a3(size),stat=statall)
   if ( statall /= 0 ) call allocate_error("a3","allreal",size)

   if(.not. present(a4)) return
   allocate(a4(size),stat=statall)
   if ( statall /= 0 ) call allocate_error("a4","allreal",size)

   if(.not. present(a5)) return
   allocate(a5(size),stat=statall)
   if ( statall /= 0 ) call allocate_error("a5","allreal",size)

   if(.not. present(a6)) return
   allocate(a6(size),stat=statall)
   if ( statall /= 0 ) call allocate_error("a6","allreal",size)

   if(.not. present(a7)) return
   allocate(a7(size),stat=statall)
   if ( statall /= 0 ) call allocate_error("a7","allreal",size)

   if(.not. present(a8)) return
   allocate(a8(size),stat=statall)
   if ( statall /= 0 ) call allocate_error("a8","allreal",size)

   return
end subroutine allreal
     

!-----------------------------------------------------------------
!--------- DEALLOCATE EVERYTHING -----------------------------------
!-----------------------------------------------------------------
subroutine deallocate_everything
   use cmap, only: deallocate_cmap
   implicit none
   integer :: ierr
   logical :: v

   v = verbose_debug

   if(v) write(outu,'(a)')"Deallocating atom_label"
   if(allocated(atom_label)) then
     deallocate(atom_label,stat=ierr)
     if (ierr /= 0) then
       write(outu,'(a)') "ERROR Deallocating lres."
       call mexit(outu,1)
     end if
   end if

   if(v) write(outu,'(a)')"Deallocating atom_type"
   if(allocated(atom_type)) then
     deallocate(atom_type,stat=ierr)
     if (ierr /= 0) then
       write(outu,'(a)') "ERROR Deallocating lres."
       call mexit(outu,1)
     end if
   end if

   if(v) write(outu,'(a)')"Deallocating lres"
   if(allocated(lres)) then
     deallocate(lres,stat=ierr)
     if (ierr /= 0) then
       write(outu,'(a)') "ERROR Deallocating lres."
       call mexit(outu,1)
     end if
   end if

   if(v) write(outu,'(a)')"Deallocating lresat"
   if(allocated(lresat)) then
     deallocate(lresat,stat=ierr)
     if (ierr /= 0) then
       write(outu,'(a)') "ERROR Deallocating lresat."
       call mexit(outu,1)
     end if
   end if

   if(v) write(outu,'(a)')"Deallocating attype_name"
   if(allocated(attype_name)) then
     deallocate(attype_name,stat=ierr)
     if (ierr /= 0) then
       write(outu,'(a)') "ERROR Deallocating attype_name."
       call mexit(outu,1)
     end if
   end if

   if(v) write(outu,'(a)')"Deallocating element"
   if(allocated(element)) then
     deallocate(element,stat=ierr)
     if (ierr /= 0) then
       write(outu,'(a)') "ERROR Deallocating element."
       call mexit(outu,1)
     end if
   end if

   if(v) write(outu,'(a)')"Deallocating psf_flags"
   if(allocated(psf_flags)) then
     deallocate(psf_flags,stat=ierr)
     if (ierr /= 0) then
       write(outu,'(a)') "ERROR Deallocating psf_flags."
       call mexit(outu,1)
     end if
   end if
   call deallocate_int(iac,imove)
   call deallocate_real(cg,amass)

   !Bonds
   call deallocate_int(icph,icpa)
   call deallocate_int(icth,icta,bndtypes,ibtyp,jbtyp,icb)
   call deallocate_int(ibh,jbh,iba,jba,attypes,attypeh,locattype)
   call deallocate_int(ib,jb,it,jt,kt)
   call deallocate_int(ibh,jbh,iba,jba,icbh,icba)

   !Angles
   call deallocate_int(angtypes,ith,jth,kth,ita,jta,kta)
   call deallocate_real(tk,teq)
   call deallocate_int(ittyp,jttyp,kttyp)
   call deallocate_int(ict)
   call deallocate_int(iph,jph,kph,lph,ipa,jpa,kpa,lpa)
   call deallocate_int(iph2,jph2,kph2,lph2,icph2)
   call deallocate_int(ipa2,jpa2,kpa2,lpa2,icpa2)


   !Dihedrals
   call deallocate_int(dihtypes,icp)
   call deallocate_int(ip,jp,kp,lp,im,jm,km,lm)
   call deallocate_int(idtyp,jdtyp,kdtyp,ldtyp)
   call deallocate_real(pk,pn,phase)

   call deallocate_int(idon,ihd1,iacc,iac1,inb,iblo)
   call deallocate_int(igpbs,igptyp,imoveg,iresid,ipres)

   !Impropers
   call deallocate_int(imp)
   call deallocate_int(impropertypes,imtyp,jmtyp,kmtyp,lmtyp)
   call deallocate_real(pk_impr,phase_impr)

   !UB
   call deallocate_int(ub_atm_i,ub_atm_k,ub_idx)
   call deallocate_real(rub,kub)

   !CMAP
   call deallocate_cmap

   !Non bonded
   call deallocate_real(cn1,cn2,rk,req)
   call deallocate_real(cn114,cn214)
   call deallocate_real(radii,screen)
   call deallocate_real(scnb_scale_factor,scee_scale_factor)


end subroutine deallocate_everything
     
!-----------------------------------------------------------------
!--------- DEALLOCATE INT -----------------------------------
!-----------------------------------------------------------------
subroutine deallocate_int(a0,a1,a2,a3,a4,a5,a6,a7,a8)
   implicit none

   integer, optional, dimension(:), pointer :: &
                     a0,a1,a2,a3,a4,a5,a6,a7,a8
   integer statall

   if(.not. present(a0)) return
   if(associated(a0)) then
      deallocate(a0,stat=statall)
      if ( statall /= 0 ) call deallocate_error("a0","deallocate_int")
   endif

   if(.not. present(a1)) return
   if(associated(a1)) then
      deallocate(a1,stat=statall)
      if ( statall /= 0 ) call deallocate_error("a1","deallocate_int")
   endif

   if(.not. present(a2)) return
   if(associated(a2)) then
      deallocate(a2,stat=statall)
      if ( statall /= 0 ) call deallocate_error("a2","deallocate_int")
   endif

   if(.not. present(a3)) return
   if(associated(a3)) then
      deallocate(a3,stat=statall)
      if ( statall /= 0 ) call deallocate_error("a3","deallocate_int")
   endif

   if(.not. present(a4)) return
   if(associated(a4)) then
      deallocate(a4,stat=statall)
      if ( statall /= 0 ) call deallocate_error("a4","deallocate_int")
   endif

   if(.not. present(a5)) return
   if(associated(a5)) then
      deallocate(a5,stat=statall)
      if ( statall /= 0 ) call deallocate_error("a5","deallocate_int")
   endif

   if(.not. present(a6)) return
   if(associated(a6)) then
      deallocate(a6,stat=statall)
      if ( statall /= 0 ) call deallocate_error("a6","deallocate_int")
   endif

   if(.not. present(a7)) return
   if(associated(a7)) then
      deallocate(a7,stat=statall)
      if ( statall /= 0 ) call deallocate_error("a7","deallocate_int")
   endif

   if(.not. present(a8)) return
   if(associated(a8)) then
      deallocate(a8,stat=statall)
      if ( statall /= 0 ) call deallocate_error("a8","deallocate_int")
   endif


end subroutine deallocate_int

subroutine deallocate_real(a0,a1,a2,a3,a4,a5,a6,a7,a8)
   implicit none
   real(kind=8), optional, dimension(:), pointer :: &
                     a0,a1,a2,a3,a4,a5,a6,a7,a8
   integer statall

   if(.not. present(a0)) return
   if(associated(a0)) then
      deallocate(a0,stat=statall)
      if ( statall /= 0 ) call deallocate_error("a0","deallocate_int")
   endif

   if(.not. present(a1)) return
   if(associated(a1)) then
      deallocate(a1,stat=statall)
      if ( statall /= 0 ) call deallocate_error("a1","deallocate_int")
   endif

   if(.not. present(a2)) return
   if(associated(a2)) then
      deallocate(a2,stat=statall)
      if ( statall /= 0 ) call deallocate_error("a2","deallocate_int")
   endif

   if(.not. present(a3)) return
   if(associated(a3)) then
      deallocate(a3,stat=statall)
      if ( statall /= 0 ) call deallocate_error("a3","deallocate_int")
   endif

   if(.not. present(a4)) return
   if(associated(a4)) then
      deallocate(a4,stat=statall)
      if ( statall /= 0 ) call deallocate_error("a4","deallocate_int")
   endif

   if(.not. present(a5)) return
   if(associated(a5)) then
      deallocate(a5,stat=statall)
      if ( statall /= 0 ) call deallocate_error("a5","deallocate_int")
   endif

   if(.not. present(a6)) return
   if(associated(a6)) then
      deallocate(a6,stat=statall)
      if ( statall /= 0 ) call deallocate_error("a6","deallocate_int")
   endif

   if(.not. present(a7)) return
   if(associated(a7)) then
      deallocate(a7,stat=statall)
      if ( statall /= 0 ) call deallocate_error("a7","deallocate_int")
   endif

   if(.not. present(a8)) return
   if(associated(a8)) then
      deallocate(a8,stat=statall)
      if ( statall /= 0 ) call deallocate_error("a8","deallocate_int")
   endif

end subroutine deallocate_real

!-----------------------------------------------------------------
!         ERR_EXIT
!-----------------------------------------------------------------
subroutine err_exit(string)
  implicit none
  character(len=*),intent(in)::string
  write(outu,*)string
  call exit(1)
  stop
end subroutine err_exit

!-----------------------------------------------------------------
!--------------------- ALLOCATE_ERROR ------------------------------
!-----------------------------------------------------------------
subroutine allocate_error(a0,a1,size)
   implicit none
   character(len=*), intent(in) :: a0,a1
   integer, intent(in) :: size

   write(outu,"(a,a,a,a,i8)")"ERROR allocating variable ",a0, &
                                               " in routine ",a1,size
   call err_exit("Cannot continue")
end subroutine allocate_error
subroutine deallocate_error(a0,a1)
   implicit none
   character(len=*), intent(in) :: a0,a1

   write(outu,"(a,a,a,a,i8)")"ERROR deallocating variable ",a0, &
                                               " in routine ",a1
   call err_exit("Cannot continue")
end subroutine deallocate_error

!-----------------------------------------------------------------
!--------------------- TAG_14_DUPLICATES- ------------------------
!-----------------------------------------------------------------
subroutine tag_14_duplicates(numdih,i_dihed,k_dihed,l_dihed, &
                             numbond,i_bond,j_bond,&
                             numangle,i_angle,k_angle)
 
 !This subroutine is responsible for identifying and tagging all
 !duplicate 1-4 interactions. In CHARMM 1-4 EEL and VDW interactions
 !are not scaled and so they are done as part of the standard EEL and
 !VDW code. In AMBER 1-4's are done separately by looping over the 
 !unique set of the dihedral list. This is to allow for the SCEE and
 !SCNB scaling which are not strictly needed for the CHARMM force field
 !however it makes sense to still do it this way since CHARMM uses
 !different A and B coefficients for the VDW terms of 1-4's and it
 !would be tricky to deal with this inside the regular 1-4 code in AMBER.
 !
 !To avoid multiple counting of 1-4 non-bond interactions LEaP flags additional
 !terms with a -ve value in either the 3rd or 4th integer of the DIHEDRAL list. 
 !Atoms K and L. A -ve value on the 4th integer implies an improper dihedral
 !whereas a negative value on the 3rd integer implies a duplicate or part of a
 !ring system.
 !
 ! Example duplicates:
 ! 
 ! CT-CT-CT-CT   1.0  5.0   180.0   1
 ! CT-CT-CT-CT   1.0 10.0     0.0   2
 !
 ! This would have the atom K number for the second dihedral set to -ve so it 
 ! doesn't have the 1-4 done twice. 
 !
 ! Duplicates in rings:
 !
 !   1 --- 2
 !   /      \
 !  /        \
 ! 5          3
 !  \---4----/ 
 !
 ! Here the dihedral 1-2-3-4 gives the 1-4 NB interaction between atoms 1 and 4.
 ! however 1 and 4 are currently part of an angle term and so should not be counted
 ! as a 1-4 term and need the 3rd atom (3 here) set to -ve.
 !
 ! Note atom j of the dihedral is not actually needed anywhere in this subroutine, hence
 ! the absence of the j_dihed array.

 implicit none

 !Passed in
 integer, intent(in) :: numdih
 integer, intent(inout) :: i_dihed(numdih), k_dihed(numdih) ,l_dihed(numdih)
 integer, intent(in) :: numbond
 integer, intent(in) :: i_bond(numbond),j_bond(numbond)
 integer, intent(in) :: numangle
 integer, intent(in) :: i_angle(numangle), k_angle(numangle)

 !Local
 integer :: i, j, current_i, current_l

 ! Step 1 - Check for duplicates
 !          We check to see if the current list of dihedrals
 !          contains multiple copies of i -- x -- x -- l
 !          In other words we just check the first and last atom numbers that
 !          make up the 1-4. This deals with multiple phase angles as well as other
 !          duplicates due to various ring structures etc.
  
 ! Loop through each dihedral in turn and check if any other dihedrals match the
 ! same i and l values.
 ! Note there is some duplication in this approach due to us processing already
 ! processed dihedrals (e.g. if dihedral 1 matches 5 and 7 we set 5 and 7 -ve when we
 ! process dihedral 1 but subsequently we set 7 negative again when we process dihedral 5.
 ! This could be improved but is done this way at present for readability.

 do i = 1, numdih
   current_i = i_dihed(i)
   current_l = l_dihed(i)

   do j = i+1,numdih
     !Process all subsequent dihedrals.
     !Check to see if i and l match this dihedral.
     if ( (current_i == i_dihed(j) .and. current_l == l_dihed(j)) &
          .or. &
          (current_i == l_dihed(j) .and. current_l == i_dihed(j)) ) then

        k_dihed(j) = sign(k_dihed(j),-1) !We have a duplicate; set k sign -ve
     end if  
     !else we just skip this and leave it set +ve.
   end do
 end do


 ! Step 2 - Ring systems.
 ! Dealing with 1-4's in ring systems is tricky because atoms can have dihedrals and thus be
 ! part of a 1-4 but also be bonded or parts of angle terms.
 ! E.g.
 !
 ! 1---2
 ! |   |
 ! 4---3
 !
 ! Here we have the following dihedrals:
 ! 1-2-3-4
 ! 2-3-4-1
 ! 3-4-1-2
 ! 4-1-2-3
 !
 ! In all cases the i and l atoms are bonded to each other so need to be skipped in the 1-4's.
 ! Additionally you need to be careful of angles.
 ! E.g. for Proline:
 !
 ! N1--C1
 ! |   |
 ! C4  C2
 !  \ /
 !  C3
 !
 ! Here the dihedral N1-C1-C2-C3 should not include the 1-4 because N1 and C3 are the end points of an
 ! angle term.
 !
 ! Loop over dihedrals and check against the bond list and angle list.
 
 do i = 1, numdih
   current_i = i_dihed(i)
   current_l = l_dihed(i)
   !Check if these two atoms are bonded.
   do j = 1, numbond
     if ( (current_i == i_bond(j) .and. current_l == j_bond(j)) &
        .or. &
         (current_i == j_bond(j) .and. current_l == i_bond(j)) ) then
        k_dihed(i) = sign(k_dihed(i),-1)
     !else we just skip this and leave it set +ve.
     end if
   end do
   !Check if these two atoms are the end points of an angle term.
   do j = 1, numangle
     if ( (current_i == i_angle(j) .and. current_l == k_angle(j)) &
        .or. &
         (current_i == k_angle(j) .and. current_l == i_angle(j)) ) then
        k_dihed(i) = sign(k_dihed(i),-1)
     !else we just skip this and leave it set +ve.
     end if
   end do
 end do 

 return 

end subroutine tag_14_duplicates


integer function chemicalSymbolToZ(symbol) result(Z)
  !Given a two character element symbol, it will return its Z number
  !This is generally used to populate the psfprm::element() array

  use psf_strings, only: strUpCase  !Convert everything to uppercase
  implicit none

  character(2),intent(in) :: symbol
!  integer                 :: Z

  
  select case ( StrUpCase(symbol) )
    case ("H ")
      Z = 1
    case ("HE")
      Z = 2
    case ("LI")
      Z = 3
    case ("BE")
      Z = 4
    case ("B ")
      Z = 5
    case ("C ")
      Z = 6
    case ("N ")
      Z = 7
    case ("O ")
      Z = 8
    case ("F ")
      Z = 9
    case ("NE")
      Z = 10
    case ("NA")
      Z = 11
    case ("MG")
      Z = 12
    case ("AL")
      Z = 13
    case ("SI")
      Z = 14
    case ("P ")
      Z = 15
    case ("S ")
      Z = 16
    case ("CL")
      Z = 17
    case ("AR")
      Z = 18
    case ("K ")
      Z = 19
    case ("CA")
      Z = 20
    case ("ZN")
      Z = 30
    case ("CU")
      Z = 29




    case default
      Z = 0

  end select
  if(verbose_debug)write(outu,*)"Symbol ",symbol,"  number ", Z
  if (Z == 0) then
! gfortran is warning with -Wall here:
!
! _psfprm.f: In function ‘chemicalsymboltoz’:
! _psfprm.f:3637: warning: ‘z’ may be used uninitialized in this function
!
! is bogus:
! http://gcc.gnu.org/bugzilla/show_bug.cgi?id=46979

    write(outu,*) ""
    write(outu,*) "chemicalSymbolToZ() cannot map symbol: ",symbol,", to a Z value"
    write(outu,*) ""
    write(outu,*) "NOTE: If the symbol reported above does not appear to be a valid"
    write(outu,*) "      element then this most likely means your topology file does not"
    write(outu,*) "      include element information. The code expects MASS lines in"
    write(outu,*) "      the topology file to have the following format:"
    write(outu,*) ""
    write(outu,*) "MASS     4 HT     1.00800  H ! TIP3P water hydrogen"
    write(outu,*) ""
    write(outu,*) "      your file most likely has"
    write(outu,*) ""
    write(outu,*) "MASS     4 HT     1.00800"
    write(outu,*) ""
    write(outu,*) "      you should check this before proceeding."

    call err_exit("chemicalSymbolToZ() failed to assign a Z value to an element")
  endif

end function chemicalSymbolToZ

subroutine reorder_dihedrals_to_remove_0_k(numdih,i_dihed,j_dihed,k_dihed,l_dihed)

!This subroutine reverses the four integers that make up the torsion index
!so that chamber would never write a negative zero on the k'th index:
!
!
!Within the sections:
!
!        DIHEDRALS_INC_HYDROGEN
!        DIHEDRALS_WITHOUT_HYDROGEN
!
!in a prmtop file, dihedrals are listed for a system being described. Each
!dihedral is described using five numbers, the first four relate to the
!indexes of the atoms in the dihedral and the fifth is a "pointer" into the
!DIHEDRAL_FORCE_CONSTANT, DIHEDRAL_PERIODICITY and DIHEDRAL_PHASE lists
!giving the associated parameters for that specific dihedral.

!The first four numbers are not (as one would expect) the atom index (given
!by the order in ATOM_NAME) of each atom in the dihedral, but instead the
!atom index number minus one, multiplied by three.
!
!For example, the dihedral listed in a prmtop as follows:
!
!%FLAG DIHEDRALS_WITHOUT_HYDROGEN
!%FORMAT(10i8)
!       0      12      18      24       2

!would correspond to torsion 1-5-7-9 and point to number two in the
!associated DIHEDRAL_* parameter lists.
!
!AMBER denotes an IMPROPER torsion in its dihedral indexes by prefixing the
!4th atom (l) with a negative sign. It does the same with the 3rd atom (k)
!to indicate that the end group interactions (i.e. any 1-4 interactions) of
!this torsion should be ignored.

! A problem can occur when k is atom one in the system and it needs
!to be flagged to ignore its 1-4 contribution, e.g. :
!
!Given a dihedral of atom index: 4-3-1-7
!
!it becomes the prmtop'ed form of:  
!9           6           0          18
!
!
!Here in lies the problem, you cannot put a sign on a zero:
!
!9           6          -0          18

!I (MJW) had initially looked (for a while) into using IEEE 754 to place 
!a sign on the integer zero, but there seemed to be a inconsistent 
!approach to formatted negative zero printing in current compilers.

!This subroutine reverses the four integers that make up the torsion index 
!so that chamber would never write a negative zero on the k'th index:

!Hence a problematic torsion:
!
!9           6          -0          18           1
!
!would be reversed and written to the prmtop as be:
!
!18          0          -6          9            1 


 implicit none

 !Passed in
 integer, intent(in)    :: numdih
 integer, intent(inout) :: i_dihed(numdih) ,j_dihed(numdih),&
                           k_dihed(numdih) ,l_dihed(numdih)
 !local
 integer                :: i
 integer                :: tmp_i, tmp_j, tmp_k, tmp_l

 do i = 1, numdih
        if(k_dihed(i) == 0) then
           write(6,'(a72)')"Possible issue with assigning a -ve sign to zero valued k dihedral index"
           write(6,'(a35)')"Hence, reversing dihedral order"

           write(6,'(a20,i8)') "Dihedral number ",i
           write(6,'(a20,4i8)') "Old internal index",i_dihed(i),j_dihed(i),k_dihed(i),l_dihed(i)
           write(6,'(a20,4i8)') "Old i,j,k,l",(i_dihed(i)/3)+1,(j_dihed(i)/3)+1,(k_dihed(i)/3)+1,(l_dihed(i)/3)+1

           !Reverse dihedral to avoid negative zero on k
           tmp_i = i_dihed(i)
           tmp_j = j_dihed(i)
           tmp_k = k_dihed(i)
           tmp_l = l_dihed(i)

           i_dihed(i) = tmp_l
           j_dihed(i) = tmp_k
           k_dihed(i) = tmp_j
           l_dihed(i) = tmp_i

           write(6,'(a1)') ""
           write(6,'(a20,4i8)') "New internal index",i_dihed(i),j_dihed(i),k_dihed(i),l_dihed(i)
           write(6,'(a20,4i8)') "New i,j,k,l",(i_dihed(i)/3)+1,(j_dihed(i)/3)+1,(k_dihed(i)/3)+1,(l_dihed(i)/3)+1

           write(6,'(a1)') ""
           write(6,'(a1)') ""
        endif
  enddo

end subroutine reorder_dihedrals_to_remove_0_k

end module psfprm
