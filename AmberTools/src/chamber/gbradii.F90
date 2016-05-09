! <compile=optimized>
#include "copyright.h"

!--------------------------------------------------------------
! CHAMBER: GB RADII ROUTINES
!
! Description: Routines for processing the Generalized Born
!              Radii.
!
! Authors: Mark Williamson (SDSC, 2009)
!          Michael Crowley (NREL, 2009)
!          Ross Walker (SDSC, 2009)
!
!--------------------------------------------------------------

subroutine get_gbradii(rs)
  !This routine assigns the set of screen(natom) and radii(natom)
  !parameters to all atoms within the system and this assignment 
  !is a function of the radius_set (rs) passed to this subroutine.
  !
  !The the logic here is from the code in
  ! leap/src/leap/unitio.c:5183 onwards
  use psfprm, only     : natom,radii,screen,outu,verbose_debug,&
                         attype_name,locattype,err_exit,&
                         element, radius_set_name,&
                         nbond, ib,jb

  use molnt, only      : bonded_list
  implicit none
  integer, intent(in) :: rs !selected_radius_set
  integer             :: i, currently_bonded_atom
  
  radii=0.d0
  screen=0.d0

  !Set the identifier string that will go into the prmtop
  select case (rs)
    case (0)
      radius_set_name = "Bondi radii (bondi)"
    case (1)
      radius_set_name = "amber6 modified Bondi radii (amber6)"
    case (2)
      radius_set_name = "modified Bondi radii (mbondi)"
      !http://dx.doi.org/10.1002/1097-0282(2000)56:4<275::AID-BIP10024>3.0.CO;2-E

    !Not supported
    !(3) Huo and Kollman optimized radii (pbamber)
    !(4) Jayaram et al. 'GB'
    !(5) Jayaram et al. 'MGB'

    case (6)
      radius_set_name = "H(N)-modified Bondi radii (mbondi2)"
      !http://dx.doi.org/10.1002/prot.20033
  end select

  if(verbose_debug) then 
    write(outu,*) "==========GET_GBRADII =============================="
    write(outu,*) "Called with radius set: ",rs,":",radius_set_name
  endif


  do i=1,natom

     select case( element( locattype(i) ) ) !Get the element Z value of each natom
     case (1) ! Atom is hydrogen
        radii(i)=1.2d0
        screen(i)=0.85d0

        ! The radii for hydrogen is a function of what is it bonded to,
        ! hence check this by looking at the bonded_list for that atom

        ! Look *only* at atom bonded to i, 
        currently_bonded_atom = element ( locattype( bonded_list(1, i) ) ) 
        ! e.g.        bonded_list(1:5, 4) is 1,2,3,0,0
        ! indicates that atom 4 is bonded to atoms 1, 2, and 3.

        if (rs == 1 .or. rs == 2) then

          select case ( currently_bonded_atom )
            case(6)
              !Carbon
              radii(i)=1.3d0
            case(7)
              !Nitrogen 
              radii(i)=1.3d0
            case(8)
              !Oxygen
              radii(i)=0.8d0
            case(16)
              !Sulphur
              radii(i)=0.8d0
          end select


        else if(rs == 6) then

          select case ( currently_bonded_atom )
            case(7)
              !Nitrogen 
              radii(i)=1.3d0
          end select

        endif !(rs == 1 .or. rs == 2)

     case (6) ! Atom is carbon
        radii(i)=1.7d0
        screen(i)=0.72d0
     !leap's original behaviour is ignored here since CHARMM does not use united atoms:
     !! if ( strncmp(sType,"C1",2) && strncmp(sType,"C2",2) && strncmp (sType,"C3",2) )
     !!      dGBrad = 1.7;
     !!    else
     !!       dGBrad = 2.2;
     !!    break;

     case (7) !N
        radii(i)=1.55d0
        screen(i)=0.79d0
     case (8) !O
        radii(i)=1.5d0
        screen(i)=0.85d0
     case (9) !F
        radii(i)=1.5d0
        screen(i)=0.88d0
     case (15)!P
        radii(i)=1.85d0
        screen(i)=0.86d0
     case (16)!S
        radii(i)=1.8d0
        screen(i)=0.96d0
     case (17)!Cl
        radii(i)=1.7d0
        screen(i)=0.8d0

     case default
        radii(i)=1.5d0
        screen(i)=0.8d0
     end select

  enddo
   !----------------Screen code from leap ------------------------------
       ! if( GDefaults.iGBparm < 4 || GDefaults.iGBparm == 6 ){
       !     /* for now, hardwire the Bondi radii  */
       !     switch( iElement ){
       !         case 1:  dGBscreen = 0.85; break;
       !         case 6:  dGBscreen = 0.72; break;
       !         case 7:  dGBscreen = 0.79; break;
       !         case 8:  dGBscreen = 0.85; break;
       !         case 9:  dGBscreen = 0.88; break;
       !         case 15: dGBscreen = 0.86; break;
       !         case 16: dGBscreen = 0.96; break;
       !         default: dGBscreen = 0.8; break;  /* or should fail?? */
       !     }
   !----------------Radii code from leap ------------------------------
       !         case  1: dGBrad = 1.2; 
       !             /* make the modifications that hydrogen radii
       !                depend upon the atoms they are bonded to.  
       !                iGBparm=1 corresponds to Amber 6, JACS 122:2489 (2000);
       !                iGBparm=2 adds the update of Biopolymers 56: 275 (2001)   
       !              */
       !             if( iAtomCoordination(saPAtom->aAtom) > 0 ) {
       !                 /* For multiply bonded Hydrogen atoms use the first
       !                  * bond for determining modified GB radii.
       !                  * WAT contains multiply bonded Hydrogen atoms 
       !                  * so do not emit a warning.
       !                  */
       !
       !                 aAtomA = aAtomBondedNeighbor(saPAtom->aAtom, 0);
       !                 if( GDefaults.iGBparm == 1 || GDefaults.iGBparm == 2 ) {
       !                     if( sAtomType(aAtomA)[0] == 'C' ||
       !                             sAtomType(aAtomA)[0] == 'c' ) dGBrad = 1.3;
       !                     if( sAtomType(aAtomA)[0] == 'O' ||
       !                             sAtomType(aAtomA)[0] == 'o' ) dGBrad = 0.8;
       !                     if( sAtomType(aAtomA)[0] == 'S' ||
       !                             sAtomType(aAtomA)[0] == 's' ) dGBrad = 0.8;
       !                     if( (sAtomType(aAtomA)[0] == 'N' ||
       !                             sAtomType(aAtomA)[0] == 'n')  &&
       !                             GDefaults.iGBparm == 2) dGBrad = 1.3;
       !                 }
       !                 else if( GDefaults.iGBparm == 6 ) { 
       !                     /* try Alexey's scheme */
       !                     if( sAtomType(aAtomA)[0] == 'N' ||
       !                              sAtomType(aAtomA)[0] == 'n' ) dGBrad = 1.3;
       !                 }
       !             }
       !             else {
       !                 VP0(( "WARNING: Unbonded Hydrogen atom %s in %s.\n"
       !                         " Cannot determine the requested GB radius"
       !                         " for this atom.\n"
       !                         " Writing the unmodified Bondi GB radius.\n",
       !                         saPAtom->aAtom->cHeader.sName,
       !                         saPAtom->aAtom->cHeader.cContainedBy->sName ));
       !             }
       !             break;
       !         case  6: dGBrad = 1.7; break;
       !         case  7: dGBrad = 1.55; break;
       !         case  8: dGBrad = 1.5; break;
       !         case  9: dGBrad = 1.5; break;
       !         case 14: dGBrad = 2.1; break;
       !         case 15: dGBrad = 1.85; break;
       !         case 16: dGBrad = 1.8; break;
       !         case 17: dGBrad = 1.7; break;
       !         default: dGBrad = 1.5; break;

  return
end subroutine get_gbradii
