! <compile=optimized>
#include "copyright.h"
#include "is_copyright.h"
#include "dprec.h"
#include "is_def.h"

module solvent_accessibility
    
   implicit none
    
   ! SAS/SAR control variables
    
   integer              :: radiopt       ! radi option
   integer              :: maxsph        ! max no of surface dots per atom
   integer              :: maxarc        ! max no of arcs per atom
   integer              :: maxarcdot     ! max no of arc dots per atom
   integer              :: nsrfdot       ! no of sa surface dots
   integer              :: narcdot       ! no of sa arc dots
   integer              :: narc          ! no of sa arcs
   _REAL_              :: minsph        ! min no of surface dots per atom
   _REAL_              :: sprob         ! solvent probe
   _REAL_              :: arcres        ! arc dot resolution, with respect to the grid spacing
   _REAL_              :: prtsas        ! total sa surface area for printout
    
   ! Atomic radii arrays
    
   _REAL_, allocatable ::   mdsig(:)    ! atomic radii (np vdw radii)
   _REAL_, allocatable ::    radi(:)    ! atomic radii (pb cavity radii)
   _REAL_, allocatable ::   radip(:)    ! atomic radii + sprob
   _REAL_, allocatable ::  radip2(:)    ! squared radip
    
   ! SAS/SAR arrays
    
   _REAL_, allocatable ::    scrd(:,:)  ! coordinates of maxsph evenly distributed dots on the unit sphere
   integer, allocatable ::  fstsdot(:)   ! first surface dot of an atom
   integer, allocatable ::  lstsdot(:)   ! last surface dot of an atom
   integer, allocatable ::  fstadot(:)   ! first arc dot of an atom
   integer, allocatable ::  lstadot(:)   ! last arc dot of an atom
   integer, allocatable ::   fstarc(:)   ! first arc of an atom
   integer, allocatable ::   lstarc(:)   ! last arc of an atom
   integer, allocatable ::     marc(:)   ! no. of arcs of an atom
   integer, allocatable ::   m2narc(:,:) ! marc to narc index for an atom
   integer, allocatable ::   arcatm(:,:) ! two atoms that generate an arc
   integer, allocatable ::   dotarc(:)   ! arc index of an arc dot
   _REAL_, allocatable ::   savarc(:,:) ! saved arc geometries
   _REAL_, allocatable ::  savactr(:,:) ! saved arc center coordinates
   _REAL_, allocatable ::   srfcrd(:,:) ! surface dot coordinates
   _REAL_, allocatable ::   arccrd(:,:) ! arc dot coordinates

   logical , allocatable :: knockout(:)  ! auxiliary flags for checking surface/arcs
   integer , allocatable ::   spharc(:)  ! auxiliary arc index
   integer , allocatable ::  spharc1(:)  ! auxiliary arc index
   _REAL_, allocatable ::   sphcrd(:,:) ! auxiliary surface/arc dot coordinates
   _REAL_, allocatable ::  sphcrd1(:,:) ! auxiliary surface/arc dot coordinates

   ! MVDW surface setting

   _REAL_              :: radinc
   _REAL_              :: expthresh
   _REAL_, allocatable ::  radipd(:)
   _REAL_, allocatable ::    nmax(:)
   _REAL_, allocatable ::    nexp(:)
   _REAL_, allocatable :: sumnmax(:)
   _REAL_, allocatable :: sumnexp(:)
   _REAL_, allocatable ::  avnmax(:)
   _REAL_, allocatable ::  avnexp(:)

   ! Spline surface setting

   integer  :: nn                        ! desired spline stencils
   _REAL_, allocatable :: dash(:)       ! spline up to 10 stencils
   _REAL_, allocatable :: spcoef(:,:)   ! spline up to 3rd power

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Cleanup routine of solvent accessible surface and arcs
subroutine sa_free ( dosas,ndosas,dosar,ndosar )

   integer dosas, dosar, ndosas, ndosar, alloc_err(16)

   alloc_err = 0

   if ( ndosas /= 0 .and. dosas == ndosas ) then
      deallocate( fstsdot, stat = alloc_err(1 ) )
      deallocate( lstsdot, stat = alloc_err(2 ) )
      deallocate(  srfcrd, stat = alloc_err(3 ) )
   end if
   if ( ndosar /= 0 .and. dosar == ndosar ) then
      deallocate( fstadot, stat = alloc_err(4 ) )
      deallocate( lstadot, stat = alloc_err(5 ) )
      deallocate(  fstarc, stat = alloc_err(6 ) ) 
      deallocate(  lstarc, stat = alloc_err(7 ) )
      deallocate(    marc, stat = alloc_err(8 ) )
      deallocate(  m2narc, stat = alloc_err(9 ) )
      deallocate(  arcatm, stat = alloc_err(10) )
      deallocate(  savarc, stat = alloc_err(11) )
      deallocate( savactr, stat = alloc_err(12) )
      deallocate(  dotarc, stat = alloc_err(13) )
      deallocate(  arccrd, stat = alloc_err(14) )
   end if
   if ( alloc_err(1)+alloc_err(2)+alloc_err(3)+alloc_err(4)+alloc_err( 5)+&
        alloc_err(6)+alloc_err(7)+alloc_err(8)+alloc_err(9)+alloc_err(10)+&
        alloc_err(11)+alloc_err(12)+alloc_err(13)+alloc_err(14) /= 0 ) then
      write(6, *) 'SA Bomb in sa_free(): Allocation aborted', alloc_err(1:14)
      call mexit(6, 1)
   end if

end subroutine sa_free
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ NP dispersion/cavity radi setup based on nonelectrostatic solvation free energies in TIP3P
subroutine sa_init( verbose, pbprint, natom, r, rp, rp2 )

   ! Passed variables

   logical verbose, pbprint
   integer natom
   !character (len=4) :: isymbl(natom)
   _REAL_ r(natom), rp(natom), rp2(natom)

   ! Local variables
    
   integer iatm, nsatm
    
   ! for InsightII surface display only

   !do iatm = 1, natom
   !   else if ( isymbl(iatm)(1:1) == 'C' .or. isymbl(iatm)(1:1) == 'c' ) then
   !      r(iatm) = 1.55d0
   !   else if ( isymbl(iatm)(1:1) == 'H' .or. isymbl(iatm)(1:1) == 'h' ) then
   !      r(iatm) = 1.10d0
   !   else if ( isymbl(iatm)(1:1) == 'N' .or. isymbl(iatm)(1:1) == 'n' ) then
   !      r(iatm) = 1.40d0
   !   else if ( isymbl(iatm)(1:1) == 'O' .or. isymbl(iatm)(1:1) == 'o' ) then
   !      r(iatm) = 1.35d0
   !   else if ( isymbl(iatm)(1:1) == 'P' .or. isymbl(iatm)(1:1) == 'p' ) then
   !      r(iatm) = 1.88d0
   !   else if ( isymbl(iatm)(1:1) == 'S' .or. isymbl(iatm)(1:1) == 's' ) then
   !      r(iatm) = 1.81d0
   !   else
   !      write(6, *) 'SA Bomb in sa_init(): No radius assigned for atom', iatm, isymbl(iatm)
   !      call mexit(6, 1)
   !   end if
   !end do
    
   if ( verbose .and. pbprint ) write(6, *) ' SA surface: setting up working radii'
   nsatm = 0
   do iatm = 1, natom
      if ( r(iatm) == ZERO ) then
         rp(iatm) = ZERO
         rp2(iatm) = ZERO
      else
         rp(iatm) = r(iatm) + sprob
         rp2(iatm) = rp(iatm)**2
         nsatm = nsatm + 1
      endif
   end do
   if ( verbose .and. pbprint ) write(6, *) ' SA surface: found nonzero radii', nsatm

end subroutine sa_init
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ VDW radi setup for Pymol surface potential using Pymol vdw radii
subroutine phi_aaradi( natom, isymbl, radi )

   ! Passed variables

   integer natom
   character (len=4) :: isymbl(natom)
   _REAL_ radi(natom)

   ! Local varialbes

   integer iatm

   do iatm = 1, natom
      if      (isymbl(iatm)(1:2) == 'Br' .or. isymbl(iatm)(1:2) == 'br') then
         radi(iatm) = 1.85d0
      else if ( isymbl(iatm)(1:2) == 'C0' ) then
         radi(iatm) = 1.80d0
      else if (isymbl(iatm)(1:2) == 'Cl' .or. isymbl(iatm)(1:2) == 'cl') then
         radi(iatm) = 1.75d0
      else if ( isymbl(iatm)(1:2) == 'CU' ) then
         radi(iatm) = 1.40d0
      else if ( isymbl(iatm)(1:2) == 'Cs' ) then
         radi(iatm) = 1.80d0
      else if ( isymbl(iatm)(1:1) == 'C' .or. isymbl(iatm)(1:1) == 'c' ) then
         radi(iatm) = 1.70d0
      else if ( isymbl(iatm)(1:2) == 'FE' ) then
         radi(iatm) = 1.80d0
      else if ( isymbl(iatm)(1:1) == 'F' .or. isymbl(iatm)(1:1) == 'f' ) then
         radi(iatm) = 1.47d0
      else if ( isymbl(iatm)(1:1) == 'H' .or. isymbl(iatm)(1:1) == 'h' ) then
         radi(iatm) = 1.20d0
      else if ( isymbl(iatm)(1:2) == 'IB' ) then
         radi(iatm) = 2.27d0
      else if ( isymbl(iatm)(1:2) == 'IM' ) then
         radi(iatm) = 1.75d0
      else if ( isymbl(iatm)(1:2) == 'IP' ) then
         radi(iatm) = 2.27d0
      else if ( isymbl(iatm)(1:1) == 'I' .or. isymbl(iatm)(1:1) == 'i' ) then
         radi(iatm) = 1.98d0
      else if ( isymbl(iatm)(1:2) == 'K ' ) then
         radi(iatm) = 2.75d0
      else if ( isymbl(iatm)(1:2) == 'Li' ) then
         radi(iatm) = 1.80d0
      else if ( isymbl(iatm)(1:2) == 'LP' ) then
         radi(iatm) = 1.20d0
      else if ( isymbl(iatm)(1:2) == 'Mg' ) then
         radi(iatm) = 1.73d0
      else if ( isymbl(iatm)(1:2) == 'Mn' ) then
         radi(iatm) = 1.73d0
      else if ( isymbl(iatm)(1:2) == 'Na' ) then
         radi(iatm) = 2.27d0
      else if ( isymbl(iatm)(1:1) == 'N' .or. isymbl(iatm)(1:1) == 'n' ) then
         radi(iatm) = 1.55d0
      else if ( isymbl(iatm)(1:1) == 'O' .or. isymbl(iatm)(1:1) == 'o' ) then
         radi(iatm) = 1.52d0
      else if ( isymbl(iatm)(1:1) == 'P' .or. isymbl(iatm)(1:1) == 'p' ) then
         radi(iatm) = 1.80d0
      else if ( isymbl(iatm)(1:2) == 'Rb' ) then
         radi(iatm) = 1.80d0
      else if ( isymbl(iatm)(1:1) == 'S' .or. isymbl(iatm)(1:1) == 's' ) then
         radi(iatm) = 1.80d0
      else if ( isymbl(iatm)(1:2) == 'Zn' ) then
         radi(iatm) = 1.39d0
      else
         write(6, *) 'SA Bomb in phi_aaradi(): No radius assigned for atom', iatm, isymbl(iatm)
         call mexit(6, 1)
      end if
   end do

   radi = radi -0.10d0

end subroutine phi_aaradi
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ PB cavity radi setup based on electrostatic solvation free energies in TIP3P
subroutine pb_aaradi( natom, nbonh, ibh, jbh, radi, acrg, ucrgh, ucrga, resid, igraph, isymbl, rin )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Authors:
   ! Chuck Tan, Luo Research Group, UC-Irvine
   !
   ! The assigned cavity radii were opitimzed from reproducing TIP3P electrostatic
   ! solvation free energies of all Amber standard residues/fragments (42) and single
   ! ions in the database files.
   !
   ! However, small organic molecules generated by Antechamber use parm file radii!
   ! 
   ! Dielectric assignment use a solvent probe of 1.6 A and solvent excluded surface.
   ! fdpb used a grid spacing of 0.5 A with the harmonic smoothing method to
   ! define the dielectric boundary. fdpb energy was averaged over 27 different
   ! grid positions. See Tan, Yang and Luo, In Prep. 2006 for more information.
   ! 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Passed variables

   integer natom, nbonh
   integer ibh(*), jbh(*)
   character (len=4) :: resid(*), igraph(*), isymbl(*)
   _REAL_ radi(*), acrg(*), ucrgh(*), ucrga(*), rin(*)

   ! Local variables
    
   integer nhn(natom)
   integer iatm, jatm, idum
    
   ! accumulate no. of H attached to N

   nhn(1:natom) = 0

   do idum = 1, nbonh

      iatm = ibh(idum)/3 + 1
      jatm = jbh(idum)/3 + 1

      if ( isymbl(iatm)(1:1) == 'N' )  then
         nhn(iatm) = nhn(iatm) + 1
      endif

      if ( isymbl(jatm)(1:1) == 'N' )  then
         nhn(jatm) = nhn(jatm) + 1
      endif

   enddo

   ! ::::: radii are zero by default :::::
   ! this takes care of All zero H's
    
   radi(1:natom) = ZERO
    
   ! ::::: radii for heavy atoms without H :::::

   do iatm = 1, natom

      if     ( isymbl(iatm)(1:1) == 'H' ) cycle

      ! ----- radii of CT without H -----
      if     ( isymbl(iatm)(1:2) == 'CT' ) then
         radi(iatm) = 2.00d0
      ! ----- radii of aromatic C without H -----
      elseif ( isymbl(iatm)(1:2) == 'CA' .or. &
               isymbl(iatm)(1:2) == 'CB' .or. &
               isymbl(iatm)(1:2) == 'CC' .or. &
               isymbl(iatm)(1:2) == 'CK' .or. &
               isymbl(iatm)(1:2) == 'CM' .or. &
               isymbl(iatm)(1:2) == 'CN' .or. &
               isymbl(iatm)(1:2) == 'CQ' .or. &
               isymbl(iatm)(1:2) == 'CR' .or. &
               isymbl(iatm)(1:2) == 'CV' .or. &
               isymbl(iatm)(1:2) == 'CW' .or. &
               isymbl(iatm)(1:2) == 'C*' ) then
         if ( ucrga(iatm) > 0.5d0 ) then
            radi(iatm) = 1.85d0
         else
            radi(iatm) = 1.75d0
         endif
      ! ----- radii of C  -----
      elseif ( isymbl(iatm)(1:2) == 'C ' ) then
         if ( ucrga(iatm) < -0.5d0 ) then
            radi(iatm) = 1.65d0
         else
            radi(iatm) = 2.00d0
         endif

      ! ----- radii of NB, NC, and N* without H -----
      ! two kinds of charge environments
      elseif ( isymbl(iatm)(1:2) == 'NB' .or. &
               isymbl(iatm)(1:2) == 'NC' .or. &
               isymbl(iatm)(1:2) == 'N*' ) then
         if ( ucrga(iatm) > 0.5d0 ) then
            radi(iatm) = 1.80d0
         else
            radi(iatm) = 1.70d0
         endif
      ! ----- radii of N without H -----
      elseif ( isymbl(iatm)(1:2) == 'N ' ) then
         radi(iatm) = 1.70d0

      ! ----- radii of O, O2, OS without H -----
      elseif ( isymbl(iatm)(1:2) == 'O ' ) then
         radi(iatm) = 1.57d0
      elseif ( isymbl(iatm)(1:2) == 'O2' ) then
         radi(iatm) = 1.25d0
      elseif ( isymbl(iatm)(1:2) == 'OS' ) then
         radi(iatm) = 1.30d0

      ! ----- radii of S atoms without H -----
      elseif ( isymbl(iatm)(1:1) == 'S' ) then
         if ( ucrga(iatm) < -0.5d0 ) then
            radi(iatm) = 1.45d0
         else
            radi(iatm) = 2.00d0
         endif

      ! ----- radii of P atoms without H -----
      elseif ( isymbl(iatm)(1:1) == 'P' ) then
         radi(iatm) = 1.95d0

      ! ----- radii of Ions -----
      elseif ( isymbl(iatm)(1:2) == 'Li' ) then
         radi(iatm) = 1.481d0
      elseif ( isymbl(iatm)(1:2) == 'Na' ) then
         radi(iatm) = 1.875d0
      elseif ( isymbl(iatm)(1:2) == 'IP' ) then
         radi(iatm) = 1.875d0
      elseif ( isymbl(iatm)(1:2) == 'K ' ) then
         radi(iatm) = 2.288d0
      elseif ( isymbl(iatm)(1:2) == 'Rb' ) then
         radi(iatm) = 2.448d0
      elseif ( isymbl(iatm)(1:2) == 'Cs' ) then
         radi(iatm) = 2.712d0
      elseif ( isymbl(iatm)(1:2) == 'F ' ) then
         radi(iatm) = 1.223d0
      elseif ( isymbl(iatm)(1:2) == 'Cl' ) then
         radi(iatm) = 1.516d0
      elseif ( isymbl(iatm)(1:2) == 'IM' ) then
         radi(iatm) = 1.815d0
      elseif ( isymbl(iatm)(1:2) == 'Br' ) then
         radi(iatm) = 1.745d0
      elseif ( isymbl(iatm)(1:2) == 'I ' ) then
         radi(iatm) = 1.870d0
      elseif ( isymbl(iatm)(1:2) == 'MG' ) then
         radi(iatm) = 1.515d0
      elseif ( isymbl(iatm)(1:2) == 'C0' ) then
         radi(iatm) = 2.126d0
      elseif ( isymbl(iatm)(1:2) == 'Zn' ) then
         radi(iatm) = 1.469d0

      endif

   end do

   ! ::::: radii for heavy atoms with H :::::

   do iatm = 1, natom

      if ( isymbl(iatm)(1:1) == 'H' ) cycle

      ! ----- radii of CT/HC_n, CT/H1_n, and CT/HP_n -----

      if     ( isymbl(iatm)(1:2) == 'CT' ) then
         if ( abs(ucrgh(iatm) - acrg(iatm)) .gt. 1.0d-4 ) then
            if ( ucrga(iatm) > 0.9d0 ) then      ! LYS
               radi(iatm) = 2.50d0
            elseif ( ucrga(iatm) < -0.9d0 ) then ! CYM
               radi(iatm) = 2.10d0
            else
               radi(iatm) = 2.25d0
            endif
         endif

     ! ----- radii of aromatic C/H -----

      elseif ( isymbl(iatm)(1:2) == 'CA' .or. &
               isymbl(iatm)(1:2) == 'CB' .or. &
               isymbl(iatm)(1:2) == 'CC' .or. &
               isymbl(iatm)(1:2) == 'CK' .or. &
               isymbl(iatm)(1:2) == 'CM' .or. &
               isymbl(iatm)(1:2) == 'CN' .or. &
               isymbl(iatm)(1:2) == 'CQ' .or. &
               isymbl(iatm)(1:2) == 'CR' .or. &
               isymbl(iatm)(1:2) == 'CV' .or. &
               isymbl(iatm)(1:2) == 'CW' .or. &
               isymbl(iatm)(1:2) == 'C*' ) then
         if ( abs(ucrgh(iatm) - acrg(iatm)) .gt. 1.0d-4 ) radi(iatm) = 2.15d0

      ! ----- radii of NA always with H -----

      elseif ( isymbl(iatm)(1:2) == 'NA' ) then
         if ( ucrga(iatm) > 0.5d0 ) then
            radi(iatm) = 2.15d0
         else
            radi(iatm) = 2.10d0
         endif

      ! ----- radii of N, N2, and N3 with H -----

      elseif ( nhn(iatm) /= 0 .and. &
             ( isymbl(iatm)(1:2) == 'N ' .or. &
               isymbl(iatm)(1:2) == 'N2' .or. &
               isymbl(iatm)(1:2) == 'N3' ) ) then

         if ( ucrga(iatm) < 0.00d0 ) then
            radi(iatm) = 1.80d0 ! LYN only
         elseif ( ucrga(iatm) > 0.80d0 ) then
            radi(iatm) = 2.50d0 ! LYS and ARG only
            if     ( isymbl(iatm)(1:2) == 'N2' ) then
               radi(iatm) = 2.56d0
            elseif ( isymbl(iatm)(1:2) == 'N3' ) then
               radi(iatm) = 2.52d0
            endif
         else
            if     ( nhn(iatm) == 2 .or. nhn(iatm) == 3 ) then
               radi(iatm) = 2.10d0
            elseif ( nhn(iatm) == 1 ) then
               radi(iatm) = 1.95d0
            endif
         endif

     ! ----- radii of OH always with H -----
     ! OH on rings, ASH/GLH behave like OH on rings

      elseif ( isymbl(iatm)(1:2) == 'OH' ) then
         radi(iatm) = 1.925d0

      ! ----- radii of OW always with H -----
      ! in TIP3P water, TIP4P/TIP5P not checked

      elseif ( isymbl(iatm)(1:2) == 'OW' ) then
         radi(iatm) = 1.85d0

      ! ----- radii of SH with H -----

      elseif ( isymbl(iatm)(1:2) == 'SH' ) then
         if ( abs(ucrgh(iatm) - acrg(iatm)) .gt. 1.0d-4 ) radi(iatm) = 2.30d0

      endif

   enddo

   ! exceptions:

   do iatm = 1, natom

      if ( isymbl(iatm)(1:1) == 'H' ) cycle

      ! Charged aromatic C/N: HIP only

      if ( resid(iatm)(1:3) == 'HIP' ) then
         if ( isymbl(iatm)(1:2) == 'CW' .or. &
              isymbl(iatm)(1:2) == 'CR' .or. &
              isymbl(iatm)(1:2) == 'NA' ) radi(iatm) = 2.65d0
      endif

      ! this is not exception
      ! OH not on rings, i.e. OH on linear alphatic CT

      if ( isymbl(iatm)(1:2) == 'OH' ) then
         if ( resid(iatm)(1:3) == 'SER' .or. &
              resid(iatm)(1:3) == 'THR' .or. &
              igraph(iatm)(1:2) == 'O5' ) radi(iatm) = 1.79d0
      endif

      ! this is not exception
      ! OS not on rings, they behave like O2

      if ( isymbl(iatm)(1:2) == 'OS' ) then
         if ( igraph(iatm)(1:2) == 'O4' ) radi(iatm) = 1.25d0 
      endif

   enddo

   ! for Antechamber ligands with atom types of lower cases
   ! use radii from the parm file

   do iatm = 1, natom
      if ( isymbl(iatm)(1:1) == 'h' .or. &
           isymbl(iatm)(1:1) == 'c' .or. &
           isymbl(iatm)(1:1) == 'n' .or. &
           isymbl(iatm)(1:1) == 'o' .or. &
           isymbl(iatm)(1:1) == 'p' .or. &
           isymbl(iatm)(1:1) == 's' ) radi(iatm) = rin(iatm)
   end do

   ! safe-guarding, very important!

   do iatm = 1, natom

      if ( isymbl(iatm)(1:1) == 'H' ) cycle

      if ( radi(iatm) == ZERO ) then
         write(6, *) 'PB Bomb in pb_aaradi(): No radius assigned for atom', iatm, igraph(iatm), isymbl(iatm)
         call mexit(6, 1)
      endif

   enddo

end subroutine pb_aaradi
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Initialize solvent exposed dots of a unit sphere
subroutine sa_sphere( maxsph,minsph,scrd )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Generates a set of about maxsph roughly evenly spaced points on the surface of a
   ! unit sphere, for use in generating probe-accessible surface. Number of points
   ! generated is maxsph, approximately equal to requested number. Modified from UHBD
   ! (Comp. Phys. Comm. 91:57-95, 1995) routine sphere() by Michael Gilson.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use constants, only : PI, TWOPI

   ! Passed variables

   integer maxsph
   _REAL_ minsph
   _REAL_ scrd(3, *)

   ! Local variables

   integer ntheta,npsi,npsimax,i,nt,np
   _REAL_ thtstp,theta,stheta,ctheta,psistp,psi,cpsi,spsi

   ! begin code

   ntheta = sqrt(PI*maxsph/FOUR)
   npsimax = TWO*ntheta
   thtstp = PI/ntheta

   i = 1
   do nt = 1, ntheta
      theta = thtstp*nt
      stheta = sin(theta)
      ctheta = cos(theta)
      npsi = nint(stheta*npsimax)
      if (npsi == 0) cycle
      psistp = TWOPI/npsi
      do np = 1, npsi
         psi = np*psistp
         cpsi = cos(psi)
         spsi = sin(psi)
         scrd(1,i) = cpsi * stheta
         scrd(2,i) = spsi * stheta
         scrd(3,i) = ctheta
         i = i + 1
      enddo
   enddo
   maxsph = i - 1
   minsph = 0.05*maxsph

end subroutine sa_sphere
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Driver of solvent accessible surface and arcs
subroutine sa_driver( verbose,pbprint,natom,dosas,ndosas,dosar,ndosar,npbstep,&
                      nsaslag,acrd,iar1pb,iprshrt,nex,iex )
    
   use constants, only: PI, TWOPI, FOURPI
    
   ! Passed variables
    
   logical verbose, pbprint
   integer natom, dosas, ndosas, dosar, ndosar, npbstep, nsaslag
   integer iar1pb(6,0:natom), iprshrt(*), nex(natom), iex(32,natom)
   _REAL_ acrd(3,1:natom)
    
   ! Local variables
    
   integer iatm
   _REAL_ wf0, wf1, exposure, increase, smallsas
    
   ! compute sas every "ndosas" FDPB calculations
    
   if ( ndosas > 0 .and. dosas == ndosas ) then
      dosas = 1
      call sa_srf( verbose, pbprint, natom )
   end if
    
   ! compute sar every "ndosas" FDPB calculations
    
   if ( ndosar > 0 .and. dosar == ndosar ) then
      dosar = 1
      call sa_arc( verbose, pbprint, natom )
   end if
    
   ! compute time-averaged atomic solvent accessible surfaces over
   ! the last nsaslag steps
    
   if ( npbstep <= nsaslag ) then
      do iatm = 1, natom
         if ( radi(iatm) == ZERO ) cycle
         sumnmax(iatm) = sumnmax(iatm) + nmax(iatm)
         avnmax(iatm) = sumnmax(iatm)/REAL(npbstep)
         sumnexp(iatm) = sumnexp(iatm) + nexp(iatm)
         avnexp(iatm) = sumnexp(iatm)/REAL(npbstep)
      end do
   else
      wf0 = ONE/REAL(nsaslag)
      wf1 = ONE - wf0
      do iatm = 1, natom
         if ( radi(iatm) == ZERO ) cycle
         avnmax(iatm) = wf1*avnmax(iatm) + wf0*nmax(iatm)
         avnexp(iatm) = wf1*avnexp(iatm) + wf0*nexp(iatm)
      end do
   end if
    
   ! compute effective dynamic cavity radii based on solvent accessibilities
    
   do iatm = 1, natom
       
      if ( radi(iatm) == ZERO ) cycle
       
      ! compute atomic solvent accessibilities, i.e. exposure below. a
      ! gaussian-based smoothing function is used for small avnmax(), i.e.
      ! for mostly buried atoms.
       
      ! buried by covalent-bonded atoms only, no need to increase radius
       
      if ( avnmax(iatm) <= minsph ) then
         exposure = ONE
       
      ! gradually tune down the large addition from one to zero with a
      ! gaussian smoothing function with sigma = minsph:
      ! exp ( - (x - minsph)^2/(2 minsph^2) )
       
      else if ( avnmax(iatm) <= 5*minsph ) then
         exposure = avnexp(iatm)/avnmax(iatm) + exp( -( (avnmax(iatm)-minsph)/minsph )**2/2 )
       
      ! normal exposed atoms
       
      else
         exposure = avnexp(iatm)/avnmax(iatm)
      end if
       
      ! compute dynamic atomic radii
       
      if ( exposure > expthresh ) then
         increase = ZERO
      else
         if ( expthresh == ZERO ) then
            increase = ZERO
         else
            increase = ONE - exposure/expthresh
         end if
      end if
      radipd(iatm) = radi(iatm) + increase*radinc
       
   end do
    
contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Solvent accessible surface calculation
subroutine sa_srf( verbose,pbprint,natom )
    
   ! Passed variables
    
   logical verbose, pbprint
   integer natom
    
   ! Local variables
    
   integer alloc_err(4)
   integer iatm, jatm, jfirst, jlast, ilast, isph, jsph, jp
   integer nsrf, nsurf
   integer fstsph(natom), lstsph(natom)
   _REAL_ sx, sy, sz, xi, yi, zi, dx, dy, dz, dxij, dyij, dzij, d2
   _REAL_, parameter :: small = 0.001d0
    
   allocate(knockout(  maxsph*natom), stat = alloc_err(1) )
   allocate(  sphcrd(3,maxsph*natom), stat = alloc_err(2) )
   if ( alloc_err( 1)+alloc_err( 2) /= 0 ) then
      write(6, *) 'SA Bomb in sa_srf(): Allocation aborted', alloc_err(1:2)
      call mexit(6, 1)
   end if
    
   ! generate solvent accessible points using excluded atoms only
   ! then save these points per atoms for later.
   ! this should be the maximum sas area of this molecule

   nsrf = 0
   do iatm = 1, natom
      if ( radip(iatm) == ZERO ) cycle
      fstsph(iatm) = nsrf + 1
      xi = acrd(1,iatm); yi = acrd(2,iatm); zi = acrd(3,iatm)
      jfirst = 1; jlast  = nex(iatm)
        
      ! working on iatm's points
        
      do isph = 1, maxsph
         sx = radip(iatm)*scrd(1,isph) + xi
         sy = radip(iatm)*scrd(2,isph) + yi
         sz = radip(iatm)*scrd(3,isph) + zi
            
         ! loop over exclusion pairs:
            
         do jp = jfirst, jlast
            jatm = iex(jp,iatm)
            if ( radip(jatm) == ZERO )  cycle
            dx = acrd(1,jatm) - sx
            dy = acrd(2,jatm) - sy
            dz = acrd(3,jatm) - sz
            d2 = dx*dx + dy*dy + dz*dz + small
            if (d2 < radip2(jatm)) goto 10
         end do

         ! if sphere point is not "knocked out" by all its neighbors
         ! update global srf point counter

         nsrf = nsrf + 1

         ! save its coord, wrt the ctr of iatm

         sphcrd(1,nsrf) = sx - xi
         sphcrd(2,nsrf) = sy - yi
         sphcrd(3,nsrf) = sz - zi

10    continue
      end do

      ! remember the last srf point of the atom in global srf list

      lstsph(iatm) = nsrf
   end do  ! iatm = 1, natom
 
   knockout(1:nsrf) = .false.

   allocate( fstsdot(         natom), stat = alloc_err(1) )
   allocate( lstsdot(         natom), stat = alloc_err(2) )
   allocate(  srfcrd(3,       nsrf ), stat = alloc_err(3) )
   if ( alloc_err( 1)+alloc_err( 2)+alloc_err( 3) /= 0 ) then
      write(6, *) 'SA Bomb in sa_srf(): Allocation aborted', alloc_err(1:3)
      call mexit(6, 1)
   end if

   ! loop over nblist with cutsa2

   ilast = natom - 1
   do iatm = 1, ilast
      if ( radip(iatm) == ZERO ) cycle
      xi = acrd(1,iatm); yi = acrd(2,iatm); zi = acrd(3,iatm)
      jfirst = iar1pb(1,iatm) + 1; jlast  = iar1pb(3,iatm)
      do jp = jfirst, jlast
         jatm = iprshrt(jp)
         if ( radip(jatm) == ZERO ) cycle

         dxij = xi - acrd(1,jatm)
         dyij = yi - acrd(2,jatm)
         dzij = zi - acrd(3,jatm)
         d2 = dxij**2 + dyij**2 + dzij**2
         if (d2 > (radip(iatm) + radip(jatm))**2) cycle
          
         ! working on iatm's points
          
         do isph = fstsph(iatm), lstsph(iatm)
            if ( knockout(isph) ) cycle
            dx = sphcrd(1,isph) + dxij
            dy = sphcrd(2,isph) + dyij
            dz = sphcrd(3,isph) + dzij
            d2 = dx*dx + dy*dy + dz*dz + small
            if ( d2 < radip2(jatm) ) knockout(isph) = .true.
         enddo
          
         ! working on jatm's points
          
         do jsph = fstsph(jatm), lstsph(jatm)
            if ( knockout(jsph) ) cycle
            dx = sphcrd(1,jsph) - dxij
            dy = sphcrd(2,jsph) - dyij
            dz = sphcrd(3,jsph) - dzij
            d2 = dx*dx + dy*dy + dz*dz + small
            if ( d2 < radip2(iatm) ) knockout(jsph) = .true.
         enddo
          
      enddo  ! jatm = iprshrt(ip), jp = jfirst, jlast
       
   enddo  ! iatm = 1, ilast
    
   ! determine the exposed atoms and buried atoms
   ! save surface dots for nonpolar solvation calculations
    
   lstsdot = 0
   nsrfdot = 0 ! sa dots for the whole molecule
   do iatm = 1, natom
      fstsdot(iatm) = nsrfdot+1
      if ( radip(iatm) == ZERO ) cycle

      nsurf = 0 ! sa dots for this atom
      do isph = fstsph(iatm), lstsph(iatm)
         if ( knockout(isph) ) cycle
         nsurf = nsurf + 1
         nsrfdot = nsrfdot + 1
         srfcrd(1,nsrfdot) = sphcrd(1,isph)! + acrd(1,iatm) ! this is for display
         srfcrd(2,nsrfdot) = sphcrd(2,isph)! + acrd(2,iatm) ! this is for display
         srfcrd(3,nsrfdot) = sphcrd(3,isph)! + acrd(3,iatm) ! this is for display
      end do
      lstsdot(iatm) = nsrfdot

      ! get maximum no. surface dots of an atom
      ! the actual no. surface dots of an atom

      if ( lstsph(iatm) <= fstsph(iatm) ) then
         nmax(iatm) = ZERO
         nexp(iatm) = REAL(nsurf)
      else
         nmax(iatm) = REAL(lstsph(iatm) - fstsph(iatm) + 1)
         nexp(iatm) = REAL(nsurf)
      end if
   enddo
    
   prtsas = ZERO
   do iatm = 1, natom
      prtsas = prtsas + FOURPI*radip2(iatm)*nexp(iatm)/REAL(maxsph)
   end do
   if ( verbose .and. pbprint ) then
      write(6, *)
      write(6, *) 'Atomic solvent accessible surface area:'
      do iatm = 1, natom
         write(6, *) iatm, FOURPI*radip2(iatm)*nexp(iatm)/REAL(maxsph)
      end do
      write(6,*) 'Number of SA srf points exposed', nsrfdot
   end if
    
   deallocate(knockout, stat = alloc_err(1) )
   deallocate(  sphcrd, stat = alloc_err(2) )
   if ( alloc_err( 1)+alloc_err( 2) /= 0 ) then
      write(6, *) 'SA Bomb in sa_srf(): Deallocation aborted', alloc_err(1:2)
      call mexit(6, 1)
   end if
    
   ! for InsightII/Sybyl display
   !print *, 'Number of SA srf points exposed', nsrfdot
   !open (unit=55, file='sasrf.pdb')
   !do isph = 1, nsrfdot
   !   write (55,'("ATOM  ",x,i4,2x,"O   ALA",3x,i3,4x,3(f8.3))') isph,1,srfcrd(1:3,isph)
   !   write (55, '("TER")')
   !end do
   !close(55)
   !open (unit=55, file='sasrf.dot')
   !write (55, '("DOTS")')
   !do isph = 1, nsrfdot
   !   write (55,'(4(f8.3,2x))') srfcrd(1,isph), srfcrd(2,isph), srfcrd(3,isph), 300.
   !end do
   !close(55)

end subroutine sa_srf
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Solvent accessible arc calculation
subroutine sa_arc( verbose,pbprint,natom )
    
   logical verbose, pbprint
   integer natom
    
   integer alloc_err(8)
   integer ksrf, lsrf, msrf, isph
   integer fstsph(natom), lstsph(natom), fstsph1(natom), lstsph1(natom)

   allocate( fstadot(            natom), stat = alloc_err(1) )
   allocate( lstadot(            natom), stat = alloc_err(2) )
   allocate(  fstarc(            natom), stat = alloc_err(3) ) 
   allocate(  lstarc(            natom), stat = alloc_err(4) )
   allocate(    marc(            natom), stat = alloc_err(5) )
   if ( alloc_err(1)+alloc_err(2)+alloc_err(3)+alloc_err(4)+alloc_err(5) /= 0 ) then
      write(6, *) 'SA Bomb in sa_arc(): Allocation aborted', alloc_err(1:5)
      call mexit(6, 1)
   end if
    
   fstarc = 0; lstarc = 0; marc = 0; fstsph = 0; lstsph = 0 
    
   ! get all possible arcs and arc dots

   allocate(  m2narc(  maxarc   ,natom), stat = alloc_err(1) )
   allocate(  arcatm(2,maxarc   *natom), stat = alloc_err(2) )
   allocate(  savarc(3,maxarc   *natom), stat = alloc_err(3) )
   allocate( savactr(3,maxarc   *natom), stat = alloc_err(4) )
   if ( alloc_err(1)+alloc_err(2)+alloc_err(3)+alloc_err(4) /= 0 ) then
      write(6, *) 'SA Bomb in sa_arc(): Allocation aborted', alloc_err(1:4)
      call mexit(6, 1)
   end if
    
   allocate(knockout(  maxarcdot*natom), stat = alloc_err(1) )
   allocate(  sphcrd(3,maxarcdot*natom), stat = alloc_err(2) )
   allocate( sphcrd1(3,maxarcdot*natom), stat = alloc_err(3) )
   allocate(  spharc(  maxarcdot*natom), stat = alloc_err(4) )
   allocate( spharc1(  maxarcdot*natom), stat = alloc_err(5) )
   if ( alloc_err(1)+alloc_err(2)+alloc_err(3)+alloc_err(4)+alloc_err(5) /= 0 ) then
      write(6, *) 'SA Bomb in sa_arc(): Allocation aborted', alloc_err(1:5)
      call mexit(6, 1)
   end if
    
   call circle( natom, ksrf, arcres, fstsph, lstsph )
    
   ! knock out points overlapped by other atoms
   ! step 1: loop over exclusion list

   call exclud( natom, 4, 1, ksrf, lsrf, fstsph, lstsph, fstsph1, lstsph1, spharc1, spharc, sphcrd, sphcrd1 )
    
   ! step 2: loop over short nblist with cutsa1 of 5 A
    
   call exclud( natom, 1, 2, lsrf, msrf, fstsph1, lstsph1, fstsph, lstsph, spharc, spharc1, sphcrd1, sphcrd )
    
   ! step 3: loop over long nblist with cutsa2 9 A
   ! note that narcdot, fstadot, lstadot, dotarc, arccrd are returned and used elsewhere
    
   allocate(  dotarc(   msrf), stat = alloc_err(1) )
   allocate(  arccrd(3, msrf), stat = alloc_err(2) )
   if ( alloc_err(1)+alloc_err(2) /= 0 ) then
      write(6, *) 'SA Bomb in sa_arc(): Allocation aborted', alloc_err(1:2)
      call mexit(6, 1)
   end if

   call exclud( natom, 2, 3, msrf, narcdot, fstsph, lstsph, fstadot, lstadot, spharc1, dotarc, sphcrd, arccrd )

   deallocate(knockout, stat = alloc_err(1) )
   deallocate(  sphcrd, stat = alloc_err(2) )
   deallocate( sphcrd1, stat = alloc_err(3) )
   deallocate(  spharc, stat = alloc_err(4) )
   deallocate( spharc1, stat = alloc_err(5) )
   if ( alloc_err(1)+alloc_err(2)+alloc_err(3)+alloc_err(4)+alloc_err(5) /= 0 ) then
      write(6, *) 'SA Bomb in sa_arc(): Deallocation aborted', alloc_err(1:5)
      call mexit(6, 1)
   end if
   
   if ( verbose .and. pbprint ) then
      write(6,*) 'Number of SA arcs generated', narc
      write(6,*) 'Number of SA arc points exposed', narcdot, ' with resolution (A)', arcres
   end if

   ! for InsightII/Sybyl display
   !print *, 'Number of SA arcs generated', narc
   !print *, 'Number of SA arc points exposed', narcdot, ' with resolution (A)', arcres
   !open (unit=55, file='saarc.pdb')
   !do isph = 1, narcdot
   !   write (55,'("ATOM  ",x,i4,2x,"O   ALA",3x,i3,4x,3(f8.3))') isph,1,arccrd(1:3,isph)
   !   write (55, '("TER")')
   !enddo
   !close(55)
   !open(unit=55,file='saarc.dot')
   !write (55, '("DOTS")')
   !do isph = 1, narcdot
   !   write(55,'(4(f8.3,2x))') arccrd(1:3,isph), 300.
   !enddo
   !close(55)

end subroutine sa_arc
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Generate maximum solvent-accessible circles
subroutine circle( natom,msrf,cstep,fstsph,lstsph )
    
   ! Passed variables
    
   integer natom, msrf, fstsph(*), lstsph(*)
   _REAL_ cstep
    
   ! Local variables
    
   integer  jp, iatm, ilast, jatm, jfirst, jlast, isph, nspha
   _REAL_ ri, xi, yi, zi, radius, aphi, rj
   _REAL_ x, y, dxij, dyij, dzij, d2, d, rd
   _REAL_ sintheta, costheta, cosaij, cosaji
   _REAL_ arcctr(3), frc_ctr
   _REAL_ psi, cosphi, sinphi
    
   _REAL_, parameter :: small = 0.001d0

   narc = 0; msrf = 0
   ilast = natom - 1
   do iatm = 1, ilast
       
      ri = radip(iatm)
      if ( ri == ZERO ) then
         fstsph(iatm) = lstsph(iatm) + 1
         cycle
      end if
      fstarc(iatm) = narc + 1
      fstsph(iatm) = msrf + 1
      xi = acrd(1,iatm); yi = acrd(2,iatm); zi = acrd(3,iatm)
       
      jfirst = iar1pb(4,iatm-1) + 1; jlast  = iar1pb(3,iatm)
      do jp = jfirst, jlast
          
         jatm = iprshrt(jp)
         rj = radip(jatm)
         if ( rj == ZERO ) cycle
          
         dxij = acrd(1,jatm) - xi
         dyij = acrd(2,jatm) - yi
         dzij = acrd(3,jatm) - zi
         d2 = dxij**2 + dyij**2 + dzij**2
         if ( d2 >= (ri + rj)**2 .or. d2 <= (ri - rj)**2 ) cycle
          
         ! setting up indexes ...
          
         narc = narc + 1
         marc(iatm) = marc(iatm) + 1; marc(jatm) = marc(jatm) + 1
         m2narc(marc(iatm),iatm) = narc
         m2narc(marc(jatm),jatm) = narc
          
         ! getting the center of the arc circle
          
         d = sqrt(d2); rd = ONE/d
         cosaij = (d2 + radip2(iatm) - radip2(jatm))/(TWO*d*ri)
         cosaji = (d2 + radip2(jatm) - radip2(iatm))/(TWO*d*rj)
         arcatm(1,narc) = jatm; arcatm(2,narc) = iatm
         savarc(1,narc) = cosaij; savarc(2,narc) = cosaji; savarc(3,narc) = rd
          
         radius = ri*sqrt(ONE - cosaij**2)
         frc_ctr = ri*cosaij*rd
         arcctr(1) = xi + frc_ctr*dxij
         arcctr(2) = yi + frc_ctr*dyij
         arcctr(3) = zi + frc_ctr*dzij
         savactr(1:3,narc) = arcctr(1:3)
          
         ! determine the orientation of the arc w.r.t. the z-axis
          
         if ( abs(dyij) < small .and. abs(dxij) < small ) then
            ! the vector is aligned with the z-axis
            cosphi = ONE
            sinphi = ZERO

         else if ( abs(dyij) < small ) then
            ! the vector is aligned with the x-axis
            cosphi = ZERO
            sinphi = ONE

            ! otherwise ...
         else if ( dyij > ZERO ) then
            aphi = atan(dxij/dyij)
            cosphi = cos(aphi)
            sinphi = sin(aphi)
         else
            aphi = PI + atan(dxij/dyij)
            cosphi = cos(aphi)
            sinphi = sin(aphi)
         end if
         costheta = dzij/d
         sintheta = sqrt(ONE - costheta**2)
          
         ! generating dots on the arc
          
         nspha = TWOPI*radius/cstep
         psi = ZERO
         do isph = 1, nspha
            msrf = msrf + 1
            x = radius*cos(psi)
            y = radius*sin(psi)
            sphcrd(1,msrf) = arcctr(1) + cosphi*(x) + sinphi*(   costheta*y)
            sphcrd(2,msrf) = arcctr(2) - sinphi*(x) + cosphi*(   costheta*y)
            sphcrd(3,msrf) = arcctr(3) +                                        ( - sintheta*y)
            spharc1(msrf) = narc
            psi = psi + TWOPI/nspha
         end do
          
      end do  ! jatm = iprshrt(ip), jp =  jfirst, jlast
      lstarc(iatm) = narc
      lstsph(iatm) = msrf
       
      if ( narc > maxarc*natom ) then
         write(6,*) 'SA Bomb in circle(): Stored surface arcs over limit', iatm, narc
         call mexit(6,1)
      end if
      if ( marc(iatm) > maxarc ) then
         write(6,*) 'SA Bomb in circle(): Stored surface arcs over limit', iatm, marc(iatm)
         call mexit(6,1)
      end if
      if ( msrf > maxarcdot*natom ) then
         write(6,*) 'SA Bomb in circle(): Stored surface points over limit', msrf
         call mexit(6,1)
      end if
        
   end do  ! iatm = 1, ilast
   fstarc(natom) = narc + 1
   lstarc(natom) = narc + 1
   fstsph(natom) = lstsph(natom) + 1
    
end subroutine circle
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Knock out dots on the circles that are overlapped by close pairs
subroutine exclud( natom,istart,istop,lsrf,msrf,fstsph,lstsph,fstsph1,lstsph1,spharc,spharc1,sphcrd,sphcrd1 )
    
   integer natom, istart, istop, lsrf, msrf
   integer fstsph(*), lstsph(*), fstsph1(*), lstsph1(*), spharc(*), spharc1(*)
   _REAL_ sphcrd(3,*), sphcrd1(3,*)
    
   integer jp, iatm, jatm, ilast, jfirst, jlast, isph, jsph
   _REAL_ xi, yi, zi, xj, yj, zj
   _REAL_ ri2, rj2, dx, dy, dz, d2
    
   _REAL_, parameter :: small = 0.001d0
    
   knockout(1:lsrf) = .false.
    
   ilast = natom - 1
   do iatm = 1, ilast
       
      ri2 = radip2(iatm) 
      if ( ri2 == ZERO ) cycle
      xi = acrd(1,iatm); yi = acrd(2,iatm); zi = acrd(3,iatm)
       
      if ( istart == 4 ) then
         jfirst = iar1pb(istart,iatm-1) + 1
      else
         jfirst = iar1pb(istart,iatm) + 1
      end if
      jlast  = iar1pb(istop,iatm)
      do jp = jfirst, jlast
          
         jatm = iprshrt(jp); rj2 = radip2(jatm)
         if ( rj2 == ZERO ) cycle
         xj = acrd(1,jatm); yj = acrd(2,jatm); zj = acrd(3,jatm)
           
         ! working on iatm's points
          
         do isph = fstsph(iatm), lstsph(iatm)
            if ( knockout(isph) ) cycle
            dx = sphcrd(1,isph) - xj
            dy = sphcrd(2,isph) - yj
            dz = sphcrd(3,isph) - zj
            d2 = dx*dx + dy*dy + dz*dz + small
            if ( d2 < rj2 ) knockout(isph) = .true.
         end do
          
         ! working on jatm's points
           
         do jsph = fstsph(jatm), lstsph(jatm)
            if ( knockout(jsph) ) cycle
            dx = sphcrd(1,jsph) - xi
            dy = sphcrd(2,jsph) - yi
            dz = sphcrd(3,jsph) - zi
            d2 = dx*dx + dy*dy + dz*dz + small
            if ( d2 < ri2 ) knockout(jsph) = .true.
         end do
       
      end do  ! jatm = iprshrt(ip), jp = jfirst, jlast
    
   end do  ! iatm = 1, ilast
    
   ! condense exposed points
    
   msrf = 0
   do iatm = 1, ilast
      if ( radip2(iatm) == ZERO ) then
         fstsph1(iatm) = lstsph1(iatm) + 1
         cycle
      end if
      fstsph1(iatm) = msrf + 1
      do isph = fstsph(iatm), lstsph(iatm)
         if ( knockout(isph) ) cycle
         msrf = msrf + 1
         sphcrd1(1,msrf) = sphcrd(1,isph)
         sphcrd1(2,msrf) = sphcrd(2,isph)
         sphcrd1(3,msrf) = sphcrd(3,isph)
         spharc1(msrf) = spharc(isph)
      end do
      lstsph1(iatm) = msrf
   end do
   fstsph1(natom) = lstsph1(natom) + 1
    
end subroutine exclud

end subroutine sa_driver

end module solvent_accessibility
