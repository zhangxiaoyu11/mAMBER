#include "copyright.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine locmem here]
subroutine locmem
   
   
   !     locmem:  partitions core array into storage for all
   !        the major arrays of the program.
   
#  include "dprec.h"
   implicit none
   
#  include "box.h"
#  include "memory.h"
#  include "md.h"
#  include "dynph.h"
#ifdef MPI
#  include "parallel.h"
#endif
   integer none,ntbond,ntangl,ntdih,m7,istartr, &
         istarti, iendr,iendi, istomp,ida_max
   integer r_ptr,i_ptr,h_ptr,maxpr
   integer ncpp
   _REAL_ maxpr_float,natom_float,n2_float
   
   ! dummy indices - those unused in pbsa and whose .h files have been removed
   integer lnmr01,inmr02
   
   !     --- Identification of REAL arrays ---
   
   !                  L05   ! polarization
   !                  L10   ! polarization
   !     CG      ...  L15   ! PARTIAL CHARGES FOR ATOMS
   !     AMASS   ...  LWINV ! ATOMIC MASSES (inverted in rdparm - see Lmass)
   !     XCHRG   ...  Lpol  ! atomic polarizibilities
   !     C       ...  LCRD  ! COORDINATES
   !     F       ...  Lforce! FORCE
   !     V       ...  Lvel  ! VELOCITY for MD, work space for min
   !     VOLD    ...  Lvel2 ! OLD VELOCITY for MD
   !     XR      ...  L45   ! Coords rel. to COM of each molecule
   !     CONP    ...  L50   ! BOND PARAMETER FOR SHAKE
   !     XC      ...  LCRD  ! POSITION COORDINATE FOR CONSTRAINT
   !     WEIT    ...  L60   ! WEIGHT FOR POSITION CONSTRAINT
   !                  L65   ! polarization
   !                  Lmass ! masses
   !     TMA     ...  L75   ! SUB-MOLECULAR WEIGHT ARRAY IN RUNMD
   !                  L95   ! 3*Natom Real Scratch (for pol.) or Natom (for nmr)
   !                        ! also used for SKIP array in shake (2*ntbond)
   !                  L96   ! GB "fs" array
   !                  L97   ! GB "rborn" array
   !                  L98   ! GB "reff" array
   !                 Lfrctmp! 3*Natom + 40 Real Scratch( fdist + pol.), mpi only
   !                  L105  ! NMR "xstore" variable
   !                  L110  ! NMR "fnoe" variable
   !                  L115  ! NMR "ddep" variable
   !                  L120  ! NMR "dddep" variable
   !                  L125  ! NMR "dorat" variable
   !                  L130  ! NMR "ddrat" variable
   !                  L135  ! NMR "rate" variable
   !                  L140  ! NMR "trp" variable
   !                  L145  ! NMR "dint" variable
   !                  L165  ! GB/SA TDND "vdwrad" array
   !                  L170  ! GB/SA LCPO "P1" array
   !                  L175  ! GB/SA LCPO "P2" array
   !                  L180  ! GB/SA LCPO "P3" array
   !                  L185  ! GB/SA LCPO "P4" array
   !                  L186  ! GB max of rborn array
   !                  L187  ! GB min of rborn array
   !                  L188  ! GB ave of rborn array
   !                  L189  ! GB rms fluct of rborn array
   !                  L190  ! TI dcharge array
   !                  Lcpcrg! Constant pHstate charges
   !                  Lcpene! Constant pHstate energies
   
   !     --- Identification of Hollerith arrays
   
   !     LBRES   ...  m02  ! RESIDUE LABEL (nres)
   !     IGRAPH  ...  m04  ! ATOM NAMES (natom)
   !     ISYMBL  ...  m06  ! ATOM SYMBOL ARRAY (natom)
   !     ITREE   ...  m08  ! ATOM TREE STRUCTURE ARRAY (natom)
   !     n14     ...  m12  ! 1*natom
   !     ni14    ...  m14  ! 15*natoms 1-4 indicies
   !     iarx    ...  m16  ! 1*natom  scratch for nonbond+1-4
   
   
   !     --- Identification of Integer arrays ---
   
   !     mark    ...  i01   ! principle (marker) atom for imaging
   !     IPRES   ...  I02      ITA     ...  I32
   !     IAC     ...  I04      JTA     ...  I34
   !     NO      ...  I06      KTA     ...  I36
   !     IBLO    ...  I08      ICTA    ...  I38
   !     INB     ...  I10      IPH     ...  I40
   !     IBH     ...  Iibh     JPH     ...  I42
   !     JBH     ...  Ijbh     KPH     ...  I44
   !     ICBH    ...  Iicbh    LPH     ...  I46
   !     IBA     ...  Iiba     ICPH    ...  I48
   !     JBA     ...  Ijba     IPA     ...  I50
   !     ICBA    ...  Iicba    JPA     ...  I52
   !     ITH     ...  I24      KPA     ...  I54
   !     JTH     ...  I26      LPA     ...  I56
   !     KTH     ...  I28      ICPA    ...  I58
   !     ICTH    ...  I30

   !     Cnstr IGROUP  ...  Icnstrgp
   
   !     Belly IGROUP  ...  Ibellygp
   !     JOIN    ...  I64
   !     ISTORE  ...  I65
   !     IROTAT  ...  I66 (IGRES if belly)
   !     NSP     ...  I70 ! SUBMOLECULE INDEX ARRAY
   !     IAR1    ...  I78
   !     NUMBOND ...  I80
   !     NUMNEAR ...  I82
   !     IAC_PERT...  I84  perturbed atom types, when icfe>0
   !     Constant pH state metadata ... Icpstinf
   !     Constant pH residue states ... Icpresst
   !     Constant pH state protonation levels .. Icpptcnt
   !     Dynamic protonation num titrating residues ... Icptrsct

   
   ! The IVMxx pointers are used for the following:
   
   !           Iifstwt  ...  IFSTWT (fast shake bond array)
   !           Iifstwr  .... IFSTWR (fast shake residue array)
   
   !    NMR restraints/weight changes require two storage areas. These are:
   
   !           WORKN   ...  LNMR01 (X array)
   !           IWORKN  ...  INMR02 (IX array)
   
   !     --- assign standard partition lengths ---
   
   maxdup = 2000
   none = 0
   ntbond = nbonh  + nbona + nbper
   ntangl = ntheth + ntheta + ngper
   ntdih  = nphih  + nphia + ndper + 2*maxdup
   m7     = nphih  + maxdup
   
   !-----------------------------------------------------------------------
   !     --- set pointers for real arrays ---
   !-----------------------------------------------------------------------
   
   ncpp=0
   r_ptr = 1
   call adj_mem_ptr( r_ptr, l05, 0 )
   call adj_mem_ptr( r_ptr, l10, 0 )
   call adj_mem_ptr( r_ptr, l15, natom )
   call adj_mem_ptr( r_ptr, lwinv, natom+ncpp )
   call adj_mem_ptr( r_ptr, lpol, 0 )
   lpolp = lpol
   call adj_mem_ptr( r_ptr, lcrd, 3*natom + 0+3*ncpp )
   call adj_mem_ptr( r_ptr, lforce, 3*natom + 0 + 40+3*ncpp )
   if (imin == 0) then
      call adj_mem_ptr( r_ptr, lvel,  3*natom + 0 )
      call adj_mem_ptr( r_ptr, lvel2, 3*natom + 0 )
   else
      call adj_mem_ptr( r_ptr, lvel, 6*(3*natom + 0) )
      call adj_mem_ptr( r_ptr, lvel2, 0 )
   end if
   call adj_mem_ptr( r_ptr, l45, 3*natom + 0+3*ncpp )
   call adj_mem_ptr( r_ptr, l50, ntbond+ncpp )
   
   ! positional restraints or carlos added targeted MD
   
   call adj_mem_ptr( r_ptr, lcrdr, 0 )
   call adj_mem_ptr( r_ptr, l60, 0 )
   call adj_mem_ptr( r_ptr, l65, 0 )
   
   !     --- real array NMR restraints/weight changes:
   
   call adj_mem_ptr( r_ptr, lmass, natom )
   call adj_mem_ptr( r_ptr, lnmr01, 0)
   
   call adj_mem_ptr( r_ptr, l75, natom )
   call adj_mem_ptr( r_ptr, l95, 2*ntbond )
   if( igb /= 0 ) then
      call adj_mem_ptr( r_ptr, l96, natom )
      call adj_mem_ptr( r_ptr, l97, natom )
      call adj_mem_ptr( r_ptr, l98, natom )
      if ( gbsa > 0 ) then
         call adj_mem_ptr( r_ptr, l165, natom )
         call adj_mem_ptr( r_ptr, l170, natom )
         call adj_mem_ptr( r_ptr, l175, natom )
         call adj_mem_ptr( r_ptr, l180, natom )
         call adj_mem_ptr( r_ptr, l185, natom )
      else
         call adj_mem_ptr( r_ptr, l165, 0 )
         call adj_mem_ptr( r_ptr, l170, 0 )
         call adj_mem_ptr( r_ptr, l175, 0 )
         call adj_mem_ptr( r_ptr, l180, 0 )
         call adj_mem_ptr( r_ptr, l185, 0 )
      end if
      if ( rbornstat == 1 ) then
         call adj_mem_ptr( r_ptr, l186, natom )
         call adj_mem_ptr( r_ptr, l187, natom )
         call adj_mem_ptr( r_ptr, l188, natom )
         call adj_mem_ptr( r_ptr, l189, natom )
      else
         call adj_mem_ptr( r_ptr, l186, 0 )
         call adj_mem_ptr( r_ptr, l187, 0 )
         call adj_mem_ptr( r_ptr, l188, 0 )
         call adj_mem_ptr( r_ptr, l189, 0 )
      end if
   else
      call adj_mem_ptr( r_ptr, l96,  0 )
      call adj_mem_ptr( r_ptr, l97,  0 )
      call adj_mem_ptr( r_ptr, l98,  0 )
      call adj_mem_ptr( r_ptr, l165,  0 )
      call adj_mem_ptr( r_ptr, l170,  0 )
      call adj_mem_ptr( r_ptr, l175,  0 )
      call adj_mem_ptr( r_ptr, l180,  0 )
      call adj_mem_ptr( r_ptr, l185,  0 )
      call adj_mem_ptr( r_ptr, l186,  0 )
      call adj_mem_ptr( r_ptr, l187,  0 )
      call adj_mem_ptr( r_ptr, l188,  0 )
      call adj_mem_ptr( r_ptr, l189,  0 )
   end if  ! ( igb /= 0 )

   if ( icfe /= 0 ) then
      call adj_mem_ptr( r_ptr, l190, natom )
   else
      call adj_mem_ptr( r_ptr, l190, 0 )
   end if
#ifdef MPI
   call adj_mem_ptr( r_ptr, lfrctmp, 3*natom + 40 )
#endif

   call adj_mem_ptr( r_ptr, l110, 0 )
   call adj_mem_ptr( r_ptr, l150, 0)

   if ( icnstph /= 0 .and. icfe == 0) then
      call adj_mem_ptr( r_ptr, l190, natom )
      call adj_mem_ptr( r_ptr, lcpcrg, ATOM_CHRG_C)
      call adj_mem_ptr( r_ptr, lcpene, TITR_STATES_C)
   else
      call adj_mem_ptr( r_ptr, lcpcrg, 0)
      call adj_mem_ptr( r_ptr, lcpene, 0)
   end if

   lastr = r_ptr
   
   !-----------------------------------------------------------------------
   !     --- Allocate Hollerith Space ---
   !-----------------------------------------------------------------------
   
   h_ptr = 1
   call adj_mem_ptr( h_ptr, m02, nres + 1 )
   call adj_mem_ptr( h_ptr, m04, natom )
   call adj_mem_ptr( h_ptr, m06, natom+ncpp )
   call adj_mem_ptr( h_ptr, m08, natom+ncpp )
   call adj_mem_ptr( h_ptr, m12, natom+ncpp )
   call adj_mem_ptr( h_ptr, m16, 2*natom+2*ncpp )
   lasth = h_ptr
   
   !-----------------------------------------------------------------------
   !     --- Static Integer Arrays ---
   !-----------------------------------------------------------------------
   
   i_ptr = 1
   call adj_mem_ptr( i_ptr, i01, natom+ncpp )
   call adj_mem_ptr( i_ptr, i02, nres + 1 )
   call adj_mem_ptr( i_ptr, i04, natom+ncpp )
   call adj_mem_ptr( i_ptr, i06, ntypes*ntypes )
   call adj_mem_ptr( i_ptr, i08, natom+ncpp )
   call adj_mem_ptr( i_ptr, i10, 2*nnb )
   iibh = i_ptr
   
   !     ----- BOND ARRAYS -----
   
   ijbh  = iibh  + ntbond
   iicbh = ijbh  + ntbond
   iiba  = iibh  + nbonh
   ijba  = ijbh  + nbonh
   iicba = iicbh + nbonh
   i24   = iicbh  + ntbond + nbper
   
   !     ----- ANGLE ARRAYS -----
   
   i26 = i24 + ntangl
   i28 = i26 + ntangl
   i30 = i28 + ntangl
   i32 = i24 + ntheth
   i34 = i26 + ntheth
   i36 = i28 + ntheth
   i38 = i30 + ntheth
   i40 = i30 + ntangl + ngper
   
   !     ----- DIHEDRAL ARRAYS -----
   
   i42 = i40 + ntdih
   i44 = i42 + ntdih
   i46 = i44 + ntdih
   i48 = i46 + ntdih
   i50 = i40 + m7
   i52 = i42 + m7
   i54 = i44 + m7
   i56 = i46 + m7
   i58 = i48 + m7
   icnstrgp = i48 + ntdih + ndper
   
   ibellygp = icnstrgp + natom + ncpp
   i64 = ibellygp + natom + ncpp
   i_ptr = i64 + natom + ncpp
   call adj_mem_ptr( i_ptr, i65, 0 )
   call adj_mem_ptr( i_ptr, i66, natom )
   call adj_mem_ptr( i_ptr, i70, natom + 1 )
   
   ! Allocate memory for NMR restraints/weight changes:
   
   call adj_mem_ptr( i_ptr, inmr02, 0 )
   
   ! Allocate the IVMxx array:
   
   call adj_mem_ptr( i_ptr, iifstwt, ntbond )
   call adj_mem_ptr( i_ptr, iifstwr, nres + 1 )
   
   !  --- right now, i78 (iar1) is not being created, needs no space:
   
   call adj_mem_ptr( i_ptr, i78, 0 )
   
   !  --- allocate array for numbond for surface area calculation:
   
   if( gbsa == 1 )then
      call adj_mem_ptr( i_ptr, i80, natom )
      call adj_mem_ptr( i_ptr, i82, 40*natom )
   else if( gbsa == 2 )then
      call adj_mem_ptr( i_ptr, i80, natom )
      call adj_mem_ptr( i_ptr, i82, 80*natom )
   else
   !  call adj_mem_ptr( i_ptr, i80, 0 )
      call adj_mem_ptr( i_ptr, i82, 0 )
   end if
   
   !  --- perturbed atom types, when icfe>0
   
   if ( icfe /= 0 ) then
      call adj_mem_ptr( i_ptr, i84, natom )
   else
      call adj_mem_ptr( i_ptr, i84, 0 )
   end if
   
   if(igb /= 0) then
      call adj_mem_ptr( i_ptr, i86, natom )
   else
      call adj_mem_ptr( i_ptr, i86, 0 )
   end if

   if (icnstph /= 0) then
      call adj_mem_ptr( i_ptr, icpstinf, TITR_RES_C*STATEINF_FLD_C)
      call adj_mem_ptr( i_ptr, icpresst, TITR_RES_C)
      call adj_mem_ptr( i_ptr, icpptcnt, TITR_STATES_C)
      call adj_mem_ptr( i_ptr, icptrsct, 1)
   else
      call adj_mem_ptr( i_ptr, icpstinf, 0)
      call adj_mem_ptr( i_ptr, icpresst, 0)
      call adj_mem_ptr( i_ptr, icpptcnt, 0)
      call adj_mem_ptr( i_ptr, icptrsct, 0)
   end if
   lasti = i_ptr
   
   !     --- crude (but useful?) estimate for MAXPR:
   
   if( igb /= 0 ) then
      maxpr = 1
   end if

   lastpr = maxpr
   
   return
end subroutine locmem 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine adj_mem_ptr here]
subroutine adj_mem_ptr(mem_ptr,assign_ptr,size)
   implicit none
   integer mem_ptr,assign_ptr,size

   assign_ptr = mem_ptr
   mem_ptr = mem_ptr + size
   return
end subroutine adj_mem_ptr 
