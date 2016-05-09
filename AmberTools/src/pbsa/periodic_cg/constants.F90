! <compile=optimized>
#include "copyright.h"
#include "dprec.h"

!+++++++++++++++++++++++++++++++++++++++
!This module contains various parameters
!and constants used by the different 
!routines that make up sander.
!
!If you want to use one of the constants
!in your routine you should include the
!line:
!
!use constants, only : xxx, yyy, zzz
!
!where xxx,yyy,zzz are the constants you plan
!to use in your routine.
!This line needs to go before the
!implicit none declaration.
!
! Written by: Ross Walker (TSRI, 2005)
!++++++++++++++++++++++++++++++++++++++++

module constants

!------------------------------------------------------------
! Generic Floating Point Constants
_REAL_, parameter :: TEN_TO_MINUS2  = 1.0d-2
_REAL_, parameter :: TEN_TO_MINUS3  = 1.0d-3
_REAL_, parameter :: TEN_TO_MINUS4  = 1.0d-4
_REAL_, parameter :: TEN_TO_MINUS5  = 1.0d-5
_REAL_, parameter :: TEN_TO_MINUS6  = 1.0d-6
_REAL_, parameter :: TEN_TO_MINUS10 = 1.0d-10
_REAL_, parameter :: TEN_TO_PLUS3   = 1.0d+3
_REAL_, parameter :: TEN_TO_PLUS10  = 1.0d+10
                                                                                                                                                             
_REAL_, parameter :: zero      = 0.0d0
_REAL_, parameter :: one       = 1.0d0
_REAL_, parameter :: two       = 2.0d0
_REAL_, parameter :: three     = 3.0d0
_REAL_, parameter :: four      = 4.0d0
_REAL_, parameter :: five      = 5.0d0
_REAL_, parameter :: six       = 6.0d0
_REAL_, parameter :: seven     = 7.0d0
_REAL_, parameter :: eight     = 8.0d0
_REAL_, parameter :: nine      = 9.0d0
_REAL_, parameter :: ten       = 10.0d0
_REAL_, parameter :: eleven    = 11.0d0
_REAL_, parameter :: twelve    = 12.0d0
_REAL_, parameter :: sixteen   = 16.0d0
_REAL_, parameter :: thirtytwo = 32.0d0
_REAL_, parameter :: sixtyfour = 64.0d0
                                                                                                                                                             
_REAL_, parameter :: half         = one/two
_REAL_, parameter :: third        = one/three
_REAL_, parameter :: fourth       = one/four
_REAL_, parameter :: fifth        = one/five
_REAL_, parameter :: sixth        = one/six
_REAL_, parameter :: seventh      = one/seven
_REAL_, parameter :: eighth       = one/eight
_REAL_, parameter :: ninth        = one/nine
_REAL_, parameter :: tenth        = one/ten
_REAL_, parameter :: eleventh     = one/eleven
_REAL_, parameter :: twelfth      = one/twelve
_REAL_, parameter :: sixteenth    = one/sixteen
_REAL_, parameter :: thirtysecond = one/thirtytwo
_REAL_, parameter :: sixtyfourth  = one/sixtyfour

_REAL_, parameter :: thirtieth    = one/30.0d0
                                                                                                                                                             
!------------------------------------------------------------


!------------------------------------------------------------
! Physical Constants
_REAL_, parameter :: AMBER_ELECTROSTATIC = 18.2223d0
_REAL_, parameter :: AMBER_ELECTROSTATIC2 = AMBER_ELECTROSTATIC * AMBER_ELECTROSTATIC
!Ratio by which to scale amber charges to get electron charges - amberchg * oneqscale = electron charges
! = 1.0 / 18.2223d0
_REAL_, parameter :: INV_AMBER_ELECTROSTATIC = 1.0d0/AMBER_ELECTROSTATIC
_REAL_, parameter :: INV_AMBER_ELECTROSTATIC2 = 1.0d0/AMBER_ELECTROSTATIC2


_REAL_, parameter :: BOHRS_TO_A = 0.529177249D0   ! Bohrs * this = angstroms - Same constants as used in dynamo v2.
!_REAL_, parameter :: BOHRS_TO_A = 0.52917706D0   ! Bohrs * this = angstroms - Same constants as used in Gaussian 98
!_REAL_, parameter :: BOHRS_TO_A = 0.529177D0     ! Bohrs * this = angstroms - Same constants as used in Mopac6 hcore.f
!_REAL_, parameter :: BOHRS_TO_A = 0.529167D0     !                            as used in Mopac6 repp.f
_REAL_, parameter :: A_TO_BOHRS = 1.0d0 / BOHRS_TO_A
!_REAL_, parameter :: A_TO_BOHRS = 1.88976D0      !Same constants as used in Mopac6 gover.f
_REAL_, parameter :: A2_TO_BOHRS2 = A_TO_BOHRS * A_TO_BOHRS !Mopac6 uses 3.5711928576D0 in gover.f for this.
_REAL_, parameter :: A3_TO_BOHRS3 = A2_TO_BOHRS2 * A_TO_BOHRS
_REAL_, parameter :: A4_TO_BOHRS4 = A2_TO_BOHRS2 * A2_TO_BOHRS2
                                                                                                                                                             
_REAL_, parameter :: AU_TO_EV = 27.21d0 !Conversion from AU to EV - Same as dynamo v2 uses and Gaussian 98
                                        !Note (RCW+MC): more precise would be: 1 a.u. 27.211396 eV
                                        !Mopac6 uses 27.21D0 in calpar.f, delri.f and repp.f but in
                                        !ffhpol.f it uses 27.2107 and in the manual it quotes 27.211
_REAL_, parameter :: HALF_AU_TO_EV = AU_TO_EV * half
_REAL_, parameter :: FOURTH_AU_TO_EV = AU_TO_EV * fourth
_REAL_, parameter :: EIGHTH_AU_TO_EV = AU_TO_EV * eighth
_REAL_, parameter :: SXNTH_AU_TO_EV = EIGHTH_AU_TO_EV*half
_REAL_, parameter :: A2_TO_BOHRS2xAU_TO_EV = A2_TO_BOHRS2*AU_TO_EV

!_REAL_, parameter :: EV_TO_KCAL = 23.060362D0      !Conversion from EV to KCAL/MOL
!Dynamo parameter
_REAL_, parameter :: EV_TO_KCAL = 23.061d0  !Dynamo's conversion
                                            !Mopac6 uses 23.061 in ffhpol.f analyt.f compfg.f datin.f dcart.f
                                            !                      delri1.f delri2.f deritr.f interp.f iter.f
                                            !                      moldat.f mopac.f
_REAL_, parameter :: KCAL_TO_EV = one / EV_TO_KCAL
                                                                                                                                                             
_REAL_, parameter :: AU_TO_KCAL = AU_TO_EV*EV_TO_KCAL

!------------------------------------------------------------
!Numeric Constants
_REAL_, parameter :: PI      = 3.1415926535897932384626433832795d0
_REAL_, parameter :: PI2     = PI*PI
_REAL_, parameter :: HALFPI = PI * 0.5d0
_REAL_, parameter :: TWOPI  = 2.0d0 * PI
_REAL_, parameter :: FOURPI = 4.0d0 * PI
_REAL_, parameter :: INVPI  = 1.0d0 / PI
_REAL_, parameter :: SQRTPI = 1.77245385090551602729816748334d0 !sqrt(PI)
_REAL_, parameter :: INVSQRTPI = 1.0d0 / SQRTPI
_REAL_, parameter :: DEG_TO_RAD = PI / 180.0d0
_REAL_, parameter :: RAD_TO_DEG = 180.0d0 / PI 

_REAL_, parameter :: SQRT2     = 1.4142135623730950488016887242097d0
_REAL_, parameter :: INVSQRT2  = 1.0d0 / SQRT2

!------------------------------------------------------------
!Generalised Born Constants
_REAL_, parameter :: alpb_alpha = 0.571412d0 !Alpha prefactor for alpb_alpha

!------------------------------------------------------------
! Unusual Constants
integer, parameter :: RETIRED_INPUT_OPTION = -10301 ! first 5 digit palindromic prime

end module constants

