!---------------------- constants.h --------------------

!-------------------------------------------------------
! Physical Constants

! This is based on old values of the electric constants.
_REAL_      AMBER_ELECTROSTATIC
parameter ( AMBER_ELECTROSTATIC = 18.2223d0 )

! 1998 value of the Bohr radius, physics.nist.gov/constants.
! a0 is the standard abbreviation, but perhaps that is too cryptic here.
! The previous value was 0.5291771 in Amber 7.
_REAL_      BOHR_RADIUS
parameter ( BOHR_RADIUS = 0.5291772083d0 )


!-------------------------------------------------------
! Numeric Constants

_REAL_      PI
parameter ( PI = 3.1415926535897932384626433832795d0 )

_REAL_      SQRT2
parameter ( SQRT2 = 1.4142135623730950488016887242097d0 )


!-------------------------------------------------------
! Unusual Constants

integer     RETIRED_INPUT_OPTION
parameter ( RETIRED_INPUT_OPTION = -10301 ) ! first 5 digit palindromic prime


!-------------------------------------------------------
! Generic Floating Point Constants

_REAL_      TEN_TO_MINUS2
parameter ( TEN_TO_MINUS2  = 1.0d-2 )
_REAL_      TEN_TO_MINUS3
parameter ( TEN_TO_MINUS3  = 1.0d-3 )
_REAL_      TEN_TO_MINUS4
parameter ( TEN_TO_MINUS4  = 1.0d-4 )
_REAL_      TEN_TO_MINUS5
parameter ( TEN_TO_MINUS5  = 1.0d-5 )
_REAL_      TEN_TO_MINUS6  
parameter ( TEN_TO_MINUS6  = 1.0d-6 )
_REAL_      TEN_TO_MINUS10
parameter ( TEN_TO_MINUS10 = 1.0d-10 )
_REAL_      TEN_TO_PLUS3
parameter ( TEN_TO_PLUS3  = 1.0d+3 )
_REAL_      TENTOPLUS10
parameter ( TENTOPLUS10 = 1.0d+10 )

_REAL_      zero
parameter ( zero   = 0.0d0 )
_REAL_      one
parameter ( one    = 1.0d0 )
_REAL_      two
parameter ( two    = 2.0d0 )
_REAL_      three
parameter ( three  = 3.0d0 )
_REAL_      four
parameter ( four   = 4.0d0 )
_REAL_      five
parameter ( five   = 5.0d0 )
_REAL_      six
parameter ( six    = 6.0d0 )
_REAL_      seven
parameter ( seven  = 7.0d0 )
_REAL_      eight
parameter ( eight  = 8.0d0 )
_REAL_      nine
parameter ( nine   = 9.0d0 )
_REAL_      ten
parameter ( ten    = 10.0d0 )
_REAL_      eleven
parameter ( eleven = 11.0d0 )
_REAL_      twelve
parameter ( twelve = 12.0d0 )
_REAL_      half
parameter ( half        = one/two )
_REAL_      third
parameter ( third       = one/three )
_REAL_      fourth
parameter ( fourth      = one/four )
_REAL_      fifth
parameter ( fifth       = one/five )
_REAL_      sixth
parameter ( sixth       = one/six )
_REAL_      seventh
parameter ( seventh     = one/seven )
_REAL_      eighth
parameter ( eighth      = one/eight )
_REAL_      ninth
parameter ( ninth       = one/nine )
_REAL_      tenth
parameter ( tenth       = one/ten )
_REAL_      eleventh
parameter ( eleventh    = one/eleven )
_REAL_      twelfth
parameter ( twelfth     = one/twelve )

!------------------ End constants.h --------------------
