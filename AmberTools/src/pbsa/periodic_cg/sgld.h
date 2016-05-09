!     Head file for the self-guided Langevin Dynamics simulation  
!
!      variables for SGLD simulation
!
      integer BC_SGLDI,BC_SGLDR,BC_SGLDL
      parameter( BC_SGLDI=5,BC_SGLDR=39,BC_SGLDL=4)
!
! ... integers:
!
!
!    SGLD applying rang
!      ISGSTA      Begining atom index applying SGLD
!      ISGEND      Ending atom index applying SGLD
!
!
!      variables for Isotropic Periodic Sum (IPS) calculation
!
!      IPS         IPS options: 1--for both ele and l-j
!                               2--for ele only
!                               3--for l-j only
!      NNBIPS      Number of nonbonded atom pairs
!      NNBIPST     Provious Number of nonbonded atom pairs
!
      integer ISGSTA,ISGEND,IPS,NNBIPST,NNBIPS
      COMMON/DTSGI/ISGSTA,ISGEND,IPS,NNBIPST,NNBIPS
!
! ... floats:
!
!
!    SGMD VARIABLES
!     TSGAVG  !  Local average time, ps
!     SGAVG0  !  Local average remains
!     SGAVG1  !  Local average factor, SGAVG1=1-SGAVG0
!     SGFT    !  Guiding factor 
!     TEMPSG  !  Guiding temperature, K
!     GAMMAS  !  friction coefficient
!     GAMMAT  !  Guiding temperature constant
!     AVGGG   !  Local Average of guiding effect
!
!    IPS parameters
!    BIPSE*    Electrostatic
!    BIPSVA*   Lennard-Jones repulsion
!    BIPSVC*   Lennard-Jones dispersion
!    RIPS*     Radius of IPS local region 
!    PIPS*0    Self IPS pair energies 
!    PIPS*C    IPS system energie components 
!    EIPSS*C   IPS system energies 
!    VIRIPS*C  IPS system virials 
!    VBOXIPS   IPS local region volume
!

      _REAL_  TSGAVG,SGAVG0,SGAVG1,SGFT,TEMPSG,TEMPSGI, &
       GAMMAS,GAMMAT,AVGGG, &
       AIPSE,BIPSE0,BIPSE1,BIPSE2,BIPSE3, & 
       AIPSVA,BIPSVA0,BIPSVA1,BIPSVA2,BIPSVA3, &
       AIPSVC,BIPSVC0,BIPSVC1,BIPSVC2,BIPSVC3, &
       RIPS,RIPS2,RIPS4,oneRIPS6,oneRIPS12,  &
       PIPSE0,PIPSVA0,PIPSVC0,PIPSEC,PIPSVAC,PIPSVCC,  &
       VBOXIPS,EIPSSNB,EIPSSEL,VIRIPS,ripsinv,rips2inv
      COMMON/DTSGR/TSGAVG,SGAVG0,SGAVG1,SGFT,TEMPSG,TEMPSGI, &
       GAMMAS,GAMMAT,AVGGG, &
       AIPSE,BIPSE0,BIPSE1,BIPSE2,BIPSE3, & 
       AIPSVA,BIPSVA0,BIPSVA1,BIPSVA2,BIPSVA3, &
       AIPSVC,BIPSVC0,BIPSVC1,BIPSVC2,BIPSVC3, &
       RIPS,RIPS2,RIPS4,oneRIPS6,oneRIPS12,  &
       PIPSE0,PIPSVA0,PIPSVC0,PIPSEC,PIPSVAC,PIPSVCC,  &
       VBOXIPS,EIPSSNB,EIPSSEL,VIRIPS,ripsinv,rips2inv

       
!
! ... flags:
!
!
!     TSGLD   ! Perform SGLD
!     TLANGV  ! The simulation is a Langevin dynamics simulation
!
!     TEIPS   ! Perform IPS for electrostatic interaction
!     TVIPS   ! Perform IPS for Lennard-Jones interaction
!
	LOGICAL TSGLD,TLANGV,TEIPS,TVIPS
      COMMON/DTSGL/TSGLD,TLANGV,TEIPS,TVIPS

