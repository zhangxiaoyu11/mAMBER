#include "copyright.h"
#include "dprec.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Open input files and read cntrl namelist.
subroutine mdread1()
   
   implicit none
#  include "box.h"
#  include "constants.h"
#  include "def_time.h"
#  include "extra.h"
#  include "files.h"
#  include "md.h"
#  include "memory.h"
#  include "parms.h"
   _REAL_      temp0les
   character(len=4) watdef(4),watnam,owtnm,hwtnm1,hwtnm2

   _REAL_      at1
   _REAL_      at2
   _REAL_      beta3
   _REAL_      dele
   _REAL_      gamma3
   _REAL_      tgtmdfrc,tgtrmsd,pencut,scalm,tausw,taumet,omega
   integer     ierr
   integer     ifcrst
   integer     ifind
   integer     imcdo
   integer     itotst
   integer     jn
   logical     mdin_cntrl  ! true if this namelist exists in mdin
   integer     mxgrp
   integer     itgtmd,ipol,iesp,n3b,nion,acon,nmropt,iscale,ipnlty,noeskp,intreq,irlreq
   character(len=8) date
   character(len=10) time
   _REAL_      timlim ! retired 

   
   namelist /cntrl/ irest,ibelly, &
         ntx,ntxo,ntcx,ig,tempi, &
         ntb,ntt,temp0,tautp, &
         ntp,pres0,comp,taup,npscal, &
         nscm,nstlim,t,dt, &
         ntc,ntcc,nconp,tol,ntf,ntn,nsnb, &
         cut,scnb,scee,dielc, &
         ntpr,ntwx,ntwv,ntwe,ntave,ntpp,ioutfm, &
         ntr,nrc,ntrx,taur,nmropt, &
         ivcap,fcap,imin,drms,dele,dx0, &
         pencut,ipnlty,iscale,scalm,noeskp, &
         maxcyc,ncyc,ntmin,vlimit, &
         mxsub,ipol,jfastw,watnam,owtnm,hwtnm1,hwtnm2, &
         iesp, &
         ntwprt,n3b,nion,at1,at2,acon,beta3,gamma3,tausw, &
         ntwr,plevel,iyammp,imcdo, &
         igb,rgbmax,saltcon,offset,gbsa,vrand, &
         surften,iwrap,nrespa,nrespai,gamma_ln,extdiel,intdiel, &
         cut_inner,icfe,clambda,klambda, &
         rbornstat,lastrst,lastist,itgtmd,tgtrmsd,tgtmdfrc, &
         idecomp,temp0les,restraintmask,restraint_wt,bellymask, &
         rdt,icnstph,solvph,ntcnstph &
         ,timlim  ! all retired 


   ! Define default water residue name and the names of water oxygen & hydrogens
   
   data watdef/'WAT ','O   ','H1  ','H2  '/
   
   !     ----- READ THE CONTROL DATA AND OPEN DIFFERENT FILES -----
   
   if (mdout /= "stdout" ) &
         call amopen(6,mdout,owrite,'F','W')
   call amopen(5,mdin,'O','F','R')
   write(6,9308)
   call date_and_time( DATE=date, TIME=time )
   write(6,'(12(a))') '| Run on ', date(5:6), '/', date(7:8), '/',  &
        date(1:4), ' at ', time(1:2), ':', time(3:4), ':', time(5:6)
   if (owrite /= 'N') write(6, '(2x,a)') '[-O]verwriting output'
   
   ! Echo the file assignments to the user:
   
   write(6,9700) 'MDIN'   ,mdin(1:70)  , 'MDOUT' ,mdout(1:70) , &
         'INPCRD' ,inpcrd(1:70), 'PARM'  ,parm(1:70)  , &
         'RESTRT',restrt(1:70) , 'REFC'  ,refc(1:70)  , &
         'MDVEL' ,mdvel(1:70)  , 'MDEN'   ,mden(1:70) , &
         'MDCRD' ,mdcrd(1:70)  , 'MDINFO' ,mdinfo(1:70), &
         'INPDIP', inpdip(1:70), 'RSTDIP', rstdip(1:70)
   
   ! Echo the input file to the user:
   
   call echoin(5,6)
   
   
   !     ----- READ DATA CHARACTERIZING THE MD-RUN -----
   
   read(5,'(20a4)') title
   
   !       ----read input in namelist format, first setting up defaults
   
   timlim = RETIRED_INPUT_OPTION
   irest = 0
   ibelly = 0
   ipol = 0
   iesp = 0
   ntx = 1
   ntxo = 1
   ig = 71277
   tempi = ZERO
   ntb = 1
   ntt = 0
   temp0 = 300.0d0

   tautp = ONE
   ntp = 0
   pres0 = ONE
   comp = 44.6d0
   taup = ONE
   npscal = 1
   nscm = 1000
   nstlim = 1
   t = ZERO
   dt = 0.001d0
   ntc = 1
   tol = 0.00001
   ntf = 1
   nsnb = 25
   cut =  EIGHT
   scnb = TWO
   scee = 1.2d0
   dielc = ONE
   ntpr = 50
   ntwr = 500
   ntwx = 0
   ntwv = 0
   ntwe = 0
   ntave = 0
   ioutfm = 0
   ntr = 0
   ntrx = 1
   ivcap = 0
   fcap = 1.5d0
   
   ! carlos targeted MD, like ntr
   
   itgtmd=0
   tgtrmsd=0.
   tgtmdfrc=0.

   pencut = 0.1d0
   taumet = 0.0001d0
   omega = 500.0d0
   ipnlty = 1
   scalm = 100.0d0
   iscale = 0
   noeskp = 1
   nmropt = 0
   tausw = 0.1d0
   imin = 0
   isftrp = 0
   rwell = ONE
   maxcyc = 1
   ncyc = 10
   ntmin = 1
   dx0 = 0.01d0
   drms = 1.0d-4
   vlimit = 20.0d0
   mxsub = 1
   jfastw = 0
   watnam = '    '
   owtnm =  '    '
   hwtnm1 = '    '
   hwtnm2 = '    '
   ntwprt = 0
   plevel = 1
   igb = 10
   rgbmax = 25.d0
   saltcon = ZERO
   offset = 0.09d0
   iyammp = 0
   imcdo = -1
   gbsa = 0
   vrand=1000
   surften = 0.005d0
   iwrap = 0
   nrespa = 1
   nrespai = 1
   irespa = 1
   gamma_ln = ZERO
   extdiel = 78.5d0
   intdiel = ONE
   gbgamma = ZERO
   gbbeta = ZERO
   gbalpha = ONE
   iconstreff = 0
   cut_inner = EIGHT
   icfe = 0
   clambda = ZERO
   klambda = 1
   rbornstat = 0
   idecomp = 0
   lastrst = 2000000
   lastist = 2000000
   restraintmask=''
   restraint_wt = ZERO
   bellymask=''

   icnstph = 0
   solvph = SEVEN
   ntcnstph = 10
   
   !     Check to see if "cntrl" namelist has been defined.
   
   mdin_cntrl=.false.
   mdin_ewald=.false.
   mdin_pb=.false.
   call nmlsrc('cntrl',5,ifind)
   if (ifind /= 0) mdin_cntrl=.true.
   call nmlsrc('ewald',5,ifind)
   if (ifind /= 0) mdin_ewald=.true.
   call nmlsrc('pb',5,ifind)
   if (ifind /= 0) mdin_pb=.true.

   rewind 5
   if ( mdin_cntrl ) then
      read(5,nml=cntrl)
   else
      write(6, '(1x,a,/)') 'Could not find cntrl namelist'
      call mexit(6,1)
   end if
   
   !     --- vars have been read ---
   
   write(6,9309)

   ! emit warnings for retired cntrl namelist variables

   if ( timlim /= RETIRED_INPUT_OPTION ) then
      write(6,'(/,a,/,a,/,a)') 'Warning: timlim has been retired.', &
            '  Check the Retired Namelist Variables Appendix in the manual.'
   end if

   call printflags()
   
   ifcrst= 0
   
   ishake = 0
   if (ntc > 1) ishake = 1
   
   !--------------------------------------------------------------------
   ! Set up some parameters for GB simulations:
   !--------------------------------------------------------------------
   
   if( igb == 2 ) then
      
      !       --- use our best guesses for Onufriev/Case GB  (GB^OBC I)
      
      gbgamma = 2.909125d0
      gbbeta = ZERO
      gbalpha = 0.8d0
   end if

   if( igb == 5 ) then
      
      !       --- use our second best guesses for Onufriev/Case GB (GB^OBC II)
      
      gbgamma = 4.851d0
      gbbeta = 0.8d0
      gbalpha = ONE
   end if
   
   !--------------------------------------------------------------------
   ! If user has requested PB electrostatics, read some more input
   !--------------------------------------------------------------------

   if ( igb == 10 .or. igb == 11 .or. igb == 12 ) then
      call pb_read( igb )
   end if
   
   ! -------------------------------------------------------------------
   ! If the user has requested NMR restraints, do a cursory read of the
   ! restraints file(s) now to determine the amount of memory necessary
   ! for these restraints:
   ! -------------------------------------------------------------------
   
   intreq = 0
   irlreq = 0
   
   ! Set the definition of the water molecule. The default definition is in
   ! WATDEF(4).
   
   read(watdef(1),'(A4)') iwtnm
   read(watdef(2),'(A4)') iowtnm
   read(watdef(3),'(A4)') ihwtnm(1)
   read(watdef(4),'(A4)') ihwtnm(2)
   if (watnam /= '    ') read(watnam,'(A4)') iwtnm
   if (owtnm /= '    ') read(owtnm, '(A4)') iowtnm
   if (hwtnm1 /= '    ') read(hwtnm1,'(A4)') ihwtnm(1)
   if (hwtnm2 /= '    ') read(hwtnm2,'(A4)') ihwtnm(2)
   
   return
   
   ! --- input file polar opts read err trapping:
   
   1155 write(6,*) ' ** EOF reading N3B,NION for 3-body option'
   call mexit(6,1)
   1156 write(6,*) ' ** EOF reading triplets for 3-body option'
   call mexit(6,1)
   
   9308 format(/10x,55('-'),/10x, &
         'Amber 9  PBSA                   Scripps/UCSF 2006', &
         /10x,55('-')/)
   9309 format(/80('-')/'   1.  RESOURCE   USE: ',/80('-')/)
   9700 format(/,'File Assignments:',/,12('|',a6,': ',a,/))
end subroutine mdread1 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Read constant pH file and initialize.
subroutine cnstphread(stateinf,resstate,protcnt,trescnt,&
      statene,chrgdat,charge)
   implicit none
#  include "dynph.h"
#  include "files.h"
#  include "random.h"
#  include "constants.h"

   type (const_ph_info) :: stateinf(0:TITR_RES_C)
   integer :: resstate(0:TITR_RES_C-1), protcnt(0:TITR_STATES_C-1)
   integer, intent(out) :: trescnt
   _REAL_ ::  statene(0:TITR_STATES_C-1), chrgdat(0:ATOM_CHRG_C-1),charge(1:*)
   integer itres, iselres, istat, iatom
   integer icumstat, icumchrg
   character(len=40) :: resname(0:TITR_STATES_C)
   
   common /cnstphresname/ resname
   namelist /cnstph/ stateinf, resstate, protcnt, chrgdat, statene,trescnt,resname

   icumstat = 0
   icumchrg = 0

   write(6,'(a,a)') 'reading charge increments from file: ',cpin
   call amopen(CNSTPH_UNIT,cpin,'O','F','R')
   read(18, nml=cnstph)
   do iatom=0, ATOM_CHRG_C-1
      chrgdat(iatom) = chrgdat(iatom) * AMBER_ELECTROSTATIC
   end do
   close(CNSTPH_UNIT)

   !     Alter charges to match specified initial states
   do itres=0,trescnt-1
      do iatom=0, stateinf(itres)%num_atoms - 1 ! For each atom in selected residue
         charge(iatom+stateinf(itres)%first_atom) & !set atom charge to
               = chrgdat(stateinf(itres)%first_charge + iatom & !corresponding atom charge value
               + resstate(itres) * stateinf(itres)%num_atoms & !from selected new state
               )
      end do
   end do
   !     Error (overflow) checking
   if (trescnt > TITR_RES_C) then
      write(6,*) 'Too many titrating residues', &
            'alter dynph.h, recompile'
      stop
   end if
   do itres=0,trescnt-1         !Find res with reference to highest # state
      if (stateinf(itres)%first_state > icumstat) then
         icumstat = stateinf(itres)%first_state
         iselres = itres
      end if
   end do
   icumstat = stateinf(iselres)%first_state + stateinf(iselres)%num_states
   icumchrg = stateinf(iselres)%first_charge + &
         stateinf(iselres)%num_atoms * stateinf(iselres)%num_states
   
   if (icumstat > TITR_STATES_C) then
      write(6,*) 'Too many titrating states', &
            ' alter dynph.h, recompile'
      stop
   end if
   if (icumchrg > ATOM_CHRG_C) then
      write(6,*) 'Too much charge data', &
            'alter dynph.h, recompile'
      stop
   end if
end subroutine cnstphread 



!======================================================================
!          MDREAD2
!======================================================================

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Initialize to defaults and print the inputable variables.
subroutine mdread2(x,ix,ih,ipairs,r_stack,i_stack)

   use decomp, only : jgroup, index
   use findmask

   implicit none
   _REAL_ x(*),r_stack(*)
   _REAL_      pencut,scalm,tausw,taumet,omega
   integer ix(*),ipairs(*),i_stack(*)
   character(len=4) ih(*)
   integer nbond
   integer atom1,atom2
   integer ntmp
   logical belly,konst
   character(len=1) atsymb,atsymb2
   character(len=2) atype
   integer ngrp,inerr,nr,iaci,ic,ir,i,mxresat,j
   integer ipol,iesp,nmropt,iscale,ipnlty,noeskp,intreq,irlreq
   _REAL_ dummy,rvdw,dcharge
#ifdef MPI
   !     =========================== AMBER/MPI ===========================
#ifdef MPI_DOUBLE_PRECISION
#  undef MPI_DOUBLE_PRECISION
#endif
#  include "mpif.h"
   !     ========================= END AMBER/MPI =========================
#endif
#  include "constants.h"
#  include "files.h"
#  include "md.h"
#  include "box.h"
#  include "memory.h"
#  include "parms.h"
#  include "extra.h"
#  include "def_time.h"
   
   close(unit=8)
   
   ! -------------------------------------------------------------------
   !     ----- SET THE DEFAULT VALUES FOR SOME VARIABLES -----
   ! -------------------------------------------------------------------
   
   nrp = natom
   
   if (ifbox == 1) write(6, '(/5x,''BOX TYPE: RECTILINEAR'')')
   if (ifbox == 2) write(6, '(/5x,''BOX TYPE: TRUNCATED OCTAHEDRON'')')
   if (ifbox == 3) write(6, '(/5x,''BOX TYPE: GENERAL'')')
   
   nsolut =  nrp
   if ( nscm > 0 .and. ntb == 0 ) then
      ndfmin = 6   ! both translation and rotation com motion removed
      if (nsolut == 1) ndfmin = 3
      if (nsolut == 2) ndfmin = 5
   else if ( nscm > 0 ) then
      ndfmin = 3    ! just translation com will be removed
   else
      ndfmin = 0
   end if
   if (ibelly > 0) then   ! No COM Motion Removal, ever.
      nscm = 9999999
      ndfmin = 0
   end if
   if(nscm <= 0) nscm = 9999999
   init = 3
   if (irest > 0) init = 4
   if (scnb == ZERO ) scnb = TWO
   if (dielc <= ZERO ) dielc = ONE
   if (tautp <= ZERO ) tautp = 0.2d0
   if (taup <= ZERO ) taup = 0.2d0
   
   !     ----- RESET THE CAP IF NEEDED -----
   
   if(ivcap == 2) ifcap = 0
   
   ! -------------------------------------------------------------------
   !     ----- PRINT DATA CHARACTERIZING THE RUN -----
   ! -------------------------------------------------------------------
   
   nr = nrp
   nmropt = 0
   iesp = 0
   ipol = 0
   write(6,9328)
   write(6,9008) title
   write(6,'(/a)') 'General flags:'
   write(6,'(5x,2(a,i8))') 'imin    =',imin,', nmropt  =',nmropt

   write(6,'(/a)') 'Nature and format of input:'
   write(6,'(5x,4(a,i8))') 'ntx     =',ntx,', irest   =',irest, &
         ', ntrx    =',ntrx

   write(6,'(/a)') 'Nature and format of output:'
   write(6,'(5x,4(a,i8))') 'ntxo    =',ntxo,', ntpr    =',ntpr, &
         ', ntrx    =',ntrx,', ntwr    =',ntwr
   write(6,'(5x,4(a,i8))') 'iwrap   =',iwrap,', ntwx    =',ntwx, &
         ', ntwv    =',ntwv,', ntwe    =',ntwe
   write(6,'(5x,3(a,i8),a,i7)') 'ioutfm  =',ioutfm, &
         ', ntwprt  =',ntwprt, &
         ', idecomp =',idecomp,', rbornstat=',rbornstat

   write(6,'(/a)') 'Potential function:'
   write(6,'(5x,5(a,i8))') 'ntf     =',ntf,', ntb     =',ntb, &
         ', igb     =',igb,', nsnb    =',nsnb
   write(6,'(5x,4(a,i8))') 'ipol    =',ipol,', gbsa    =',gbsa, &
         ', iesp    =',iesp
   write(6,'(5x,3(a,f10.5))') 'dielc   =',dielc, &
         ', cut     =',cut,', intdiel =',intdiel
   
   if( igb /= 0 .and. igb < 10 ) then
      write(6,'(5x,3(a,f10.5))') 'saltcon =',saltcon, &
            ', offset  =',offset,', gbalpha= ',gbalpha
      write(6,'(5x,3(a,f10.5))') 'gbbeta  =',gbbeta, &
            ', gbgamma =',gbgamma,', surften =',surften
      write(6,'(5x,3(a,f10.5))') 'rdt     =',rdt, &
            ', rgbmax  =',rgbmax
   end if

   write(6,'(5x,3(a,f10.5))') 'scnb    =',scnb, &
         ', scee    =',scee
   
   write(6,'(/a)') 'Frozen or restrained atoms:'
   write(6,'(5x,4(a,i8))') 'ibelly  =',ibelly,', ntr     =',ntr

   if( imin /= 0 ) then
      write(6,'(/a)') 'Energy minimization:'
      ! print inputable variables applicable to all minimization methods.
      write(6,'(5x,4(a,i8))') 'maxcyc  =',maxcyc,', ncyc    =',ncyc, &
            ', ntmin   =',ntmin
      write(6,'(5x,2(a,f10.5))') 'dx0     =',dx0, ', drms    =',drms

      ! Input flag ntmin determines the method of minimization
      select case ( ntmin )
      case ( 0, 1, 2 )
         ! no specific output
      case default
         ! invalid ntmin
         write(6,'(/2x,a,i3,a)') 'Error: Invalid NTMIN (',ntmin,').'
         stop
      end select

   else
      write(6,'(/a)') 'Molecular dynamics:'
      write(6,'(5x,4(a,i8))') 'nstlim  =',nstlim,', nscm    =',nscm, &
            ', nrespa  =',nrespa
      write(6,'(5x,3(a,f10.5))') 't       =',t, &
            ', dt      =',dt,', vlimit  =',vlimit
      
      if( ntt == 1 ) then
         write(6,'(/a)') 'Berendsen (weak-coupling) temperature regulation:'
         write(6,'(5x,3(a,f10.5))') 'temp0   =',temp0, &
               ', tempi   =',tempi,', tautp   =', tautp
      else if( ntt == 2 ) then
         write(6,'(/a)') 'Anderson (strong collision) temperature regulation:'
         write(6,'(5x,4(a,i8))') 'ig      =',ig, ', vrand   =',vrand
         write(6,'(5x,3(a,f10.5))') 'temp0   =',temp0, ', tempi   =',tempi
      else if( ntt == 3 ) then
         write(6,'(/a)') 'Langevin dynamics temperature regulation:'
         write(6,'(5x,4(a,i8))') 'ig      =',ig
         write(6,'(5x,3(a,f10.5))') 'temp0   =',temp0, &
               ', tempi   =',tempi,', gamma_ln=', gamma_ln
      end if

   end if

   if( ntc /= 1 ) then
      write(6,'(/a)') 'SHAKE:'
      write(6,'(5x,4(a,i8))') 'ntc     =',ntc,', jfastw  =',jfastw
      write(6,'(5x,3(a,f10.5))') 'tol     =',tol
   end if

   if( ifcap /= 0 ) then
      write(6,'(/a)') 'Water cap:'
      write(6,'(5x,4(a,i8))') 'ivcap   =',ivcap,', natcap  =',natcap
      write(6,'(5x,2(a,f10.5))') 'fcap    =',fcap, ', cutcap  =',cutcap
      write(6,'(5x,3(a,f10.5))') 'xcap    =',xcap, ', ycap    =',ycap,    &
                                 ', zcap    =',zcap
   end if

   if( icfe /= 0 ) then
      write(6,'(/a)') 'Free energy options:'
      write(6,'(5x,4(a,i8))') 'klambda =',klambda
      write(6,'(5x,3(a,f10.5))') 'clambda =',clambda
   end if

   if( icnstph /= 0) then
      write(6, '(/a)') 'Constant pH options:'
      write(6, '(5x,a,i8)') 'ntcnstph =', ntcnstph
      write(6, '(5x,a,f10.5)') 'solvph =', solvph
   end if
   
   cut = cut*cut
   cut_inner = cut_inner*cut_inner
   
   
   !------------------------------------------------------------------------
   ! If user has requested generalized born electrostatics, set up variables
   !------------------------------------------------------------------------
   
   if( gbsa == 2 .and. &
       ((imin == 0 .and. nstlim > 1) .or. &
        (imin == 1 .and. maxcyc > 1)) ) then
      write(0,*) 'GBSA=2 only works for single point energy calc'
      call mexit( 6,1 )
   end if
   if( igb /= 0 .and. igb < 10) then
      
      !       put fs(i)*(rborn(i) - offset) into the "fs" array
      
      fsmax = 0.d0
      do i=1,natom
         x(l96-1+i) = x(l96-1+i)*( x(l97-1+i) - offset )
         fsmax = max( fsmax, x(l96-1+i) )
         if (rbornstat == 1) then
            x(l186-1+i) = 0.d0
            x(l187-1+i) = 999.d0
            x(l188-1+i) = 0.d0
            x(l189-1+i) = 0.d0
         end if
      end do
      
      !     ---------------------------------------------------------------------
      !       ---get Debye-Huckel kappa (A**-1) from salt concentration (M), assuming:
      !         T = 298.15, epsext=78.5,
      
      kappa = sqrt( 0.10806d0 * saltcon )
      
      !       ---scale kappa by 0.73 to account(?) for lack of ion exlcusions:
      
      kappa = 0.73d0* kappa
      
      !     ---------------------------------------------------------------------
      
      if ( gbsa == 1 ) then
         
         !     --- assign parameters for calculating SASA according to the
         !         LCPO method ---
         
         do i=1,natom
            ix(i80+i-1)=0
         end do
         
         !         --- get the number of bonded neighbors for each atom:
         
         do i=1,nbona
            atom1=ix(iiba-1+i)/3+1
            atom2=ix(ijba-1+i)/3+1
            ix(i80+atom1-1)=ix(i80+atom1-1)+1
            ix(i80+atom2-1)=ix(i80+atom2-1)+1
         end do
         
         !         --- construct parameters for SA calculation; note that the
         !             radii stored in L165 are augmented by 1.4 Ang.
         
         do i=1,natom
            write(atype,'(a2)') ih(m06+i-1)
            nbond=ix(i80+i-1)
            if (atype == 'CT') then
               if (nbond == 1) then
                  x(l165-1+i) = 1.70d0 + 1.4d0
                  x(l170-1+i) = 0.7887d0
                  x(l175-1+i) = -0.28063d0
                  x(l180-1+i) = -0.0012968d0
                  x(l185-1+i) = 0.00039328d0
               else if (nbond == 2) then
                  x(l165-1+i) = 1.70d0 + 1.4d0
                  x(l170-1+i) = 0.56482d0
                  x(l175-1+i) = -0.19608d0
                  x(l180-1+i) = -0.0010219d0
                  x(l185-1+i) = 0.0002658d0
               else if (nbond == 3) then
                  x(l165-1+i) = 1.70d0 + 1.4d0
                  x(l170-1+i) = 0.23348d0
                  x(l175-1+i) = -0.072627d0
                  x(l180-1+i) = -0.00020079d0
                  x(l185-1+i) = 0.00007967d0
               else if (nbond == 4) then
                  x(l165-1+i) = 1.70d0 + 1.4d0
                  x(l170-1+i) = 0.00000d0
                  x(l175-1+i) = 0.00000d0
                  x(l180-1+i) = 0.00000d0
                  x(l185-1+i) = 0.00000d0
               end if
            else if (atype(1:1) == 'C' .or. atype(1:1) == 'c') then
               if (nbond == 2) then
                  x(l165-1+i) = 1.70d0 + 1.4d0
                  x(l170-1+i) = 0.51245d0
                  x(l175-1+i) = -0.15966d0
                  x(l180-1+i) = -0.00019781d0
                  x(l185-1+i) = 0.00016392d0
               else if (nbond == 3) then
                  x(l165-1+i) = 1.70d0 + 1.4d0
                  x(l170-1+i) = 0.070344d0
                  x(l175-1+i) = -0.019015d0
                  x(l180-1+i) = -0.000022009d0
                  x(l185-1+i) = 0.000016875d0
               end if
            else if (atype == 'O ') then
               x(l165-1+i) = 1.60d0 + 1.4d0
               x(l170-1+i) = 0.68563d0
               x(l175-1+i) = -0.1868d0
               x(l180-1+i) = -0.00135573d0
               x(l185-1+i) = 0.00023743d0
            else if (atype == 'O2') then
               x(l165-1+i) = 1.60d0 + 1.4d0
               x(l170-1+i) = 0.88857d0
               x(l175-1+i) = -0.33421d0
               x(l180-1+i) = -0.0018683d0
               x(l185-1+i) = 0.00049372d0
            else if (atype(1:1) == 'O' .or. atype(1:1) == 'o') then
               if (nbond == 1) then
                  x(l165-1+i) = 1.60d0 + 1.4d0
                  x(l170-1+i) = 0.77914d0
                  x(l175-1+i) = -0.25262d0
                  x(l180-1+i) = -0.0016056d0
                  x(l185-1+i) = 0.00035071d0
               else if (nbond == 2) then
                  x(l165-1+i) = 1.60d0 + 1.4d0
                  x(l170-1+i) = 0.49392d0
                  x(l175-1+i) = -0.16038d0
                  x(l180-1+i) = -0.00015512d0
                  x(l185-1+i) = 0.00016453d0
               end if
            else if (atype == 'N3') then
               if (nbond == 1) then
                  x(l165-1+i) = 1.65d0 + 1.4d0
                  x(l170-1+i) = 0.078602d0
                  x(l175-1+i) = -0.29198d0
                  x(l180-1+i) = -0.0006537d0
                  x(l185-1+i) = 0.00036247d0
               else if (nbond == 2) then
                  x(l165-1+i) = 1.65d0 + 1.4d0
                  x(l170-1+i) = 0.22599d0
                  x(l175-1+i) = -0.036648d0
                  x(l180-1+i) = -0.0012297d0
                  x(l185-1+i) = 0.000080038d0
               else if (nbond == 3) then
                  x(l165-1+i) = 1.65d0 + 1.4d0
                  x(l170-1+i) = 0.051481d0
                  x(l175-1+i) = -0.012603d0
                  x(l180-1+i) = -0.00032006d0
                  x(l185-1+i) = 0.000024774d0
               end if
            else if (atype(1:1) == 'N' .or. atype(1:1) == 'n') then
               if (nbond == 1) then
                  x(l165-1+i) = 1.65d0 + 1.4d0
                  x(l170-1+i) = 0.73511d0
                  x(l175-1+i) = -0.22116d0
                  x(l180-1+i) = -0.00089148d0
                  x(l185-1+i) = 0.0002523d0
               else if (nbond == 2) then
                  x(l165-1+i) = 1.65d0 + 1.4d0
                  x(l170-1+i) = 0.41102d0
                  x(l175-1+i) = -0.12254d0
                  x(l180-1+i) = -0.000075448d0
                  x(l185-1+i) = 0.00011804d0
               else if (nbond == 3) then
                  x(l165-1+i) = 1.65d0 + 1.4d0
                  x(l170-1+i) = 0.062577d0
                  x(l175-1+i) = -0.017874d0
                  x(l180-1+i) = -0.00008312d0
                  x(l185-1+i) = 0.000019849d0
               end if
            else if (atype == 'SH') then
               x(l165-1+i) = 1.90d0 + 1.4d0
               x(l170-1+i) = 0.7722d0
               x(l175-1+i) = -0.26393d0
               x(l180-1+i) = 0.0010629d0
               x(l185-1+i) = 0.0002179d0
            else if (atype(1:1) == 'S' .or. atype(1:1) == 's') then
               x(l165-1+i) = 1.90d0 + 1.4d0
               x(l170-1+i) = 0.54581d0
               x(l175-1+i) = -0.19477d0
               x(l180-1+i) = -0.0012873d0
               x(l185-1+i) = 0.00029247d0
            else if (atype(1:1) == 'P' .or. atype(1:1) == 'p') then
               if (nbond == 3) then
                  x(l165-1+i) = 1.90d0 + 1.4d0
                  x(l170-1+i) = 0.3865d0
                  x(l175-1+i) = -0.18249d0
                  x(l180-1+i) = -0.0036598d0
                  x(l185-1+i) = 0.0004264d0
               else if (nbond == 4) then
                  x(l165-1+i) = 1.90d0 + 1.4d0
                  x(l170-1+i) = 0.03873d0
                  x(l175-1+i) = -0.0089339d0
                  x(l180-1+i) = 0.0000083582d0
                  x(l185-1+i) = 0.0000030381d0
               end if
            else if (atype(1:1) == 'Z') then
               x(l165-1+i) = 0.00000d0 + 1.4d0
               x(l170-1+i) = 0.00000d0
               x(l175-1+i) = 0.00000d0
               x(l180-1+i) = 0.00000d0
               x(l185-1+i) = 0.00000d0
            else if (atype(1:1) == 'H' .or. atype(1:1) == 'h') then
               x(l165-1+i) = 0.00000d0 + 1.4d0
               x(l170-1+i) = 0.00000d0
               x(l175-1+i) = 0.00000d0
               x(l180-1+i) = 0.00000d0
               x(l185-1+i) = 0.00000d0
            else if (atype == 'MG') then
               !             Mg radius = 0.99A: ref. 21 in J. Chem. Phys. 1997, 107, 5422
               !             Mg radius = 1.18A: ref. 30 in J. Chem. Phys. 1997, 107, 5422
               !             Mg radius = 1.45A: Aqvist 1992
               x(l165-1+i) = 1.18d0 + 1.4d0
               !             The following values were taken from O.sp3 with two bonded neighbors
               !             -> O has the smallest van der Waals radius compared to all other
               !             elements which had been parametrized
               x(l170-1+i) = 0.49392d0
               x(l175-1+i) = -0.16038d0
               x(l180-1+i) = -0.00015512d0
               x(l185-1+i) = 0.00016453d0
            else
               write( 0,* ) 'bad atom type: ',atype
               call mexit( 6,1 )
            end if  ! (atype == 'CT')
            !           write(6,*) i,x(L165-1+i),x(L170-1+i),x(L175-1+i)
            !           write(6,*) x(L180-1+i),x(L185-1+i)
         end do  !  i=1,natom
         !
      else if ( gbsa == 2 ) then

         !     --- assign parameters for calculating SASA according to the
         !         ICOSA method; the radii are augmented by 1.4 A ---

         do i=1,natom
            write(atype,'(a2)') ih(m06+i-1)
            if (atype(1:1) == 'N' .or. atype(1:1) == 'n') then
               x(L165-1+i) = 1.55d0 + 1.4d0
            else if (atype(1:1) == 'C' .or. atype(1:1) == 'c') then
               x(L165-1+i) = 1.70d0 + 1.4d0
            else if (atype(1:1) == 'H' .or. atype(1:1) == 'h' .or. &
                     atype == '1H' .or. &
                     atype == '2H' .or. &
                     atype == '3H') then
               x(L165-1+i) = 1.20d0 + 1.4d0
            else if (atype(1:1) == 'O' .or. atype(1:1) == 'o') then
               x(L165-1+i) = 1.50d0 + 1.4d0
            else if (atype(1:1) == 'P' .or. atype(1:1) == 'p') then
               x(L165-1+i) = 1.80d0 + 1.4d0
            else if (atype(1:1) == 'S' .or. atype(1:1) == 's') then
               x(L165-1+i) = 1.80d0 + 1.4d0
            else if (atype == 'MG' .or. atype == 'mg') then
               !             Mg radius = 0.99A: ref. 21 in J. Chem. Phys. 1997, 107, 5422
               !             Mg radius = 1.18A: ref. 30 in J. Chem. Phys. 1997, 107, 5422
               !             Mg radius = 1.45A: Aqvist 1992
               x(L165-1+i) = 1.18d0 + 1.4d0
            else
               write( 0,* ) 'bad atom type: ',atype
               call mexit( 6,1 )
            end if  !  atype(1:1) == 'N'
            x(L170-1+i) = 0.0d0
            x(L175-1+i) = 0.0d0
            x(L180-1+i) = 0.0d0
            x(L185-1+i) = 0.0d0
            !  write(6,*) i,' ',atype,x(L165-1+i)
         end do  !  i=1,natom

      end if ! ( gbsa == 1 )
      
   end if  ! ( igb /= 0 .and. igb < 10)
   
   !------------------------------------------------------------------------
   ! If user has requested Poisson-Boltzmann electrostatics, set up variables
   !------------------------------------------------------------------------
    
   if ( igb >= 10 ) then
      call pb_init(ifcap,natom,nres,ntypes,nbonh,nbona,ix(i02),ix(i04),ix(i06),ix(i08),ix(i10),&
                   ix(iibh),ix(ijbh),ix(iiba),ix(ijba),ih(m02),ih(m04),ih(m06),x(l15),x(l97))
   end if  ! ( igb >= 10 ) 
 
   if (icnstph /= 0) then
      !     Read charge data and alter current charges accordingly
      call cnstphread(ix(icpstinf),ix(icpresst),ix(icpptcnt), &
            ix(icptrsct),x(lcpene),x(lcpcrg),x(l15))

      !     Fill proposed charges array from current charges
      do i=1,natom
         x(l190-1+i) = x(l15-1+i)
      end do
   end if

   if( iyammp /= 0 ) &
         write( 6, '(a)' ) '  Using yammp non-bonded potential'
   
   ! -------------------------------------------------------------------
   !     --- checks on bogus data ---
   ! -------------------------------------------------------------------
   
   inerr = 0
   
   if (igb /= 1 .and. igb /= 2 .and. igb /= 5 .and. igb > 12) then
      write(6,'(/2x,a,i3,a)') 'IGB (',igb,') must be 1,2,5 or 10.'
      inerr = 1
   end if
   if (irest /= 0 .and. irest /= 1) then
      write(6,'(/2x,a,i3,a)') 'IREST (',irest,') must be 0 or 1.'
      inerr = 1
   end if
   if (ibelly /= 0 .and. ibelly /= 1) then
      write(6,'(/2x,a,i3,a)') 'IBELLY (',ibelly,') must be 0 or 1.'
      inerr = 1
   end if
   if (imin < 0) then
      write(6,'(/2x,a,i3,a)') 'IMIN (',imin,') must be >= 0.'
      inerr = 1
   end if

   if (ntx < 1 .or. ntx > 7) then
      write(6,'(/2x,a,i3,a)') 'NTX (',ntx,') must be in 1..7'
      inerr = 1
   end if
   if (ntxo /= 0 .and. ntxo /= 1) then
      write(6,'(/2x,a,i3,a)') 'NTXO (',ntxo,') must be 1 or 0.'
      inerr = 1
   end if
   
   if (ntb /= 0 .and. ntb /= 1 .and. ntb /= 2) then
      write(6,'(/2x,a,i3,a)') 'NTB (',ntb,') must be 0, 1 or 2.'
      inerr = 1
   end if
   
   if (ntt < 0 .or. ntt > 3) then
      write(6,'(/2x,a,i3,a)') 'NTT (',ntt,') must be between 0 and 3.'
      inerr = 1
   end if
   if (ntt == 1 .and. tautp < dt) then
      write(6, '(/2x,a,f6.2,a)') 'TAUTP (',tautp,') < DT (step size)'
      inerr = 1
   end if
   if( ntt /= 3 .and. gamma_ln > 0.d0 ) then
      write(6,'(a)') 'ntt must be 3 if gamma_ln > 0'
      inerr = 1
   end if
   
   if (ntp /= 0 .and. ntp /= 1 .and. ntp /= 2) then
      write(6,'(/2x,a,i3,a)') 'NTP (',ntp,') must be 0, 1 or 2.'
      inerr = 1
   end if
   if (ntp > 0 .and. taup < dt) then
      write(6, '(/2x,a,f6.2,a)') 'TAUP (',taup,') < DT (step size)'
      inerr = 1
   end if
   if (npscal /= 1) then
      write(6,'(/2x,a,i3,a)') 'NPSCAL (',npscal,') must be 1.'
      inerr = 1
   end if
   
   if (ntc < 1 .or. ntc > 4) then
      write(6,'(/2x,a,i3,a)') 'NTC (',ntc,') must be 1,2,3 or 4.'
      inerr = 1
   end if
   if (jfastw < 0 .or. jfastw > 4) then
      write(6,'(/2x,a,i3,a)') 'JFASTW (',jfastw,') must be 0->4.'
      inerr = 1
   end if
   
   if (ntf < 1 .or. ntf > 8) then
      write(6,'(/2x,a,i3,a)') 'NTF (',ntf,') must be in 1..8.'
      inerr = 1
   end if
   
   if (scee == 0.0d0) then
      write(6,'(/2x,a)') 'SCEE must be set explicitly'
      inerr = 1
   end if
   
   if (ioutfm /= 0 .and. ioutfm /= 1) then
      write(6,'(/2x,a,i3,a)') 'IOUTFM (',ioutfm,') must be 0 or 1.'
      inerr = 1
   end if
   
   if (ntpr < 0) then
      write(6,'(/2x,a,i3,a)') 'NTPR (',ntpr,') must be >= 0.'
      inerr = 1
   end if
   if (ntwx < 0) then
      write(6,'(/2x,a,i3,a)') 'NTWX (',ntwx,') must be >= 0.'
      inerr = 1
   end if
   if (ntwv < 0) then
      write(6,'(/2x,a,i3,a)') 'NTWV (',ntwv,') must be >= 0.'
      inerr = 1
   end if
   if (ntwe < 0) then
      write(6,'(/2x,a,i3,a)') 'NTWE (',ntwe,') must be >= 0.'
      inerr = 1
   end if
   if (ntave < 0) then
      write(6,'(/2x,a,i3,a)') 'NTAVE (',ntave,') must be >= 0.'
      inerr = 1
   end if
   if (ntr /= 0 .and. ntr /= 1) then
      write(6,'(/2x,a,i3,a)') 'NTR (',ntr,') must be 0 or 1.'
      inerr = 1
   end if
   if (ntrx /= 0 .and. ntrx /= 1) then
      write(6,'(/2x,a,i3,a)') 'NTRX (',ntrx,') must be 1 or 0.'
      inerr = 1
   end if
   
   if (idecomp < 0 .or. idecomp > 4) then
      write(6,'(/2x,a)') 'IDECOMP must be 0..4'
      inerr = 1
   end if
   if (idecomp > 0 .and. &
       (imin == 0 .or. &
        (imin == 1 .and. maxcyc > 1))) then
      write(6,*) 'IDECOMP>1 only works for IMIN=1 and MAXCYC=1'
      inerr = 1
   end if
   
   !     -- consistency checking
   
   if (imin > 0.and.nrespa > 1)  then
      write(6,'(/2x,a)') 'For minimization, set nrespa,nrespai=1'
      inerr = 1
   end if
   if (ntp > 0 .and. nrespa > 1) then
      write(6,'(/2x,a)') 'nrespa must be 1 if ntp>0'
      inerr = 1
   end if
   if  (ntx < 4.and.init /= 3)  then
      write(6,'(/2x,a)') 'NTX / IREST inconsistency'
      inerr = 1
   end if
   if (ntb == 2 .and. ntp == 0) then
      write(6,'(/2x,a)') 'NTB set but no NTP option (must be 1 or 2)'
      inerr = 1
   end if
   if (ntp /= 0 .and. ntb /= 2) then
      write(6,'(/,a,a)')' NTP > 0 but not constant pressure P.B.C.', &
            ' (NTB = 2) must be used'
      inerr = 1
   end if
   if (ntb /= 0 .and. ifbox == 0 .and. ntp /= 0) then
      write(6,'(/,a)') ' (NTB /= 0 && NTP /= 0) but IFBOX == 0'
      write(6,'(/,a)') ' This combination is not supported'
      inerr = 1
   end if
   if (ntb /= 0 .and. &
         ( box(1) < 1.d0  .or. &
         box(2) < 1.d0  .or. &
         box(3) < 1.d0 ) ) then
      write(6,'(/,a,3f10.3)') ' BOX is too small: ',box(1),box(2),box(3)
      inerr = 1
   else if (ntb /= 0 .and. &
         (sqrt(cut) >= box(1)*0.5d0 .or. &
         sqrt(cut) >= box(2)*0.5d0 .or. &
         sqrt(cut) >= box(3)*0.5d0) ) then
      write(6,'(/,a)') ' CUT must be < half smallest box dimension'
      inerr = 1
   end if
   if (ntb /= 0 .and. igb > 0 ) then
      write(6,'(/,a)') ' igb>0 is only compatible with ntb=0'
      inerr = 1
   end if
   if ( ntb == 0 .and. sqrt(cut) < 8.05 .and. igb < 10 ) then
      write(6,'(/,a,f8.2)') ' unreasonably small cut for non-periodic run: ', &
         sqrt(cut)
      inerr = 1
   end if
   if ( rgbmax < 5.d0*fsmax ) then
      write(6,'(/,a,f8.2)') ' rgbmax must be at least ', 5.d0*fsmax
      inerr = 1
   end if
   if (icfe /= 0 .and. ibelly /= 0 ) then
      write(6,'(/,a)') ' ibelly cannot be used with icfe'
      inerr = 1
   end if
   if (klambda < 1 .or. klambda > 6 ) then
      write(6,'(/,a)') ' klambda must be between 1 and 6'
      inerr = 1
   end if
   if (clambda < 0.d0 .or. clambda > 1.d0 ) then
      write(6,'(/,a)') ' clambda must be between 0 and 1'
      inerr = 1
   end if
   
   if (idecomp > 0 .and. (igb >= 10)) then
      write(6,'(/,a)') 'IDECOMP is not compatible with EWALD or PB'
      inerr = 1
   end if
   if (idecomp > 0 .and. (ntr > 0 .or. ibelly > 0)) then
      write(6,'(/,a)') 'IDECOMP is not compatible with NTR or IBELLY'
      inerr = 1
   end if
   if (icnstph /= 0) then
      if (igb >= 10) then
         write(6, '(/,a)') 'Constant pH requires GB implicit solvent'
         inerr = 1
      end if
      if (icfe /= 0) then
         write(6, '(/,a)') &
         'Constant pH and thermodynamic integration are incompatable'
         inerr = 1
      end if
   end if
   
   !------------------------------------------------------------------------------
   !        ---WARNINGS:
   
#ifdef MPI
   
   !     -- Parallelization level
   
   if (plevel < 0 .or. plevel > 3) plevel = 1
   if (plevel == 0) then
      write (6, '(a)') &
            '| PLEVEL = 0: force parallelization only'
   else if (plevel == 1) then
      write (6, '(a)') &
            '| PLEVEL = 1: runmd parallelization, no EKCMR'
   else if (plevel == 2) then
      write (6, '(a)') &
            '| PLEVEL = 2: runmd parallelization with EKCMR'
   end if
#endif
   if (inerr == 1) then
      write(6, '(/,a)') ' *** input error(s)'
      call mexit(6,1)
   end if
   
   !     ----- LOAD THE CONSTRAINED ATOMS OR THE BELLY
   !           ATOMS. THESE ARE READ AS GROUPS -----
   
   ! used for RESTRAINED MD (NTR=1) AND ALSO TARGETED MD (ITGTMD=1)
   
   konst = ntr > 0
   belly = .false.
   natc = 0
   ngrp = 0
   natbel = 0
   nrc = 0
   if(konst) then
      if (ntrx <= 0) then
         call amopen(10,refc,'O','U','R')
      else
         call amopen(10,refc,'O','F','R')
      end if
      if (konst) write(6,9408)

      call rdrest(natom,ntx,ntrx,x(lcrdr))
      close(10)
      
         if( len_trim(restraintmask) == 0 ) then
            call rgroup(natom,natc,nres,ngrp,ix(i02),ih(m02), &
                  ih(m04),ih(m06),ih(m08),ix(icnstrgp),jgroup,index,npdec, &
                  x(l60),x(lcrdr),konst,.false.,belly,idecomp,5,.true.)
         else
            call atommask( natom, nres, 0, ih(m04), ih(m06), &
               ix(i02), ih(m02), x(lcrd), restraintmask, ix(icnstrgp) )

            ! for now, emulate the "GATHER ALL THE CONSTRAINED ATOMS TOGETHER"
            ! section of rgroup(); later, the various masks should be done
            ! differently, i.e. without the "gather", as in the following:
            !     x(l60:l60+natom-1) = restraint_wt
            !     natc = sum(ix(icnstrgp:icnstrgp+natom-1))

            natc = 0
            do i=1,natom
              if( ix(icnstrgp-1+i) <= 0 ) cycle
              natc = natc + 1
              ix(icnstrgp-1+natc) = i
              x(l60-1+natc) = restraint_wt
            end do
            write(6,'(a,a,a,i5,a)') '     Mask ', &
            restraintmask(1:len_trim(restraintmask)), ' matches ',natc,' atoms'
         end if
      nrc = natc
   end if  ! (konst)
   
   konst = .false.
   belly = ibelly > 0
   ngrp = 0
   if(belly) then
      write(6,9418)
      if( len_trim(bellymask) == 0 ) then
         call rgroup(natom,natbel,nres,ngrp,ix(i02),ih(m02), &
            ih(m04),ih(m06),ih(m08),ix(ibellygp), &
            jgroup,index,npdec, &
            x(l60),x(lcrdr),konst,.false.,belly,idecomp,5,.true.)
      else
         call atommask( natom, nres, 0, ih(m04), ih(m06), &
            ix(i02), ih(m02), x(lcrd), bellymask, ix(ibellygp) )
         natbel = sum(ix(ibellygp:ibellygp+natom-1))
         write(6,'(a,a,a,i5,a)') '     Mask ', &
            bellymask(1:len_trim(bellymask)), ' matches ',natbel,' atoms'
      end if
   end if

   call setvar(ix,belly)
   
   konst = .false.
   belly = .false.
   if(idecomp > 0) then
      write(6,9428)
      call rgroup(natom,ntmp,nres,ngrp,ix(i02),ih(m02), &
            ih(m04),ih(m06),ih(m08),ix(ibellygp), &
            jgroup,index,npdec, &
            x(l60),x(lcrdr),konst,.false.,belly,idecomp,5,.true.)
   end if
   
   if( ibelly > 0 .and. igb > 0 ) then
      
      !          ---here, the only allowable belly has just the first
      !             NATBEL atoms in the moving part.  Check to see that this
      !             requirement is satisfied:
      
      do i=natbel+1,natom
         if( ix(ibellygp+i-1) /= 0 ) then
            write(6,*) 'When igb>0, the moving part must be at the'
            write(6,*) '   start of the molecule.  This does not seem'
            write(6,*) '   to be the case here.'
            write(6,*) 'natbel,i,igroup(i) = ' &
                  ,natbel,i,ix(ibellygp+i-1)
            call mexit(6,1)
         end if
      end do
   end if
   
   !     ----- CALCULATE THE SQUARE OF THE BOND PARAMETERS FOR SHAKE
   !           THE PARAMETERS ARE PUT SEQUENTIALLY IN THE ARRAY CONP -----
   
   do i=1,nbonh + nbona + nbper
      j = ix(iicbh+i-1)
      x(l50+i-1) = req(j)**2
   end do
   
   !     rewind (5)

   !  DEBUG input; force checking
   call load_debug(5)

   return
   ! -------------------------------------------------------------------------
   ! Standard format statements:
   
   9328 format(/80('-')/,'   2.  CONTROL  DATA  FOR  THE  RUN',/80('-')/)
   9408 format(/4x,'LOADING THE CONSTRAINED ATOMS AS GROUPS',/)
   9409 format(/4x,'LOADING THE TARGETED MD ATOMS AS GROUPS',/)
   9418 format(/4x,'LOADING THE BELLY ATOMS AS GROUPS',/)
   9428 format(/4x,'LOADING THE DECOMP ATOMS AS GROUPS',/)
   9008 format(a80)
end subroutine mdread2 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Emit defined preprocessor names, ie, flags.
subroutine printflags()

   implicit none
   integer     max_line_length
   parameter ( max_line_length = 80 )

   character(len=max_line_length) line  ! output string of active flags
   integer n                            ! len(line)
   
   line = '| Flags:'
   n = 8
   
#ifdef ISTAR2
   call printflags2(' ISTAR2',7,n,line,.false.)
#endif
#ifdef SGIFFT
   call printflags2(' SGIFFT',7,n,line,.false.)
#endif
#ifdef MPI
   call printflags2(' MPI',4,n,line,.false.)
#endif
#ifdef noBTREE
   call printflags2(' noBTREE',8,n,line,.false.)
#endif
#ifdef NMODE
   call printflags2(' NMODE',6,n,line,.false.)
#endif
#ifdef HAS_10_12
   call printflags2(' HAS_10_12',10,n,line,.false.)
#endif
#ifdef NO_SIGN
   call printflags2(' NO_SIGN',8,n,line,.false.)
#endif
#ifdef CHARMM
   call printflags2(' CHARMM',7,n,line,.false.)
#endif
#ifdef DNA_SHIFT
   call printflags2(' DNA_SHIFT',10,n,line,.false.)
#endif
#ifdef CHARGE_MIXING
   call printflags2(' CHARGE_MIXING',14,n,line,.false.)
#endif
#ifdef NO_EWGRPS
   call printflags2(' NO_EWGRPS',10,n,line,.false.)
#endif
   
   call printflags2(' ',1,n,line,.true.)
   return
end subroutine printflags 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Primitive pre-Fortran90 implementation of printflags.
subroutine printflags2(flag,flag_len,line_len,line,last)

   implicit none
   integer     max_line_length
   parameter ( max_line_length = 80 )

   character(*) flag                ! flag name with blank prefix, intent(in)
   integer flag_len                 ! len(flag), intent(in)
   integer line_len                 ! len(line), intent(inout)
   character(len=max_line_length) line ! intent(inout)
   logical last                     ! is this the last flag ?, intent(in)

   if (line_len + flag_len > max_line_length) then
      write( 6,'(a)') line
      ! begin another line
      line = '| Flags:'
      line_len=8
   end if
   line=line(1:line_len) // flag(1:flag_len)
   line_len=line_len+flag_len
   if(last)write( 6,'(a)') line
   return
end subroutine printflags2 
!-------------------------------------------------


!     --- FLOAT_LEGAL_RANGE ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Check the range of a float; abort on illegal values.
subroutine float_legal_range(string,param,lo,hi)
   implicit none
   _REAL_ param,lo,hi
   character(len=*)string

   if ( param < lo .or. param > hi )then
      write(6,59)
      write(6,60)string,param
      write(6,61)
      write(6,62)lo,hi
      write(6,63)
      call mexit(6,1)
   end if
   59 format(/,1x,'Ewald PARAMETER RANGE CHECKING: ')
   60 format(1x,'parameter ',a,' has value ',e12.5)
   61 format(1x,'This is outside the legal range')
   62 format(1x,'Lower limit: ',e12.5,' Upper limit: ',e12.5)
   63 format(1x,'Check ew_legal.h')
   return
end subroutine float_legal_range 
!-------------------------------------------------


!     --- INT_LEGAL_RANGE ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Check the range of an integer; abort on illegal values.
subroutine int_legal_range(string,param,lo,hi)
   implicit none
   integer param,lo,hi
   character(len=*)string

   if ( param < lo .or. param > hi )then
      write(6,59)
      write(6,60)string,param
      write(6,61)
      write(6,62)lo,hi
      write(6,63)
      call mexit(6,1)
   end if
   59 format(/,1x,'PARAMETER RANGE CHECKING: ')
   60 format(1x,'parameter ',a,' has value ',i8)
   61 format(1x,'This is outside the legal range')
   62 format(1x,'Lower limit: ',i8,' Upper limit: ',i8)
   63 format(1x,'The limits may be adjustable; search in the .h files ')
   return
end subroutine int_legal_range 
!-------------------------------------------------


!     --- OPT_LEGAL_RANGE ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Check the range of an integer option; abort on illegal values.
subroutine opt_legal_range(string,param,lo,hi)
   implicit none
   integer param,lo,hi
   character(len=*)string

   if ( param < lo .or. param > hi )then
      write(6,59)
      write(6,60)string,param
      write(6,61)
      write(6,62)lo,hi
      write(6,63)
      call mexit(6,1)
   end if
   59 format(/,1x,'Ewald OPTION CHECKING: ')
   60 format(1x,'option ',a,' has value ',i5)
   61 format(1x,'This is outside the legal range')
   62 format(1x,'Lower limit: ',i5,' Upper limit: ',i5)
   63 format(1x,'Check the manual')
   return
end subroutine opt_legal_range 
!-------------------------------------------------


!     --- SANDER_BOMB ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Print an error message and quit
subroutine sander_bomb(routine,string1,string2)
   implicit none
   character(len=*) routine,string1,string2

   write(6, '(1x,2a)') &
         'SANDER BOMB in subroutine ', routine
   write(6, '(1x,a)') string1
   write(6, '(1x,a)') string2
   call mexit(6,1)
end subroutine sander_bomb
!-------------------------------------------------

