#include "copyright.h"
#include "dprec.h"


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Emit coordinates or velocities, r(istart:n), to unit nf.
subroutine corpac(r,istart,n,nf,loutfm)
   
   implicit none
   integer istart,n,nf
   _REAL_  r(n)
   logical loutfm  ! true for formatted output

   integer i,j
   logical three_digit_fractional_part
   _REAL_  rmax,rmin
   parameter (rmax = 9999.99d0)
   parameter (rmin = -999.99d0)
   
   if (istart > n) return
   
   if (.not.loutfm) then

      ! Unformatted writes:

      write(nf) (r(j),j=istart,n)

   else
      
      ! Formatted writes:
      
      three_digit_fractional_part = .true.
      do i = istart,n
         if (r(i) > rmax .or. r(i) < rmin) then
            three_digit_fractional_part = .false.
            exit
         end if
      end do
      
      if (three_digit_fractional_part) then
         write(nf,'(10f8.3)') (r(j),j=istart,n)
      else
         write(nf,'(10f8.2)') (r(j),j=istart,n)
      end if
   end if
   
   return
end subroutine corpac 
!-----------------------------------------------------------------------

#ifdef MPI

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine ekcmr here]
subroutine ekcmr(nspm,nsp,tma,ekcmt,xr,v,amass,istart,iend)
#else

subroutine ekcmr(nspm,nsp,tma,ekcmt,xr,v,amass)
#endif
   
   implicit _REAL_ (a-h,o-z)
   
   !     ----- ROUTINE TO CALCULATE THE TOTAL KINETIC ENERGY OF THE
   !           CENTER OF MASS OF THE SUB-MOLECULES AND ALSO THE
   !           COORDINATES OF THE MOLECULES RELATIVE TO THE CENTER OF
   !           MASS -----
   
#  include "box.h"
   
   dimension nsp(*),tma(*),ekcmt(*),xr(*),v(*),amass(*)
   dimension xcm(3),vcm(3),ekcml(3)
   data zero, one /0.0d0, 1.0d0/
   
   i3 = 0
   iat = 0

   ekcml(1) = zero
   ekcml(2) = zero
   ekcml(3) = zero
   
   
   do 900 n = 1,nspm
      nn = nsp(n)
      tmn = one/tma(n)
      
      !       ----- CALCULATE THE CENTER OF MASS AND THEN MOVE EACH
      !             SUB-MOLECULE TO ITS CENTER OF MASS -----
      
      j3 = i3
      xcm(1) = zero
      xcm(2) = zero
      xcm(3) = zero
      vcm(1) = zero
      vcm(2) = zero
      vcm(3) = zero
      do 220 j = 1,nn
         aamass = amass(iat+j)
         xcm(1) = xcm(1)+xr(j3+1)*aamass
         xcm(2) = xcm(2)+xr(j3+2)*aamass
         xcm(3) = xcm(3)+xr(j3+3)*aamass
#ifdef MPI
         if (iat+j >= istart .and. iat+j <= iend) then
#endif
            vcm(1) = vcm(1)+v(j3+1)*aamass
            vcm(2) = vcm(2)+v(j3+2)*aamass
            vcm(3) = vcm(3)+v(j3+3)*aamass
#ifdef MPI
         end if
#endif

         j3 = j3+3
      220 continue
      xcm(1) = xcm(1)*tmn
      xcm(2) = xcm(2)*tmn
      xcm(3) = xcm(3)*tmn
      
      j3 = i3
      do 240 j = 1,nn
         xr(j3+1) = xr(j3+1)-xcm(1)
         xr(j3+2) = xr(j3+2)-xcm(2)
         xr(j3+3) = xr(j3+3)-xcm(3)
         j3 = j3+3
      240 continue
      
      ekcml(1) = ekcml(1)+tmn*vcm(1)*vcm(1)
      ekcml(2) = ekcml(2)+tmn*vcm(2)*vcm(2)
      ekcml(3) = ekcml(3)+tmn*vcm(3)*vcm(3)
      
      !       ----- END OF CALCULATION FOR EACH SUB-MOLECULE -----
      
      i3 = j3
      iat = iat+nn
   900 continue
   
   ekcmt(1) = ekcml(1)
   ekcmt(2) = ekcml(2)
   ekcmt(3) = ekcml(3)
   
   return
end subroutine ekcmr 
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Open the coordinate, velocity, energy and constant pH output files.
subroutine open_dump_files

   implicit none
#  include "files.h"
#  include "md.h"
#  include "extra.h"
   
   !     subr amopen(lun,fname,fstat,fform,facc)

   if ( master ) then
      if (ioutfm <= 0) then
         
         !     ----- FORMATTED DUMPING -----
         
         if (ntwx > 0) then
            if (imin == 5) then
               call amopen(mdcrd_unit,mdcrd,'O','F','R')
               read(mdcrd_unit,1000) title
               write (6,1000) title
            else
               call amopen(mdcrd_unit,mdcrd,owrite,'F','W')
               write(mdcrd_unit,1000) title
            end if
         end if
         if (ntwv > 0) then
            call amopen(mdvel_unit,mdvel,owrite,'F','W')
            write(mdvel_unit,1000) title
         end if
         if (ntwe > 0) then
            call amopen(mden_unit,mden,owrite,'F','W')
            ! Probably the line below should not be commented, but now,
            ! since at least amber 6, it is a feature.
            !               write(MDEN_UNIT,1000) title
         end if
      else
         
         !     ----- UNFORMATTED DUMPING -----
         
         if (ntwx > 0) then
            call amopen(mdcrd_unit,mdcrd,owrite,'U','W')
            write(mdcrd_unit) title
         end if
         if (ntwv > 0) then
            call amopen(mdvel_unit,mdvel,owrite,'U','W')
            write(mdvel_unit) title
         end if
         if (ntwe > 0) then
            call amopen(mden_unit,mden,owrite,'U','W')
            write(mden_unit) title
         end if
      end if  ! (ioutfm <= 0)
      if (icnstph /= 0) then
         call amopen(CPOUT_UNIT, cpout, owrite, 'F', 'W')
      end if
   end if  ! ( master )
   1000 format(20a4)
   return
end subroutine open_dump_files 
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Close the coordinate, velocity, energy and constant pH output files.
subroutine close_dump_files

   implicit none
#  include "files.h"
#  include "extra.h"
#  include "md.h"   
   
   if ( master ) then
      if ( ntwx > 0 ) close( mdcrd_unit )
      if ( ntwv > 0 ) close( mdvel_unit )
      if ( ntwe > 0 ) close( mden_unit )
      if ( icnstph /= 0 ) close ( cpout_unit )
   end if
   return
end subroutine close_dump_files 
!----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Energy output for md, in human-readable form.
subroutine prntmd(nstep,nitp,nits,time,ener,fac,iout7,rms)
   
   implicit _REAL_ (a-h,o-z)
   logical rms   !  true if this is a print of RMS fluctuations

#  include "files.h"
#  include "md.h"
   dimension ener(*),fac(*)
   
   !     ----- DEFINE VARIOUS TERMS -----
   
   etot   = ener(1)
   ektot  = ener(2)
   temp   = ektot/fac(1)
   eksolt = ener(3)/fac(2)
   
   if(ntt == 5) then
      eksolv = ener(4)/fac(3)
   else
      rms_pbs = ener(4)
   end if
   
   scaltp = ener(5)
   
   boxx   = ener(7)
   boxy   = ener(8)
   boxz   = ener(9)
   volume = ener(10)
   densit = ener(42)
   
   presx  = ener(11)
   presy  = ener(12)
   presz  = ener(13)
   press  = ener(14)
   
   ekcmx  = ener(15)
   ekcmy  = ener(16)
   ekcmz  = ener(17)
   ekcmt  = ener(18)
   
   virx   = ener(19)
   viry   = ener(20)
   virz   = ener(21)
   virt   = ener(22)
   
   epot   = ener(23)
   enonb  = ener(24)
   eel    = ener(25)
   ehbond = ener(26)
   ebond  = ener(27)
   eangle = ener(28)
   edihed = ener(29)
   enb14  = ener(30)
   eel14  = ener(31)
   econst = ener(32)
   epol   = ener(33)
   aveper   = ener(34)
   aveind   = ener(35)
   avetot   = ener(36)
   esurf = ener(37)
   e3bod = ener(38)
   dvdl = ener(39)
   edisp = ener(40)
   virvsene = ener(43)
   diprms = ener(44)
   dipiter = ener(45)
   dipole_temp = ener(46)
   
   write(6,9018) nstep,time,temp,press
   write(6,9028) etot,ektot,epot
   write(6,9038) ebond,eangle,edihed
   write(6,9048) enb14,eel14,enonb
   if (igb < 10) then
      write(6,9059) eel,ehbond,econst
   else
      write(6,9060) eel,ehbond,econst
   endif
   if (gbsa > 0) write(6,9077) esurf
   if (igb >= 10) write(6,9074) esurf,edisp
   if (econst /= 0.0) write(6,9076) epot-econst
   if ( dvdl /= 0.d0) write(6,9089) dvdl
   if (rms .and. ntt==0 ) write(6,9075) rms_pbs
   if (volume /= 0.0) write(6,9078) ekcmt,virt,volume

   if (epol /= 0.0 .or. e3bod /= 0.0) &
         write(6,9070) epol,e3bod
   !     if (induced.gt.0) WRITE(6,9088) aveper,aveIND,avetot
   if (induced > 0.and.indmeth < 3) write(6,9090)diprms,dipiter
   if (induced > 0.and.indmeth == 3) write(6,9091)diprms, &
         dipole_temp
   if (volume /= 0.0) write(6,9079) densit

   write(6,8088)
   
   !     --- flush i/o buffer ---
   
   call amflsh(6)
   if (iout7 == 0) return
   
   !       ----- OUTPUT THE INFO FILE if requested -----
   
   write(7,9018) nstep,time,temp,press
   write(7,9028) etot,ektot,epot
   write(7,9038) ebond,eangle,edihed
   write(7,9048) enb14,eel14,enonb
   write(7,9059) eel,ehbond,econst
   if (gbsa > 0) write(7,9077) esurf
   if (igb >= 10) write(7,9074) esurf,edisp
   if (econst /= 0.0) write(7,9076) epot-econst
   if ( dvdl /= 0.d0) write(7,9089) dvdl
   if (rms .and. ntt==0 ) write(6,9075) rms_pbs
   if (volume /= 0.0) write(7,9078) ekcmt,virt,volume
   if (epol /= 0.0 .or. e3bod /= 0.0) &
         write(7,9070) epol,e3bod
   !     if (induced.gt.0) WRITE(7,9088) aveper,aveIND,avetot
   if (induced > 0.and.indmeth < 3) write(7,9090)diprms,dipiter
   if (induced > 0.and.indmeth == 3) write(7,9091)diprms, &
         dipole_temp
   if (volume /= 0.0) write(7,9079) densit
   
   8088 format(t2,78('-'),/)
   9018 format(/1x, 'NSTEP =',i9,3x,'TIME(PS) =',f12.3,2x, &
         'TEMP(K) =',f9.2,2x,'PRESS =',f8.1)
   9028 format (1x,'Etot   = ',f14.4,2x,'EKtot   = ',f14.4,2x, &
         'EPtot      = ',f14.4)
   9038 format (1x,'BOND   = ',f14.4,2x,'ANGLE   = ',f14.4,2x, &
         'DIHED      = ',f14.4)
   9048 format (1x,'1-4 NB = ',f14.4,2x,'1-4 EEL = ',f14.4,2x, &
         'VDWAALS    = ',f14.4)
   9058 format (1x,'EELEC  = ',f14.4,2x,'EHBOND  = ',f14.4,2x, &
         'RESTRAINT  = ',f14.4)
   9059 format (1x,'EELEC  = ',f14.4,2x,'EGB     = ',f14.4,2x, &
         'RESTRAINT  = ',f14.4)
   9060 format (1x,'EELEC  = ',f14.4,2x,'EPB     = ',f14.4,2x, &
         'RESTRAINT  = ',f14.4)
   9075 format ('|E(PBS) = ',f14.4)
   9076 format (1x,'EAMBER (non-restraint)  = ',f14.4)
   9077 format (1x,'ESURF= ',f14.4)
   9074 format (1x,'ECAVITY= ',f14.4,2x,'EDISPER = ',f14.4)
   9078 format (1x,'EKCMT  = ',f14.4,2x,'VIRIAL  = ',f14.4,2x, &
         'VOLUME     = ',f14.4)
   9079 format (52x,'Density    = ',f14.4)
   
   9070 format (1x,'EPOLZ  = ',f14.4,2x,'E3BODY  = ',f14.4)
   
   9088 format (1x,'DIPOLE MOMENTS/RESIDUE:',/ &
         ,1x,'PERMEN = ',f8.3,2x,'INDUCED = ',f8.3, &
         2x,'VECTOR SUM =',f8.3)
   9089 format (1x,'DV/DL  = ',f14.4)
   9090 format(1x,'Dipole convergence: rms = ',e10.3,' iters = ',f6.2)
   9091 format(1x,'Dipole convergence: rms = ',e10.3, &
         ' temperature = ',f6.2)
   9188 format (1x,'Ewald error estimate: ', e12.4)
   return
end subroutine prntmd 
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Assign velocities from a Maxwellian distribution.
subroutine setvel(nr,v,winv,tempi,init,iscale,scalm)
   
   implicit none
   _REAL_ v(*),winv(*),tempi,scalm
   integer nr,init,iscale

   integer i,j,nr3
   _REAL_ boltz,sd
   
   nr3 = 3*nr
   
   !     ----- Assign velocities from a Maxwellian distribution
   
   if (tempi < 1.d-6) then
      v(1:nr3+iscale) = 0.d0
      return
   end if
   
   boltz = 8.31441d-3*tempi/4.184d0
   i = 0
   do j=1,nr
      sd =  sqrt(boltz*winv(j))
      call gauss(0.d0,sd,v(i+1))
      call gauss(0.d0,sd,v(i+2))
      call gauss(0.d0,sd,v(i+3))
      i = i+3
   end do
   if (iscale > 0) then
      sd =  sqrt(boltz/scalm)
      do j=1,iscale
         call gauss(0.d0,sd,v(i+j))
      end do
   end if
   
   return
end subroutine setvel 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine mdeng here]
subroutine mdeng(nf,nstep,time,ener,fac,ntp)
   
   implicit none
#  include "box.h"
   integer nf,nstep,ntp,i
   _REAL_ ener(*),fac(*),time
   _REAL_ zero
   logical first
   character(len=16) labs(41)
   save first
   save labs
   data zero/ 0.d0 /
   data first/.true./
   data labs/'Nsteps  ','time(ps)  ','Etot  ','EKinetic  ', &
         'Temp  ','T_solute ','T_solv  ','Pres_scal_solu ', &
         'Pres_scal_solv ','BoxX  ','BoxY  ','BoxZ  ', &
         'volume  ','pres_X  ','pres_Y  ','pres_Z  ', &
         'Pressure ','EKCoM_x ','EKCoM_y ','EKCoM_z', &
         'EKComTot ','VIRIAL_x ','VIRIAL_y ','VIRIAL_z ', &
         'VIRIAL_tot ','E_pot  ','E_vdw  ','E_el  ', &
         'E_hbon  ','E_bon  ','E_angle  ','E_dih  ', &
         'E_14vdw  ','E_14el  ','E_const  ','E_pol  ', &
         'AV_permMoment ','AV_indMoment ','AV_totMoment ', &
         'Density', 'dV/dlambda'/
   
   !     ----- DEFINE VARIOUS TERMS -----
   
   if (first) then
      !       -- up to Ekinetic:
      write(nf,1) 'L0 ', (labs(i),i=1,4)
      !       -- up to Pres_scal_solu:
      write(nf,1) 'L1 ', (labs(i),i=5,8)
      !       -- up to boxZ:
      write(nf,1) 'L2 ', (labs(i),i=9,12)
      !       -- up to pres_Z:
      write(nf,1) 'L3 ', (labs(i),i=13,16)
      !       -- up to EKCoM_z:
      write(nf,1) 'L4 ', (labs(i),i=17,20)
      !       -- up to VIRIAL_z:
      write(nf,1) 'L5 ', (labs(i),i=21,24)
      !       -- up to E_el:
      write(nf,1) 'L6 ', (labs(i),i=25,28)
      !       -- up to E_dih:
      write(nf,1) 'L7 ', (labs(i),i=29,32)
      !       -- up to E_pol:
      write(nf,1) 'L8 ', (labs(i),i=33,36)
      !       -- up to Density:
      write(nf,1) 'L9 ', (labs(i),i=37,41)
      1 format(a,10(1x,a))
      first = .false.
   end if
   
   !     ----- write values for this step -----
   
   !     -- up to Ekinetic:
   write(nf, 2) 'L0 ', nstep, time, ener(1), ener(2)
   !     -- up to Pres_scal_solu:
   write(nf, 3) 'L1 ', ener(2)/fac(1), ener(3)/fac(2), &
         ener(4)/fac(3), ener(5)
   !     -- up to boxZ:
   write(nf, 3) 'L2 ', ener(6), box(1), box(2), box(3)
   !     -- up to pres_Z:
   write(nf, 3) 'L3 ', (ener(i), i=10,13)
   !     -- up to EKCoM_z:
   write(nf, 3) 'L4 ', (ener(i), i=14,17)
   !     -- up to VIRIAL_z:
   if( ntp > 0 ) then
      write(nf, 3) 'L5 ', (ener(i), i=18,21)
   else
      write(nf, 3) 'L5 ', ener(18),zero,zero,zero
   end if
   !     -- up to E_el:
   write(nf, 3) 'L6 ', (ener(i), i=22,25)
   !     -- up to E_dih:
   write(nf, 3) 'L7 ', (ener(i), i=26,29)
   !      -- up to E_pol:
   write(nf, 3) 'L8 ', (ener(i), i=30,33)
   !     -- up to Density:
   write(nf, 3) 'L9 ', (ener(i), i=34,36), ener(42), ener(39)
   
   2 format(a, i8, 20(2x,e16.10))
   3 format(a, 20(e16.10,2x))
   return
end subroutine mdeng 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine mden2 here]
subroutine mden2(nf,nstep,time,ener,fac)
   
   implicit _REAL_ (a-h,o-z)
   dimension ener(*),fac(*)
   
   !     ----- DEFINE VARIOUS TERMS -----
   
   write(nf,'(''Nsteps = '',i10)')nstep
   write(nf,'(''Time               = '',f20.10,'' p.s. '')')time
   write(nf,'(''Etotal             = '',f20.10)')ener(1)
   write(nf,'(''EKinetic           = '',f20.10)')ener(2)
   write(nf,'(''Temperature        = '',f20.10)')ener(2)/fac(1)
   write(nf,'(''Temp_solute        = '',f20.10)')ener(3)/fac(2)
   write(nf,'(''Temp_solvent       = '',f20.10)')ener(4)/fac(3)
   write(nf,'(''Press_SCALE_solute = '',f20.10)')ener(5)
   write(nf,'(''Press_SCALE_solvent= '',f20.10)')ener(6)
   write(nf,'(''BOX_x              = '',f20.10)')ener(7)
   write(nf,'(''BOX_y              = '',f20.10)')ener(8)
   write(nf,'(''BOX_z              = '',f20.10)')ener(9)
   write(nf,'(''VOLUME             = '',f20.10)')ener(10)
   write(nf,'(''PRES_x             = '',f20.10)')ener(11)
   write(nf,'(''PRES_y             = '',f20.10)')ener(12)
   write(nf,'(''PRES_z             = '',f20.10)')ener(13)
   write(nf,'(''PRESSURE           = '',f20.10)')ener(14)
   write(nf,'(''EKCoM_x            = '',f20.10)')ener(15)
   write(nf,'(''EKCoM_y            = '',f20.10)')ener(16)
   write(nf,'(''EKCoM_z            = '',f20.10)')ener(17)
   write(nf,'(''EKCoMTotal         = '',f20.10)')ener(18)
   write(nf,'(''VIRIAL_x           = '',f20.10)')ener(19)
   write(nf,'(''VIRIAL_y           = '',f20.10)')ener(20)
   write(nf,'(''VIRIAL_z           = '',f20.10)')ener(21)
   write(nf,'(''VIRIAL_Total       = '',f20.10)')ener(22)
   write(nf,'(''Epotential         = '',f20.10)')ener(23)
   write(nf,'(''vanderwaals        = '',f20.10)')ener(24)
   write(nf,'(''electrostatic      = '',f20.10)')ener(25)
   write(nf,'(''h-bond             = '',f20.10)')ener(26)
   write(nf,'(''bond               = '',f20.10)')ener(27)
   write(nf,'(''angle              = '',f20.10)')ener(28)
   write(nf,'(''dihedral           = '',f20.10)')ener(29)
   write(nf,'(''1-4v.d.w.          = '',f20.10)')ener(30)
   write(nf,'(''1-4electrostatic   = '',f20.10)')ener(31)
   write(nf,'(''restraint          = '',f20.10)')ener(32)
   write(nf,'(''polarization       = '',f20.10)')ener(33)
   write(nf,'(''ave perm moment    = '',f20.10)')ener(34)
   write(nf,'(''ave ind  moment    = '',f20.10)')ener(35)
   write(nf,'(''ave total moment   = '',f20.10)')ener(36)
   write(nf,'(''density            = '',f20.10)')ener(42)
   write(nf,'(''-----------------------------------------'')')
   
   return
end subroutine mden2 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine cenmas here]
subroutine cenmas(natom,x,v,tmass,tmassinv,amass, &
      ekcm,xcm,vcm,acm,ekrot,ocm,icm)
   implicit _REAL_ (a-h,o-z)
#  include "extra.h"
   
   !     ----- ROUTINE TO CALCULATE THE TRANSLATIONAL AND ROTATIONAL
   !           KINETIC ENERGIES AND VELOCITIES -----
   
   !     icm=0   return, do nothing
   !     icm=1   just return com coords in XCM
   !     icm=2   also calculate energy of com motion
   !     icm=3   also calcualte the angular momemtum
   !     icm=4   also calculate the inertia tensor
   
   dimension x(*),v(*),amass(*),xcm(*),vcm(*),acm(*),ocm(*)
   dimension tcm(3,3),lh(3),mh(3)
   
   data crit/1.0d-06/
   
   if (icm == 0) return
   
   i0 = 3*natom
   
   !     ----- CALCULATE THE CENTER OF MASS COORDINATES -----
   
   xcm(1) = 0.d0
   xcm(2) = 0.d0
   xcm(3) = 0.d0
   
   i = 0
   do j = 1,natom
      aamass = amass(j)
      xcm(1) = xcm(1) + x(i+1)*aamass
      xcm(2) = xcm(2) + x(i+2)*aamass
      xcm(3) = xcm(3) + x(i+3)*aamass
      i = i + 3
   end do
   xcm(1) = xcm(1) * tmassinv
   xcm(2) = xcm(2) * tmassinv
   xcm(3) = xcm(3) * tmassinv
   
   if (iabs(icm) < 2) return
   
   !     ----- CALCULATE THE VELOCITY AND THE TRANSLATIONAL
   !           KINETIC ENERGY OF THE CENTRE OF MASS -----
   
   ekcm = 0.d0
   vcm(1) = 0.0d0
   vcm(2) = 0.0d0
   vcm(3) = 0.0d0
   
   i = 0
   do j = 1,natom
      aamass = amass(j)
      vcm(1) = vcm(1) + v(i+1)*aamass
      vcm(2) = vcm(2) + v(i+2)*aamass
      vcm(3) = vcm(3) + v(i+3)*aamass
      i = i + 3
   end do
   do i=1,3
      vcm(i) = vcm(i) * tmassinv
      ekcm = ekcm + vcm(i)*vcm(i)
   end do
   ekcm = ekcm * tmass * 0.5d0
   comvel = sqrt(vcm(1)*vcm(1)+vcm(2)*vcm(2)+vcm(3)*vcm(3))
   
   if(iabs(icm) < 3) return
   
   !     ----- CALCULATE THE ANGULAR MOMENTUM ABOUT THE
   !           CENTER OF MASS ----
   
   acm(1) = 0.0d0
   acm(2) = 0.0d0
   acm(3) = 0.0d0

   i = 0
   do j = 1,natom
      aamass = amass(j)
      acm(1) = acm(1) + (x(i+2)*v(i+3)-x(i+3)*v(i+2)) * aamass
      acm(2) = acm(2) + (x(i+3)*v(i+1)-x(i+1)*v(i+3)) * aamass
      acm(3) = acm(3) + (x(i+1)*v(i+2)-x(i+2)*v(i+1)) * aamass
      i = i+3
   end do

   acm(1) = acm(1) - (xcm(2)*vcm(3)-xcm(3)*vcm(2)) * tmass
   acm(2) = acm(2) - (xcm(3)*vcm(1)-xcm(1)*vcm(3)) * tmass
   acm(3) = acm(3) - (xcm(1)*vcm(2)-xcm(2)*vcm(1)) * tmass

   if (iabs(icm) < 4) return
   
   !     ----- CALCULATE THE INERTIA TENSOR -----
   
   xx = 0.d0
   xy = 0.d0
   xz = 0.d0
   yy = 0.d0
   yz = 0.d0
   zz = 0.d0

   i = 0
   do j = 1,natom
      x1 = x(i+1)-xcm(1)
      x2 = x(i+2)-xcm(2)
      x3 = x(i+3)-xcm(3)
      aamass = amass(j)
      xx = xx+x1*x1*aamass
      xy = xy+x1*x2*aamass
      xz = xz+x1*x3*aamass
      yy = yy+x2*x2*aamass
      yz = yz+x2*x3*aamass
      zz = zz+x3*x3*aamass
      i = i+3
   end do
   tcm(1,1) = yy+zz
   tcm(2,1) = -xy
   tcm(3,1) = -xz
   tcm(1,2) = -xy
   tcm(2,2) = xx+zz
   tcm(3,2) = -yz
   tcm(1,3) = -xz
   tcm(2,3) = -yz
   tcm(3,3) = xx+yy
   
   !     ----- INVERT THE INERTIA TENSOR -----
   
   call matinv(tcm,3,d,lh,mh)
   if(abs(d) <= crit) then
      write(6,307)
      307 format(/5x,'%CENMAS-F-INERTIA_TENSOR,  determinant', &
            ' is zero ... stop')
      call mexit(6,1)
   end if
   
   !     ----- CALCULATE THE ANGULAR VELOCITY ABOUT THE CENTER OF
   !           MASS AND THE ROTATIONAL KINETIC ENERGY -----
   
   ekrot = 0.d0
   do m = 1,3
      ocm(m) = 0.d0
      do n = 1,3
        ocm(m) = ocm(m)+tcm(m,n)*acm(n)
      end do
      ekrot = ekrot+ocm(m)*acm(m)
   end do
   ekrot = ekrot * 0.5d0

   if (icm < 0) return
   if (master) then
      write(6,'(/3x,a,f11.4,3x,a,f11.4,3x,a,f12.6)') 'KE Trans =', &
            ekcm, 'KE Rot =', ekrot,'C.O.M. Vel =',comvel
      call amflsh(6)
   end if
   return
end subroutine cenmas 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Remove Center of Mass transrotational motion.
subroutine stopcm(nr,x,v,xcm,vcm,ocm)
   
   implicit none
   integer nr
   _REAL_  x(*),v(*),xcm(*),vcm(*),ocm(*)

#ifdef MPI
#  include "parallel.h"
#endif
#  include "extra.h"
   
   integer i, j, m
   _REAL_  x1, x2, x3

   !     ----- STOP THE CENTER OF MASS TRANSLATION -----
   
   i = 0
   do j = 1,nr
      do m = 1,3
         i = i+1
         v(i) = v(i)-vcm(m)
      end do
   end do
   
   !     ----- STOP THE ROTATION ABOUT THE CENTER OF MASS -----
   
   i = 0
   do j = 1,nr
      x1 = x(i+1)-xcm(1)
      x2 = x(i+2)-xcm(2)
      x3 = x(i+3)-xcm(3)
      v(i+1) = v(i+1)-ocm(2)*x3+ocm(3)*x2
      v(i+2) = v(i+2)-ocm(3)*x1+ocm(1)*x3
      v(i+3) = v(i+3)-ocm(1)*x2+ocm(2)*x1
      i = i+3
   end do
   if (master) write(6,9008)
   9008 format(/3x,'Translational and rotational motion removed')
   return
end subroutine stopcm 
