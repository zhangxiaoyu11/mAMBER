! <compile=optimized>
#include "copyright.h"
#include "dprec.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine bond here]
subroutine bond(nbin,ib,jb,icb,x,xx,ix,f,eb,nocrst)
   

   !     ----- ROUTINE TO GET BOND ENERGY AND FORCES FOR THE POTENTIAL
   !           OF CB*(B-B0)**2
   !------------------------------------------------------------------------

   use decomp, only : decpair
   implicit none
#ifdef MPI
#  include "parallel.h"
#endif
#  include "md.h"
   
   !-------------passed-in variables  ---------------------------
   integer nbin
   integer ib(*),jb(*),icb(*),ix(*)
   _REAL_ x(*),xx(*),f(*),eb
   logical nocrst
   
#  include "box.h"
#  include "parms.h"
#  include "flocntrl.h"
   
   !------------- local variables  ---------------------------
   integer max190
   parameter (max190=190)
   integer nb,ist,maxlen,ksum,i3,j3,ic,jn,ii,jj
   _REAL_ xij(max190),yij(max190),zij(max190),eaw(max190), &
         rij(max190),dfw(max190)
   _REAL_ da,df,xa,ya,za,rij0,ebl
   logical skip
   integer piece,start,end,newnb
#ifndef MPI
   integer numtasks,mytaskid
   
   numtasks=1
   mytaskid=0
#endif
   
   if ( do_bond == 0 .or. nbin == 0 ) return
   
   if(numtasks > 1)then
      piece = nbin/numtasks
      start = mytaskid*piece+1
      end   = mytaskid*piece+piece
      if(mytaskid == (numtasks-1)) end = nbin
   else
      start=1
      end=nbin
   end if
   nb = end
   ist = start -1
   
   ebl = 0.0d+00
#ifdef PSANDER
   ist=0
   nb=nbin
#endif
   
   !     ----- GRAND LOOP FOR THE bond STUFF -----
   
   4200 continue
   maxlen = max190
   skip = (ist+maxlen) > nb
   if (skip) maxlen = nb-ist
   if (maxlen > 0) then
      
      do jn = 1,maxlen
         i3 = ib(jn+ist)
         j3 = jb(jn+ist)
         
         !           ----- CALCULATION OF THE bond vector -----
         
         xij(jn) = x(i3+1)-x(j3+1)
         yij(jn) = x(i3+2)-x(j3+2)
         zij(jn) = x(i3+3)-x(j3+3)
      end do
      
      do jn = 1,maxlen
         rij0 = xij(jn)*xij(jn)+yij(jn)*yij(jn)+zij(jn)*zij(jn)
         rij(jn) = sqrt(rij0)
      end do
      
      !         ----- CALCULATION OF THE ENERGY AND DER -----
      
      do jn = 1,maxlen
         ic = icb(jn+ist)
         rij0 = rij(jn)
         da = rij0-req(ic)
         !                                 for rms deviation from ideal bonds:
         df = rk(ic)*da
         eaw(jn) = df*da
         if(idecomp == 1 .or. idecomp == 2) then
            ii = (ib(jn+ist) + 3)/3
            jj = (jb(jn+ist) + 3)/3
            call decpair(4,ii,jj,eaw(jn))
         end if
         dfw(jn) = (df+df)/rij0
      end do
      
      !         ----- CALCULATION OF THE FORCE -----
      
      do jn = 1,maxlen
         i3 = ib(jn+ist)
         j3 = jb(jn+ist)
         df = dfw(jn)
         xa = df*xij(jn)
         ya = df*yij(jn)
         za = df*zij(jn)
         f(i3+1) = f(i3+1)-xa
         f(i3+2) = f(i3+2)-ya
         f(i3+3) = f(i3+3)-za
         f(j3+1) = f(j3+1)+xa
         f(j3+2) = f(j3+2)+ya
         f(j3+3) = f(j3+3)+za
      end do
      do ksum = 1, maxlen
         ebl = ebl + eaw(ksum)
      end do
      ist = ist+maxlen
   end if  ! (maxlen > 0)
   if(.not.skip) goto 4200
   
   !     ----- ALL DONE -----
   
   eb = ebl
   return
end subroutine bond 
!-----------------------------------------------------------------------

subroutine bondcp(nbin,ib,jb,icb,x,xx,ix,f,eb,nocrst,iqmshk)
   
   !************************************************************************
   !                              AMBER                                   **
   !                                                                      **
   !               Copyright (c) 1986, 1991, 1995, 1997, 1999             **
   !                Regents of the University of California               **
   !                       All Rights Reserved.                           **
   !                                                                      **
   !  This software provided pursuant to a license agreement containing   **
   !  restrictions on its disclosure, duplication, and use. This software **
   !  contains confidential and proprietary information, and may not be   **
   !  extracted or distributed, in whole or in part, for any purpose      **
   !  whatsoever, without the express written permission of the authors.  **
   !  This notice, and the associated author list, must be attached to    **
   !  all copies, or extracts, of this software. Any additional           **
   !  restrictions set forth in the license agreement also apply to this  **
   !  software.                                                           **
   !************************************************************************
   
   !     ----- ROUTINE TO GET BOND ENERGY AND FORCES FOR THE POTENTIAL
   !           OF CB*(B-B0)**2
   !------------------------------------------------------------------------

   use decomp, only : decpair
   implicit none
#ifdef MPI
#  include "parallel.h"
#endif
#  include "md.h"

   !-------------passed-in variables  ---------------------------
   integer nbin
   integer ib(*),jb(*),icb(*),ix(*)
   integer iqmshk(*)
   _REAL_ x(*),xx(*),f(*),eb
   logical nocrst

#  include "box.h"
#  include "parms.h"
#  include "flocntrl.h"

   !------------- local variables  ---------------------------
   integer max190
   parameter (max190=190)
   integer nb,ist,maxlen,ksum,i3,j3,ic,jn,ii,jj
   _REAL_ xij(max190),yij(max190),zij(max190),eaw(max190), &
         rij(max190),dfw(max190)
   _REAL_ da,df,xa,ya,za,rij0,ebl
   logical skip
   integer piece,start,end,newnb
#ifndef MPI
   integer numtasks,mytaskid

   numtasks=1
   mytaskid=0
#endif

   if ( do_bond == 0 .or. nbin == 0 ) return

   if(numtasks > 1)then
      piece = nbin/numtasks
      start = mytaskid*piece+1
      end   = mytaskid*piece+piece
      if(mytaskid == (numtasks-1)) end = nbin
   else
      start=1
      end=nbin
   end if
   nb = end
   ist = start -1

   ebl = 0.0d+00
#ifdef PSANDER
   ist=0
   nb=nbin
#endif

   !     ----- GRAND LOOP FOR THE bond STUFF -----

   4200 continue
   maxlen = max190
   skip = (ist+maxlen) > nb
   if (skip) maxlen = nb-ist
   if (maxlen > 0) then

      do jn = 1,maxlen
         if(iqmshk(jn+ist).eq.1) exit
         i3 = ib(jn+ist)
         j3 = jb(jn+ist)

         !           ----- CALCULATION OF THE bond vector -----

         xij(jn) = x(i3+1)-x(j3+1)
         yij(jn) = x(i3+2)-x(j3+2)
         zij(jn) = x(i3+3)-x(j3+3)
      end do

      do jn = 1,maxlen
         if(iqmshk(jn+ist).eq.1) exit
         rij0 = xij(jn)*xij(jn)+yij(jn)*yij(jn)+zij(jn)*zij(jn)
         rij(jn) = sqrt(rij0)
      end do

      !         ----- CALCULATION OF THE ENERGY AND DER -----

      do jn = 1,maxlen
         if(iqmshk(jn+ist).eq.1) exit
         ic = icb(jn+ist)
         rij0 = rij(jn)
         da = rij0-req(ic)
         !                                 for rms deviation from ideal bonds:
         df = rk(ic)*da
         eaw(jn) = df*da
         if(idecomp == 1 .or. idecomp == 2) then
            ii = (ib(jn+ist) + 3)/3
            jj = (jb(jn+ist) + 3)/3
            call decpair(4,ii,jj,eaw(jn))
         end if
         dfw(jn) = (df+df)/rij0
      end do

      !         ----- CALCULATION OF THE FORCE -----

      do jn = 1,maxlen
         if(iqmshk(jn+ist).eq.1) exit
         i3 = ib(jn+ist)
         j3 = jb(jn+ist)
         df = dfw(jn)
         xa = df*xij(jn)
         ya = df*yij(jn)
         za = df*zij(jn)
         f(i3+1) = f(i3+1)-xa
         f(i3+2) = f(i3+2)-ya
         f(i3+3) = f(i3+3)-za
         f(j3+1) = f(j3+1)+xa
         f(j3+2) = f(j3+2)+ya
         f(j3+3) = f(j3+3)+za
      end do
      do ksum = 1, maxlen
         if(iqmshk(ksum+ist).eq.1) exit
         ebl = ebl + eaw(ksum)
      end do
      ist = ist+maxlen
   end if  ! (maxlen > 0)
   if(.not.skip) goto 4200

   !     ----- ALL DONE -----

   eb = ebl
   return
end subroutine bondcp 
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine angl here]
subroutine angl(nbain,it,jt,kt,ict,x,xx,ix,f,eba)

   use decomp, only : decpair, decangle
   implicit none

#ifdef MPI
#  include "parallel.h"
#endif
   !-------------passed-in variables  ---------------------------
   logical nocrst
   _REAL_ x(*),xx(*),f(*),eba
   integer ix(*),it(*),jt(*),kt(*),ict(*),nbain
   !------------- local variables  ---------------------------
   logical skip

   !     ----- ROUTINE TO GET THE ANGLE ENERGIES AND FORCES FOR THE
   !           POTENTIAL OF THE TYPE CT*(T-T0)**2

#  include "box.h"
#  include "parms.h"
#  include "md.h"
#  include "flocntrl.h"
#ifdef CHARMM
   dimension ebw(190)
#endif
   integer max190
   parameter (max190=190)
   _REAL_ dt1,dt2,dt3,dt4,dt5,dt6,dt7,dt8,dt9
   _REAL_ sth,cik,cii,ckk,da,df,ant0,ct0,ct1,ct2,st
   _REAL_ rij0,rik0,rkj0,ebal,ecnl,eub,ecn
   integer nba,maxlen,lenc,ksum,i3,j3,k3,ic,istc,jn,ist
   _REAL_ xij(190),yij(190),zij(190),xkj(190),ykj(190), &
         zkj(190),cst(190),eaw(190),rij(190),rkj(190),rik(190), &
         dfw(190),ant(190)

   ! previously declared, see passed-in variables.  srb
   !     DIMENSION IT(*),JT(*),KT(*),ICT(*),X(*),XX(*),IX(*),F(*)
   integer piece,start,end,nbtmp,newnb
   integer ii,jj,kk
   _REAL_ pt999
   data pt999 /0.9990d0/

   if ( do_angle == 0 .or. nbain == 0 )return

   nba=nbain
#ifdef MPI
   piece = nba/numtasks
   start = mytaskid*piece+1
   end   = mytaskid*piece+piece
   if(mytaskid == (numtasks-1)) end = nba
   nbtmp = nba
   !     JV Use actual count this PE will do, reset at end of routine
   nba = end
   ist = start -1
#else
   ! FLOW CONTROL FLAG (debug)
   ist = 0
#endif

   ebal = 0.0d0
   ecnl = 0.0d0
   eub  = 0.0d0
#ifdef PSANDER
   ist=0
   nba=nbain
#endif

   !     ----- GRAND LOOP FOR THE angle STUFF -----

   4200 continue
   maxlen = max190
   skip = (ist+maxlen) > nba
   if (skip) maxlen = nba-ist
   if (maxlen <= 0) goto 220

   do jn = 1,maxlen
      i3 = it(jn+ist)
      j3 = jt(jn+ist)
      k3 = kt(jn+ist)

      !           ----- CALCULATION OF THE angle -----

      xij(jn) = x(i3+1)-x(j3+1)
      yij(jn) = x(i3+2)-x(j3+2)
      zij(jn) = x(i3+3)-x(j3+3)
      xkj(jn) = x(k3+1)-x(j3+1)
      ykj(jn) = x(k3+2)-x(j3+2)
      zkj(jn) = x(k3+3)-x(j3+3)
   end do

   do jn = 1,maxlen
      rij0 = xij(jn)*xij(jn)+yij(jn)*yij(jn)+zij(jn)*zij(jn)
      rkj0 = xkj(jn)*xkj(jn)+ykj(jn)*ykj(jn)+zkj(jn)*zkj(jn)
      rik0 = sqrt(rij0*rkj0)
      ct0 = (xij(jn)*xkj(jn)+yij(jn)*ykj(jn)+zij(jn)*zkj(jn))/rik0
      ct1 = max(-pt999,ct0)
      ct2 = min(pt999,ct1)
      cst(jn) = ct2
      ant(jn) = acos(ct2)
      rij(jn) = rij0
      rkj(jn) = rkj0
      rik(jn) = rik0
   end do

   !         ----- CALCULATION OF THE ENERGY AND DER -----

   do jn = 1,maxlen
      ic = ict(jn+ist)
      ant0 = ant(jn)
      da = ant0-teq(ic)
      !                                 for rms deviation from ideal angles:
      df = tk(ic)*da
      eaw(jn) = df*da
      if(idecomp == 1 .or. idecomp == 2) then
         ii = (it(jn+ist) + 3)/3
         jj = (jt(jn+ist) + 3)/3
         kk = (kt(jn+ist) + 3)/3
         call decangle(ii,jj,kk,eaw(jn))
      end if
      dfw(jn) = -(df+df)/sin(ant0)
   end do

   !         ----- CALCULATION OF THE FORCE -----

   do jn = 1,maxlen
      i3 = it(jn+ist)
      j3 = jt(jn+ist)
      k3 = kt(jn+ist)

      st = dfw(jn)
      sth = st*cst(jn)
      cik = st/rik(jn)
      cii = sth/rij(jn)
      ckk = sth/rkj(jn)
      dt1 = cik*xkj(jn)-cii*xij(jn)
      dt2 = cik*ykj(jn)-cii*yij(jn)
      dt3 = cik*zkj(jn)-cii*zij(jn)
      dt7 = cik*xij(jn)-ckk*xkj(jn)
      dt8 = cik*yij(jn)-ckk*ykj(jn)
      dt9 = cik*zij(jn)-ckk*zkj(jn)
      dt4 = -dt1-dt7
      dt5 = -dt2-dt8
      dt6 = -dt3-dt9

      f(i3+1) = f(i3+1)-dt1
      f(i3+2) = f(i3+2)-dt2
      f(i3+3) = f(i3+3)-dt3
      f(j3+1) = f(j3+1)-dt4
      f(j3+2) = f(j3+2)-dt5
      f(j3+3) = f(j3+3)-dt6
      f(k3+1) = f(k3+1)-dt7
      f(k3+2) = f(k3+2)-dt8
      f(k3+3) = f(k3+3)-dt9
   end do
#ifdef CHARMM

   !         --- include Urey-Bradley terms:

   do jn = 1,maxlen
      i3 = it(jn+ist)
      j3 = kt(jn+ist)
      xij(jn) = x(i3+1)-x(j3+1)
      yij(jn) = x(i3+2)-x(j3+2)
      zij(jn) = x(i3+3)-x(j3+3)
   end do

   do jn = 1,maxlen
      rij0 = xij(jn)*xij(jn)+yij(jn)*yij(jn)+zij(jn)*zij(jn)
      rij(jn) = sqrt(rij0)
   end do

   !         ----- CALCULATION OF THE ENERGY AND DER -----

   do jn = 1,maxlen
      ic = ict(jn+ist)
      rij0 = rij(jn)
      da = rij0-rub(ic)
      df = rkub(ic)*da
      ebw(jn) = df*da
      if(idecomp == 1 .or. idecomp == 2) then
         ii = (it(jn+ist) + 3) / 3
         kk = (kt(jn+ist) + 3) / 3
         call decpair(4,ii,kk,ebw(jn))
      end if
      dfw(jn) = (df+df)/rij0
   end do

   !         ----- CALCULATION OF THE FORCE -----

   do jn = 1,maxlen
      i3 = it(jn+ist)
      j3 = kt(jn+ist)
      df = dfw(jn)
      xa = df*xij(jn)
      ya = df*yij(jn)
      za = df*zij(jn)
      f(i3+1) = f(i3+1)-xa
      f(i3+2) = f(i3+2)-ya
      f(i3+3) = f(i3+3)-za
      f(j3+1) = f(j3+1)+xa
      f(j3+2) = f(j3+2)+ya
      f(j3+3) = f(j3+3)+za
   end do
   do ksum = 1, maxlen
      eub = eub + ebw(ksum)
   end do
#endif
   do ksum=1, maxlen
      ebal = ebal + eaw(ksum)
   end do
   ist = ist+maxlen

   220 if(.not.skip) goto 4200
   eba = ebal
#ifdef CHARMM
   !     write(6,*) 'Urey-Bradley energy: ',eub
   eba = eba + eub
#endif
#ifdef MPI
   nba = nbtmp
#endif
   return
end subroutine angl
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine ephi here]
subroutine ephi(nphiin,ip,jp,kp,lp,icp,cg,iac,x,xx,ix,f,dvdl, &
      ep,enbp,eelp,dcharge)

   use decomp, only : decpair, decphi
   implicit none

#ifdef MPI
#  include "parallel.h"
#else
   integer numtasks,mytaskid
   parameter (numtasks=1,mytaskid=0)
#endif
   _REAL_ lfac
#  include "parms.h"
#  include "box.h"
#  include "md.h"
#  include "flocntrl.h"

   !-------------passed-in variables  ---------------------------
   logical nocrst
   integer ix(*),ip(*),jp(*),kp(*),lp(*),icp(*),iac(*)
   _REAL_   ep,enbp,eelp,ecn
   _REAL_ cg(*),dcharge(*),dvdl
   _REAL_ x(*),xx(*),f(*),eba
   integer nphiin,mba,mphi

   !------------- local variables  ---------------------------
   logical skip

   _REAL_ intdieli
   integer max190,maxlen
   parameter (max190=190)
   _REAL_ xij(max190),yij(max190),zij(max190),xkj(max190), &
         ykj(max190),zkj(max190),xkl(max190),ykl(max190), &
         zkl(max190),dx(max190),dy(max190),dz(max190),gx(max190), &
         gy(max190),gz(max190),ct(max190),cphi(max190),sphi(max190), &
         z1(max190),z2(max190),fxi(max190),fyi(max190),fzi(max190), &
         fxj(max190),fyj(max190),fzj(max190),fxk(max190),fyk(max190), &
         fzk(max190),fxl(max190),fyl(max190),fzl(max190),df(max190)

   _REAL_ xa,ya,za,f1,f2,r1,r2,r6,r12,dfn,g,fmuln
   _REAL_ dr1,dr2,dr3,dr4,dr5,dr6,drx,dry,drz
   _REAL_ dc1,dc2,dc3,dc4,dc5,dc6
   _REAL_ df0,df1,dflim,dums,cosnp,sinnp
   _REAL_ z10,z20,z12,z11,z22,ap0,ap1,ct0,ct1,s,ftem
   _REAL_ epl,ecnl,enbpl,eelpl,scnb0,scee0
   integer ksum,lenc,ia1,ia2,ibig,isml,ii,jj,kk,ll
   integer ic,inc,ic0,iduml,idumi,kdiv
   integer i3,j3,k3,l3,k3t,l3t
   integer jn,istc,ist,nphi
   integer eedmeth


   _REAL_ epw(max190)
   _REAL_ gmul(10)


   data gmul/0.0d+00,2.0d+00,0.0d+00,4.0d+00,0.0d+00,6.0d+00, &
         0.0d+00,8.0d+00,0.0d+00,10.0d+00/
   _REAL_ tm24,tm06,tenm3
   data tm24,tm06,tenm3/1.0d-18,1.0d-06,1.0d-03/
   _REAL_ pi
   data pi/3.141592653589793d+00/
   _REAL_ zero,one,two,four,six,twelve
   data zero,one,two,four,six,twelve/0.d0,1.d0,2.d0,4.d0,6.d0,12.d0/

   !     ---- ARRAYS GAMC = PK*COS(PHASE) AND GAMS = PK*SIN(PHASE) ----

   integer piece,start,end,nbtmp,newnb

   eedmeth = 4
   if ( do_ephi == 0 .or. nphiin == 0 )return
   nbtmp=nphiin
   if(numtasks > 1)then
      piece = nphiin/numtasks
      start = mytaskid*piece+1
      end   = mytaskid*piece+piece
      if(mytaskid == (numtasks-1)) end = nphiin
   else
      start=1
      end=nphiin
   end if
   nphi = end
   ist = start -1
   intdieli = 1.d0/intdiel
   epl = zero
   ecnl = zero
   enbpl = zero
   eelpl = zero
   scnb0 = one/scnb
   scee0 = one/scee

   !     ----- GRAND LOOP FOR THE DIHEDRAL STUFF -----

   4200 continue
   maxlen = 190
   skip = (ist+maxlen) > nphi
   if(skip) maxlen = nphi-ist
   if(maxlen <= 0) goto 820

   do jn = 1,maxlen
      i3 = ip(jn+ist)
      j3 = jp(jn+ist)
      k3t = kp(jn+ist)
      l3t = lp(jn+ist)
      k3 = iabs(k3t)
      l3 = iabs(l3t)

      !           ----- CALCULATION OF ij, kj, kl VECTORS -----

      xij(jn) = x(i3+1)-x(j3+1)
      yij(jn) = x(i3+2)-x(j3+2)
      zij(jn) = x(i3+3)-x(j3+3)
      xkj(jn) = x(k3+1)-x(j3+1)
      ykj(jn) = x(k3+2)-x(j3+2)
      zkj(jn) = x(k3+3)-x(j3+3)
      xkl(jn) = x(k3+1)-x(l3+1)
      ykl(jn) = x(k3+2)-x(l3+2)
      zkl(jn) = x(k3+3)-x(l3+3)
   end do

   !         ----- GET THE NORMAL VECTOR -----

   do jn = 1,maxlen
      dx(jn) = yij(jn)*zkj(jn)-zij(jn)*ykj(jn)
      dy(jn) = zij(jn)*xkj(jn)-xij(jn)*zkj(jn)
      dz(jn) = xij(jn)*ykj(jn)-yij(jn)*xkj(jn)
      gx(jn) = zkj(jn)*ykl(jn)-ykj(jn)*zkl(jn)
      gy(jn) = xkj(jn)*zkl(jn)-zkj(jn)*xkl(jn)
      gz(jn) = ykj(jn)*xkl(jn)-xkj(jn)*ykl(jn)
   end do

   do jn = 1,maxlen
      fxi(jn) = sqrt(dx(jn)*dx(jn) &
            +dy(jn)*dy(jn) &
            +dz(jn)*dz(jn)+tm24)
      fyi(jn) = sqrt(gx(jn)*gx(jn) &
            +gy(jn)*gy(jn) &
            +gz(jn)*gz(jn)+tm24)
      ct(jn) = dx(jn)*gx(jn)+dy(jn)*gy(jn)+dz(jn)*gz(jn)
   end do

   !         ----- BRANCH IF LINEAR DIHEDRAL -----

   do jn = 1,maxlen
#ifdef CRAY_PVP
      bit = one/fxi(jn)
      bik = one/fyi(jn)
      z10 = cvmgt(zero,bit,tenm3 > fxi(jn))
      z20 = cvmgt(zero,bik,tenm3 > fyi(jn))
#else
      z10 = one/fxi(jn)
      z20 = one/fyi(jn)
      if (tenm3 > fxi(jn)) z10 = zero
      if (tenm3 > fyi(jn)) z20 = zero
#endif
      z12 = z10*z20
      z1(jn) = z10
      z2(jn) = z20
#ifdef CRAY_PVP
      ftem = cvmgz(zero,one,z12)
#else
      ftem = zero
      if (z12 /= zero) ftem = one
#endif
      fzi(jn) = ftem
      ct0 = min(one,ct(jn)*z12)
      ct1 = max(-one,ct0)
      s = xkj(jn)*(dz(jn)*gy(jn)-dy(jn)*gz(jn))+ &
            ykj(jn)*(dx(jn)*gz(jn)-dz(jn)*gx(jn))+ &
            zkj(jn)*(dy(jn)*gx(jn)-dx(jn)*gy(jn))
      ap0 = acos(ct1)
      ap1 = pi-sign(ap0,s)
      ct(jn) = ap1
      cphi(jn) = cos(ap1)
      sphi(jn) = sin(ap1)
   end do

   !         ----- CALCULATE THE ENERGY AND THE DERIVATIVES WITH RESPECT TO
   !               COSPHI -----

   do jn = 1,maxlen
      ic = icp(jn+ist)
      inc = ipn(ic)
      ct0 = pn(ic)*ct(jn)
      cosnp = cos(ct0)
      sinnp = sin(ct0)
      epw(jn) = (pk(ic)+cosnp*gamc(ic)+sinnp*gams(ic))*fzi(jn)
      if(idecomp == 1 .or. idecomp == 2) then
         ii = (ip(jn+ist) + 3)/3
         jj = (jp(jn+ist) + 3)/3
         kk = (iabs(kp(jn+ist)) + 3)/3
         ll = (iabs(lp(jn+ist)) + 3)/3
         call decphi(ii,jj,kk,ll,epw(jn))
      end if
      df0 = pn(ic)*(gamc(ic)*sinnp-gams(ic)*cosnp)
      dums = sphi(jn)+sign(tm24,sphi(jn))
      dflim = gamc(ic)*(pn(ic)-gmul(inc)+gmul(inc)*cphi(jn))
#ifdef CRAY_PVP
      df1 = cvmgt(dflim,df0/dums,tm06 > abs(dums))
#else
      df1 = df0/dums
      if(tm06 > abs(dums)) df1 = dflim
#endif
      df(jn) = df1*fzi(jn)
   end do

   !         ----- NOW DO TORSIONAL FIRST DERIVATIVES -----

   do jn = 1,maxlen

      !           ----- NOW, SET UP ARRAY DC = FIRST DER. OF COSPHI W/RESPECT
      !                 TO THE CARTESIAN DIFFERENCES T -----

      z11 = z1(jn)*z1(jn)
      z12 = z1(jn)*z2(jn)
      z22 = z2(jn)*z2(jn)
      dc1 = -gx(jn)*z12-cphi(jn)*dx(jn)*z11
      dc2 = -gy(jn)*z12-cphi(jn)*dy(jn)*z11
      dc3 = -gz(jn)*z12-cphi(jn)*dz(jn)*z11
      dc4 =  dx(jn)*z12+cphi(jn)*gx(jn)*z22
      dc5 =  dy(jn)*z12+cphi(jn)*gy(jn)*z22
      dc6 =  dz(jn)*z12+cphi(jn)*gz(jn)*z22

      !           ----- UPDATE THE FIRST DERIVATIVE ARRAY -----

      dr1 = df(jn)*(dc3*ykj(jn) - dc2*zkj(jn))
      dr2 = df(jn)*(dc1*zkj(jn) - dc3*xkj(jn))
      dr3 = df(jn)*(dc2*xkj(jn) - dc1*ykj(jn))
      dr4 = df(jn)*(dc6*ykj(jn) - dc5*zkj(jn))
      dr5 = df(jn)*(dc4*zkj(jn) - dc6*xkj(jn))
      dr6 = df(jn)*(dc5*xkj(jn) - dc4*ykj(jn))
      drx = df(jn)*(-dc2*zij(jn) + dc3*yij(jn) + &
            dc5*zkl(jn) - dc6*ykl(jn))
      dry = df(jn)*( dc1*zij(jn) - dc3*xij(jn) - &
            dc4*zkl(jn) + dc6*xkl(jn))
      drz = df(jn)*(-dc1*yij(jn) + dc2*xij(jn) + &
            dc4*ykl(jn) - dc5*xkl(jn))
      fxi(jn) = - dr1
      fyi(jn) = - dr2
      fzi(jn) = - dr3
      fxj(jn) = - drx + dr1
      fyj(jn) = - dry + dr2
      fzj(jn) = - drz + dr3
      fxk(jn) = + drx + dr4
      fyk(jn) = + dry + dr5
      fzk(jn) = + drz + dr6
      fxl(jn) = - dr4
      fyl(jn) = - dr5
      fzl(jn) = - dr6
   end do  !  jn = 1,maxlen

   !         ----- END OF A STRIP OF DIHEDRALS AND START OF 1-4 NB -----

   !         --- for ewald calcs, 1-4 interactions are done in get_14_cg()
   !             or get_14_dipole().

   if( igb == 0 ) goto 701

   do jn = 1,maxlen
      i3 = ip(jn+ist)
      l3t = lp(jn+ist)
      l3 = iabs(l3t)
      xij(jn) = x(i3+1)-x(l3+1)
      yij(jn) = x(i3+2)-x(l3+2)
      zij(jn) = x(i3+3)-x(l3+3)
   end do

   !             ----- NOW LOOP OVER ALL THE DIHEDRALS (CONSTANT DIEL) -----

   do jn = 1,maxlen
      ct(jn) = xij(jn)*xij(jn)+yij(jn)*yij(jn)+zij(jn)*zij(jn)
   end do

   do jn = 1,maxlen
      i3 = ip(jn+ist)
      k3t = kp(jn+ist)
      l3t = lp(jn+ist)
      ic0 = icp(jn+ist)
      idumi = isign(1,k3t)
      iduml = isign(1,l3t)
      kdiv = (2+idumi+iduml)/4
      l3 = iabs(l3t)
      fmuln = kdiv*fmn(ic0)

      ii = (i3+3)/3
      jj = (l3+3)/3
      ia1 = iac(ii)
      ia2 = iac(jj)
      ibig = max0(ia1,ia2)
      isml = min0(ia1,ia2)
      ic = ibig*(ibig-1)/2+isml

      !             ----- CALCULATE THE 14-EEL ENERGY -----

      r2 = fmuln/ct(jn)
      r1 = sqrt(r2)
      lfac = 1.d0
      if( eedmeth == 5 ) then
         g = cg(ii)*cg(jj)*r2*lfac
      else
         g = cg(ii)*cg(jj)*r1*lfac*intdieli
      end if
      if (icnstph /= 0 .and. mod(irespa,ntcnstph) == 0) then
         dvdl = dvdl + scee0*(dcharge(ii)* &
               dcharge(jj)*r1*lfac*intdieli-g)
      end if
      sphi(jn) = g
      if(idecomp > 0) then
         ii = (i3+3)/3
         jj = (l3+3)/3
         if(idecomp == 1) then
            call decpair(4,ii,jj,sphi(jn)*scee0)
         else if(idecomp == 2) then
            call decpair(2,ii,jj,sphi(jn)*scee0)
            !             else if(idecomp.eq.3) then
            !               --- not considered since
            !                     no pairwise decomp for internal energies
         else if(idecomp == 4) then
            call decpair(-2,ii,jj,sphi(jn)*scee0)
         end if
      end if
      r6 = r2*r2*r2
      r12 = r6*r6
#ifdef CHARMM
      f1 = cn114(ic)*r12*lfac
      f2 = cn214(ic)*r6*lfac
#else
      f1 = cn1(ic)*r12*lfac
      f2 = cn2(ic)*r6*lfac
#endif
      cphi(jn) = f1-f2
      if(idecomp > 0) then
         ii = (i3+3)/3
         jj = (l3+3)/3
         if(idecomp == 1) then
            call decpair(4,ii,jj,cphi(jn)*scnb0)
         else if(idecomp == 2) then
            call decpair(3,ii,jj,cphi(jn)*scnb0)
            !             else if(idecomp.eq.3) then
            !               --- not considered since
            !                     no pairwise decomp for internal energies
         else if(idecomp == 4) then
            call decpair(-3,ii,jj,cphi(jn)*scnb0)
         end if
      end if

      if( eedmeth == 5 ) then
         dfn =((-twelve*f1+six*f2)*scnb0-two*g*scee0)*r2
      else
         dfn =((-twelve*f1+six*f2)*scnb0-g*scee0)*r2
      end if

      xa = xij(jn)*dfn
      ya = yij(jn)*dfn
      za = zij(jn)*dfn
      fxi(jn) = fxi(jn)-xa
      fyi(jn) = fyi(jn)-ya
      fzi(jn) = fzi(jn)-za
      fxl(jn) = fxl(jn)+xa
      fyl(jn) = fyl(jn)+ya
      fzl(jn) = fzl(jn)+za
   end do  !  jn = 1,maxlen

   !         ----- SUMUP ALL THE GRADIENTS -----

   701 do jn = 1,maxlen
      i3 = ip(jn+ist)
      j3 = jp(jn+ist)
      k3 = iabs(kp(jn+ist))
      l3 = iabs(lp(jn+ist))

      f(i3+1) = f(i3+1) + fxi(jn)
      f(i3+2) = f(i3+2) + fyi(jn)
      f(i3+3) = f(i3+3) + fzi(jn)
      f(j3+1) = f(j3+1) + fxj(jn)
      f(j3+2) = f(j3+2) + fyj(jn)
      f(j3+3) = f(j3+3) + fzj(jn)
      f(k3+1) = f(k3+1) + fxk(jn)
      f(k3+2) = f(k3+2) + fyk(jn)
      f(k3+3) = f(k3+3) + fzk(jn)
      f(l3+1) = f(l3+1) + fxl(jn)
      f(l3+2) = f(l3+2) + fyl(jn)
      f(l3+3) = f(l3+3) + fzl(jn)
   end do
   do ksum = 1,maxlen
      enbpl = enbpl+cphi(ksum)
      eelpl = eelpl+sphi(ksum)
      epl   = epl  +epw(ksum)
   end do

   ist = ist+maxlen
   820 if(.not.skip) goto 4200

   !     ---- ALL DONE -----

   enbp = enbpl*scnb0
   eelp = eelpl*scee0
   ep   = epl
   ecn = ecnl
#ifdef MPI
   nphiin = nbtmp
#endif
   return
end subroutine ephi
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine capwat here]
subroutine capwat(nat,x,f)

   implicit _REAL_ (a-h,o-z)
#ifdef MPI
#  include "parallel.h"
#endif

   !     ----- ROUTINE TO CALCULATE THE CAP FORCE -----

#  include "box.h"
#  include "flocntrl.h"
   dimension x(3,*),f(3,*)
   data tm34,zero/1.0d-34,0.0d0/

   ! FLOW CONTROL FLAG (debug)
   if ( do_cap == 0 )return
#ifdef MPI
   do i = natcap+1+mytaskid,nat,numtasks
#else
   do i = natcap+1,nat
#endif

      xa = xcap-x(1,i)
      ya = ycap-x(2,i)
      za = zcap-x(3,i)
      da = sqrt(xa*xa+ya*ya+za*za+tm34)
      df = fcap*max(zero,da-cutcap)/da
      f(1,i) = f(1,i)+df*xa
      f(2,i) = f(2,i)+df*ya
      f(3,i) = f(3,i)+df*za
   end do
   return
end subroutine capwat
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ 
subroutine xconst(natc,econ,igroup,x,f,xc,weit &
#ifdef PSANDER
      ,spmylist,nmyatm,natom &
#endif
      )

   implicit _REAL_ (a-h,o-z)
#ifdef MPI
#  include "parallel.h"
#endif
#  include "flocntrl.h"
#ifdef PSANDER
   integer natom,nmyatm,spmylist(nmyatm)
   logical myatm(natom)
#endif

   !     ----- ROUTINE TO PUT HARMONIC CONSTRAINTS FOR POSITION -----

   dimension igroup(*),x(*),f(*),xc(*),weit(*)

   ! FLOW CONTROL FLAG (debug)
   if ( doxconst == 0 )return
   econ = 0.0d+00

#ifdef PSANDER
   myatm(1:natom)=.false.
   do i=1,nmyatm
      myatm(spmylist(i))=.true.
   enddo
   do ii = 1,natc
      i = igroup(ii)
      if(myatm(i))then
#else
# ifdef MPI
   do ii = 1+mytaskid,natc,numtasks
# else
   do ii = 1,natc
# endif
      i = igroup(ii)
#endif
      wt = weit(ii)
      i3 = 3*i-3
      ax = x(i3+1)-xc(i3+1)
      ay = x(i3+2)-xc(i3+2)
      az = x(i3+3)-xc(i3+3)
      wx = wt*ax
      wy = wt*ay
      wz = wt*az
      eadd = wx*ax+wy*ay+wz*az
      econ = econ+eadd
      f(i3+1) = f(i3+1)-(wx+wx)
      f(i3+2) = f(i3+2)-(wy+wy)
      f(i3+3) = f(i3+3)-(wz+wz)
#ifdef PSANDER
      endif
#endif
   end do
   return
end subroutine xconst
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine bellyf here]
subroutine bellyf(nat,igrp,f)

   implicit _REAL_ (a-h,o-z)
   dimension igrp(*),f(*)
   data zero/0.0d0/
   do i = 1,nat
      if(igrp(i) == 0) then
         i3 = 3*i-3
         f(i3+1) = zero
         f(i3+2) = zero
         f(i3+3) = zero
      end if
   end do
   return
end subroutine bellyf
