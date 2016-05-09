! <compile=optimized>
#include "copyright.h"
#include "is_copyright.h"
#include "dprec.h"
#include "is_def.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ dielectric assignment for fdpb
module pb_dielectrics

   implicit none
    
contains

subroutine pb_epsmap( pbverbose,dielopt,smoothopt,epsin,epsout,sprob,&
                      level,nfocus,xm,ym,zm,xmymzm,h,gox,goy,goz,&
                      gcrd,arccrd,radi,nbnd,iepsav,atmsas,insas,epsx,epsy,epsz )

   ! Passed variables

   logical pbverbose
   integer dielopt,smoothopt
   _REAL_ epsin,epsout,sprob
   integer level,nfocus,xm,ym,zm,xmymzm
   _REAL_ h,gox,goy,goz
   _REAL_ gcrd(3,*),arccrd(3,*),radi(*)
   integer nbnd,iepsav(4,*)
   integer atmsas(xmymzm),insas(xmymzm)
   _REAL_ epsx(xmymzm),epsy(xmymzm),epsz(xmymzm)
    
   ! Local variables
    
   _REAL_ rh
    
   rh = ONE/h
    
   ! initialize all grids to solvent by default
    
   epsx(1:xmymzm) = epsout; epsy(1:xmymzm) = epsout; epsz(1:xmymzm) = epsout
    
   ! save dielectric boundary edges for db energy and forces
    
   if ( level == nfocus ) then
      select case ( dielopt )
      case ( 2 )
         call epsbnd_spline( xm, ym, zm, atmsas, insas )  
      case ( 1 )
         call epsbnd_mvdw( xm, ym, zm, atmsas, insas )
      case ( 0 )
         call epsbnd_ses( xm, ym, zm, atmsas, insas )
      case default
         ! invalid surfopt
         write(6,"(a)") "PB Bomb in pb_epsmap(): invalid dielopt"
         call mexit (1)
      end select
   endif
    
   ! use the insas grid to setup epsx, epsy and epsz maps
   
   if ( epsout /= epsin ) & 
   call epsmap( xm, ym, zm, atmsas, insas, epsx, epsy, epsz )
    
contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Save dielectric boundary grid points
subroutine epsbnd_spline ( xm, ym, zm, atmsas,insas )
    
   integer xm, ym, zm
   integer atmsas(xm,ym,zm), insas(xm,ym,zm)
    
   ! Local variables
    
   return
end subroutine epsbnd_spline
subroutine epsbnd_mvdw ( xm, ym, zm, atmsas,insas )

   integer xm, ym, zm
   integer atmsas(xm,ym,zm), insas(xm,ym,zm)
    
   ! Local variables
    
   return
end subroutine epsbnd_mvdw
subroutine epsbnd_ses ( xm, ym, zm, atmsas,insas )
   
   integer xm, ym, zm 
   integer atmsas(xm,ym,zm), insas(xm,ym,zm)
    
   ! Local variables
    
   logical boundary
   integer nwarn, buffer, i, j, k, clstmp
    
   nwarn = 0
   nbnd = 0
   buffer = 1
   do k = 1+buffer, zm-buffer; do j = 1+buffer, ym-buffer; do i = 1+buffer, xm-buffer
       
      ! set up condition for a boundary grid point
       
      boundary = .false.
      if ( (insas(i,j,k) > 0) .and.&
           (insas(i-1,j,k) < 0 .or. insas(i+1,j,k) < 0 .or. &
            insas(i,j-1,k) < 0 .or. insas(i,j+1,k) < 0 .or. &
            insas(i,j,k-1) < 0 .or. insas(i,j,k+1) < 0) ) then 
            boundary = .true.
      else if ( (insas(i,j,k) < 0) .and.&
           (insas(i-1,j,k) > 0 .or. insas(i+1,j,k) > 0 .or.&
            insas(i,j-1,k) > 0 .or. insas(i,j+1,k) > 0 .or.&
            insas(i,j,k-1) > 0 .or. insas(i,j,k+1) > 0) ) then
            boundary = .true.
      end if
      if ( .not. boundary ) cycle
 
      nbnd = nbnd + 1; iepsav(1,nbnd) = i; iepsav(2,nbnd) = j; iepsav(3,nbnd) = k
 
      ! for a grid point in contact region +/- 2 or in a solvent probe, simply use
      ! the atom/probe that marks it
 
      if ( abs(insas(i,j,k)) == 2 .or. insas(i,j,k) == -1 .or. insas(i,j,k) == -3 ) then
         clstmp = atmsas(i,j,k)
         if ( clstmp == 0 ) then
            if ( pbverbose ) write(6, '(a,4i4)') &
            'PB Warning in epsbnd_ses(): No neighbor found for exposed contact grid', &
            i, j, k, insas(i,j,k)
            nwarn = nwarn + 1
         end if
 
      ! for a buried reentry grid point, find the atom that marked its neighoring exposed
      ! reentry grid points. Note that this may not be possible when grid spacing is large
 
      else if ( insas(i,j,k) == 1 ) then
         clstmp = fndcls( i, j, k, xm, ym, zm, insas, atmsas )
         if ( clstmp == 0 ) then
            if ( pbverbose ) write(6, '(a,4i4)') &
            'PB Warning in epsbnd_ses(): No neighbor found for exposed reentry grid', &
            i, j, k, insas(i,j,k)
            nwarn = nwarn + 1
         end if
      end if
 
      iepsav(4,nbnd) = clstmp
   end do; end do; end do
   if ( nwarn > 0 ) then
      if ( pbverbose ) write(6, '(a,i4)') &
      'PB Warning in epsbnd_ses(): No neighbor found for boundary grids total:', nwarn
   end if

end subroutine epsbnd_ses
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Find the closest reentry probe for a reentry boundary grid
function fndcls( i,j,k,xm,ym,zm,insas,atmsas )
    
   implicit none
    
   ! Passed variables
    
   integer fndcls, i, j, k
   integer xm, ym, zm
   integer insas(xm,ym,zm), atmsas(xm,ym,zm)
    
   ! Local variables
    
   integer iatm, l, lp, ip, jp, kp, iip(6), jjp(6), kkp(6), clsatm(6)
   _REAL_ xg, yg, zg
   _REAL_ dx, dy, dz, d, clsdst, clscrd(3,6)

   ! first stack these candidates into a 1-d list
    
   iip(1)=i-1; iip(2)=i+1; jjp(1:2)=j; kkp(1:2)=k
   iip(3:4)=i; jjp(3)=j-1; jjp(4)=j+1; kkp(3:4)=k
   iip(5:6)=i; jjp(5:6)=j; kkp(5)=k-1; kkp(6)=k+1
   lp = 0
   do l = 1, 6
      ip = iip(l); jp = jjp(l); kp = kkp(l)
      if ( atmsas(ip,jp,kp) == 0 .or. insas(ip,jp,kp) /= -1 ) cycle
      lp = lp + 1; iatm = atmsas(ip,jp,kp); clsatm(lp) = iatm
      clscrd(1,lp) = arccrd(1,iatm)
      clscrd(2,lp) = arccrd(2,iatm)
      clscrd(3,lp) = arccrd(3,iatm)
   end do
 
   ! now find the closest
 
   xg = gox + i*h; yg = goy + j*h; zg = goz + k*h
   clsdst = 999.d0
   fndcls = 0
   do ip = 1, lp
      dx = clscrd(1,ip) - xg; dy = clscrd(2,ip) - yg; dz = clscrd(3,ip) - zg
      d = abs(sqrt(dx**2 + dy**2 + dz**2) - sprob)
      if ( d >= clsdst ) cycle
      clsdst = d
      fndcls = clsatm(ip)
   end do
 
end function fndcls
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Map insas into epsmap
subroutine epsmap( xm,ym,zm,atmsas,insas,epsx,epsy,epsz )

   integer xm, ym, zm
   integer atmsas(xm,ym,zm), insas(xm,ym,zm)
   _REAL_ epsx(xm,ym,zm), epsy(xm,ym,zm), epsz(xm,ym,zm)

   integer i, j, k, i1, j1, k1, a, b, c, d, a1, b1, c1, d1
   _REAL_ epsint

   !integer total_volume, reentry_volume, radiin_volume
   !total_volume=0
   !reentry_volume=0
   !radiin_volume=0 
   !open(10,file='volume.dat')
   !open(11,file='reentry1.dat')
   !open(12,file='reentry3.dat')

   epsint = TWO*epsin*epsout/(epsin+epsout)
   do k = 1, zm; do j = 1, ym; do i = 1, xm

      ! for PBC
      ! for NPBC, 1 and max are both should be in solvent

      i1 = i+1; if ( i1 > xm ) i1 = 1
      j1 = j+1; if ( j1 > ym ) j1 = 1
      k1 = k+1; if ( k1 > zm ) k1 = 1

      !if(insas(i,j,k).ge.1)total_volume=total_volume+1 
      !if(insas(i,j,k)==1) then
      !reentry_volume=reentry_volume+1
      !write(011,*)i,j,k
      !endif 
      !if(insas(i,j,k)==3) then
      !radiin_volume=radiin_volume+1
      !write(012,*)i,j,k
      !endif

      a = insas(i,j,k)
      b = insas(i1,j,k)
      c = insas(i,j1,k)
      d = insas(i,j,k1)
      a1 = atmsas(i,j,k)
      b1 = atmsas(i1,j,k)
      c1 = atmsas(i,j1,k)
      d1 = atmsas(i,j,k1)
      ! x-edges
      if ( (a > 0 .and. b > 0) ) then
         epsx(i,j,k) = epsin
      end if
      if ( (a > 0 .and. b <= 0) .or. (a <= 0 .and. b > 0) ) then 
         if ( smoothopt == 1 ) call epsfracx(i,j,k,a,b,a1,b1,rh,epsint,epsin,epsout)
         epsx(i,j,k) = epsint
      end if 
      ! y-edges
      if ( (a > 0 .and. c > 0) ) then
         epsy(i,j,k) = epsin
      end if
      if ( (a > 0 .and. c <= 0) .or. (a <= 0 .and. c > 0) ) then
         if ( smoothopt == 1 ) call epsfracy(i,j,k,a,c,a1,c1,rh,epsint,epsin,epsout)
         epsy(i,j,k) = epsint
      end if
      ! z-edges
      if ( (a > 0 .and. d > 0) ) then
         epsz(i,j,k) = epsin
      end if
      if ( (a > 0 .and. d <= 0) .or. (a <= 0 .and. d > 0) ) then
         if ( smoothopt == 1 ) call epsfracz(i,j,k,a,d,a1,d1,rh,epsint,epsin,epsout)
         epsz(i,j,k) = epsint
      end if
   end do; end do; end do

   ! write(10,*)total_volume,reentry_volume,radiin_volume
   ! close(12)
   ! close(11)
   ! close(10)
!  do k = 1, zm
!      write(20, *) 'plane', k
!   do j = 1, ym
!      write(20, '(100f6.1)') epsx(1:xm,j,k)/eps0
!   end do
!   end do
!   do k = 1, zm
!      write(21, *) 'plane', k
!   do i = 1, xm
!      write(21, '(100f6.1)') epsy(i,1:ym,k)/eps0
!   end do
!   end do
!   do j = 1, ym
!      write(22, *) 'plane', j
!   do i = 1, xm
!      write(22, '(100f6.1)') epsz(i,j,1:zm)/eps0
!   end do
!   end do

end subroutine epsmap
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Assign fractional eps values for x-edges
subroutine epsfracx( i, j, k, a, b, a1, b1, rh, epsint, depsin, depsout )

   implicit none
   integer  i, j, k, a, b, a1, b1
   _REAL_ rh, epsint
   _REAL_ depsin, depsout

   integer iatm
   _REAL_ range1, range3, xi, yi, zi, aa

   ! locate the atom that is crossing this x-edge, note the order iatm is assigned

   if ( a == 2 ) then
      iatm = a1
   else if ( b == 2 ) then
      iatm = b1
   else if ( a == -1 ) then
      iatm = a1
   else if ( b == -1 ) then
      iatm = b1
   end if

   ! obtain the position and radius of the atom (probe) in grid unit

   if ( a == 2 .or. b == 2 ) then
      range1 = radi(iatm)*rh
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      range3 = sqrt(range1**2-(zi-k)**2-(yi-j)**2)
      if ( range3 > ZERO ) then
         if ( b == 2 ) then
            aa = range3 - xi + REAL(i+1)
         else
            aa = range3 + xi - REAL(i)
         end if
         epsint = (depsout*depsin)/(depsin*(ONE-aa) + depsout*aa)
      else
         epsint = depsout
      end if
   else if ( a == -1 .or. b == -1 ) then
      range1 = sprob*rh
      xi = (arccrd(1,iatm) - gox)*rh; yi = (arccrd(2,iatm) - goy)*rh; zi = (arccrd(3,iatm) - goz)*rh
      range3 = sqrt(range1**2-(zi-k)**2-(yi-j)**2)
      if ( range3 > ZERO ) then
         if ( b == -1 ) then
            aa = range3 - xi + REAL(i+1)
         else
            aa = range3 + xi - REAL(i)
         end if
         epsint = (depsin*depsout)/(depsout*(ONE-aa) + depsin*aa)
      else
         epsint = depsin
      end if
   end if

   ! other situations will not be considered and use default epsint value

end subroutine epsfracx
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Assign fractional eps values for y-edges
subroutine epsfracy( i, j, k, a, b, a1, b1, rh, epsint, depsin, depsout )

   implicit none
   integer  i, j, k, a, b, a1, b1
   _REAL_ rh, epsint
   _REAL_ depsin, depsout

   integer iatm
   _REAL_ range1, range3, xi, yi, zi, aa

   ! locate the atom that is crossing this y-edge, note the order iatm is assigned

   if ( a == 2 ) then
      iatm = a1
   else if ( b == 2 ) then
      iatm = b1
   else if ( a == -1 ) then
      iatm = a1
   else if ( b == -1 ) then
      iatm = b1
   end if

   ! obtain the position and radius of the atom (probe) in grid unit

   if ( a == 2 .or. b == 2 ) then
      range1 = radi(iatm)*rh
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      range3 = sqrt(range1**2-(zi-k)**2-(xi-i)**2)
      if ( range3 > ZERO ) then
         if ( b == 2 ) then
            aa = range3 - yi + REAL(j+1)
         else
            aa = range3 + yi - REAL(j)
         end if
         epsint = (depsout*depsin)/(depsin*(ONE-aa) + depsout*aa)
      else
         epsint = depsout
      end if
   else if ( a == -1 .or. b == -1 ) then
      range1 = sprob*rh
      xi = (arccrd(1,iatm) - gox)*rh; yi = (arccrd(2,iatm) - goy)*rh; zi = (arccrd(3,iatm) - goz)*rh
      range3 = sqrt(range1**2-(zi-k)**2-(xi-i)**2)
      if ( range3 > ZERO ) then
         if ( b == -1 ) then
            aa = range3 - yi + REAL(j+1)
         else
            aa = range3 + yi - REAL(j)
         end if
         epsint = (depsin*depsout)/(depsout*(ONE-aa) + depsin*aa)
      else
         epsint = depsin
      end if
   end if

   ! other situations will not be considered and use default epsint value

end subroutine epsfracy
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Assign fractional eps values for z-edges
subroutine epsfracz( i, j, k, a, b, a1, b1, rh, epsint, depsin, depsout )

   implicit none
   integer  i, j, k, a, b, a1, b1
   _REAL_ rh, epsint
   _REAL_ depsin, depsout

   integer iatm
   _REAL_ range1, range3, xi, yi, zi, aa

   ! locate the atom that is crossing this z-edge, note the order iatm is assigned

   if ( a == 2 ) then
      iatm = a1
   else if ( b == 2 ) then
      iatm = b1
   else if ( a == -1 ) then
      iatm = a1
   else if ( b == -1 ) then
      iatm = b1
   end if

   if ( a == 2 .or. b == 2 ) then
      range1 = radi(iatm)*rh
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      range3 = sqrt(range1**2-(yi-j)**2-(xi-i)**2)
      if ( range3 > ZERO ) then
         if ( b == 2 ) then
            aa = range3 - zi + REAL(k+1)
         else
            aa = range3 + zi - REAL(k)
         end if
         epsint = (depsout*depsin)/(depsin*(ONE-aa) + depsout*aa)
      else
         epsint = depsout
      end if
   else if ( a == -1 .or. b == -1 ) then
      range1 = sprob*rh
      xi = (arccrd(1,iatm) - gox)*rh; yi = (arccrd(2,iatm) - goy)*rh; zi = (arccrd(3,iatm) - goz)*rh
      range3 = sqrt(range1**2-(yi-j)**2-(xi-i)**2)
      if ( range3 > ZERO ) then
         if ( b == -1 ) then
            aa = range3 - zi + REAL(k+1)
         else
            aa = range3 + zi - REAL(k)
         end if
         epsint = (depsin*depsout)/(depsout*(ONE-aa) + depsin*aa)
      else
         epsint = depsin
      end if
   end if

   ! other situations will not be considered and use default epsint value

end subroutine epsfracz

end subroutine pb_epsmap
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Exlusion of ions from protein interior, insas holds flags of sas surf.
subroutine pb_ionmap ( xm,ym,zm,insas,pv )
    
   ! Passed variables
    
   integer xm, ym, zm
   integer insas(xm,ym,zm)
   _REAL_ pv(xm,ym,zm)
    
   ! Local variables
    
   integer i, j, k, buffer
   _REAL_ exclusion
    
   ! for InsightII display
   !_REAL_ g(3)
   ! 
   !open (unit=55, file='ions.dot')
   !write (55, '("DOTS")')

   buffer = 1
   do k = 1+buffer, zm-buffer; do j = 1+buffer, ym-buffer; do i = 1+buffer, xm-buffer
      exclusion = ZERO
      if ( insas(i-1,j,k) == 1 .or. insas(i  ,j,k) == 1 .or. insas(i+1,j,k) == 1 .or. &
           insas(i,j-1,k) == 1 .or. insas(i,j+1,k) == 1 .or.&
           insas(i,j,k-1) == 1 .or. insas(i,j,k+1) == 1 ) then
         exclusion = SIX
         ! for InsightII display
         !g(1) = gox + h*i; g(2) = goy + h*j; g(3) = goz + h*k
         !write (55,'(4(f8.3,2x))') g(1:3), 300.
      end if
      pv(i,j,k) = exclusion
   end do; end do; end do

   ! for InsightII display
   !close(55)
   !stop

end subroutine pb_ionmap

end module pb_dielectrics
