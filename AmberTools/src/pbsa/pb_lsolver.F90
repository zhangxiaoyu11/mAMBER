#include "copyright.h"
#include "../include/dprec.fh"

module pb_lsolver

   implicit none

#  include "pb_constants.h"

   integer l_xm, l_ym, l_zm, l_xmym, l_xmymzm
   _REAL_ l_fmiccg,l_wsor,l_fmiccg2
   _REAL_ l_norm, l_inorm, l_accept, l_epsout, l_pbkappa, l_h
   integer l_itn,l_maxitn,l_bcopt
   integer mg_nlevel,ncyc_before,ncyc_after

!  All
   _REAL_,allocatable :: l_am1(:), l_am2(:), l_am3(:), l_ad(:)
   _REAL_,allocatable :: l_am4(:), l_am5(:), l_am6(:)
   _REAL_,allocatable :: l_bv(:)
   _REAL_,allocatable :: l_pv(:), l_tv(:), l_zv(:), l_rd(:)
!PICCG
   _REAL_,allocatable :: pl_am1(:,:,:), pl_am2(:,:,:), pl_am3(:,:,:)
   _REAL_,allocatable :: pl_am4(:,:,:), pl_am5(:,:,:), pl_am6(:,:,:)
   _REAL_,allocatable :: pl_bv(:,:,:), pl_ad(:,:,:), pl_rd(:,:,:)
   _REAL_,allocatable :: pl_pv(:,:,:), pl_tv(:,:,:), pl_zv(:,:,:)
   _REAL_,allocatable :: pl_xs(:,:,:), pl_phi(:,:,:)
   integer*4,allocatable :: pl_ind(:,:)
!  MG
   integer,allocatable ::  mg_index(:), mg_index_ext(:),mg_x_idx(:),mg_size(:,:)
   _REAL_,allocatable ::  mg_onorm(:)
   _REAL_,allocatable ::  l_rv(:), l_iv(:), l_bz(:), l_xv(:)

contains

!===========================================================================

subroutine init_param(nx,ny,nz,nxny,nxnynz,p_maxitn,p_bcopt,p_fmiccg,p_fmiccg2,p_accept,p_pbkappa,p_epsout,p_h,p_wsor)

   implicit none

   integer nx,ny,nz,nxny,nxnynz,p_maxitn,p_bcopt
   _REAL_ p_fmiccg,p_accept,p_epsout,p_pbkappa,p_wsor,p_h,p_fmiccg2

   l_xm = nx
   l_ym = ny
   l_zm = nz
   l_xmym = nxny
   l_xmymzm = nxnynz
   l_maxitn = p_maxitn
   l_bcopt = p_bcopt
   l_accept = p_accept

!  ICCG
   l_fmiccg = p_fmiccg
   l_fmiccg2= p_fmiccg2
!  MG
   mg_nlevel = 4
   ncyc_before = 10
   ncyc_after = 10
   l_pbkappa = p_pbkappa
   l_epsout = p_epsout
   l_h       = p_h
!  SOR
   l_wsor = p_wsor

end subroutine

!===========================================================================

subroutine allocate_array(solvopt)

   implicit none
   integer solvopt

   integer l,m,n,i,lb,ub

   select case (solvopt)
   case (1)
      if (l_bcopt /= 10 ) then
         allocate(l_ad(1:l_xmymzm+l_xmym),l_am1(1-l_xmym:l_xmymzm+l_xmym))
         allocate(l_am2(1-l_xmym:l_xmymzm+l_xmym), l_am3(1-l_xmym:l_xmymzm+l_xmym))
         allocate(l_rd(1-l_xmym:l_xmymzm),l_bv(1-l_xmym:l_xmymzm))
         allocate(l_tv(1-l_xmym:l_xmymzm+l_xmym),l_zv(1-l_xmym:l_xmymzm+l_xmym))
         allocate(l_pv(1-l_xmym:l_xmymzm+l_xmym))
      else
         lb=1
         ub=(1+l_xm)*(1+l_ym)*(1+l_zm)
         allocate( l_ad(lb:ub),l_am1(lb:ub))
         allocate(l_am2(lb:ub),l_am3(lb:ub))
         allocate(l_am4(lb:ub),l_am5(lb:ub))
         allocate(l_am6(lb:ub))
         allocate( l_rd(lb:ub), l_bv(lb:ub))
         allocate( l_tv(lb:ub), l_zv(lb:ub))
         allocate( l_pv(lb:ub))

         allocate( pl_ad(0:l_xm,0:l_ym,0:l_zm),pl_rd(0:l_xm,0:l_ym,0:l_zm) )
         allocate( pl_bv(0:l_xm,0:l_ym,0:l_zm),pl_tv(1:l_xm+1,1:l_ym+1,1:l_zm+1) )
         allocate( pl_zv(0:l_xm,0:l_ym,0:l_zm),pl_pv(0:l_xm,0:l_ym,0:l_zm) )
         allocate( pl_am1(0:l_xm,0:l_ym,0:l_zm),pl_am2(0:l_xm,0:l_ym,0:l_zm) )
         allocate( pl_am3(0:l_xm,0:l_ym,0:l_zm),pl_am4(0:l_xm,0:l_ym,0:l_zm) )
         allocate( pl_am5(0:l_xm,0:l_ym,0:l_zm),pl_am6(0:l_xm,0:l_ym,0:l_zm) )
         allocate( pl_ind(1:l_xmymzm,3), pl_xs(0:l_xm,0:l_ym,0:l_zm))
         allocate( pl_phi(0:l_xm,0:l_ym,0:l_zm) )
      end if
   case (3)
      allocate(l_ad(1:l_xmymzm),l_am1(1-l_xmym:l_xmymzm),l_am2(1-l_xmym:l_xmymzm),l_am3(1-l_xmym:l_xmymzm))
      allocate(l_am4(1-l_xmym:l_xmymzm),l_am5(1-l_xmym:l_xmymzm),l_am6(1-l_xmym:l_xmymzm))
      allocate(l_bv(1:l_xmymzm),l_pv(1-l_xmym:l_xmymzm+l_xmym),l_zv(1:l_xmymzm))
   case (4)
      allocate(l_ad(1:l_xmymzm),l_am1(1-l_xmym:l_xmymzm),l_am2(1-l_xmym:l_xmymzm),l_am3(1-l_xmym:l_xmymzm))
      allocate(l_bv(1:l_xmymzm),l_zv(1:l_xmymzm))
   case (5)
      allocate(l_ad(1:l_xmymzm),l_am1(1-l_xmym:l_xmymzm),l_am2(1-l_xmym:l_xmymzm),l_am3(1-l_xmym:l_xmymzm))
      allocate(l_am4(1-l_xmym:l_xmymzm),l_am5(1-l_xmym:l_xmymzm),l_am6(1-l_xmym:l_xmymzm))
      allocate(l_bv(1:l_xmymzm),l_pv(1-l_xmym:l_xmymzm+l_xmym),l_zv(1:l_xmymzm))
   case (8)
         lb=1
         ub=(1+l_xm)*(1+l_ym)*(1+l_zm)
         allocate( l_ad(lb:ub),l_am1(lb:ub))
         allocate(l_am2(lb:ub),l_am3(lb:ub))
         allocate(l_am4(lb:ub),l_am5(lb:ub))
         allocate(l_am6(lb:ub))
         allocate( l_rd(lb:ub), l_bv(lb:ub))
         allocate( l_tv(lb:ub), l_zv(lb:ub))
         allocate( l_pv(lb:ub))

         allocate( pl_ad(0:l_xm,0:l_ym,0:l_zm),pl_rd(0:l_xm,0:l_ym,0:l_zm) )
         allocate( pl_bv(0:l_xm,0:l_ym,0:l_zm),pl_tv(1:l_xm+1,1:l_ym+1,1:l_zm+1) )
         allocate( pl_zv(0:l_xm,0:l_ym,0:l_zm),pl_pv(0:l_xm,0:l_ym,0:l_zm) )
         allocate( pl_am1(0:l_xm,0:l_ym,0:l_zm),pl_am2(0:l_xm,0:l_ym,0:l_zm) )
         allocate( pl_am3(0:l_xm,0:l_ym,0:l_zm),pl_am4(0:l_xm,0:l_ym,0:l_zm) )
         allocate( pl_am5(0:l_xm,0:l_ym,0:l_zm),pl_am6(0:l_xm,0:l_ym,0:l_zm) )
         allocate( pl_ind(1:l_xmymzm,3), pl_xs(0:l_xm,0:l_ym,0:l_zm))
         allocate( pl_phi(0:l_xm,0:l_ym,0:l_zm) )
   case (2)
      allocate ( mg_index(1:mg_nlevel+1), mg_index_ext(1:mg_nlevel+1))
      allocate ( mg_x_idx(1:mg_nlevel+1),mg_size(1:3,1:mg_nlevel) )
      allocate ( mg_onorm(1:mg_nlevel) )

      mg_index_ext(1) = 1
      mg_index(1) = 1
      mg_x_idx(1) = 1
      mg_size(1,1) = l_xm
      mg_size(2,1) = l_ym
      mg_size(3,1) = l_zm
      m = l_xmymzm
      l = m + l_xmym
      n = l + l_xmym
      allocate( l_zv(1:m) )
      do i = 2, mg_nlevel
         mg_index_ext(i) = 1 + l
         mg_index(i) = 1 + m
         mg_x_idx(i) = 1 + n
         mg_size(1:3,i) = ( mg_size(1:3,i-1) - 1 ) / 2
         m = m + mg_size(1,i) * mg_size(2,i) * mg_size(3,i)
         l = l + mg_size(1,i) * mg_size(2,i) * mg_size(3,i) + mg_size(1,i) * mg_size(2,i)
         n = n + mg_size(1,i) * mg_size(2,i) * mg_size(3,i) + 2 * mg_size(1,i) * mg_size(2,i)
      end do
      mg_index_ext(i) = 1 + l
      mg_index(i) = 1 + m
      mg_x_idx(i) = 1 + n
      allocate ( l_ad(1:m), l_bv(1:m), l_rv(1:m), l_iv(1:m), l_bz(1:m) )
      allocate ( l_am1(1:l), l_am2(1:l), l_am3(1:l) )
      allocate ( l_xv(1:n) )
   end select

end subroutine allocate_array

!===========================================================================

subroutine deallocate_array(solvopt)

   implicit none
   integer solvopt

   select case (solvopt)
   case (1)
      if ( l_bcopt /= 10 ) then
         deallocate(l_ad,l_am1)
         deallocate(l_am2,l_am3)
         deallocate(l_rd,l_bv)
         deallocate(l_tv,l_zv)
         deallocate(l_pv)
      else
         deallocate( l_ad,l_am1)
         deallocate(l_am2,l_am3)
         deallocate(l_am4,l_am5)
         deallocate(l_am6)
         deallocate( l_rd, l_bv)
         deallocate( l_tv, l_zv)
         deallocate( l_pv)

         deallocate( pl_ad,pl_rd )
         deallocate( pl_bv,pl_tv )
         deallocate( pl_zv,pl_pv )
         deallocate( pl_am1,pl_am2 )
         deallocate( pl_am3,pl_am4 )
         deallocate( pl_am5,pl_am6 )
         deallocate( pl_ind,pl_xs, pl_phi )
      end if
   case (3)
      deallocate(l_ad,l_am1,l_am2,l_am3)
      deallocate(l_am4,l_am5,l_am6)
      deallocate(l_bv,l_pv,l_zv)
   case (4)
      deallocate(l_ad,l_am1,l_am2,l_am3)
      deallocate(l_bv,l_zv)
   case (5)
      deallocate(l_ad,l_am1,l_am2,l_am3)
      deallocate(l_am4,l_am5,l_am6)
      deallocate(l_bv,l_pv,l_zv)
   case (8)
         deallocate( l_ad,l_am1)
         deallocate(l_am2,l_am3)
         deallocate(l_am4,l_am5)
         deallocate(l_am6)
         deallocate( l_rd, l_bv)
         deallocate( l_tv, l_zv)
         deallocate( l_pv)

         deallocate( pl_ad,pl_rd )
         deallocate( pl_bv,pl_tv )
         deallocate( pl_zv,pl_pv )
         deallocate( pl_am1,pl_am2 )
         deallocate( pl_am3,pl_am4 )
         deallocate( pl_am5,pl_am6 )
         deallocate( pl_ind,pl_xs, pl_phi )
   case (2)
      deallocate( mg_index, mg_index_ext )
      deallocate( mg_x_idx,mg_size )
      deallocate( mg_onorm )
      deallocate( l_zv )
      deallocate( l_ad, l_bv, l_rv, l_iv, l_bz )
      deallocate( l_am1, l_am2, l_am3 )
      deallocate( l_xv )
   end select

end subroutine deallocate_array

!==============================================================================

!===========================================================================

subroutine init_array( solvopt,epsx,epsy,epsz,p_bv,p_iv,p_xs )
   
   implicit none

   integer solvopt
   _REAL_ epsx(*),epsy(*),epsz(*)
   _REAL_ p_bv(1:l_xmymzm),p_iv(1:l_xmymzm)
   _REAL_ p_xs(1-l_xmym:l_xmymzm+l_xmym)

   integer lxmym,l,m,n,i,j,k,ii
   _REAL_,allocatable :: lepsx(:), lepsy(:), lepsz(:) 
   _REAL_ lfactor

   integer lb,ub

   select case (solvopt)
   case (1)
      if ( l_bcopt /= 10 ) then
         call feedepsintoam(l_xm, l_ym, l_zm, l_am1(1:l_xmymzm),  &
                         l_am2(1:l_xmymzm), l_am3(1:l_xmymzm), &
                                              epsx, epsy, epsz )
         l_am1(l_xmymzm+1:l_xmymzm+l_xmym) = ZERO
         l_am2(l_xmymzm+1:l_xmymzm+l_xmym) = ZERO
         l_am3(l_xmymzm+1:l_xmymzm+l_xmym) = ZERO
         if (l_pbkappa == ZERO) then
            call feedepsintoad(l_xm, l_ym, l_zm, l_ad, epsx, epsy, epsz)
         else
            lfactor = l_epsout*(l_h*l_pbkappa)**2
            call feedepsintoad(l_xm, l_ym, l_zm, l_ad, epsx, epsy, epsz)
            l_ad(1:l_xmymzm) = lfactor*p_iv(1:l_xmymzm) + l_ad(1:l_xmymzm)
         end if
         l_am1(1-l_xmym:0) = ZERO
         l_am2(1-l_xmym:0) = ZERO
         l_am3(1-l_xmym:0) = ZERO
         call pb_setupper(l_am1(1),l_am2(1),l_am3(1))
         l_ad (l_xmymzm+1:l_xmymzm+l_xmym) = ZERO
         l_bv (1-l_xmym:0) = ZERO
         l_bv (1:l_xmymzm) = p_bv(1:l_xmymzm)
         l_pv (1-l_xmym:0) = ZERO
         l_pv (l_xmymzm+1:l_xmymzm+l_xmym) = ZERO
         l_tv = ZERO
         l_zv = ZERO
         l_rd = ONE
      else
         lb=1

         ub=(1+l_xm)*(1+l_ym)*(1+l_zm)
         !initialize all to zero first
         pl_ad=ZERO
         pl_am1=ZERO; pl_am2=ZERO; pl_am3=ZERO 
         pl_am4=ZERO; pl_am5=ZERO; pl_am6=ZERO
         pl_bv=ZERO
         pl_tv=ZERO
         pl_zv=ZERO
         pl_pv=ZERO
         pl_rd=ZERO

         !set up the 1D to 3D indexing array, pl_ind
         do k=1,l_zm;do j=1,l_ym;do i=1,l_xm
            ii=i+l_xm*(j-1+l_ym*(k-1))
            pl_ind(ii,1) = i
            pl_ind(ii,2) = j
            pl_ind(ii,3) = k
         end do;end do; end do

         !set l_bv
         call pb_piccg_loadbv(pl_bv,p_bv(1:l_xmymzm))

         call feedepsintoam_piccg(  l_xm,l_ym,l_zm,pl_am1,  &
                                 pl_am2,pl_am3,     &
                                 epsx,epsy,epsz )
         !set central diagonal array
         call feedepsintoad_piccg(l_xm, l_ym, l_zm, pl_ad, epsx, epsy, epsz,&
                               pl_pv)
         call pb_setoutter_piccg2(pl_am1,pl_am2,pl_am3, &
                                pl_am4,pl_am5,pl_am6)
      end if
   case (3)
      call feedepsintoam(l_xm, l_ym, l_zm, l_am1(1:l_xmymzm),  &
                         l_am2(1:l_xmymzm), l_am3(1:l_xmymzm), &
                                              epsx, epsy, epsz )
      if (l_pbkappa == ZERO) then
         call feedepsintoad(l_xm, l_ym, l_zm, l_ad, epsx, epsy, epsz)
      else
         lfactor = l_epsout*(l_h*l_pbkappa)**2
         call feedepsintoad(l_xm, l_ym, l_zm, l_ad, epsx, epsy, epsz)
         l_ad(1:l_xmymzm) = lfactor*p_iv(1:l_xmymzm) + l_ad(1:l_xmymzm)
      end if
      l_am1(1-l_xmym:0) = ZERO
      l_am2(1-l_xmym:0) = ZERO
      l_am3(1-l_xmym:0) = ZERO
      if (l_bcopt == 10 ) then
         call pb_setupper2(l_am1(1),l_am2(1),l_am3(1), &
                           l_am4(1),l_am5(1),l_am6(1))
      else
         call pb_setupper(l_am1(1),l_am2(1),l_am3(1))
      end if
      l_bv (1:l_xmymzm) = p_bv(1:l_xmymzm)
      l_pv (1-l_xmym:0) = ZERO
      l_pv (l_xmymzm+1:l_xmymzm+l_xmym) = ZERO
   case (4)
      call feedepsintoam(l_xm, l_ym, l_zm, l_am1(1:l_xmymzm),  &
                         l_am2(1:l_xmymzm), l_am3(1:l_xmymzm), &
                                              epsx, epsy, epsz )
      if (l_pbkappa == ZERO) then
         call feedepsintoad(l_xm, l_ym, l_zm, l_ad, epsx, epsy, epsz)
      else
         lfactor = l_epsout*(l_h*l_pbkappa)**2
         call feedepsintoad(l_xm, l_ym, l_zm, l_ad, epsx, epsy, epsz)
         l_ad(1:l_xmymzm) = lfactor*p_iv(1:l_xmymzm) + l_ad(1:l_xmymzm)
      end if
      l_am1(1-l_xmym:0) = ZERO
      l_am2(1-l_xmym:0) = ZERO
      l_am3(1-l_xmym:0) = ZERO
      call pb_setupper(l_am1(1),l_am2(1),l_am3(1))
      l_bv (1:l_xmymzm) = p_bv(1:l_xmymzm)
   case (5)
      call feedepsintoam(l_xm, l_ym, l_zm, l_am1(1:l_xmymzm),  &
                         l_am2(1:l_xmymzm), l_am3(1:l_xmymzm), &
                                              epsx, epsy, epsz )
      if (l_pbkappa == ZERO) then
         call feedepsintoad(l_xm, l_ym, l_zm, l_ad(1:l_xmymzm), epsx, epsy, epsz)
      else
         lfactor = l_epsout*(l_h*l_pbkappa)**2
         call feedepsintoad(l_xm, l_ym, l_zm, l_ad(1:l_xmymzm), epsx, epsy, epsz)
         l_ad(1:l_xmymzm) = lfactor*p_iv(1:l_xmymzm) + l_ad(1:l_xmymzm)
      end if
      l_am1(1-l_xmym:0) = ZERO
      l_am2(1-l_xmym:0) = ZERO
      l_am3(1-l_xmym:0) = ZERO
      !call pb_setupper(l_am1(1),l_am2(1),l_am3(1))
      if ( l_bcopt == 10 .or. l_bcopt == 1 ) then
!         : init_array: calling setupper2';flush(6)
        call pb_setupper2(l_am1(1:l_xmymzm),l_am2(1:l_xmymzm),l_am3(1:l_xmymzm), &
                                           l_am4(1:l_xmymzm),l_am5(1:l_xmymzm),l_am6(1:l_xmymzm))
      else
        call pb_setupper(l_am1(1:l_xmymzm),l_am2(1:l_xmymzm),l_am3(1:l_xmymzm))
      end if
      l_bv (1:l_xmymzm) = p_bv(1:l_xmymzm)
      l_pv (1-l_xmym:0) = ZERO
      l_pv (l_xmymzm+1:l_xmymzm+l_xmym) = ZERO
!WMBS - Will use case 6 for periodic iccg
   case (8)
!     ! : init_array: calling feedepsintoam';flush(6)
      !array upper and lower dimension bounds  
      lb=1-l_xmymzm+l_xmym
      ub=l_xmymzm+l_xmymzm-l_xmym

      !initialize all to zero first
      l_ad( lb:ub)=ZERO
      l_am1(lb:ub)=ZERO; l_am2(lb:ub)=ZERO; l_am3(lb:ub)=ZERO 
      l_am4(lb:ub)=ZERO; l_am5(lb:ub)=ZERO; l_am6(lb:ub)=ZERO
      l_bv( lb:ub)=ZERO
      l_tv( lb:ub)=ZERO
      l_zv( lb:ub)=ZERO
      l_pv( lb:ub)=ZERO
      l_rd( lb:ub)=ZERO; l_rd(1:l_xmymzm)=ONE

      !set l_bv
      l_bv(1:l_xmymzm)=p_bv(1:l_xmymzm)

      !set off diagonal arrays
      call feedepsintoam(l_xm, l_ym, l_zm, l_am1(1:l_xmymzm),  &
                         l_am2(1:l_xmymzm), l_am3(1:l_xmymzm), &
                                              epsx, epsy, epsz )
      !set central diagonal array
      if (l_pbkappa == ZERO) then
!        ! : init_array: calling feedepsintoad';flush(6)
         call feedepsintoad(l_xm, l_ym, l_zm, l_ad(1:l_xmymzm), epsx, epsy, epsz)
      else
         lfactor = l_epsout*(l_h*l_pbkappa)**2
!!        ! : init_array: calling feedepsintoad';flush(6)
         call feedepsintoad(l_xm, l_ym, l_zm, l_ad(1:l_xmymzm), epsx, epsy, epsz)
         l_ad(1:l_xmymzm) = lfactor*p_iv(1:l_xmymzm) + l_ad(1:l_xmymzm)
      end if
 
      !set l_am1-3 upper bounds to zero and set up periodic faces
!     ! : init_array: calling pb_setupper2';flush(6)
      if ( l_bcopt == 10 ) then
        call pb_setupper2(l_am1(1:l_xmymzm),l_am2(1:l_xmymzm),l_am3(1:l_xmymzm), &
                               l_am4(1:l_xmymzm),l_am5(1:l_xmymzm),l_am6(1:l_xmymzm))
      else
        call pb_setupper(l_am1(1:l_xmymzm),l_am2(1:l_xmymzm),l_am3(1:l_xmymzm))
      end if

      !Ensure that array paddings are zeroed
      l_ad( lb:0)=ZERO; l_ad( l_xmymzm+1:ub)=ZERO
      l_am1(lb:0)=ZERO; l_am1(l_xmymzm+1:ub)=ZERO 
      l_am2(lb:0)=ZERO; l_am2(l_xmymzm+1:ub)=ZERO 
      l_am3(lb:0)=ZERO; l_am3(l_xmymzm+1:ub)=ZERO
      l_am4(lb:0)=ZERO; l_am4(l_xmymzm+1:ub)=ZERO
      l_am5(lb:0)=ZERO; l_am5(l_xmymzm+1:ub)=ZERO
      l_am6(lb:0)=ZERO; l_am6(l_xmymzm+1:ub)=ZERO
      l_bv( lb:0)=ZERO; l_bv( l_xmymzm+1:ub)=ZERO 
      l_tv( lb:0)=ZERO; l_tv( l_xmymzm+1:ub)=ZERO
      l_zv( lb:0)=ZERO; l_zv( l_xmymzm+1:ub)=ZERO
      l_pv( lb:0)=ZERO; l_pv( l_xmymzm+1:ub)=ZERO
      l_rd( lb:0)=ZERO; l_rd( l_xmymzm+1:ub)=ZERO      
 
!     ! : init_array: array initialization complete';flush(6)
   case (2)
      l_ad  = ZERO
      l_am1 = ZERO
      l_am2 = ZERO
      l_am3 = ZERO
      l_iv  = ZERO
      l_bv  = ZERO
      l_rv  = ZERO
      l_zv  = ZERO
      l_xv  = ZERO
      l_xv(1+l_xmym:l_xmymzm+l_xmym) = p_xs(1:l_xmymzm)
      l_bv(1:l_xmymzm) = p_bv(1:l_xmymzm)

      m = 0
      do i = 1, mg_nlevel
         m = m + (mg_size(1,i)+1) * (mg_size(2,i)+1) * (mg_size(3,i)+1)
      end do
      allocate ( lepsx(1:m), lepsy(1:m), lepsz(1:m) )
      call feedepsintoam(l_xm, l_ym, l_zm, lepsx(1:l_xmymzm),  &
                         lepsy(1:l_xmymzm), lepsz(1:l_xmymzm), &
                                              epsx, epsy, epsz )
      ! This is not accurate but I haven't figured out how to solve this problem.

      lfactor = l_epsout*(l_h*l_pbkappa)**2
      l_iv(1:l_xmymzm) = p_iv(1:l_xmymzm)

      i = 1
      m = mg_index(i)
      n = mg_index_ext(i)
      lxmym = mg_size(1,i)*mg_size(2,i)
      call set_am_ad(lepsx(m),lepsy(m),lepsz(m),l_iv(m), &
                     l_am1(n+lxmym), l_am2(n+lxmym), l_am3(n+lxmym), &
                     l_ad(m), l_bz(m), &
                     mg_size(1,i),mg_size(2,i),mg_size(3,i),lfactor,l_epsout)

      do i = 2, mg_nlevel
         l = mg_index(i-1)
         m = mg_index(i)
         n = mg_index_ext(i)
         call restrict_eps_map(lepsx(l),lepsy(l),lepsz(l),mg_size(1,i-1),mg_size(2,i-1),mg_size(3,i-1), &
                               lepsx(m),lepsy(m),lepsz(m),mg_size(1,i),mg_size(2,i),mg_size(3,i))
         call restrict_iv(l_iv(l),mg_size(1,i-1),mg_size(2,i-1),mg_size(3,i-1),&
                          l_iv(m),mg_size(1,i),mg_size(2,i),mg_size(3,i) )
         lfactor = lfactor * 4
         lxmym = mg_size(1,i)*mg_size(2,i)

         call set_am_ad(lepsx(m),lepsy(m),lepsz(m),l_iv(m), &
                        l_am1(n+lxmym), l_am2(n+lxmym), l_am3(n+lxmym), &
                        l_ad(m), l_bz(m), &
                        mg_size(1,i),mg_size(2,i),mg_size(3,i),lfactor,l_epsout)
      end do
      deallocate(lepsx,lepsy,lepsz)
   end select

contains

subroutine feedepsintoam(xm, ym, zm, am1, am2, am3, eps1, eps2, eps3)
    implicit none
    integer xm, ym, zm
    _REAL_ am1(1:xm,1:ym,1:zm)
    _REAL_ am2(1:xm,1:ym,1:zm)
    _REAL_ am3(1:xm,1:ym,1:zm)
    _REAL_ eps1(0:xm,1:ym,1:zm)
    _REAL_ eps2(1:xm,0:ym,1:zm)
    _REAL_ eps3(1:xm,1:ym,0:zm)
    am1(1:xm,1:ym,1:zm) = eps1(1:xm,1:ym,1:zm)
    am2(1:xm,1:ym,1:zm) = eps2(1:xm,1:ym,1:zm)
    am3(1:xm,1:ym,1:zm) = eps3(1:xm,1:ym,1:zm)
end subroutine feedepsintoam

subroutine feedepsintoad(xm, ym, zm, ad, eps1, eps2, eps3)
    implicit none
    integer xm, ym, zm
    _REAL_ ad(1:xm,1:ym,1:zm)
    _REAL_ eps1(0:xm,1:ym,1:zm)
    _REAL_ eps2(1:xm,0:ym,1:zm)
    _REAL_ eps3(1:xm,1:ym,0:zm)
    ad(1:xm,1:ym,1:zm) =                    eps1(1:xm,  1:ym,  1:zm)
    ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps1(0:xm-1,1:ym,  1:zm)
    ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps2(1:xm,  1:ym,  1:zm)
    ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps2(1:xm,  0:ym-1,1:zm)
    ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps3(1:xm,  1:ym,  1:zm)
    ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps3(1:xm,  1:ym,  0:zm-1)
end subroutine feedepsintoad

subroutine feedepsintoam_piccg(xm, ym, zm, am1, am2, am3, eps1, eps2, eps3)
   implicit none
   integer xm,ym,zm
   _REAL_ am1(0:xm,0:ym,0:zm)
   _REAL_ am2(0:xm,0:ym,0:zm)
   _REAL_ am3(0:xm,0:ym,0:zm)
   _REAL_ eps1(0:xm,1:ym,1:zm)
   _REAL_ eps2(1:xm,0:ym,1:zm)
   _REAL_ eps3(1:xm,1:ym,0:zm)
   am1(1:xm,1:ym,1:zm) = eps1(1:xm,1:ym,1:zm)
   am2(1:xm,1:ym,1:zm) = eps2(1:xm,1:ym,1:zm)
   am3(1:xm,1:ym,1:zm) = eps3(1:xm,1:ym,1:zm)
end subroutine feedepsintoam_piccg

subroutine feedepsintoad_piccg(xm, ym, zm, ad, eps1, eps2, eps3,piv)
   implicit none
   integer xm,ym,zm
   _REAL_ ad(0:xm,0:ym,0:zm)
   _REAL_ eps1(0:xm,1:ym,1:zm)
   _REAL_ eps2(1:xm,0:ym,1:zm)
   _REAL_ eps3(1:xm,1:ym,0:zm)
   _REAL_ piv(1:xm,1:ym,1:zm)

   _REAL_ lfactor

    ad(1:xm,1:ym,1:zm) =                    eps1(1:xm,  1:ym,  1:zm)
    ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps1(0:xm-1,1:ym,  1:zm)
    ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps2(1:xm,  1:ym,  1:zm)
    ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps2(1:xm,  0:ym-1,1:zm)
    ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps3(1:xm,  1:ym,  1:zm)
    ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps3(1:xm,  1:ym,  0:zm-1)
    if (.not. (l_pbkappa == ZERO)) then
      lfactor = l_epsout*(l_h*l_pbkappa)**2
      ad(1:xm,1:ym,1:zm) = piv(1:xm,1:ym,1:zm)*lfactor + ad(1:xm,1:ym,1:zm)
    end if  
end subroutine feedepsintoad_piccg

end subroutine init_array

!===========================================================================

subroutine pb_iccg(phi,xs)

!  use poisson_boltzmann, only: level
   implicit none

! Passed variables

   _REAL_ phi(l_xmymzm), xs(1-l_xmym:l_xmymzm+l_xmym)

   !phi is solution (output) array
   !xs seems to be an initial guess / storage for the solution
   !it has additional "padding" above and below 1:l_xmymzm
   !to accomodate matrix multiplication operations

! Local variables

   logical uconvg
   integer i, j
   _REAL_ alpha, beta, pdotz, bdotb1, bdotb2

   ! initialization

   !WMBS- periodicity loops for l_rd(i) go here 
   do i = 1, l_xmymzm
   !WMBS -  This seems to modify the diagonal row (d(i)^-1) using:
   !        the formula from Luo R., David L., and Gilson M. K. J.C.P. 2002
   !        see page 1246
   !        l_rd seems to store the reciprical of the diagonal matrix D
   !        from the M=(L+D)D^-1(D+U) preconditioner splitting
      l_rd(i) = ONE/( l_ad(i) - &
         l_am1(i-1   )*(       l_am1(i-1   )+l_fmiccg*l_am2(i-1   )+l_fmiccg*l_am3(i-1   ))*l_rd(i-1   ) - &
         l_am2(i-l_xm  )*(l_fmiccg*l_am1(i-l_xm  )+       l_am2(i-l_xm  )+l_fmiccg*l_am3(i-l_xm  ))*l_rd(i-l_xm  ) - &
         l_am3(i-l_xmym)*(l_fmiccg*l_am1(i-l_xmym)+l_fmiccg*l_am2(i-l_xmym)+       l_am3(i-l_xmym))*l_rd(i-l_xmym) )
   end do

   !WMBS -  periodicity loops for 
   !        l_ad(i) and l_am arrays
   !        go here
   do i = 1, l_xmymzm
      !use l_rd to modify the l_ad terms. In the unmodified this would have been
      !skipped since the diagonal would be unity.
      l_ad(i) = l_ad(i)*l_rd(i)
      l_rd(i) = sqrt(l_rd(i))
      l_bv(i) = l_bv(i)*l_rd(i) !first step to precondition initial guess in l_bv (diagonal mult.)
      !WMBS- L matrix sub-diagonal term building
      l_am1(i-1   ) = l_am1(i-1   )*l_rd(i)*l_rd(i-1   )       !L sub-diagonal of previous x terms
      l_am2(i-l_xm  ) = l_am2(i-l_xm  )*l_rd(i)*l_rd(i-l_xm  ) !L sub-diagonal of previous y terms
      l_am3(i-l_xmym) = l_am3(i-l_xmym)*l_rd(i)*l_rd(i-l_xmym) !L sub-diagonal of previous z terms
   end do


   !WMBS- periodicity loops for l_bv go here
   l_inorm = ZERO
   do i = 1, l_xmymzm 
   !WMBS- this loop could probably be combined w/ above to a single loop to save
   !overhead
      l_ad(i) = l_ad(i) - TWO !this would yield -1 in unmodified version

      !WMBS-apply L matrix to l_bv. Will need to apply periodic face part here
      l_bv(i) = l_bv(i) + l_am1(i-1   )*l_bv(i-1   ) &
                    + l_am2(i-l_xm  )*l_bv(i-l_xm  ) &
                    + l_am3(i-l_xmym)*l_bv(i-l_xmym)
      l_inorm = l_inorm + abs(l_bv(i))
   end do

!   write(10,*) 'l_inorm, accept*l_inorm: ',l_inorm,l_accept*l_inorm

   !WMBS- periodicity loops for l_tv go here
   do i = l_xmymzm, 1, -1
      !WMBS - l_tv looks like the transpose of the lower triangular
      !       this loop is applying (D+U) part of M to initial guess xs
      !       l_tv(i) = 0 when i > l_xmymzm or i < 1 here 
      l_tv(i) = xs(i) + l_am1(i     )*l_tv(i+1   ) &
                    + l_am2(i     )*l_tv(i+l_xm  ) &
                    + l_am3(i     )*l_tv(i+l_xmym)
   end do

   !WMBS- periodicity loops for l_zv go here
   do i = 1, l_xmymzm
      l_zv(i) = xs(i) + l_ad (i     )*l_tv(i     ) &
                    + l_am1(i-1   )*l_zv(i-1   ) &
                    + l_am2(i-l_xm  )*l_zv(i-l_xm  ) &
                    + l_am3(i-l_xmym)*l_zv(i-l_xmym)
   end do

   !WMBS-   need to add periodicity for l_tv. This could be a bit
   !        if we use additional l_am4 - l_am6 with nonzero terms
   !        only on boundaries, we can simply extend l_tv setting
   !        line
   bdotb1 = ZERO
   do i = l_xmymzm, 1, -1
      l_zv(i) = l_zv(i) + l_tv(i)
      l_bv(i) = l_bv(i) - l_zv(i)

      ! iteration 0.

      bdotb1 = bdotb1 + l_bv(i)*l_bv(i)
      l_pv(i)  = l_bv(i)

      ! first step of the matrix vector multiplication, see below
      !WMBS- apply (D+U) part of M 
      l_tv(i) = l_pv(i) + l_am1(i     )*l_tv(i+1   ) &
                    + l_am2(i     )*l_tv(i+l_xm  ) &
                    + l_am3(i     )*l_tv(i+l_xmym)
   end do

   l_itn = 0
   uconvg = .true.

   ! the main loop of iccg solver
   !write(10,*) 'l_itn, l_norm,sum(l_zv),sum(l_pv),sum(l_tv)'; flush(10)

   do while ( uconvg )

      ! second and third steps of the matrix vector multiplication

      pdotz = ZERO
      do i = 1, l_xmymzm+l_xmym
      !WMBS - Apply (L+D)D^1 part of m to 
      !       goes to l_xmymzm+l_xmym to accomodate
      !       (D+U).p matrix vector multiplication op
         l_zv(i) = l_pv(i) + l_ad (i     )*l_tv(i     ) &
                       + l_am1(i-1   )*l_zv(i-1   ) &
                       + l_am2(i-l_xm  )*l_zv(i-l_xm  ) &
                       + l_am3(i-l_xmym)*l_zv(i-l_xmym)

         j = i - l_xmym
         l_zv(j) = l_zv(j) + l_tv(j)

         pdotz = pdotz + l_pv(j)*l_zv(j)
      end do
      alpha = bdotb1/pdotz

      l_norm = ZERO
      bdotb2 = ZERO
      l_itn = l_itn + 1
      do i = 1, l_xmymzm
         xs(i) = xs(i) + alpha*l_pv(i)
         l_bv(i) = l_bv(i) - alpha*l_zv(i)
         l_norm  = l_norm  + abs(l_bv(i))

         bdotb2= bdotb2+ l_bv(i)*l_bv(i)
      end do

      ! check convergence

!      write(10,*) l_itn,sngl(l_norm),sngl(sum(abs(l_zv))),sngl(sum(abs(l_pv))),sngl(sum(abs(l_tv)));flush(10)

      if ( l_itn >= l_maxitn .or. l_norm <= l_accept*l_inorm ) then

         uconvg = .false.
         if ( l_itn >= l_maxitn ) then
            write(6, *) 'PB warning in pb_miccg(): CG l_maxitn exceeded!'
         end if

      else

         beta = bdotb2/bdotb1
         bdotb1 = bdotb2

         ! first step of the matrix vector multiplication

         do i = l_xmymzm, 1, -1
            l_pv(i) = l_bv(i) + beta*l_pv(i)

            l_tv(i) = l_pv(i) + l_am1(i)*l_tv(i+1   ) &
                          + l_am2(i)*l_tv(i+l_xm  ) &
                          + l_am3(i)*l_tv(i+l_xmym)
         end do
      end if
   end do  !  while ( uconvg ), end of the main iccg loop

   ! back scaling of the solution

   do i = l_xmymzm, 1, -1
      l_tv(i)  = xs(i) + l_am1(i)*l_tv(i+1   ) &
                     + l_am2(i)*l_tv(i+l_xm  ) &
                     + l_am3(i)*l_tv(i+l_xmym)

      phi(i) = l_tv(i)*l_rd(i)
   end do

end subroutine pb_iccg

!===========================================================================

subroutine pb_piccg ( phi, xs)

!  use poisson_boltzmann, only: level
!   use poisson_boltzmann!, only: fixed_zbp,zubp,zlbp
      !importing fixed z potential directly from pb_force.
      !this is a bit inelegant, and should be fixed later
      !once testing is done
   implicit none

! Passed variables

   _REAL_ phi(l_xmymzm), xs(1-l_xmym:l_xmymzm+l_xmym)

   !phi is solution (output) array
   !xs seems to be an initial guess / storage for the solution
   !it has additional "padding" above and below 1:l_xmymzm
   !to accomodate matrix multiplication operations

! Local variables

   logical uconvg
   integer i, j, k, ii, ip, jp, kp
   integer xi, yi, zi, ub, l_itn
   _REAL_ alpha, beta, pdotz, bdotb1, bdotb2,neutcg
   ! initialization

!   _REAL_ l_zubp, l_zlbp

   ub = (1+l_xm)*(1+l_ym)*(1+l_zm)  

 
   ip = l_xm-1          !x=1 => x=xm
   jp = l_xmym-l_xm     !y=1 => y=ym
   kp = l_xmymzm-l_xmym !z=1 => z=zm

  !WMBS- load phi and xs into working arrays 
   call pb_piccg_load_xs(pl_xs,xs(1:l_xmymzm))

!Code for piccg using mixed boundary conditions commented out for now.
!   l_zubp=zubp 
!   l_zlbp=zlbp

!   if (fixed_zbp) then 


!      call pb_piccg_initialize_fixz(pl_ind,pl_rd,pl_ad,pl_am1,pl_am2,pl_am3,  &
!            pl_am4,pl_am5,pl_am6,pl_bv,pl_tv,pl_pv,pl_zv,            &
!            alpha,bdotb1,bdotb2,pdotz,beta,l_inorm,uconvg,           &
!            pl_xs,l_zubp,l_zlbp)


!      uconvg = .true.

!      call pb_piccg_main_loop_fixz(pl_ind,pl_ad,pl_tv,pl_zv,pl_pv,pl_bv,pl_am1,pl_am2,   &
!                           pl_am3,pl_am4,pl_am5,pl_am6,alpha,bdotb1,bdotb2,pdotz,&
!                           pl_xs,beta,l_norm,l_inorm,l_accept,l_itn,    &
!                           l_maxitn,uconvg,l_zubp,l_zlbp)

!      call pb_piccg_backscale_fixz(pl_ind,pl_xs,pl_phi,pl_tv,pl_rd, &
!                           pl_am1,pl_am2,pl_am3,pl_am4,pl_am5,pl_am6,l_zubp,l_zlbp)

   
!   else

      call pb_piccg_initialize(pl_ind,pl_rd,pl_ad,pl_am1,pl_am2,pl_am3,  &
            pl_am4,pl_am5,pl_am6,pl_bv,pl_tv,pl_pv,pl_zv,            &
            alpha,bdotb1,bdotb2,pdotz,beta,l_inorm,uconvg,           &
            pl_xs)

      uconvg = .true.
      call pb_piccg_main_loop(pl_ind,pl_ad,pl_tv,pl_zv,pl_pv,pl_bv,pl_am1,pl_am2,   &
                           pl_am3,pl_am4,pl_am5,pl_am6,alpha,bdotb1,bdotb2,pdotz,&
                           pl_xs,beta,l_norm,l_inorm,l_accept,l_itn,    &
                           l_maxitn,uconvg)

   ! back scaling of the solution
      call pb_piccg_backscale(pl_ind,pl_xs,pl_phi,pl_tv,pl_rd, &
                           pl_am1,pl_am2,pl_am3,pl_am4,pl_am5,pl_am6)


!   end if

   !copy updated working array copies of phi and xs back to originals
   call pb_piccg_return_xs(pl_xs,xs(1:l_xmymzm))
   call pb_piccg_return_phi(pl_phi,phi(1:l_xmymzm))

end subroutine pb_piccg
!   contains
subroutine pb_piccg_load_xs(l_xs,xs)

   implicit none
   _REAL_ l_xs(0:l_xm,0:l_ym,0:l_zm)
   _REAL_ xs(1:l_xm,1:l_ym,1:l_zm)

   l_xs= ZERO
   l_xs(1:l_xm,1:l_ym,1:l_zm)=xs(1:l_xm,1:l_ym,1:l_zm)

end subroutine

subroutine pb_piccg_initialize(llind,llrd,llad,llam1,llam2,llam3,   &
            llam4,llam5,llam6,llbv,lltv,llpv,llzv,              &
            alpha,bdotb1,bdotb2,pdotz,beta,inorm,unconv, &
            lxs)
!subroutine pb_piccg_initialize(alpha,bdotb1,bdotb2,pdotz,beta,inorm,unconv, &
!            lxs)

   integer*4 llind(1:l_xmymzm,3)
   _REAL_ alpha, bdotb1, bdotb2, pdotz, beta, inorm
   _REAL_ llrd (0:l_xm  ,0:l_ym  ,0:l_zm),llad (0:l_xm  ,0:l_ym  ,0:l_zm)
   _REAL_ llam1(0:l_xm  ,0:l_ym  ,0:l_zm),llam2(0:l_xm  ,0:l_ym  ,0:l_zm)
   _REAL_ llam3(0:l_xm  ,0:l_ym  ,0:l_zm),llam4(0:l_xm  ,0:l_ym  ,0:l_zm)
   _REAL_ llam5(0:l_xm  ,0:l_ym  ,0:l_zm),llam6(0:l_xm  ,0:l_ym  ,0:l_zm)
   _REAL_ llbv (0:l_xm  ,0:l_ym  ,0:l_zm),llpv (0:l_xm  ,0:l_ym  ,0:l_zm)
   _REAL_ llzv (0:l_xm  ,0:l_ym  ,0:l_zm),lltv (1:l_xm+1,1:l_ym+1,1:l_zm+1)
   _REAL_ lxs  (0:l_xm,0:l_ym,0:l_zm)
 
   character*14 larrayname, ldataname, lfilename
   character(len=19) arrayname,dataname 
   _REAL_ tempreal   

   logical unconv
   integer xi,yi,zi,xp,yp,zp,i,j,k,ii

   llrd(1:l_xm,1:l_ym,1:l_zm)=ONE

   !first step in setup of l_rd working array
   !WMBS -  This seems to modify the diagonal row (d(i)^-1) using:
   !        the formula from Luo R., David L., and Gilson M. K. J.C.P. 2002
   !        see page 1246
   !        l_rd seems to store the reciprical of the diagonal matrix D
   !        from the M=(L+D)D^-1(D+U) preconditioner splitting
   if ( ABS(SUM(llbv(1:l_xm,1:l_ym,1:l_zm))) > 1.0d0 ) then
      llbv(1:l_xm,1:l_ym,1:l_zm) = llbv(1:l_xm,1:l_ym,1:l_zm) -              &
                            SUM(llbv(1:l_xm,1:l_ym,1:l_zm)/(l_xmymzm))
   end if 

   do ii=1,l_xmymzm
      xi = llind(ii,1); yi = llind(ii,2); zi = llind(ii,3)
      xp=Mod(xi,l_xm)+1;yp=Mod(yi,l_ym)+1;zp=Mod(zi,l_zm)+1 
      llrd(xi,yi,zi)= ONE/( llad(xi,yi,zi)-                                     &  
                     llam1(xi-1,yi,zi)*(llam1(xi-1,yi,zi)     +l_fmiccg   *(    & 
                        llam2(xi-1,yi,zi)+ llam3(xi-1,yi,zi)) +l_fmiccg2  *(    &
                        llam4(xi-1,yi,zi)+ llam5(xi-1,yi,zi)  +llam6(xi-1,yi,zi) &
                        ))*llrd(xi-1,yi,zi) -                                  &
                     llam2(xi,yi-1,zi)*(llam2(xi,yi-1,zi)     +l_fmiccg   *(    &
                        llam1(xi,yi-1,zi)+ llam3(xi,yi-1,zi)) +l_fmiccg2  *(    &
                        llam4(xi,yi-1,zi)+ llam5(xi,yi-1,zi)  +llam6(xi,yi-1,zi) &
                        ))*llrd(xi,yi-1,zi) -                                  &
                     llam3(xi,yi,zi-1)*(llam3(xi,yi,zi-1)     +l_fmiccg   *(    &
                        llam1(xi,yi,zi-1)+ llam2(xi,yi,zi-1)) +l_fmiccg2  *(    &
                        llam4(xi,yi,zi-1)+ llam5(xi,yi,zi-1)  +llam6(xi,yi,zi-1) &
                        ))*llrd(xi,yi,zi-1) -                                  &
                     llam4(xp,yi,zi)*(llam4(xp,yi,zi)     +l_fmiccg2  *(        &
                        llam5(xp,yi,zi)+ llam6(xp,yi,zi)) +l_fmiccg   *(    &
                        llam1(xp,yi,zi)+ llam2(xp,yi,zi)  +llam3(xp,yi,zi) &
                        ))*llrd(xp,yi,zi) -                                  &
                     llam5(xi,yp,zi)*(llam5(xi,yp,zi)     +l_fmiccg2  *(        &
                        llam4(xi,yp,zi)+ llam6(xi,yp,zi)) +l_fmiccg   *(    &
                        llam1(xi,yp,zi)+ llam2(xi,yp,zi)  +llam3(xi,yp,zi) &
                        ))*llrd(xi,yp,zi) -                                  &
                     llam6(xi,yi,zp)*(llam6(xi,yi,zp)     +l_fmiccg2  *(        &
                        llam4(xi,yi,zp)+ llam5(xi,yi,zp)) +l_fmiccg   *(    &
                        llam1(xi,yi,zp)+ llam2(xi,yi,zp)  +llam3(xi,yi,zp) &
                        ))*llrd(xi,yi,zp)                                    &
                     )
      if (llrd(xi,yi,zi) .le. 0.0d0 ) then
      write(6,*) 'PB Bomb in pb_piccg_initialize: llrd corrupted at (xi,yi,zi):',xi,yi,zi
       call mexit(6,1)
      end if
   end do

   do ii=1,l_xmymzm
      xi = llind(ii,1); yi = llind(ii,2); zi = llind(ii,3)
      tempreal=llrd(xi,yi,zi)
      if ((tempreal .ne. 0.0d0) .and. (.not. ((tempreal .ge. 0.0d0) .or. (tempreal .le. 0.0d0) )).or.(tempreal.lt.0)) then
         write(6,*) 'PB Bomb in pb_piccg_initialize: llrd corrupted during first precondition pass'
         call mexit(6,1)
      end if

      llad(xi,yi,zi) = llad(xi,yi,zi)*llrd(xi,yi,zi)
      llrd(xi,yi,zi) = sqrt(llrd(xi,yi,zi))
      llbv(xi,yi,zi) = llbv(xi,yi,zi)*llrd(xi,yi,zi)!preconditions inital guess
      !WMBS - L matrix sub-diagonal term building
      llam1(xi-1,yi,zi) = llam1(xi-1,yi,zi)*llrd(xi,yi,zi)*llrd(xi-1,yi,zi)
      llam2(xi,yi-1,zi) = llam2(xi,yi-1,zi)*llrd(xi,yi,zi)*llrd(xi,yi-1,zi)
      llam3(xi,yi,zi-1) = llam3(xi,yi,zi-1)*llrd(xi,yi,zi)*llrd(xi,yi,zi-1)
   
      tempreal=llbv(xi,yi,zi)
      if ((tempreal .ne. 0.0d0) .and. (.not. ((tempreal .ge. 0.0d0) .or. (tempreal .le. 0.0d0) ))) then
         write(6,*) 'PB Bomb in pb_piccg_initialize: l_bv corrupted during first precondition pass'
         call mexit(6,1)
      end if
   end do

   do yi = 1, l_ym; do zi = 1, l_zm
      llam4(1,yi,zi)   = llam4(1,yi,zi)*llrd(1,yi,zi)*llrd(l_xm,yi,zi)
   end do; end do
   do xi = 1, l_xm; do zi = 1, l_zm
      llam5(xi,1,zi)   = llam5(xi,1,zi)*llrd(xi,1,zi)*llrd(xi,l_ym,zi)
   end do; end do
   do xi = 1, l_xm; do yi = 1, l_ym
      llam6(xi,yi,1)   = llam6(xi,yi,1)*llrd(xi,yi,1)*llrd(xi,yi,l_zm)
   end do; end do

  !Continue preconditioning of initial guess - l_bv
  !WMBS- periodicity loops for l_bv go here
   inorm = ZERO

   do ii=1,l_xmymzm
      xi = llind(ii,1); yi = llind(ii,2); zi = llind(ii,3)
      xp=MOD(xi,l_xm)+1;yp=MOD(yi,l_ym)+1;zp=MOD(zi,l_zm)+1
      llad(xi,yi,zi) = llad(xi,yi,zi) - TWO
      llbv(xi,yi,zi) = llbv(xi,yi,zi)                         &
                        +  llam1(xi-1,yi,zi)*llbv(xi-1,yi,zi) &
                        +  llam2(xi,yi-1,zi)*llbv(xi,yi-1,zi) &
                        +  llam3(xi,yi,zi-1)*llbv(xi,yi,zi-1) &
                        +  llam4(xp,yi,zi  )*llbv(xp,yi,zi  ) &
                        +  llam5(xi,yp,zi  )*llbv(xi,yp,zi  ) &
                        +  llam6(xi,yi,zp  )*llbv(xi,yi,zp  ) 
      inorm = inorm + abs(llbv(xi,yi,zi))
   end do

   !U*p explicit conditioning
   !Interior grid loops
   do ii=l_xmymzm,1,-1
      xi = llind(ii,1); yi = llind(ii,2); zi = llind(ii,3)
      xp = MOD(xi+l_xm-2,l_xm)+1
      yp = MOD(yi+l_ym-2,l_ym)+1
      zp = MOD(zi+l_zm-2,l_zm)+1
      lltv(xi,yi,zi) = lxs(xi,yi,zi)  + llam1(xi,yi,zi)*lltv(xi+1,yi,zi) &
                                    + llam2(xi,yi,zi)*lltv(xi,yi+1,zi) &
                                    + llam3(xi,yi,zi)*lltv(xi,yi,zi+1) &
                                    + llam4(xi,yi,zi)*lltv(xp,yi,zi  ) &
                                    + llam5(xi,yi,zi)*lltv(xi,yp,zi  ) &
                                    + llam6(xi,yi,zi)*lltv(xi,yi,zp  ) 
   end do   
   !Boundary grid loops
   !WMBS- loops for applying the L and D matrices
   !L*p explicit conditioning
   !Boundary Grids

   !Interior Grids
   do ii=1,l_xmymzm
      xi = llind(ii,1); yi = llind(ii,2); zi = llind(ii,3)
      xp=MOD(xi,l_xm)+1;yp=MOD(yi,l_ym)+1;zp=MOD(zi,l_zm)+1
     llzv(xi,yi,zi) =lxs(xi,yi,zi)+ llad(xi,yi,zi   )*lltv(xi,yi,zi  )  &
                                 + llam1(xi-1,yi,zi)*llzv(xi-1,yi,zi)  &
                                 + llam2(xi,yi-1,zi)*llzv(xi,yi-1,zi)  &
                                 + llam3(xi,yi,zi-1)*llzv(xi,yi,zi-1)  &
                                 + llam4(xp,yi,zi  )*llzv(xp,yi,zi  )  &
                                 + llam5(xi,yp,zi  )*llzv(xi,yp,zi  )  &
                                 + llam6(xi,yi,zp  )*llzv(xi,yi,zp  )      
   end do

   !WMBS-   need to add periodicity for l_tv. This could be a bit
   !        if we use additional l_am4 - l_am6 with nonzero terms
   !        only on boundaries, we can simply extend l_tv setting
   !        line
   bdotb1 = ZERO
   do ii=l_xmymzm,1,-1
      xi = llind(ii,1); yi = llind(ii,2); zi = llind(ii,3)
      xp=Mod(xi+l_xm-2,l_xm)+1 !wrap x=1 to x=l_xm
      yp=Mod(yi+l_ym-2,l_ym)+1 !wrap y=1 to y=l_ym
      zp=Mod(zi+l_zm-2,l_zm)+1 !wrap z=1 to z=l_zm

      llzv(xi,yi,zi)  = llzv(xi,yi,zi) + lltv(xi,yi,zi)
      llbv(xi,yi,zi)  = llbv(xi,yi,zi) - llzv(xi,yi,zi)

      bdotb1 = bdotb1 + llbv(xi,yi,zi)**2
      llpv(xi,yi,zi) = llbv(xi,yi,zi)

      !WMBS - apply (D+U) part of M
      !am4 through 6 store periodic dielectric on 0 indices
      !lltv needs to wrap to element 1 for periodicity
      lltv(xi,yi,zi) = llpv(xi,yi,zi) +llam1(xi,yi,zi) *lltv(xi+1,yi,zi) &
                                    +llam2(xi,yi,zi) *lltv(xi,yi+1,zi) &
                                    +llam3(xi,yi,zi) *lltv(xi,yi,zi+1) &
                                    +llam4(xi,yi,zi) *lltv(xp,yi,zi)   &
                                    +llam5(xi,yi,zi) *lltv(xi,yp,zi)   &
                                    +llam6(xi,yi,zi) *lltv(xi,yi,zp)
                      
   end do
   l_itn = 0
   unconv = .true.


end subroutine pb_piccg_initialize

subroutine pb_piccg_main_loop(lind,lad,ltv,lzv,lpv,lbv,           &
                              lam1,lam2,lam3,lam4,lam5,lam6,      &
                              alpha,bdotb1,bdotb2,pdotz,lxs,beta, &
                              norm,inorm,accept,itn,maxitn,unconv)
!subroutine pb_piccg_main_loop(alpha,bdotb1,bdotb2,pdotz,lxs,beta, &
!                              norm,inorm,accept,itn,maxitn,unconv)



   integer*4 lind(1:l_xmymzm,3)
   integer itn, maxitn, i, j, k, ii 
   integer xi, yi, zi, xj, yj, zj, xp, yp, zp, zpp
   logical unconv
   _REAL_ norm, accept, inorm 
   _REAL_ alpha, beta, bdotb1, bdotb2, pdotz 
   _REAL_ lad(0:l_xm,0:l_ym,0:l_zm),ltv(1:l_xm+1,1:l_ym+1,1:l_zm+1)
   _REAL_ lzv(0:l_xm,0:l_ym,0:l_zm),lpv(0:l_xm,0:l_ym,0:l_zm)
   _REAL_ lbv(0:l_xm,0:l_ym,0:l_zm)
   _REAL_ lxs(0:l_xm,0:l_ym,0:l_zm)
   _REAL_ lam1(0:l_xm,0:l_ym,0:l_zm),lam2(0:l_xm,0:l_ym,0:l_zm)
   _REAL_ lam3(0:l_xm,0:l_ym,0:l_zm),lam4(0:l_xm,0:l_ym,0:l_zm)
   _REAL_ lam5(0:l_xm,0:l_ym,0:l_zm),lam6(0:l_xm,0:l_ym,0:l_zm)

   itn=0

   do while ( unconv )
      ! second and third steps of the matrix vector multiplication
      pdotz = ZERO
     do ii=1,l_xmymzm-l_xmym
         xi = lind(ii,1); yi = lind(ii,2); zi = lind(ii,3)
         xp=MOD(xi,l_xm)+1;yp=MOD(yi,l_ym)+1;zp=MOD(zi,l_zm)+1
         lzv(xi,yi,zi) = lpv(xi,yi,zi)                        &
                           + lad(xi,yi,zi)    *ltv(xi,yi,zi  )&
                           + lam1(xi-1,yi,zi) *lzv(xi-1,yi,zi)&
                           + lam2(xi,yi-1,zi) *lzv(xi,yi-1,zi)&
                           + lam3(xi,yi,zi-1) *lzv(xi,yi,zi-1)&
                           + lam4(xp,yi,zi  ) *lzv(xp,yi,zi  )&
                           + lam5(xi,yp,zi  ) *lzv(xi,yp,zi  )&
                           + lam6(xi,yi,zp  ) *lzv(xi,yi,zp  )    
      end do

      zi=l_zm
      do xi = 1, l_xm; do yi = 1, l_ym
         xp=MOD(xi,l_xm)+1;yp=MOD(yi,l_ym)+1;zp=MOD(zi,l_zm)+1
         lzv(xi,yi,zi) = lpv(xi,yi,zi)                        &
                           + lad(xi,yi,zi)    *ltv(xi,yi,zi  )&
                           + lam1(xi-1,yi,zi) *lzv(xi-1,yi,zi)&
                           + lam2(xi,yi-1,zi) *lzv(xi,yi-1,zi)&
                           + lam3(xi,yi,zi-1) *lzv(xi,yi,zi-1)&
                           + lam4(xp,yi,zi  ) *lzv(xp,yi,zi  )&
                           + lam5(xi,yp,zi  ) *lzv(xi,yp,zi  )&
                           + lam6(xi,yi,zp  ) *lzv(xi,yi,zp  )
         lzv(xi,yi,1) = lzv(xi,yi,1) + ltv(xi,yi,1)
         pdotz = pdotz + lpv(xi,yi,1) * lzv(xi,yi,1)
      end do; end do

      lzv(1:l_xm,1:l_ym,2:l_zm) = lzv(1:l_xm,1:l_ym,2:l_zm) + &
                                  ltv(1:l_xm,1:l_ym,2:l_zm)
      pdotz = pdotz + SUM(lpv(1:l_xm,1:l_ym,2:l_zm)*lzv(1:l_xm,1:l_ym,2:l_zm))

      alpha = bdotb1/pdotz
      norm = ZERO
      bdotb2 = ZERO
      itn = itn + 1
      l_itn=itn

      do ii=1,l_xmymzm
         xi = lind(ii,1); yi = lind(ii,2); zi = lind(ii,3)
         lxs(xi,yi,zi) = lxs(xi,yi,zi) + alpha*lpv(xi,yi,zi)
         lbv(xi,yi,zi) = lbv(xi,yi,zi) - alpha*lzv(xi,yi,zi)
         norm = norm + abs(lbv(xi,yi,zi))

         bdotb2 = bdotb2 + lbv(xi,yi,zi)**2
      end do
     ! check convergence

      if ( itn >= maxitn .or. norm <= accept*inorm ) then
         !WMBS- either we converged or ran out of iterations
         unconv = .false.
         if ( itn >= maxitn ) then
            !WMBS- we ran out of iterations
            write(6, *) 'PB warning in pb_miccg(): CG l_maxitn exceeded!'
         end if

      else
         !WMBS-
         !haven't converged yet, and still have iterations left,
         !get ready for next loop
         beta = bdotb2/bdotb1
         bdotb1 = bdotb2

         ! first step of the matrix vector multiplication
         lpv(1:l_xm,1:l_ym,1:l_zm)=lbv(1:l_xm,1:l_ym,1:l_zm)+&
                                    beta*lpv(1:l_xm,1:l_ym,1:l_zm)
         do ii=l_xmymzm,1,-1
            xi = lind(ii,1); yi = lind(ii,2); zi = lind(ii,3)
            !WMBS-   need to update tv for next loop_
            !        as with startup, periodicity can't be done ahead
            !        since tv needs an updated pv which is tangled in
            !        the same loop. We will soak extra zero addition
            !        steps instead of building bv ahead of time in
            !        extra loop over all grids
            xp=Mod(xi+l_xm-2,l_xm)+1 !wrap x=1 to x=l_xm
            yp=Mod(yi+l_ym-2,l_ym)+1 !wrap y=1 to y=l_ym
            zp=Mod(zi+l_zm-2,l_zm)+1 !wrap z=1 to z=l_zm
            ltv(xi,yi,zi) = lpv(xi,yi,zi)                     &
                              + lam1(xi,yi,zi)*ltv(xi+1,yi,zi)&
                              + lam2(xi,yi,zi)*ltv(xi,yi+1,zi)&
                              + lam3(xi,yi,zi)*ltv(xi,yi,zi+1)&
                              + lam4(xi,yi,zi)*ltv(xp,yi,zi  )&
                              + lam5(xi,yi,zi)*ltv(xi,yp,zi  )&
                              + lam6(xi,yi,zi)*ltv(xi,yi,zp  )

         end do
      end if

   end do  !  while ( uconvg ), end of the main iccg loop

end subroutine pb_piccg_main_loop

subroutine pb_piccg_backscale(lind,lxs,lphi,ltv,lrd,lam1,lam2,lam3,&
                           lam4,lam5,lam6)
!subroutine pb_piccg_backscale(lxs,phi)

   integer*4 lind(1:l_xmymzm,3)
   integer xi, yi, zi, ii
   _REAL_ lxs(0:l_xm,0:l_ym,0:l_zm), lphi(0:l_xm,0:l_ym,0:l_zm)
   _REAL_ ltv(1:l_xm+1,1:l_ym+1,1:l_zm+1),lrd(0:l_xm,0:l_ym,0:l_zm)
   _REAL_ lam1(0:l_xm,0:l_ym,0:l_zm),lam2(0:l_xm,0:l_ym,0:l_zm)
   _REAL_ lam3(0:l_xm,0:l_ym,0:l_zm),lam4(0:l_xm,0:l_ym,0:l_zm)
   _REAL_ lam5(0:l_xm,0:l_ym,0:l_zm),lam6(0:l_xm,0:l_ym,0:l_zm) 

  do ii=l_xmymzm,1,-1
      xi=lind(ii,1);yi=lind(ii,2);zi=lind(ii,3)   
      ltv(xi,yi,zi) = lxs(xi,yi,zi)                                          &
                        + lam1(xi,yi,zi)          *ltv(xi+1,yi,zi)          &
                        + lam2(xi,yi,zi)          *ltv(xi,yi+1,zi)          &
                        + lam3(xi,yi,zi)          *ltv(xi,yi,zi+1)          &
                        + lam4(xi,yi,zi)*ltv(Mod(xi+l_xm-2,l_xm)+1,yi,zi)    &
                        + lam5(xi,yi,zi)*ltv(xi,Mod(yi+l_ym-2,l_ym)+1,zi)    &
                        + lam6(xi,yi,zi)*ltv(xi,yi,Mod(zi+l_ym-2,l_zm)+1)

      lphi(xi,yi,zi) = ltv(xi,yi,zi)*lrd(xi,yi,zi)

   end do

end subroutine pb_piccg_backscale

subroutine pb_piccg_return_xs(l_xs,xs)

   implicit none
   _REAL_ l_xs(0:l_xm,0:l_ym,0:l_zm)
   _REAL_ xs(1:l_xm,1:l_ym,1:l_zm)

   xs(1:l_xm,1:l_ym,1:l_zm)=l_xs(1:l_xm,1:l_ym,1:l_zm)

end subroutine

subroutine pb_piccg_return_phi(l_phi,phi)

   implicit none
   _REAL_ l_phi(0:l_xm,0:l_ym,0:l_zm)
   _REAL_ phi(1:l_xm,1:l_ym,1:l_zm)

   phi(1:l_xm,1:l_ym,1:l_zm)=l_phi(1:l_xm,1:l_ym,1:l_zm)

end subroutine

subroutine pb_piccg_loadbv(lbv,pbv)
   implicit none
   _REAL_ lbv(0:l_xm,0:l_ym,0:l_zm)
   _REAL_ pbv(1:l_xm,1:l_ym,1:l_zm)

   _REAL_ neutcg
   integer i,j,k,ii

   character(len=19) arrayname, dataname

   lbv(1:l_xm,1:l_ym,1:l_zm)=pbv(1:l_xm,1:l_ym,1:l_zm)

!   do k=1,l_zm; do j=1,l_ym; do i=1,l_xm
!      ii=i+l_xm*(j-1+l_zm*(k-1))
!      lbv(i,j,k)=pbv(ii)
!   end do; end do; end do

      if (Abs(Sum(lbv)) > 1.0d-6/l_xmymzm) then
         neutcg=sum(lbv)/l_xmymzm
         do k=1,l_zm; do j=1,l_ym; do i=1,l_xm
            lbv(i,j,k) = lbv(i,j,k) - neutcg
         end do; end do; end do
      end if

end subroutine

subroutine pb_setoutter_piccg2(l_am1,l_am2,l_am3, &
                                l_am4,l_am5,l_am6)

  _REAL_  l_am1(0:l_xm,0:l_ym,0:l_zm)
  _REAL_  l_am2(0:l_xm,0:l_ym,0:l_zm)
  _REAL_  l_am3(0:l_xm,0:l_ym,0:l_zm)
  _REAL_  l_am4(0:l_xm,0:l_ym,0:l_zm)
  _REAL_  l_am5(0:l_xm,0:l_ym,0:l_zm)
  _REAL_  l_am6(0:l_xm,0:l_ym,0:l_zm)

   integer xi,yi,zi
   do yi = 1, l_ym; do zi = 1, l_zm
      l_am4(1,yi,zi)    = l_am1(l_xm,yi,zi)
!      l_am4(l_xm,yi,zi) = l_am1(l_xm,yi,zi)
      l_am1(l_xm,yi,zi) = 0
   end do; end do

   do xi = 1, l_xm; do zi = 1, l_zm
      l_am5(xi,1,zi)    = l_am2(xi,l_ym,zi)
!      l_am5(xi,l_ym,zi) = l_am2(xi,l_ym,zi)
      l_am2(xi,l_ym,zi) = 0
   end do; end do

   do xi = 1, l_xm; do yi = 1, l_ym
      l_am6(xi,yi,1)    = l_am3(xi,yi,l_zm)
!      l_am6(xi,yi,l_zm) = l_am3(xi,yi,l_zm)
      l_am3(xi,yi,l_zm) = 0
   end do; end do

end subroutine pb_setoutter_piccg2

!end subroutine pb_piccg

!===========================================================================


subroutine pb_sor(phi,xs)

! Passed variables

   _REAL_ phi(l_xmymzm), xs(1-l_xmym:l_xmymzm+l_xmym)

! Local variables

   logical uconvg
   integer ii
   _REAL_ wsor1

!  calculate initial norm

   l_inorm = sum( abs( l_bv(1:l_xmymzm) ) )
   wsor1 = ONE - l_wsor
   uconvg = .true.
   l_itn = 0

   l_zv(1:l_xmymzm) = ONE/l_ad(1:l_xmymzm)

   do while ( uconvg )

      do ii = 1, l_xmymzm
         xs(ii) = wsor1*xs(ii) + l_wsor * (l_am1(ii-1     ) * xs(ii-1   ) + &
                                           l_am1(ii       ) * xs(ii+1   ) + &
                                           l_am2(ii-l_xm  ) * xs(ii-l_xm  ) + &
                                           l_am2(ii       ) * xs(ii+l_xm  ) + &
                                           l_am3(ii-l_xmym) * xs(ii-l_xmym) + &
                                           l_am3(ii       ) * xs(ii+l_xmym) + &
                                           l_bv(ii        )  ) * l_zv(ii)
      end do

      l_itn = l_itn + 1

      do ii = 1,l_xmymzm
         phi(ii) = l_am1(ii-1) * xs(ii-1) + l_am1(ii)* xs(ii+1) + &
                   l_am2(ii-l_xm) * xs(ii-l_xm) + l_am2(ii)*xs(ii+l_xm) + &
                   l_am3(ii-l_xmym) * xs(ii-l_xmym) + l_am3(ii)*xs(ii+l_xmym) + &
                   l_bv(ii) - l_ad(ii)* xs(ii)
      end do
      l_norm = sum(abs(phi(1:l_xmymzm)))

      if ( l_itn >= l_maxitn .or. l_norm <= l_accept*l_inorm ) then
         uconvg = .false.
         if ( l_itn .ge. l_maxitn ) then
            write(6, *) 'PBMD WARNING: SOR maxitn exceeded!'
         end if
      end if

   end do

   phi(1:l_xmymzm) = xs(1:l_xmymzm)

end subroutine pb_sor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! pb_pcg
!
! CG core routine for linearized FDPB equation.  
!
! Authors:
! Ray Luo, Jun Wang
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine pb_pcg ( phi, xs )
! Passed variables
      _REAL_ phi(l_xmymzm), xs(1-l_xmym:l_xmymzm+l_xmym)
!
! Local variables
!
      logical uconvg
      integer i, j, k, ii, jj
      _REAL_ alpha, beta, pdotz, bdotb1, bdotb2

!
! begin code
!
!
! initialization
!
!   compute ||b||
!
!      : pb_pcg: setting l_inorm';flush(6)
      l_inorm = sum(ABS(l_bv(1:l_xmymzm)))

!        write(5201,*) l_am1(1:l_xmymzm)
!        write(5202,*) l_am2(1:l_xmymzm)
!        write(5203,*) l_am3(1:l_xmymzm)
!        write(5204,*) l_am4(1:l_xmymzm)
!        write(5205,*) l_am5(1:l_xmymzm)
!        write(5206,*) l_am6(1:l_xmymzm)
!        write(5207,*) l_ad(1:l_xmymzm)
!        write(5208,*) l_bv(1:l_xmymzm)
!        write(5209,*) l_tv(1:l_xmymzm)
!        write(5210,*) l_zv(1:l_xmymzm)
!        write(5211,*) l_pv(1:l_xmymzm)
!        stop

!     write(6, *)  'itn & norm ', 0, inorm
!
!   compute b - A * x(0) and save it in r(0)
!   p(0) = r(0)
!
! iteration 0:
!   compute <r(0),r(0)>
!
!      : pb_pcg: checking net charge';flush(6)
      if (Abs(Sum(l_bv)) > 1.0d-6/l_xmymzm) then
        l_bv(1:l_xmymzm) = l_bv(1:l_xmymzm) - sum(l_bv)/l_xmymzm
      end if 

      l_itn = 0
      bdotb1 = ZERO
! WJ
!      : pb_pcg: running edge periodicity loops';flush(6)
      if ( l_bcopt == 10 .or. l_bcopt == 1 ) then
         do j = 1, l_ym; do k = 1, l_zm
!            : pb_pcg: xface (y,z)',j,k;flush(6)
            ii = 1+(j-1)*l_xm+(k-1)*l_xmym
            jj = ii + l_xm - 1
            l_bv(ii)=l_bv(ii)+l_am4(ii)*xs(jj)
            l_bv(jj)=l_bv(jj)+l_am4(ii)*xs(ii)
         end do; end do
         do i = 1, l_xm; do k = 1, l_zm
            ii = i+(k-1)*l_xmym
            jj = ii + l_xmym - l_xm
            l_bv(ii)=l_bv(ii)+l_am5(ii)*xs(jj)
            l_bv(jj)=l_bv(jj)+l_am5(ii)*xs(ii)
          end do; end do
          do i = 1, l_xm; do j = 1, l_ym
             ii = i+(j-1)*l_xm
             jj = ii + l_xmymzm - l_xmym
             l_bv(ii)=l_bv(ii)+l_am6(ii)*xs(jj)
             l_bv(jj)=l_bv(jj)+l_am6(ii)*xs(ii)
          end do; end do
       end if
!
!      : pb_pcg: running central update loop';flush(6)
      do i = 1,l_xmymzm
         l_bv(i)  = l_bv(i)  + l_am3(i-l_xmym)*xs(i-l_xmym) &
                         + l_am2(i-l_xm  )*xs(i-l_xm  ) &
                         + l_am1(i-1   )*xs(i-1   ) &
                         - l_ad(i)      *xs(i     ) &
                         + l_am1(i     )*xs(i+1   ) &
                         + l_am2(i     )*xs(i+l_xm  ) &
                         + l_am3(i     )*xs(i+l_xmym)
         bdotb1 = bdotb1 + l_bv(i)*l_bv(i)
         l_pv(i)  = l_bv(i)
      end do

!      do i=1,l_xm; do j=1,l_ym; do k=1,l_zm
!      ii=i+(j-1)*l_xm+(k-1)*l_xmym
!         write(5100+l_itn,*) ii,i,j,k,l_bv(ii)
!      end do; end do; end do

!
!
! the main loop of the CG solver
!
!      : pb_pcg: starting convergence loop';flush(6)
      uconvg = .true.
      do while ( uconvg )
!
! iteration i:
!
!   compute Ap(i) = A * p(i)
!   compute alpha(i) = <r(i),r(i)>/<p(i),Ap(i)>
!
         pdotz = ZERO
! WJ
         l_zv = ZERO
         if ( l_bcopt == 10 .or. l_bcopt == 1 ) then
            do j = 1, l_ym; do k = 1, l_zm
               ii = 1+(j-1)*l_xm+(k-1)*l_xmym
               jj = ii + l_xm - 1
               l_zv(ii)=l_zv(ii)-l_am4(ii)*l_pv(jj)
               l_zv(jj)=l_zv(jj)-l_am4(ii)*l_pv(ii)
            end do; end do
            do i = 1, l_xm; do k = 1, l_zm
               ii = i+(k-1)*l_xmym
               jj = ii + l_xmym - l_xm
               l_zv(ii)=l_zv(ii)-l_am5(ii)*l_pv(jj)
               l_zv(jj)=l_zv(jj)-l_am5(ii)*l_pv(ii)
            end do; end do
            do i = 1, l_xm; do j = 1, l_ym
               ii = i+(j-1)*l_xm
               jj = ii + l_xmymzm - l_xmym
               l_zv(ii)=l_zv(ii)-l_am6(ii)*l_pv(jj)
               l_zv(jj)=l_zv(jj)-l_am6(ii)*l_pv(ii)
            end do; end do
         end if
!
         do i = 1,l_xmymzm
            l_zv(i) =   l_ad(i)      *l_pv(i)      &
                    - l_am3(i-l_xmym)*l_pv(i-l_xmym) &
                    - l_am2(i-l_xm  )*l_pv(i-l_xm  ) &
                    - l_am1(i-1   )*l_pv(i-1   ) &
                    - l_am1(i     )*l_pv(i+1   ) &
                    - l_am2(i     )*l_pv(i+l_xm  ) &
                    - l_am3(i     )*l_pv(i+l_xmym) &
! WJ
                    + l_zv(i)
!
            pdotz = pdotz + l_pv(i)*l_zv(i)
         end do

         l_itn = l_itn + 1

!         if (l_itn<10) then
!            do i=1,l_xm; do j=1,l_ym; do k=1,l_zm
!               ii=i+(j-1)*l_xm+(k-1)*l_xmym
!               write(5100+l_itn,*) ii,i,j,k,l_zv(ii)
!            end do; end do; end do
!         else
!            stop
!         end if
!
!   update x(i+1) = x(i) + alpha(i) p(i)
!          r(i+1) = r(i) - alpha(i) Ap(i)
!
         alpha  = bdotb1/pdotz
         l_norm   = ZERO
         bdotb2 = ZERO
         do i = 1,l_xmymzm
            xs(i)       = xs(i)       + alpha*l_pv(i)
            l_bv(i)       = l_bv(i)       - alpha*l_zv(i)
            l_norm        = l_norm        +   ABS(l_bv(i))
!
!   compute beta(i) = <r(i+1),r(i+1)>/<r(i),r(i)>, part one
!
            bdotb2      = bdotb2      + l_bv(i)*l_bv(i)
         end do
!        write(6, *)  'l_itn & l_norm ',l_itn, l_norm
!
!   check convergence
!
!         write(5000,*) l_itn,l_norm

!         if ( l_itn > l_inorm*10 ) then
!            write(6,*) '-PB BOMB: pb_piccg: l_norm appears to be diverging!';flush(6)
!            call mexit(6,1)
!         end if

         if ( l_itn .ge. l_maxitn .or. l_norm .le. l_accept*l_inorm ) then

            uconvg = .false.
            if ( l_itn .ge. l_maxitn ) then
               write(6, *) 'PBMD WARNING: CG maxitn exceeded!'
            endif

         else
!
!   compute beta(i) = <r(i+1),r(i+1)>/<r(i),r(i)>, part two
!
            beta   = bdotb2/bdotb1
            bdotb1 = bdotb2
!
!   update p(i+1) = r(i+1) + beta(i) p(i)
!
            l_pv(1:l_xmymzm) = l_bv(1:l_xmymzm) + beta*l_pv(1:l_xmymzm)
         endif
      enddo
!      : pb_pcg: finished convergence loop';flush(6)
!
! end of the main CG loop
!
!
      phi(1:l_xmymzm) = xs(1:l_xmymzm)
!
!
      return

end subroutine pb_pcg

!===========================================================================

subroutine pb_cg(phi,xs)

    implicit none

! Passed variables

   _REAL_ phi(l_xmymzm), xs(1-l_xmym:l_xmymzm+l_xmym)

! Local variables

   logical uconvg
   integer i, j
   _REAL_ alpha, beta, pdotz, bdotb1, bdotb2

   l_inorm = sum(abs(l_bv(1:l_xmymzm)))

!  compute b - A * x(0) and save it in r(0)
!  p(0) = r(0)
!
!  iteration 0:
!  compute <r(0),r(0)>
!
   l_itn = 0
   bdotb1 = ZERO
   do i = 1,l_xmymzm
      l_bv(i)  = l_bv(i)  + l_am3(i-l_xmym)*xs(i-l_xmym) &
                      + l_am2(i-l_xm  )*xs(i-l_xm  ) &
                      + l_am1(i-1   )*xs(i-1   ) &
                      - l_ad(i)      *xs(i     ) &
                      + l_am1(i     )*xs(i+1   ) &
                      + l_am2(i     )*xs(i+l_xm  ) &
                      + l_am3(i     )*xs(i+l_xmym)
      bdotb1 = bdotb1 + l_bv(i)*l_bv(i)
      l_pv(i)  = l_bv(i)
   end do
!
! the main loop of the CG solver
!
   uconvg = .true.
   do while ( uconvg )
!
! iteration i:
!
!   compute Ap(i) = A * p(i)
!   compute alpha(i) = <r(i),r(i)>/<p(i),Ap(i)>
!
      pdotz = ZERO
      do i = 1,l_xmymzm
         l_zv(i) = - l_am3(i-l_xmym)*l_pv(i-l_xmym) &
                 - l_am2(i-l_xm  )*l_pv(i-l_xm  ) &
                 - l_am1(i-1   )*l_pv(i-1   ) &
                 - l_am1(i     )*l_pv(i+1   ) &
                 - l_am2(i     )*l_pv(i+l_xm  ) &
                 - l_am3(i     )*l_pv(i+l_xmym) &
                 + l_ad(i)      *l_pv(i)  
         pdotz = pdotz + l_pv(i)*l_zv(i)
      end do
!
! iteration i+1:
!
      l_itn = l_itn + 1
!
!   update x(i+1) = x(i) + alpha(i) p(i)
!          r(i+1) = r(i) - alpha(i) Ap(i)
!
      alpha  = bdotb1/pdotz
      l_norm   = ZERO
      bdotb2 = ZERO
      do i = 1,l_xmymzm
         xs(i) = xs(i) + alpha*l_pv(i)
         l_bv(i) = l_bv(i) - alpha*l_zv(i)
         l_norm  = l_norm  +  abs(l_bv(i))
!
!   compute beta(i) = <r(i+1),r(i+1)>/<r(i),r(i)>, part one
!
         bdotb2 = bdotb2 + l_bv(i)*l_bv(i)
      end do
!     write(6, *)  'itn & norm ',l_itn, l_norm
!
!   check convergence
!
      if ( l_itn .ge. l_maxitn .or. l_norm .le. l_accept*l_inorm ) then

         uconvg = .false.
         if ( l_itn .ge. l_maxitn ) then
            write(6, *) 'PBMD WARNING: CG l_maxitn exceeded!'
         endif

      else
!
!   compute beta(i) = <r(i+1),r(i+1)>/<r(i),r(i)>, part two
!
         beta   = bdotb2/bdotb1
         bdotb1 = bdotb2
!
!   update p(i+1) = r(i+1) + beta(i) p(i)
!
         l_pv(1:l_xmymzm) = l_bv(1:l_xmymzm) + beta*l_pv(1:l_xmymzm)
      endif
   enddo
!
! end of the main CG loop
!
   phi(1:l_xmymzm) = xs(1:l_xmymzm)

end subroutine pb_cg

!===========================================================================

subroutine pb_mg(phi,xs)

   implicit none

! Passed variables

   _REAL_ phi(l_xmymzm), xs(1-l_xmym:l_xmymzm+l_xmym)

! Local variables
   logical uconvg
   integer j,lxly,lxlylz,p1,p2

   l_inorm = sum(abs(l_bv(1:l_xmymzm)))

   l_itn = 0
   uconvg = .true.

   do while ( uconvg )
      mg_onorm = 9.9d99
      do j = 1, mg_nlevel-1
         lxly = mg_size(1,j)*mg_size(2,j)
         lxlylz = lxly*mg_size(3,j)
         call relax(l_xv(mg_x_idx(j)),mg_size(1,j),mg_size(2,j),mg_size(3,j), &
                    l_am1(mg_index_ext(j)),l_am2(mg_index_ext(j)),l_am3(mg_index_ext(j)), &
                    l_ad(mg_index(j)),l_bv(mg_index(j)),l_rv(mg_index(j)), &
                    lxly,lxlylz,ncyc_before,l_accept,mg_onorm(j) )
         call restrict_bv(l_rv(mg_index(j)), mg_size(1,j), mg_size(2,j), mg_size(3,j), &
                          l_bv(mg_index(j+1)), mg_size(1,j+1),mg_size(2,j+1),mg_size(3,j+1) )
         l_xv(mg_x_idx(j+1):mg_x_idx(j+2)-1) = ZERO
      end do
      lxly = mg_size(1,j)*mg_size(2,j)
      lxlylz = lxly*mg_size(3,j)
      call relax(l_xv(mg_x_idx(j)),mg_size(1,j),mg_size(2,j),mg_size(3,j), &
                 l_am1(mg_index_ext(j)),l_am2(mg_index_ext(j)),l_am3(mg_index_ext(j)), &
                 l_ad(mg_index(j)),l_bv(mg_index(j)),l_rv(mg_index(j)), &
                 lxly,lxlylz,-1,l_accept,mg_onorm(j) )
      do j = mg_nlevel-1, 1, -1
         lxly = mg_size(1,j)*mg_size(2,j)
         lxlylz = lxly*mg_size(3,j)
         p1 = mg_x_idx(j+1)+mg_size(1,j+1)*mg_size(2,j+1)
         p2 = mg_x_idx(j)+lxly
         call interpolate(l_xv(p1), mg_size(1,j+1), mg_size(2,j+1), mg_size(3,j+1),&
                          l_xv(p2), mg_size(1,j) ,  mg_size(2,j), mg_size(3,j) , &
                          l_am1(mg_index_ext(j)),l_am2(mg_index_ext(j)), &
                          l_am3(mg_index_ext(j)),l_bz(mg_index(j)),l_epsout)
         call relax(l_xv(mg_x_idx(j)),mg_size(1,j),mg_size(2,j),mg_size(3,j), &
                    l_am1(mg_index_ext(j)),l_am2(mg_index_ext(j)),l_am3(mg_index_ext(j)), &
                    l_ad(mg_index(j)),l_bv(mg_index(j)), l_rv(mg_index(j)), &
                    lxly,lxlylz,ncyc_after,l_accept,mg_onorm(j) )
      end do
      l_itn = l_itn + 1
      l_norm = sum(abs(l_rv(1:l_xmymzm)))
!  
!  convergence
!  
      if ( l_itn .ge. l_maxitn .or. l_norm .le. l_inorm*l_accept ) then
         uconvg = .false.
!        write(6, *)  '   Multigrid itn & norm', litn, lnorm
         if ( l_itn .ge. l_maxitn ) then
            write(6, *) 'PBMD WARNING: Multigrid maxitn exceeded!'
         endif
      end if
   end do

   xs (1:l_xmymzm) = l_xv(l_xmym+1:l_xmym+l_xmymzm)
   phi(1:l_xmymzm) = xs (1:l_xmymzm)

end subroutine pb_mg

!===========================================================================

subroutine relax(xs,nx,ny,nz,lam1,lam2,lam3,lad,lbv,lrv,nxny,nxnynz,ncyc,accept,onorm)
 
   implicit none

   integer nx,ny,nz,nxny,nxnynz,ncyc
   _REAL_ xs(1-nxny:nxnynz+nxny), lam1(1-nxny:nxnynz), lam2(1-nxny:nxnynz), lam3(1-nxny:nxnynz)
   _REAL_ lad(1:nxnynz), lbv(1:nxnynz), lrv(1:nxnynz)
   _REAL_ accept, onorm

   logical luconvg
   integer ii,litn,itmax
   _REAL_ wsor, wsor1, linorm, lnorm
   integer itn_checknorm

   if (ncyc>0) then
      itn_checknorm = ncyc
      wsor = 1.0d0
      wsor1 = ONE - wsor
   else
      itn_checknorm = 10
      wsor = 1.9d0
      wsor1 = ONE - wsor
   end if

   linorm = sum(abs(lbv(1:nxnynz)))
!  write(6, *)  '      relax itn & norm ', litn, linorm!, onorm
   l_zv(1:nxnynz) = ONE/lad(1:nxnynz)

   litn = 0
!  lnorm = 0
   itmax = 1000
   luconvg = .true.

   do while ( luconvg )

!     start the sor iteration ...
 
      do ii = 1, nxnynz, 1
         xs(ii) = wsor1*xs(ii) + wsor * (lam1(ii-1   ) * xs(ii-1   ) + &
                                         lam1(ii     ) * xs(ii+1   ) + &
                                         lam2(ii-nx  ) * xs(ii-nx  ) + &
                                         lam2(ii     ) * xs(ii+nx  ) + &
                                         lam3(ii-nxny) * xs(ii-nxny) + &
                                         lam3(ii     ) * xs(ii+nxny) + &
                                          lbv(ii     )             ) * l_zv(ii)
      end do

      litn = litn + 1

      if ( mod(litn,itn_checknorm) == 0 ) then
         do ii = 1,nxnynz
            lrv(ii) = lam1(ii-1) * xs(ii-1) + lam1(ii)* xs(ii+1) + &
                     lam2(ii-nx) * xs(ii-nx) + lam2(ii)*xs(ii+nx) + &
                     lam3(ii-nxny) * xs(ii-nxny) + lam3(ii)*xs(ii+nxny) + &
                     lbv(ii) - lad(ii)* xs(ii)
         end do
         lnorm = sum(abs(lrv(1:nxnynz)))

!        write(6, *)  '      relax itn & norm ', litn, lnorm!, onorm
!  
!        check convergence
!  
         if ( litn .ge. itmax .or. ( ncyc .gt. 0 .and. (litn .ge. ncyc .and. lnorm < onorm ) ) &
              .or. lnorm .le. accept*linorm ) then

            luconvg = .false.
! WJ
            if ( ncyc .gt. 0 .and. litn .ge. ncyc .and. lnorm > onorm ) then
               write(6,*) "MG FAILED: ncyc, itn, norm, onorm", ncyc, litn, lnorm, onorm
               stop
            end if
!
            if ( ncyc .gt. 0 ) onorm = lnorm 
            if ( litn .ge. itmax ) then
               write(6, *) 'PBMD WARNING: SOR maxitn exceeded!'
            endif
         end if
      end if

   end do

end subroutine relax

!===========================================================================
      
subroutine set_am_ad( epsx,epsy,epsz,iv,lam1,lam2,lam3,lad,lbz, &
                      xn,yn,zn,lfactor,epsout )

    implicit none

   _REAL_  lfactor, epsout

   integer xn,yn,zn
   _REAL_  epsx(xn,yn,zn), epsy(xn,yn,zn), epsz(xn,yn,zn), iv(xn,yn,zn)
   _REAL_  lam1(xn,yn,zn),lam2(xn,yn,zn),lam3(xn,yn,zn)
   _REAL_  lad(xn,yn,zn),lbz(xn,yn,zn)

   integer  i,j,k,i1,j1,k1
   _REAL_ lam1t,lam2t,lam3t

   lam1(1:xn,1:yn,1:zn) = epsx(1:xn,1:yn,1:zn)
   lam2(1:xn,1:yn,1:zn) = epsy(1:xn,1:yn,1:zn)
   lam3(1:xn,1:yn,1:zn) = epsz(1:xn,1:yn,1:zn)

   do k = 1, zn
      k1 = k-1
      do j = 1, yn
         j1 = j-1
         do i = 1, xn
            i1 = i-1
            lam1t = epsout
            if ( i1 /= 0 ) lam1t = lam1(i1,j,k)
            lam2t = epsout
            if ( j1 /= 0 ) lam2t = lam2(i,j1,k)
            lam3t = epsout
            if ( k1 /= 0 ) lam3t = lam3(i,j,k1)
            lad(i,j,k) = lam1t + lam1(i,j,k) + lam2t + lam2(i,j,k) &
                       + lam3t + lam3(i,j,k) 
         end do
      end do
   end do

   lbz = lfactor*iv
   lad = lad + lbz

   do k = 1, zn; do j = 1, yn
      lam1(xn,j,k) = ZERO
   end do; end do
   do k = 1, zn; do i = 1, xn
      lam2(i,yn,k) = ZERO
   end do; end do
   do j = 1, yn; do i = 1, xn
      lam3(i,j,zn) = ZERO
   end do; end do

end subroutine set_am_ad

!===========================================================================

subroutine restrict_eps_map( epsx, epsy, epsz, xn, yn, zn, epsxr, epsyr, epszr, xnr, ynr, znr )

   implicit none
 
   integer xn, yn, zn, xnr, ynr, znr
   _REAL_ epsx(xn,yn,zn),epsxr(xnr,ynr,znr)
   _REAL_ epsy(xn,yn,zn),epsyr(xnr,ynr,znr)
   _REAL_ epsz(xn,yn,zn),epszr(xnr,ynr,znr)

   integer i,j,k,i2,j2,k2

   if ( xn /= xnr*2 + 1 .or. yn /= ynr*2 + 1 .or. zn /= znr*2 + 1 ) then
      write (6,*) "Restriction of eps map Failed because of incorrect dimension"
      stop
   end if
   do k = 1, znr
      k2 = 2*k
      do j = 1, ynr
         j2 = 2*j
         do i = 1, xnr
            i2 = 2*i
            epsxr(i,j,k) = hmav(epsx(i2  ,j2  ,k2  ),epsx(i2+1,j2  ,k2  ))/4.d0 + &
                          (hmav(epsx(i2  ,j2-1,k2  ),epsx(i2+1,j2-1,k2  )) + &
                           hmav(epsx(i2  ,j2+1,k2  ),epsx(i2+1,j2+1,k2  )) + &
                           hmav(epsx(i2  ,j2  ,k2-1),epsx(i2+1,j2  ,k2-1)) + &
                           hmav(epsx(i2  ,j2  ,k2+1),epsx(i2+1,j2  ,k2+1)))/8.d0 + &
                          (hmav(epsx(i2  ,j2-1,k2-1),epsx(i2+1,j2-1,k2-1)) + &
                           hmav(epsx(i2  ,j2+1,k2-1),epsx(i2+1,j2+1,k2-1)) + &
                           hmav(epsx(i2  ,j2-1,k2+1),epsx(i2+1,j2-1,k2+1)) + &
                           hmav(epsx(i2  ,j2+1,k2+1),epsx(i2+1,j2+1,k2+1)))/16.d0 

            epsyr(i,j,k) = hmav(epsy(i2  ,j2  ,k2  ),epsy(i2  ,j2+1,k2  ))/4.d0 + &
                          (hmav(epsy(i2-1,j2  ,k2  ),epsy(i2-1,j2+1,k2  )) + &
                           hmav(epsy(i2+1,j2  ,k2  ),epsy(i2+1,j2+1,k2  )) + &
                           hmav(epsy(i2  ,j2  ,k2-1),epsy(i2  ,j2+1,k2-1)) + &
                           hmav(epsy(i2  ,j2  ,k2+1),epsy(i2  ,j2+1,k2+1)))/8.d0 + &
                          (hmav(epsy(i2-1,j2  ,k2-1),epsy(i2-1,j2+1,k2-1)) + &
                           hmav(epsy(i2+1,j2  ,k2-1),epsy(i2+1,j2+1,k2-1)) + &
                           hmav(epsy(i2-1,j2  ,k2+1),epsy(i2-1,j2+1,k2+1)) + &
                           hmav(epsy(i2+1,j2  ,k2+1),epsy(i2+1,j2+1,k2+1)))/16.d0 

            epszr(i,j,k) = hmav(epsz(i2  ,j2  ,k2  ),epsz(i2  ,j2  ,k2+1))/4.d0 + &
                          (hmav(epsz(i2  ,j2-1,k2  ),epsz(i2  ,j2-1,k2+1)) + &
                           hmav(epsz(i2  ,j2+1,k2  ),epsz(i2  ,j2+1,k2+1)) + &
                           hmav(epsz(i2-1,j2  ,k2  ),epsz(i2-1,j2  ,k2+1)) + &
                           hmav(epsz(i2+1,j2  ,k2  ),epsz(i2+1,j2  ,k2+1)))/8.d0 + &
                          (hmav(epsz(i2-1,j2-1,k2  ),epsz(i2-1,j2-1,k2+1)) + &
                           hmav(epsz(i2-1,j2+1,k2  ),epsz(i2-1,j2+1,k2+1)) + &
                           hmav(epsz(i2+1,j2-1,k2  ),epsz(i2+1,j2-1,k2+1)) + &
                           hmav(epsz(i2+1,j2+1,k2  ),epsz(i2+1,j2+1,k2+1)))/16.d0 

         end do
      end do
   end do

contains 

function hmav(a,b)

   _REAL_ hmav, a, b

   hmav = 2.d0*a*b/(a+b)

end function hmav

function hmav4(a,b,c,d)

   _REAL_ hmav4, a,b,c,d

   hmav4 = 4.d0/(1.d0/a+1.d0/b+1.d0/c+1.d0/d)

end function hmav4

end subroutine restrict_eps_map

!===========================================================================

subroutine restrict_iv(ivf, xn, yn, zn, ivr, xnr, ynr, znr)

   implicit none

   integer xn, yn, zn, xnr, ynr, znr
   _REAL_ ivf(xn,yn,zn),ivr(xnr,ynr,znr)

   integer i,j,k,i2,j2,k2

   if ( xn /= xnr*2 + 1 .or. yn /= ynr*2 + 1 .or. zn /= znr*2 + 1 ) then
      write (6,*) "Restriction of iv Failed because of incorrect dimension"
      stop
   end if

   do k = 1, znr
      k2 = k * 2
      do j = 1, ynr
         j2 = j * 2
         do i = 1, xnr
            i2 = i * 2
            ivr(i,j,k) =       ( ivf(i2-1,j2-1,k2-1) + TWO * ivf(i2  ,j2-1,k2-1) + ivf(i2+1,j2-1,k2-1) ) +  &
                        TWO  * ( ivf(i2-1,j2  ,k2-1) + TWO * ivf(i2  ,j2  ,k2-1) + ivf(i2+1,j2  ,k2-1) ) +  &
                               ( ivf(i2-1,j2+1,k2-1) + TWO * ivf(i2  ,j2+1,k2-1) + ivf(i2+1,j2+1,k2-1) ) +  &
                        TWO  * ( ivf(i2-1,j2-1,k2  ) + TWO * ivf(i2  ,j2-1,k2  ) + ivf(i2+1,j2-1,k2  ) ) +  &
                        FOUR * ( ivf(i2-1,j2  ,k2  ) + TWO * ivf(i2  ,j2  ,k2  ) + ivf(i2+1,j2  ,k2  ) ) +  &
                        TWO  * ( ivf(i2-1,j2+1,k2  ) + TWO * ivf(i2  ,j2+1,k2  ) + ivf(i2+1,j2+1,k2  ) ) +  &
                               ( ivf(i2-1,j2-1,k2+1) + TWO * ivf(i2  ,j2-1,k2+1) + ivf(i2+1,j2-1,k2+1) ) +  &
                        TWO  * ( ivf(i2-1,j2  ,k2+1) + TWO * ivf(i2  ,j2  ,k2+1) + ivf(i2+1,j2  ,k2+1) ) +  &
                               ( ivf(i2-1,j2+1,k2+1) + TWO * ivf(i2  ,j2+1,k2+1) + ivf(i2+1,j2+1,k2+1) )
            ivr(i,j,k) = ivr(i,j,k) / 64.d0
         end do
      end do
   end do

end subroutine restrict_iv

!===========================================================================

subroutine restrict_bv(bvf, xn, yn, zn, bvr, xnr, ynr, znr)
 
   implicit none 

   integer xn, yn, zn, xnr, ynr, znr
   _REAL_ bvf(xn,yn,zn),bvr(xnr,ynr,znr)

   integer i,j,k,i2,j2,k2

   if ( xn /= xnr*2 + 1 .or. yn /= ynr*2 + 1 .or. zn /= znr*2 + 1 ) then
      write (6,*) "Restriction of bv Failed because of incorrect dimension "
      stop
   end if

   do k = 1, znr
      k2 = k * 2
      do j = 1, ynr
         j2 = j * 2
         do i = 1, xnr
            i2 = i * 2
            bvr(i,j,k) =       ( bvf(i2-1,j2-1,k2-1) + TWO * bvf(i2  ,j2-1,k2-1) + bvf(i2+1,j2-1,k2-1) ) +  &
                        TWO  * ( bvf(i2-1,j2  ,k2-1) + TWO * bvf(i2  ,j2  ,k2-1) + bvf(i2+1,j2  ,k2-1) ) +  &
                               ( bvf(i2-1,j2+1,k2-1) + TWO * bvf(i2  ,j2+1,k2-1) + bvf(i2+1,j2+1,k2-1) ) +  &
                        TWO  * ( bvf(i2-1,j2-1,k2  ) + TWO * bvf(i2  ,j2-1,k2  ) + bvf(i2+1,j2-1,k2  ) ) +  &
                        FOUR * ( bvf(i2-1,j2  ,k2  ) + TWO * bvf(i2  ,j2  ,k2  ) + bvf(i2+1,j2  ,k2  ) ) +  &
                        TWO  * ( bvf(i2-1,j2+1,k2  ) + TWO * bvf(i2  ,j2+1,k2  ) + bvf(i2+1,j2+1,k2  ) ) +  &
                               ( bvf(i2-1,j2-1,k2+1) + TWO * bvf(i2  ,j2-1,k2+1) + bvf(i2+1,j2-1,k2+1) ) +  &
                        TWO  * ( bvf(i2-1,j2  ,k2+1) + TWO * bvf(i2  ,j2  ,k2+1) + bvf(i2+1,j2  ,k2+1) ) +  &
                               ( bvf(i2-1,j2+1,k2+1) + TWO * bvf(i2  ,j2+1,k2+1) + bvf(i2+1,j2+1,k2+1) )
            bvr(i,j,k) = bvr(i,j,k) / 16.d0
         end do
      end do
   end do

end subroutine restrict_bv

!===========================================================================

   subroutine interpolate(v, xn, yn ,zn, vi, xni, yni, zni, lam1, lam2, lam3, lbz, epsout )


   integer xn,yn,zn,xni,yni,zni
   _REAL_ v(1:xn*yn*zn),vi(1:xni*yni*zni), epsout
   _REAL_ lam1(1-xni*yni:xni*yni*zni),lam2(1-xni*yni:xni*yni*zni),lam3(1-xni*yni:xni*yni*zni)
   _REAL_ lbz(1:xni*yni*zni)

   integer i,j,k,xniyni,xniynizni,ii,ii2

   if ( xn*2+1 /= xni .or. yn*2+1 /= yni .or. zn*2+1 /= zni ) then
      write (6,*) "Interpolation Failed because of incorrect dimension"
      stop
   end if

   xniyni = xni*yni
   xniynizni = xni*yni*zni

   lam1(1-xniyni:0) = epsout
   lam2(1-xniyni:0) = epsout
   lam3(1-xniyni:0) = epsout
   do k = 1, zni; do j = 1, yni
      lam1(j*xni+(k-1)*xniyni) = epsout
   end do; end do
   do k = 1, zni; do i = 1, xni
      lam2(i-xni+k*xniyni) = epsout
   end do; end do
   do j = 1, yni; do i = 1, xni
      lam3(i+(j-1)*xni+xniynizni-xniyni) = epsout
   end do; end do

   ii = 0
   do k = 1, zn; do j=1 ,yn; do i=1, xn
      ii = ii + 1
      ii2 = (k*2-1)*xniyni + (j*2-1)*xni + i*2
      vi(ii2) = vi(ii2) + v(ii)
      call ipl_chain(vi,xniyni,xniynizni,ii2,v(ii),lbz,lam1,-1,lam2,xni,lam3,xniyni)
      call ipl_chain(vi,xniyni,xniynizni,ii2,v(ii),lbz,lam1,+1,lam2,xni,lam3,xniyni)
      call ipl_chain(vi,xniyni,xniynizni,ii2,v(ii),lbz,lam2,-xni,lam1,1,lam3,xniyni)
      call ipl_chain(vi,xniyni,xniynizni,ii2,v(ii),lbz,lam2,+xni,lam1,1,lam3,xniyni)
      call ipl_chain(vi,xniyni,xniynizni,ii2,v(ii),lbz,lam3,-xniyni,lam2,xni,lam1,1)
      call ipl_chain(vi,xniyni,xniynizni,ii2,v(ii),lbz,lam3,+xniyni,lam2,xni,lam1,1)
   end do; end do; end do

   lam1(1-xniyni:0) = ZERO
   lam2(1-xniyni:0) = ZERO
   lam3(1-xniyni:0) = ZERO
   do k = 1, zni; do j = 1, yni
      lam1(j*xni+(k-1)*xniyni) = ZERO
   end do; end do
   do k = 1, zni; do i = 1, xni
      lam2(i-xni+k*xniyni) = ZERO
   end do; end do
   do j = 1, yni; do i = 1, xni
      lam3(i+(j-1)*xni+xniynizni-xniyni) = ZERO
   end do; end do

contains

subroutine ipl_chain(vi,xnyn,xnynzn,l,v,lbz,am_1,shift_1,am_2,shift_2,am_3,shift_3)

   integer xnyn,xnynzn,l,shift_1,shift_2,shift_3
   _REAL_ vi(1:xnynzn),v,am_1(1-xnyn:xnynzn),am_2(1-xnyn:xnynzn),am_3(1-xnyn:xnynzn)
   _REAL_ lbz(1:xnynzn)

   integer l1
   _REAL_ v1

   v1 = ipl_comp1(v,l,lbz,am_1,xnyn,xnynzn,shift_1)
   l1 = l + shift_1 
   vi(l1) = vi(l1) + v1

   call ipl_chain2(vi,xnyn,xnynzn,l1,v1,lbz,am_1,shift_1,am_2,-shift_2,am_3,shift_3)
   call ipl_chain2(vi,xnyn,xnynzn,l1,v1,lbz,am_1,shift_1,am_2,shift_2,am_3,shift_3)
   call ipl_chain2(vi,xnyn,xnynzn,l1,v1,lbz,am_1,shift_1,am_3,-shift_3,am_2,shift_2)
   call ipl_chain2(vi,xnyn,xnynzn,l1,v1,lbz,am_1,shift_1,am_3,shift_3,am_2,shift_2)

end subroutine ipl_chain

subroutine ipl_chain2(vi,xnyn,xnynzn,l1,v1,lbz,am_1,shift_1,am_2,shift_2,am_3,shift_3)

   integer xnyn,xnynzn,l1,shift_1,shift_2,shift_3
   _REAL_ vi(1:xnynzn),v1,am_1(1-xnyn:xnynzn),am_2(1-xnyn:xnynzn),am_3(1-xnyn:xnynzn)
   _REAL_ lbz(1:xnynzn)

   integer l2
   _REAL_ v2

   v2 = ipl_comp2(v1,l1,lbz,am_1,am_2,xnyn,xnynzn,shift_1,shift_2)
   l2 = l1 + shift_2
   vi(l2) = vi(l2) + v2
   vi(l2-shift_3)= vi(l2-shift_3) + ipl_comp3(v2,l2,lbz,am_1,am_2,am_3,xnyn,xnynzn,shift_1,shift_2,-shift_3)
   vi(l2+shift_3)= vi(l2+shift_3) + ipl_comp3(v2,l2,lbz,am_1,am_2,am_3,xnyn,xnynzn,shift_1,shift_2,+shift_3)

end subroutine ipl_chain2

function ipl_comp1(v,l,lbz,am_1,xnyn,xnynzn,shift_1) 

   integer l,xnyn,xnynzn,shift_1
   _REAL_ ipl_comp1,v,am_1(1-xnyn:xnynzn),lbz(1:xnynzn)
   if ( shift_1 < 0 ) then
      ipl_comp1 = v * am_1(l+shift_1) / ( lbz(l+shift_1) + am_1(l+2*shift_1) + am_1(l+shift_1) )
   else
      ipl_comp1 = v * am_1(l) / ( lbz(l+shift_1) + am_1(l) + am_1(l+shift_1) )
   end if
   return

end function ipl_comp1

function ipl_comp2(v,l,lbz,am_1,am_2,xnyn,xnynzn,shift_1,shift_2) 

   integer l,xnyn,xnynzn,shift_1,shift_2
   _REAL_ ipl_comp2,v,am_1(1-xnyn:xnynzn),am_2(1-xnyn:xnynzn)
   _REAL_ lbz(1:xnynzn)

   _REAL_ lad

   lad = am_1(l+shift_2) + am_1(l+shift_2-abs(shift_1)) + lbz(l+shift_2)
   if ( shift_2 < 0 ) then
      ipl_comp2 = v * am_2(l+shift_2) / ( am_2(l+2*shift_2) + am_2(l+shift_2) + lad)
   else
      ipl_comp2 = v * am_2(l) / ( am_2(l) + am_2(l+shift_2) + lad )
   end if
   return

end function ipl_comp2
            
function ipl_comp3(v,l,lbz,am_1,am_2,am_3,xnyn,xnynzn,shift_1,shift_2,shift_3) 

   integer l,xnyn,xnynzn,shift_1,shift_2,shift_3
   _REAL_ ipl_comp3,v,am_1(1-xnyn:xnynzn),am_2(1-xnyn:xnynzn),am_3(1-xnyn:xnynzn)
   _REAL_ lbz(1:xnynzn)

   _REAL_ lad

   lad = am_1(l+shift_3) + am_1(l+shift_3-abs(shift_1)) + am_2(l+shift_3) + &
         am_2(l+shift_3-abs(shift_2)) + lbz(l+shift_3)
   if ( shift_3 < 0 ) then
      ipl_comp3 = v * am_3(l+shift_3) / ( am_3(l+2*shift_3) + am_3(l+shift_3) + lad)
   else
      ipl_comp3 = v * am_3(l) / ( am_3(l) + am_3(l+shift_3) + lad )
   end if
   return

end function ipl_comp3

end subroutine interpolate

!===========================================================================

subroutine pb_setupper( l_am1, l_am2, l_am3 )

   implicit none

   _REAL_ l_am1(l_xm,l_ym,l_zm), l_am2(l_xm,l_ym,l_zm), l_am3(l_xm,l_ym,l_zm)

   integer i,j,k

   do j = 1, l_ym; do k = 1, l_zm
      l_am1(l_xm,j,k) = ZERO
   end do; end do
   do i = 1, l_xm; do k = 1, l_zm
      l_am2(i,l_ym,k) = ZERO
   end do; end do
   do i = 1, l_xm; do j = 1, l_ym
      l_am3(i,j,l_zm) = ZERO
   end do; end do

end subroutine pb_setupper

subroutine pb_setupper2( l_am1, l_am2, l_am3, l_am4, l_am5, l_am6 )

   _REAL_ l_am1(l_xm,l_ym,l_zm), l_am2(l_xm,l_ym,l_zm), l_am3(l_xm,l_ym,l_zm)
   _REAL_ l_am4(l_xm,l_ym,l_zm), l_am5(l_xm,l_ym,l_zm), l_am6(l_xm,l_ym,l_zm)
   integer i, j, k

!   : pb_setupper2: setting x faces';flush(6)
   do j = 1, l_ym; do k = 1, l_zm
      if ( l_bcopt == 10 .or. l_bcopt == 1 ) then 
!!        ! : pb_setupper2: setting x=1 face';flush(6)
         l_am4(1,j,k) = l_am1(l_xm,j,k)
!!        ! : pb_setupper2: setting x=xm face';flush(6)
!         l_am4(l_xm,j,k) = l_am4(1,j,k)
      end if
      l_am1(l_xm,j,k) = ZERO
   end do; end do
!   : pb_setupper2: setting y faces';flush(6)
   do i = 1, l_xm; do k = 1, l_zm
      if ( l_bcopt == 10  .or. l_bcopt == 1) then 
         l_am5(i,1,k) = l_am2(i,l_ym,k)
!         l_am5(i,l_ym,k) = l_am5(i,1,k)
        ! l_am5(i,2:l_ym,k)=ZERO
      end if
      l_am2(i,l_ym,k) = ZERO
   end do; end do
!   : pb_setupper2: setting z faces';flush(6)
   do i = 1, l_xm; do j = 1, l_ym
      if ( l_bcopt == 10  .or. l_bcopt == 1) then
         l_am6(i,j,1) = l_am3(i,j,l_zm)
!         l_am6(i,j,l_zm) = l_am6(i,j,1)
      end if
      l_am3(i,j,l_zm) = ZERO
   end do; end do
!   : pb_setupper2: done';flush(6)
end subroutine pb_setupper2

subroutine pb_setupper3( l_am1, l_am2, l_am3, l_am4, l_am5, l_am6 )

   _REAL_ l_am1(l_xm,l_ym,l_zm), l_am2(l_xm,l_ym,l_zm), l_am3(l_xm,l_ym,l_zm)
   _REAL_ l_am4(l_xm,l_ym,l_zm), l_am5(l_xm,l_ym,l_zm), l_am6(l_xm,l_ym,l_zm)
   integer i, j, k

!!  ! : pb_setupper2: setting x faces';flush(6)
   do j = 1, l_ym; do k = 1, l_zm
      if ( l_bcopt == 10 .or. l_bcopt == 1 ) then 
!!        ! : pb_setupper2: setting x=1 face';flush(6)
         l_am4(1,j,k) = l_am1(l_xm,j,k)
!!        ! : pb_setupper2: setting x=xm face';flush(6)
         l_am4(l_xm,j,k) = l_am1(1,j,k)
      end if
      l_am1(l_xm,j,k) = ZERO
      l_am1(1,j,k) = ZERO
   end do; end do
!!  ! : pb_setupper2: setting y faces';flush(6)
   do i = 1, l_xm; do k = 1, l_zm
      if ( l_bcopt == 10  .or. l_bcopt == 1) then 
         l_am5(i,1,k) = l_am2(i,l_ym,k)
         l_am5(i,l_ym,k) = l_am2(i,1,k)
      end if
      l_am2(i,l_ym,k) = ZERO
      l_am2(i,1,k) = ZERO
   end do; end do
!!  ! : pb_setupper2: setting z faces';flush(6)
   do i = 1, l_xm; do j = 1, l_ym
      if ( l_bcopt == 10  .or. l_bcopt == 1) then
         l_am6(i,j,1) = l_am3(i,j,l_zm)
         l_am6(i,j,l_zm) = l_am3(i,j,1)
      end if
      l_am3(i,j,l_zm) = ZERO
      l_am3(i,j,1) = ZERO
   end do; end do
end subroutine pb_setupper3

end module pb_lsolver
