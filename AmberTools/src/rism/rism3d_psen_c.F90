! <compile=optimized>

!The 3D-RISM-KH software found here is copyright (c) 2011-2012 by 
!Andriy Kovalenko, Tyler Luchko and David A. Case.
!
!This program is free software: you can redistribute it and/or modify it
!under the terms of the GNU General Public License as published by the Free
!Software Foundation, either version 3 of the License, or (at your option)
!any later version.
!
!This program is distributed in the hope that it will be useful, but
!WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
!or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!for more details.
!
!You should have received a copy of the GNU General Public License in the
!../../LICENSE file.  If not, see <http://www.gnu.org/licenses/>.
!
!Users of the 3D-RISM capability found here are requested to acknowledge
!use of the software in reports and publications.  Such acknowledgement
!should include the following citations:
!
!1) A. Kovalenko and F. Hirata. J. Chem. Phys., 110:10095-10112  (1999); 
!ibid. 112:10391-10417 (2000).   
!
!2) A. Kovalenko,  in:  Molecular  Theory  of  Solvation,  edited  by  
!F. Hirata  (Kluwer Academic Publishers, Dordrecht, 2003), pp.169-275.  
!
!3) T. Luchko, S. Gusarov, D.R. Roe, C. Simmerling, D.A. Case, J. Tuszynski,
!and  A. Kovalenko, J. Chem. Theory Comput., 6:607-624 (2010). 

#include "../include/dprec.fh"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Partial series expansion of order-n (PSE-n) closure class for 3D-RISM.
!!!Stefan M. Kast and Thomas Kloss. J. Chem. Phys. 129, 236101 (2008)
!!!
!!!Interpolates between the Kovalenko-Hirata (KH) and hypernetted chain equation
!!!closures (HNC).  Order-1 gives KH and order-infinity gives HNC.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  module rism3d_psen_c
    use rism3d_potential_c
    use rism3d_grid_c
    !the PSEN type
    type rism3d_psen
       type(rism3d_potential),pointer :: pot => NULL()
       !grid : points to grid in potential object
       type(rism3d_grid),pointer :: grid => NULL()
       !order    : order of the series expansion.  >= 1
       !order1   : order+1
       integer :: order=0, order1
       !order1fac: (order+1)!
       _REAL_ :: order1fac
    end type rism3d_psen

    public rism3d_psen_new, rism3d_psen_destroy!, rism3d_psen_guv
  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Initializes the PSEN closure
!!!IN:
!!!   this : PSEN object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_psen_new(this,pot, order)
      implicit none
      type(rism3d_psen), intent(inout) :: this
      type(rism3d_potential), target, intent(in) :: pot
      integer, intent(in) :: order
      integer :: i, orderfac
      this%pot => pot
      this%grid => this%pot%grid
      if(order <1)then
         call rism_report_error('(a,i4)',"PSE-n closure: order must be > 0:",order)
      end if
      this%order = order
      this%order1 = order+1
      this%order1fac = 1
      orderfac=1

      !calculate the factorial of (order+1) for future reference.  Check that we 
      !are not losing precision or overflowing the precision
      do i=2,this%order1
         orderfac = this%order1fac
         this%order1fac = this%order1fac*i
         !keep 16 significant digits
         if( abs((this%order1fac/i)/orderfac -1d0) >1d-16 )then
            call rism_report_error('(a,i4)',&
                 "PSE-n: numerical precision exceeded.  Use a smaller order. Try ",i-1)
         end if
      end do
    end subroutine rism3d_psen_new

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates Guv from Uuv, Huv, and Cuv using the PSEN closure
!!!IN:
!!!   this : the PSEN closure object
!!!   guv  : site-site pair correlation function
!!!   huv  : site-site total correlation function
!!!   cuv  : site-site direct correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_psen_guv(this,guv, huv, cuv)
      implicit none
      type(rism3d_psen), intent(in) :: this
      _REAL_, intent(out) :: guv(:,:)
      _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
      integer :: i, iv, ir, ix, iy, iz, ig, ig1
      _REAL_ :: tuv, orderfac

      do iv = 1,this%pot%solv%natom
         do iz = 1, this%grid%nr(3)
            do iy = 1, this%grid%nr(2)
               do ix = 1, this%grid%nr(1)
                  ig1 = ix + (iy-1)*this%grid%nr(1) + (iz-1)*this%grid%nr(2)*this%grid%nr(1)
#if defined(MPI)
                  ig = ix + (iy-1)*(this%grid%nr(1)+2) + (iz-1)*(this%grid%nr(1)+2)*this%grid%nr(2)
#else
                  ig = ix + (iy-1)*this%grid%nr(1) + (iz-1)*this%grid%nr(1)*this%grid%nr(2)
#endif /*defined(MPI)*/
                  tuv = -this%pot%uuv(ix,iy,iz,iv) + huv(ig,iv) - cuv(ix,iy,iz,iv)
                  if(tuv >= 0d0)then
                     guv(ig,iv) = 1d0
                     orderfac = 1d0
                     do i=1,this%order
                        orderfac = orderfac*i
                        guv(ig,iv) = guv(ig,iv) + (tuv**i)/orderfac
                     end do
                  else
                     guv(ig,iv) = exp(tuv)
                  endif
               end do
            end do
         end do
      end do
    end subroutine rism3d_psen_guv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the excess chemical potential in kT for each site
!!!IN:
!!!   this : the closure object
!!!   huv  : site-site total correlation function
!!!   cuv  : site-site direct correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_psen_exchem(this, huv, cuv) result(exchem)
      implicit none
      type(rism3d_psen), intent(inout) :: this
      _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
      _REAL_ :: exchem(this%pot%solv%natom)
      _REAL_ :: tuv, tsuv
      integer :: ix, iy, iz, iv, ig, igk
      exchem = 0.d0
      do iv=1,this%pot%solv%natom
         do iz=1,this%grid%nr(3)
            do iy=1,this%grid%nr(2)
               do ix=1,this%grid%nr(1)
                  ig = ix + (iy-1)*this%grid%nr(1) + (iz-1)*this%grid%nr(2)*this%grid%nr(1)
#if defined(MPI)
                  igk = ix + (iy-1)*(this%grid%nr(1)+2) + (iz-1)*this%grid%nr(2)*(this%grid%nr(1)+2)
#else
                  igk = ix + (iy-1)*this%grid%nr(1) + (iz-1)*this%grid%nr(2)*this%grid%nr(1)
#endif /*defined(MPI)*/
                  tuv = huv(igk,iv) - cuv(ix,iy,iz,iv)
                  exchem(iv) = exchem(iv) + 0.5d0*huv(igk,iv)*tuv - cuv(ix,iy,iz,iv)
                  if (huv(igk,iv) > 0d0)  then
                     tsuv = tuv - this%pot%uuv(ix,iy,iz,iv)
                     exchem(iv) = exchem(iv) - tsuv**(this%order1)/this%order1fac
                  endif
               end do
            end do
         end do
         exchem(iv) =  exchem(iv)*this%pot%solv%rho(iv)*this%grid%voxel
      enddo

    end function rism3d_psen_exchem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the short range part of the excess chemical potential w/ 
!!!asymptotic correction in kT for each site
!!!IN:
!!!   this : the closure object
!!!   huv  : site-site total correlation function
!!!   cuv  : site-site direct correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_psen_aexchem(this, huv, cuv) result(exchem)
      implicit none
      type(rism3d_psen), intent(inout) :: this
      _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
      _REAL_ :: exchem(this%pot%solv%natom)
      !exchemh2lr :: contribution from asymhr**2
      !exchemhclr :: contribution from asymhr*asymcr
      !exchemhnlr :: contribution from asymhr**n      
      _REAL_ :: exchemh2lr(this%pot%solv%natom),exchemhclr(this%pot%solv%natom),&
           exchemhnlr(this%pot%solv%natom)
      _REAL_ :: tuv, tsuv, huvlr, cuvlr, tuvlr
      integer :: ix, iy, iz, iv, ig, igk

      if(.not.this%pot%solu%charged)then
         exchem = rism3d_psen_exchem(this,huv,cuv)
         return
      end if


      !
      !long range part
      !
      call rism3d_potential_int_h2_hc(this%pot,exchemh2lr,exchemhclr)
      exchemhnlr = rism3d_potential_int_hn(this%pot,this%order1)
      
      do iv=1,this%pot%solv%natom
         if(this%pot%solu%totcharge*this%pot%solv%charge_sp(iv) > 0.d0) &
              exchemhnlr(iv)=0
      end do

      !
      !short range part
      !
      exchem = 0.d0
      huvlr=0d0
      do iv=1,this%pot%solv%natom
         do iz=1,this%grid%nr(3)
            do iy=1,this%grid%nr(2)
               do ix=1,this%grid%nr(1)
                  ig = ix + (iy-1)*this%grid%nr(1) + (iz-1)*this%grid%nr(2)*this%grid%nr(1)
#if defined(MPI)
                  igk = ix + (iy-1)*(this%grid%nr(1)+2) + (iz-1)*this%grid%nr(2)*(this%grid%nr(1)+2)
#else
                  igk = ix + (iy-1)*this%grid%nr(1) + (iz-1)*this%grid%nr(2)*this%grid%nr(1)
#endif /*defined(MPI)*/
                  tuv = huv(igk,iv) - cuv(ix,iy,iz,iv)
                  if(this%pot%solv%ionic)&
                       huvlr = this%pot%solv%charge_sp(iv)*this%pot%asymhr(ig)
                  cuvlr = this%pot%solv%charge(iv)*this%pot%asymcr(ig)
                  tuvlr = huvlr - cuvlr
                  exchem(iv) = exchem(iv) + 0.5d0*huv(igk,iv)*tuv - cuv(ix,iy,iz,iv)
                  exchem(iv) = exchem(iv) - (0.5d0*huvlr*tuvlr - cuvlr)
                  if (huv(igk,iv) > 0d0)  then
                     tsuv = tuv - this%pot%uuv(ix,iy,iz,iv)
                     exchem(iv) = exchem(iv) - tsuv**(this%order1)/this%order1fac
                  endif

                  if (this%pot%solu%totcharge*this%pot%solv%charge_sp(iv) <= 0.d0) then
                     !recall that cuvlr = - uuvlr so 
                     !(huvlr - cuvlr - uuvlr)**n = huvlr**n
                     exchem(iv) = exchem(iv)  + huvlr**(this%order1)/this%order1fac
                  endif
               end do
            end do
         end do
         exchem(iv) =  this%pot%solv%rho(iv)*(exchem(iv)*this%grid%voxel &
              +  (exchemh2lr(iv) - exchemhclr(iv))/2d0 - exchemhnlr(iv)/this%order1fac)
      enddo

    end function rism3d_psen_aexchem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Frees memory and resets the PSEN closure
!!!IN:
!!!   this : PSEN object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_psen_destroy(this)
      implicit none
      type(rism3d_psen), intent(inout) :: this
      nullify(this%pot)
      nullify(this%grid)
    end subroutine rism3d_psen_destroy
  end module rism3d_psen_c
