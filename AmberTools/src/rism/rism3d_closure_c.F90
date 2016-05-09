!<compile=optimized>
!The 3D-RISM-KH software found here is copyright (c) 2010-2012 by 
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Closure super class for 3D-RISM.  Closure sub-classes, (i.e., actual closure 
!!!implementations) are registered here.  Subroutine calls then call the 
!!!appropriate subroutine of the subclass. interface.  This is an explicit 
!!!implementation of class inheritance. See V. K. Decyk, C. D. Norton, 
!!!B. K. Szymanski.  How to express C++ concepts in Fortran 90. Scientific 
!!!Programming. 6, 363-390 (1997).
!!!
!!!Some closure independent properties are calculated within this class. Uvv is
!!!the site-site potential.
!!!
!!!In gerneral, this class is MPI aware only through the rism3d_grid
!!!class.  I.e., it knows about the total system size and its own
!!!piece of it.  It does not know about processes and does not perform
!!!reductions.  I.E., all thermodynamic quantities are distributed and
!!!each process has only the contribution of the local slab.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "../include/dprec.fh"

  module rism3d_closure_c
    !add new closure modules here
    use rism3d_potential_c
    use rism3d_grid_c
    use rism3d_kh_c
    use rism3d_hnc_c
    use rism3d_psen_c
    use rism_report_c
    use safemem
    implicit none

    type rism3d_closure
       character(len=4) :: type
       type(rism3d_kh), pointer :: kh=>NULL()
       type(rism3d_hnc), pointer :: hnc=>NULL()
       type(rism3d_psen), pointer :: psen=>NULL()
       type(rism3d_potential),pointer :: pot => NULL()
       !grid : points to grid in potential object
       type(rism3d_grid),pointer :: grid => NULL()
       !solv : points to solvent in potential object
       type(rism3d_solv),pointer :: solv => NULL()
       !solu : points to solute in potential object
       type(rism3d_solu),pointer :: solu => NULL()
    end type rism3d_closure

  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Creates a new closure object of the requested type.
!!!IN:
!!!   this  : the closure object
!!!   pot   : rism3d_potential object.  Must be initialized.
!!!   type : one of 'KH', 'HNC', 'PSEn' , where 'n' is the order of
!!!          the PSE-n closure
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_closure_new(this,type,pot)
      use rism_util, only : caseup
      implicit none
      type(rism3d_closure), intent(inout) :: this
      type(rism3d_potential), target, intent(in)  :: pot
      character(len=*), intent(in) :: type
      integer :: order, iostat
      this%pot => pot
      this%grid => this%pot%grid
      this%solv => this%pot%solv
      this%solu => this%pot%solu
      this%type=trim(type)
      call caseup(this%type)
      if(this%type .eq. "KH") then
         allocate(this%kh)
         call rism3d_kh_new(this%kh,this%pot)
      elseif(index(this%type,"PSE") ==1) then
         read(this%type(4:),*, iostat=iostat) order
         if(iostat/=0)&
              call rism_report_error("'"//trim(this%type)//"' not a valid closure")
         allocate(this%psen)
         call rism3d_psen_new(this%psen,this%pot,order)
      elseif(trim(this%type) .eq. "HNC") then
         allocate(this%hnc)
         call rism3d_hnc_new(this%hnc,this%pot)
      else
         call rism_report_error("'"//trim(this%type)//"' not a valid closure")
      end if
    end subroutine rism3d_closure_new

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!check if we can perform a temperature derivative calculation
!!!(i.e. does the closure support temperature derivatives)
!!!IN:
!!!   this : rism3d_closure object
!!!OUT:
!!!    .true. if we can, .false. if we can't
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism3d_closure_canCalc_DT(this) result(can_dT)
    implicit none
    type(rism3d_closure), intent(in) :: this
    logical :: can_dT
    can_dT=.false.
    if(associated(this%kh))then
       can_dT=.true.
!      elseif(associated(this%psen))then
!      elseif(associated(this%hnc))then
    end if
  end function rism3d_closure_canCalc_DT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns a identifier string for the closure type
!!!IN:
!!!   this : the closure object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_closure_type(this) result(type)
      implicit none
      type(rism3d_closure), intent(in) :: this
      character(len=4) :: type
      type=this%type
    end function rism3d_closure_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates Guv from Uuv, Huv, and Cuv using the associated closure
!!!IN:
!!!   this : the closure object
!!!   guv  : site-site pair correlation function
!!!   huv  : site-site total correlation function
!!!   cuv  : site-site direct correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_closure_guv(this,guv, huv, cuv)
      implicit none
      type(rism3d_closure), intent(inout) :: this
      _REAL_, intent(out) :: guv(:,:)
      _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
      if(associated(this%kh))then
         call rism3d_kh_guv(this%kh,guv,huv,cuv)
      elseif(associated(this%psen))then
         call rism3d_psen_guv(this%psen,guv,huv,cuv)
      elseif(associated(this%hnc))then
         call rism3d_hnc_guv(this%hnc,guv,huv,cuv)
      end if
    end subroutine rism3d_closure_guv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates Guv from Uuv, Huv, and Cuv using the associated closure
!!!IN:
!!!   this : the closure object
!!!   guv  : site-site pair correlation function
!!!   huv  : site-site total correlation function
!!!   cuv  : site-site direct correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_closure_guv_dt(this,guv_dT, huv_dT, cuv_dT, guv)
      implicit none
      type(rism3d_closure), intent(inout) :: this
      _REAL_, intent(out) :: guv_dT(:,:)
      _REAL_, intent(in) :: huv_dT(:,:),cuv_dT(:,:,:,:), guv(:,:)
      if(associated(this%kh))then
         call rism3d_kh_guv_dt(this%kh,guv_dT,huv_dT,cuv_dT,guv)
!      elseif(associated(this%psen))then
!         call rism3d_psen_guv_dT(this%psen,guv_dT,huv_dT,cuv_dT,guv)
!      elseif(associated(this%hnc))then
!         call rism3d_hnc_guv_dT(this%hnc,guv_dT,huv_dT,cuv_dT,guv)
      end if
    end subroutine rism3d_closure_guv_dt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the excess chemical potential in kT for each site
!!!IN:
!!!   this : the closure object
!!!   huv  : site-site total correlation function
!!!   cuv  : site-site direct correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_closure_exchem(this, huv, cuv) result(exchem)
      implicit none
      type(rism3d_closure), intent(inout) :: this
      _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
      _REAL_ :: exchem(this%solv%natom)
      if(associated(this%kh))then
         exchem = rism3d_kh_exchem(this%kh,huv,cuv)
      elseif(associated(this%psen))then
         exchem = rism3d_psen_exchem(this%psen,huv,cuv)
      elseif(associated(this%hnc))then
         exchem = rism3d_hnc_exchem(this%hnc,huv,cuv)
      end if
    end function rism3d_closure_exchem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the excess chemical potential w/ asymptotic correction in kT for 
!!!each site
!!!IN:
!!!   this : the closure object
!!!   huv  : site-site total correlation function
!!!   cuv  : site-site direct correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_closure_aexchem(this, huv, cuv) result(exchem)
      use rism_util, only : gaussquad_legendre, checksum
      implicit none
      type(rism3d_closure), intent(inout) :: this
      _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
      _REAL_ :: exchem(this%solv%natom)
      if(associated(this%kh))then
         exchem = rism3d_kh_aexchem(this%kh,huv,cuv)
      elseif(associated(this%psen))then
         exchem = rism3d_psen_aexchem(this%psen,huv,cuv)
      elseif(associated(this%hnc))then
         exchem = rism3d_hnc_aexchem(this%hnc,huv,cuv)
      end if
    end function rism3d_closure_aexchem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the excess chemical potential in kT for each site using the
!!!Guassian fluctuation expression.  This is closure independent.
!!!IN:
!!!   this : the closure object
!!!   huv  : site-site total correlation function
!!!   cuv  : site-site direct correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_closure_exchemGF(this, huv, cuv) result(exchem)
      implicit none
      type(rism3d_closure), intent(inout) :: this
      _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
      _REAL_ :: exchem(this%solv%natom)
      _REAL_ :: tuv
      integer :: ix, iy, iz, iv, ig, igk
      exchem = 0.d0
      do iv=1,this%solv%natom
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
                  exchem(iv) = exchem(iv) - (1.d0+0.5d0*huv(igk,iv))*cuv(ix,iy,iz,iv)
               end do
            end do
         end do
         exchem(iv) =  exchem(iv)*this%solv%rho(iv)*this%grid%voxel
      enddo
    end function rism3d_closure_exchemGF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the excess chemical potential in kT for each site using the
!!!Guassian fluctuation expression.  This is closure independent.
!!!IN:
!!!   this : the closure object
!!!   huv  : site-site total correlation function
!!!   cuv  : site-site direct correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_closure_aexchemGF(this, huv, cuv) result(exchem)
      implicit none
      type(rism3d_closure), intent(inout) :: this
      _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
      _REAL_ :: exchem(this%solv%natom)
      _REAL_ :: exchemlr(this%solv%natom),exchemh2lr(this%solv%natom),exchemhclr(this%solv%natom)
      _REAL_ :: tuv, cuvlr, huvlr
      integer :: ix, iy, iz, iv, ig, igk

      if(.not.this%solu%charged)then
         exchem = rism3d_closure_exchemGF(this,huv,cuv)
         return
      end if

      !
      !Closure independent, long-range part
      !
      call rism3d_potential_int_h2_hc(this%pot,exchemh2lr,exchemhclr)
      exchemlr = exchemh2lr+exchemhclr
      !
      !Short-range part
      !
      huvlr=0d0
      exchem = 0.d0
      do iv=1,this%solv%natom
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
                  if(this%solv%ionic)then
                     huvlr = this%solv%charge_sp(iv)*this%pot%asymhr(ig)
                  end if
                  cuvlr = this%solv%charge(iv)*this%pot%asymcr(ig)
                  exchem(iv) = exchem(iv) - (1.d0+0.5d0*huv(igk,iv))*cuv(ix,iy,iz,iv)&
                       + (1.d0+0.5d0*huvlr)*cuvlr
               end do
            end do
         end do
         exchem(iv) =  exchem(iv)*this%grid%voxel
      enddo
      exchem = this%solv%rho*(exchem+exchemlr)
    end function rism3d_closure_aexchemGF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the solvation energy, dE = dmu + dTS, from the temperature derivative w/o
!!!asymptotic correction in kT for each site
!!!IN:
!!!   this : the closure object
!!!   huv_dT  : temperature derivative site-site total correlation function
!!!   cuv_dT  : temperature derivative site-site direct correlation function
!!!   huv     : site-site total correlation function
!!!   cuv     : site-site direct correlation function
!!!OUT:
!!!    excess chemical potential with long-range correction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_closure_solvEne(this, huv_dT, cuv_dT, huv, cuv) result(solvEne)
      use rism_util, only : gaussquad_legendre, checksum
      implicit none
      type(rism3d_closure), intent(inout) :: this
      _REAL_, intent(in) :: huv_dT(:,:),cuv_dT(:,:,:,:)
      _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
      _REAL_ :: solvEne(this%solv%natom)
      
      solvEne = huge(1d0)
      if(rism3d_solv_canCalc_dT(this%solv))then
         return
      elseif(associated(this%kh))then
         solvEne = rism3d_kh_solvEne(this%kh,huv_dT,cuv_dT,huv,cuv)
!      elseif(associated(this%psen))then
!         solvEne = rism3d_psen_aexchem(this%psen,huv_dT,cuv_dT,huv,cuv)
!      elseif(associated(this%hnc))then
!         solvEne = rism3d_hnc_aexchem(this%hnc,huv_dT,cuv_dT,huv,cuv)
      end if
    end function rism3d_closure_solvEne

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the solvation energy, dE = dmu + dTS, from the temperature derivative w/
!!!asymptotic correction in kT for each site
!!!IN:
!!!   this : the closure object
!!!   huv_dT  : temperature derivative site-site total correlation function
!!!   cuv_dT  : temperature derivative site-site direct correlation function
!!!   huv     : site-site total correlation function
!!!   cuv     : site-site direct correlation function
!!!OUT:
!!!    excess chemical potential with long-range correction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_closure_asolvEne(this, huv_dT, cuv_dT, huv, cuv) result(solvEne)
      use rism_util, only : gaussquad_legendre, checksum
      implicit none
      type(rism3d_closure), intent(inout) :: this
      _REAL_, intent(in) :: huv_dT(:,:),cuv_dT(:,:,:,:)
      _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
      _REAL_ :: solvEne(this%solv%natom)

      solvEne = huge(1d0)
      if(.not.rism3d_solv_canCalc_dT(this%solv))then
         return
      elseif(associated(this%kh))then
         solvEne = rism3d_kh_asolvEne(this%kh,huv_dT,cuv_dT,huv,cuv)
!      elseif(associated(this%psen))then
!         solvEne = rism3d_psen_aexchem(this%psen,huv_dT,cuv_dT,huv,cuv)
!      elseif(associated(this%hnc))then
!         solvEne = rism3d_hnc_aexchem(this%hnc,huv_dT,cuv_dT,huv,cuv)
      end if
    end function rism3d_closure_asolvEne

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculate the total solvation interaction energy: de = rho sum g*u for
!!!each solvent site.  I.e., the direct intection potential energy of
!!!solute and solvent and not the total solvation energy (see solvEne).
!!! IN:
!!!    this :: rism3d object with computed solution
!!!    guv  :: site-site pair distribution function
!!!OUT:
!!!    the contribution of each solvent site to the total
!!!    solute-solvent potential energy [kT]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_closure_solvPotEne (this,guv) result(ene)
      implicit none
      type(rism3d_closure) :: this
      _REAL_, intent(in) :: guv(:,:)
      integer ::  igk, iv, ix,iy,iz,iatu
      _REAL_ ::  ene(this%solv%natom), ulj, uc
      ene = 0.d0
#ifdef RISM_DEBUG
      write(6,*)"EXENER"
#endif /*RISM_DEBUG*/

      !!FIX - can we avoid the loops and use BLAS array operations?
      do iv=1,this%solv%natom
         do iz=1,this%grid%nr(3)
            do iy=1,this%grid%nr(2)
               do ix=1,this%grid%nr(1)
#if defined(MPI)
                  igk = ix + (iy-1)*(this%grid%nr(1)+2) + (iz-1)*this%grid%nr(2)*(this%grid%nr(1)+2)
#else
                  igk = ix + (iy-1)*this%grid%nr(1) + (iz-1)*this%grid%nr(2)*this%grid%nr(1)
#endif /*defined(MPI)*/
                  ene(iv) = ene(iv) + guv(igk,iv)*this%pot%uuv(ix,iy,iz,iv)
               end do
            end do
         end do
      enddo
      !!endfix
      ene = ene*this%solv%rho*0.5d0*this%grid%voxel
    end function rism3d_closure_solvPotEne

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!calculates the partial molar volume.  Distributed solution.
!!!function
!!!IN:
!!!   this :: rism3d object with calculated solution    
!!!   cuv  :: (nr(1),nr(2),nr(3),solv%natom) site-site direct correlation function
!!!OUT:
!!!   the calculated partial molar volume [A^3]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_closure_pmv(this,cuv) result(pmv)
      implicit none
      type(rism3d_closure) :: this
      _REAL_, intent(in) :: cuv(:,:,:,:)
      _REAL_ :: pmv
      integer :: iv, ierr
      !........................... calculating PMV ...........................
      pmv=0
      do iv=1,this%solv%natom
         pmv = pmv + sum(cuv(:,:,:,iv))*this%solv%rho(iv)
      enddo
      pmv = -pmv*this%grid%voxel*this%solv%xikt
      if(this%grid%mpirank==0)&
           pmv= pmv + this%solv%xikt

    end function rism3d_closure_pmv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculate the excess number of particles about the solute compared to the 
!!!bulk solvation.  No attempt is made to account for excluded volume.  This is
!!!related to the Kirkwood G-factor.
!!! IN:
!!!    this :: rism3d object with computed solution
!!!    guv  :: site-site pair distribution function
!!!OUT:
!!!    The number of excess particles for each solvent site.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_closure_exNum (this,guv) result(num)
      implicit none
      type(rism3d_closure) :: this
      _REAL_, intent(in) :: guv(:,:)
      _REAL_ ::  num(this%solv%natom)
      num = this%solv%rho*rism3d_closure_kirkwoodG(this,guv)
    end function rism3d_closure_exNum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculate the excess number of particles about the solute compared to the 
!!!bulk solvation corrected with the longrange asymptotics.  No attempt is made 
!!!to account for excluded volume.  This is related to the Kirkwood G-factor.
!!! IN:
!!!    this :: rism3d object with computed solution
!!!    guv  :: site-site pair distribution function
!!!OUT:
!!!    The number of excess particles for each solvent site.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_closure_aexNum (this,guv) result(num)
      implicit none
      type(rism3d_closure) :: this
      _REAL_, intent(in) :: guv(:,:)
      _REAL_ ::  num(this%solv%natom)
      num = this%solv%rho*rism3d_closure_akirkwoodG(this,guv)
    end function rism3d_closure_aexNum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculate the Kirkwood G-factor for the solute.
!!!
!!!This expression is from Perkyns et al. J. Chem. Phys. 132, 064106 (2010) which 
!!!cites Hansen and McDonald. Theory of Simple Liquids. 2nd Ed.
!!! IN:
!!!    this :: rism3d object with computed solution
!!!    guv  :: site-site pair distribution function
!!!OUT:
!!!    Kirkwood G-factor for each solvent site
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_closure_kirkwoodG (this,guv) result(G)
      implicit none
      type(rism3d_closure) :: this
      _REAL_, intent(in) :: guv(:,:)
      _REAL_ ::  G(this%solv%natom)
      integer :: iv, ix, iy, iz, igk
      G=0
      do iv=1,this%solv%natom
         do iz=1,this%grid%nr(3)
            do iy=1,this%grid%nr(2)
               do ix=1,this%grid%nr(1)
#if defined(MPI)
                  igk = ix + (iy-1)*(this%grid%nr(1)+2) + (iz-1)*this%grid%nr(2)*(this%grid%nr(1)+2)
#else
                  igk = ix + (iy-1)*this%grid%nr(1) + (iz-1)*this%grid%nr(2)*this%grid%nr(1)
#endif /*defined(MPI)*/
                  G(iv) = G(iv) + (guv(igk,iv)-1d0)
               end do
            end do
         end do
      enddo
      !!endfix
      G = G*this%grid%voxel
    end function rism3d_closure_kirkwoodG

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculate the Kirkwood G-factor for the solute with long range asymptotic 
!!!correction.  For the correction integral, the charge smear parameter is set 
!!!to zero to ensure charge neutrality
!!!
!!!This expression is from Perkyns et al. J. Chem. Phys. 132, 064106 (2010) which 
!!!cites Hansen and McDonald. Theory of Simple Liquids. 2nd Ed.
!!! IN:
!!!    this :: rism3d object with computed solution
!!!    guv  :: site-site pair distribution function
!!!OUT:
!!!    Kirkwood G-factor for each solvent site
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_closure_akirkwoodG (this,guv) result(G)
      implicit none
      type(rism3d_closure) :: this
      _REAL_, intent(in) :: guv(:,:)
      _REAL_ ::  G(this%solv%natom)
      integer :: iv, ix, iy, iz, ig, igk

      if(.not.this%solv%ionic)then
         G=rism3d_closure_kirkwoodG(this,guv)
         return
      end if
      G=0
      do iv=1,this%solv%natom
         do iz=1,this%grid%nr(3)
            do iy=1,this%grid%nr(2)
               do ix=1,this%grid%nr(1)
                  ig = ix + (iy-1)*this%grid%nr(1) + (iz-1)*this%grid%nr(2)*this%grid%nr(1)
#if defined(MPI)
                  igk = ix + (iy-1)*(this%grid%nr(1)+2) + (iz-1)*this%grid%nr(2)*(this%grid%nr(1)+2)
#else
                  igk = ix + (iy-1)*this%grid%nr(1) + (iz-1)*this%grid%nr(2)*this%grid%nr(1)
#endif /*defined(MPI)*/
                  G(iv) = G(iv) + (guv(igk,iv)-1d0)&
                       - this%pot%solv%charge_sp(iv)*this%pot%asymhr(ig)
               end do
            end do
         end do
      enddo
      G = G*this%grid%voxel
      !use a smear of 0 to ensure electroneutrality
      if(sum(this%grid%nrOff) == 0)&
           G=G+rism3d_potential_int_h(this%pot)
    end function rism3d_closure_akirkwoodG

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!calculates the forces on the solute contributed by the solvent according
!!!to 3D-RISM.  In fact, this subroutine calls the appropriate subroutines
!!!to calculate this
!!!IN:
!!!   this :: rism3d closure object with computed solution
!!!   ff   :: 3D-RISM forces [kT/A]
!!!   guv  :: site-site pair distribution function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism3d_closure_force (this,ff,guv)
    implicit none
    type(rism3d_closure):: this
    _REAL_,intent(out) :: ff(3,this%solu%natom)
    _REAL_, intent(in) :: guv(:,:)
    !................... Calculation force fields    .......................

#ifdef RISM_DEBUG
    write(6,*) "RISM_FF"
    call flush(6)
#endif /*RISM_DEBUG*/    
    ff=0
    if(this%pot%cut2 > sum(this%grid%boxlen**2))then
       call force_brute(this,&
            ff,guv)
    else
       call force_lj(this,&
            ff,guv)
       call force_coulomb(this,&
            ff,guv)
    end if
  end subroutine rism3d_closure_force

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Frees memory and resets object state
!!!IN:
!!!   this : the closure object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_closure_destroy(this)
      use safemem
      implicit none
      type(rism3d_closure), intent(inout) :: this
      if(associated(this%kh))then
         call rism3d_kh_destroy(this%kh)
         deallocate(this%kh)
      end if
      if(associated(this%psen))then
         call rism3d_psen_destroy(this%psen)
         deallocate(this%psen)
      end if
      if(associated(this%hnc))then
         call rism3d_hnc_destroy(this%hnc)
         deallocate(this%hnc)
      end if
      nullify(this%pot)
      nullify(this%grid)
      nullify(this%solu)
      nullify(this%solv)
    end subroutine rism3d_closure_destroy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                         PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Lennard-Jones contributions to the total mean solvation forces
!!!Forces remain distributed over the nodes for MPI calculations.
!!!IN:
!!!  this :: rism3d closure object with computed solution
!!!  ff   :: force array [kT/A]
!!!  guv  : site-site pair correlation function  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine force_lj (this, ff, guv)
      implicit none

      type(rism3d_closure), intent(in):: this
      _REAL_,intent(inout) :: ff(3,this%solu%natom)
      _REAL_, intent(in) :: guv(:,:)
      integer ::  igx,igy,igz,igA(3), ig,ig2, igk, iu, iv, id, idim1,idim2
      _REAL_ ::  rx,ry,rz, dX(3), dz2,dyz2,rs2,rs6i, r5,r4,r3,r2, ra
      _REAL_ ::  rcor, rcor2,offset
      parameter (rcor=0.002d0, rcor2=rcor**2)

      _REAL_ ::  dUlj_dr,dU_dr(3)
      !linear spacing of the grid
      _REAL_ :: ff_temp(3,this%solu%natom)
      !number of gridpoints in each direction to use
      !closest point to solute atom
      integer :: grdpnts(3),cp(3),first(3),last(3),iu2,id2,ierr

#ifdef RISM_DEBUG
      write(6,*) "FORCE_LJ"; call flush(6)
#endif /*RISM_DEBUG*/
      !......................... clearing summators ..........................
      ff = 0d0

      offset = this%grid%grdspc(3)*(this%grid%nrOff(3))
      !calculate the number of grid point necessary to cover this cutoff range.
      !In case of a large cutoff, we need to protect against integer overflow
      grdpnts = nint(min(dble(huge(1)-1),sqrt(this%pot%cut2)/this%grid%grdspc))+1
      do iu=1,this%solu%natom
         do id=1,3
            cp(id) = nint(this%solu%ratu(id,iu)/this%grid%grdspc(id))
            first(id) = max(0,cp(id)-grdpnts(id))
            !Note: we have to protect against cp(id)+grdpnts(id) overflowing
            last(id) = min(this%grid%ngr(id)-1,&
                 cp(id)+min(huge(1) -cp(id),grdpnts(id)))
         end do
         cp(3) = cp(3) - this%grid%nrOff(3)
         first(3) = max(0,first(3) - this%grid%nrOff(3))
         last(3) = min(this%grid%nr(3)-1,last(3) - this%grid%nrOff(3))
         do igz=first(3),last(3)
            rz = igz*this%grid%grdspc(3)+offset
            dX(3) = this%solu%ratu(3,iu) - rz
            dz2 = dX(3)*dX(3)
            do igy=first(2),last(2)
               ry = igy*this%grid%grdspc(2)
               dX(2) = this%solu%ratu(2,iu) - ry
               dyz2=dX(2)*dX(2)+dz2
               do igx=first(1),last(1)
                  rx = igx*this%grid%grdspc(1)
#if defined(MPI)
                  ig = 1 + igx + igy*(this%grid%ngr(1)+2) + igz*this%grid%ngr(2)*(this%grid%ngr(1)+2)
#else
                  ig = 1 + igx + igy*this%grid%ngr(1) + igz*this%grid%ngr(2)*this%grid%ngr(1)
#endif /*defined(MPI)*/

                  !......... site separation subject to minimal image condition ..........
                  dX(1) = this%solu%ratu(1,iu) - rx

                  r2 = dX(1)*dX(1) + dyz2
                  if (r2 < rcor2)  r2 = rcor2
                  if(r2<this%pot%cut2)then
                     !...................... maintaining potential Ulj ......................
                     do iv=1,this%solv%natom
                        rs2 = r2/this%pot%siguv(iu,iv)**2
                        rs6i = 1d0/rs2**3
                        dUlj_dr=12.d0*this%pot%epsuv(iu,iv)*rs6i*(rs6i-1.d0)&
                             /r2*this%solv%rho(iv)*guv(ig,iv)
                        do id = 1,3 
                           ff(id,iu) = &
                                ff(id,iu) + dUlj_dr*dX(id)
                        end do
                     enddo
                  endif
               enddo
            enddo
         enddo
      enddo
      ff = ff*this%grid%voxel
    end subroutine force_lj

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Electrostatic contributions to the total mean solvation forces
!!!Forces remain distributed over the nodes for MPI calculations.
!!!IN:
!!!  this :: rism3d closure object with computed solution
!!!  ff   :: force array [kT/A]
!!!  guv  : site-site pair correlation function  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine force_coulomb(this,ff,guv)
      implicit none
      type(rism3d_closure),intent(in) :: this
      _REAL_,intent(inout) ::  ff(3,this%solu%natom)
      _REAL_, intent(in) :: guv(:,:)
      _REAL_ :: frc(3)
      integer ::  igx,igy,igz,igA(3), ig, igk, iu, iv, id, idim1,idim2
      _REAL_ ::  rx,ry,rz, dX(3),dz2,dyz2, r2, ra
      _REAL_ ::  rcor, rcor2
      parameter (rcor=0.002d0, rcor2=rcor**2)

      _REAL_ ::  cff
      _REAL_ ::  dUc_dr,dU_dr(3)
      _REAL_ :: qutemp(this%solu%natom),qvtemp(this%solv%natom)
      !linear spacing of the grid
      _REAL_ :: offset,ff_temp(3,this%solu%natom)
      !smstart :: initial grid point where the sparse and fine grids have the same value
      integer :: smstart
      integer :: ierr,irank
      !number of gridpoints in each direction to use
      !closest point to solute atom
      integer :: grdpnts(3),cp(3),first(3),last(3),iu2,id2

#ifdef RISM_DEBUG
      write(6,*) "FORCE_COULOMB",this%pot%cut2,sqrt(this%pot%cut2); call flush(6)
#endif /*RISM_DEBUG*/
      smstart = mod(this%grid%nrOff(3),2)
      qutemp(1:this%solu%natom)=this%solu%charge(1:this%solu%natom)!*AMBER_ELECTROSTATIC
      qvtemp(1:this%solv%natom)=this%solv%charge(1:this%solv%natom)!*AMBER_ELECTROSTATIC
      offset = this%grid%grdspc(3)*this%grid%nrOff(3)
      !calculate the number of grid point necessary to cover this cutoff range.
      !In case of a large cutoff, we need to protect against integer overflow
      grdpnts = nint(min(dble(huge(1)-1),sqrt(this%pot%cut2)/this%grid%grdspc))+1!make it an odd number
      do iu=1,this%solu%natom
         frc = 0
         do id=1,3
            cp(id) = nint(this%solu%ratu(id,iu)/this%grid%grdspc(id))
            first(id) = max(0,cp(id)-grdpnts(id))
            first(id) = first(id) + mod(first(id),2) !start on an even number
            last(id) = min(this%grid%ngr(id)-1,first(id)+2*grdpnts(id)-1) !end on an odd number
         end do
         cp(3) = cp(3) - this%grid%nrOff(3)
         first(3) = max(0,first(3) - this%grid%nrOff(3))
         last(3) = min(this%grid%nr(3)-1,last(3) - this%grid%nrOff(3))
         do igz=first(3),last(3)
            rz = igz*this%grid%grdspc(3)+offset
            dX(3) =  this%solu%ratu(3,iu) - rz
            dz2=dX(3)*dX(3)
            do igy=first(2),last(2)
               ry = igy*this%grid%grdspc(2)
               dX(2) =  this%solu%ratu(2,iu) - ry
               dyz2=dX(2)*dX(2)+dz2
               do igx=first(1),last(1)

#if defined(MPI)
                  ig = 1 + igx + igy*(this%grid%ngr(1)+2) + igz*this%grid%ngr(2)*(this%grid%ngr(1)+2)
#else
                  ig = 1 + igx + igy*this%grid%ngr(1) + igz*this%grid%ngr(2)*this%grid%ngr(1)
#endif /*defined(MPI)*/
                  rx = igx*this%grid%grdspc(1)

                  !......... site separation subject to minimal image condition ..........
                  dX(1) =  this%solu%ratu(1,iu) - rx
                  r2 = dX(1)*dX(1) + dyz2
                  if(r2<rcor2)r2 = rcor2
                  ra=sqrt(r2)
                  !                  ra = max(sqrt(r2),rcor)
                  do iv=1,this%solv%natom
                     dUc_dr  = qutemp(iu)*qvtemp(iv)/r2/ra&
                          *this%solv%rho(iv)*guv(ig,iv)
                     do id = 1,3 
                        frc(id) = &
                             frc(id) + dUc_dr*dX(id)
                     end do
                  end do
               end do
            enddo
         enddo
         ff(:,iu) = ff(:,iu) + frc*this%grid%voxel
         frc = 0
         !long range coarse grid
         do igz=smstart,this%grid%nr(3)-1,2
            rz = igz*this%grid%grdspc(3)+offset
            dX(3) =  this%solu%ratu(3,iu) - rz
            dz2=dX(3)*dX(3)
            do igy=0,this%grid%nr(2)-1,2
               ry = igy*this%grid%grdspc(2)
               dX(2) =  this%solu%ratu(2,iu) - ry
               dyz2=dX(2)*dX(2)+dz2
               do igx=0,this%grid%nr(1)-1,2
                  if((igx >= first(1) .and. igx <= last(1)) .and.&
                       (igy >= first(2) .and. igy <= last(2)) .and.&
                       (igz >= first(3) .and. igz <= last(3))) cycle

#if defined(MPI)
                  ig = 1 + igx + igy*(this%grid%nr(1)+2) + igz*this%grid%nr(2)*(this%grid%nr(1)+2)
#else
                  ig = 1 + igx + igy*this%grid%nr(1) + igz*this%grid%nr(2)*this%grid%nr(1)
#endif /*defined(MPI)*/
                  rx = igx*this%grid%grdspc(1)

                  !......... site separation subject to minimal image condition ..........
                  dX(1) =  this%solu%ratu(1,iu) - rx
                  r2 = dX(1)*dX(1) + dyz2
                  ra = sqrt(r2)
                  do iv=1,this%solv%natom
                     dUc_dr  = -qutemp(iu)*qvtemp(iv)/r2/ra&
                          *this%solv%rho(iv)*guv(ig,iv)
                     do id = 1,3 
                        frc(id) = &
                             frc(id) - dUc_dr*dX(id)
                     end do
                  end do
               end do
            enddo
         enddo
         ff(:,iu) = ff(:,iu) + frc*this%grid%voxel*8d0
      enddo
    end subroutine force_coulomb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates both Lennard-Jones and electrostatic solvation forces over the 
!!!entire grid without cutoffs.
!!!IN:
!!!   this : closure object
!!!   ff   : array for forces [kT/A]
!!!   guv  : site-site pair correlation function  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine force_brute (this,ff,guv)
      implicit none
      type(rism3d_closure), intent(in) :: this
      _REAL_, intent(inout) ::  ff(3,this%solu%natom)
      _REAL_, intent(in) :: guv(:,:)

      integer ::  igx,igy,igz, ig, igoffset, iu, iv, id, idim1,idim2
      _REAL_ ::  rx,ry,rz, dX(3), rs2,rs6i, r5,r4,r3,r2, ra
      _REAL_ :: ffcoef
      _REAL_ ::  rcor, rcor2
      parameter (rcor=0.002d0, rcor2=rcor**2)

      _REAL_ ::  sum, dUlj_dr, dUc_dr,dU_dr(3)
      _REAL_ :: qutemp(this%solu%natom),qvtemp(this%solv%natom)
      _REAL_ :: offset
      !......................... clearing summators ..........................
#ifdef RISM_DEBUG
      write(6,*)"rismff"
#endif /*RISM_DEBUG*/
      offset = this%grid%grdspc(3)*(this%grid%nrOff(3))

      ff = 0d0
      qutemp(1:this%solu%natom)=this%solu%charge(1:this%solu%natom)!*AMBER_ELECTROSTATIC
      qvtemp(1:this%solv%natom)=this%solv%charge(1:this%solv%natom)!*AMBER_ELECTROSTATIC
      !..................... enumerating box grid points .....................
      do igz=0,this%grid%nr(3)-1
         rz = igz*this%grid%grdspc(3)+offset
         do igy=0,this%grid%nr(2)-1
            ry = igy*this%grid%grdspc(2)
            do igx=0,this%grid%nr(1)-1
#if defined(MPI)
               ig = 1 + igx + igy*(this%grid%ngr(1)+2) + igz*this%grid%ngr(2)*(this%grid%ngr(1)+2)
#else
               ig = 1 + igx + igy*this%grid%ngr(1) + igz*this%grid%ngr(2)*this%grid%ngr(1)
#endif /*defined(MPI)*/
               rx = igx*this%grid%grdspc(1)

               !...................... summing over solute sites ......................
               do iu=1,this%solu%natom

                  !......... getting site separation ..........
                  dX(1) =  this%solu%ratu(1,iu) - rx
                  dX(2) =  this%solu%ratu(2,iu) - ry
                  dX(3) =  this%solu%ratu(3,iu) - rz

                  !                r2 = max(rcor2,(dX(1)*dX(1) + dX(2)*dX(2) + dX(3)*dX(3)))
                  r2 = dX(1)*dX(1) + dX(2)*dX(2) + dX(3)*dX(3)                
                  ra = sqrt(r2)
                  r3 = ra*r2
                  r4 = r2*r2
                  r5 = ra**5
                  !...................... maintaining potential dU_dR ....................
                  if(ra <rcor) then
                     cycle
                  end if
                  do iv=1,this%solv%natom
                     rs2 = r2/this%pot%siguv(iu,iv)**2
                     rs6i = 1.d0/rs2**3
                     dUlj_dr = -12.d0*this%pot%epsuv(iu,iv)*rs6i*(rs6i-1.d0)&
                          /ra
                     dUc_dr  = -qutemp(iu)*qvtemp(iv)/r2
                     do id = 1,3
                        dU_dR(id) = (dUlj_dr + dUc_dr)*dX(id)/ra
                        ff(id,iu) = &
                             ff(id,iu) - dU_dR(id)*this%solv%rho(iv)*guv(ig,iv)
                     end do
                  enddo
               enddo
            enddo
         enddo
      end do
      ff = ff*this%grid%voxel
      return
    end subroutine force_brute
  end module rism3d_closure_c
