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

#include "../include/dprec.fh"

module rism3d_fft_c
#ifdef MKLdisabled
  use mkl_dfti
#else /*FFTW*/
  use FFTW3
#endif /*MKLdisabled FFTW*/
  use rism3d_grid_c
#ifdef RISM3D_DEBUG
!  use rism3d_debug_c
#endif
  implicit none
#ifdef MKLdisabled
  integer, parameter :: FFT_ALIGNED=0, FFT_UNALIGNED=1
  integer, parameter :: FFT_ESTIMATE=0, FFT_MEASURE=1,&
       FFT_PATIENT=2, FFT_EXHAUSTIVE=3
#else /*FFTW*/
  integer, parameter :: FFT_ALIGNED=0, FFT_UNALIGNED=FFTW_UNALIGNED
  integer, parameter :: FFT_ESTIMATE=FFTW_ESTIMATE, FFT_MEASURE=FFTW_MEASURE,&
       FFT_PATIENT=FFTW_PATIENT, FFT_EXHAUSTIVE=FFTW_EXHAUSTIVE
#endif
  type rism3d_fft
     private
#ifdef MKLdisabled
     !planfwd :: MKL plan for forward 3D-FFT
     !planfwd :: MKL plan for backward 3D-FFT
     type(DFTI_DESCRIPTOR), pointer ::  planfwd=>NULL(),planbwd=>NULL()
#else /*FFTW*/
     !planfwd :: FFTW plan for forward 3D-FFT
     !planfwd :: FFTW plan for backward 3D-FFT
     type(C_PTR) ::  planfwd=C_NULL_PTR, planbwd=C_NULL_PTR
#endif /*MKLdisabled FFTW*/
     !aligned :: Needed for SIMD.
     !           FFTW_ALIGNED -   assume that the memory is 16-byte aligned.  
     !           FFTW_UNALIGNED - assume that memory is not aligned.
     integer :: aligned
     !localtrans :: transpose data locally.  Not used yet.
     logical :: localtrans
     !grid  :: grid object with array dimensions
     type(rism3d_grid),pointer :: grid=>NULL()
     !narray :: the number of arrays to transform
     integer :: narray
     !wrk :: work space to transform to and from Numerical Recipes memory layout
     !wrk_trans :: work space to do local transpose
     _REAL_, pointer :: wrk(:)=>NULL(), wrk_trans(:)=>NULL()
  end type rism3d_fft
  public rism3d_fft_global_init, rism3d_fft_global_finalize,&
       rism3d_fft_setgrid, rism3d_fft_new, rism3d_fft_fwd, rism3d_fft_bwd, rism3d_fft_destroy
#ifdef MKLdisabled
#else
  private r2c_pointer
#endif

  private rlft3i
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Call this from each process (MPI or serial) before calling anyother FFT 
!!!routines or after rism3d_fft_global_finalize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism3d_fft_global_init()
    implicit none
#ifdef MPI    
    call fftw_mpi_init()
#endif
  end subroutine rism3d_fft_global_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Call this from each process (MPI or serial) to release all FFT memory after 
!!!destroying all FFT objects
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism3d_fft_global_finalize()
    implicit none
#ifdef MPI    
    call fftw_mpi_cleanup()
#endif
  end subroutine rism3d_fft_global_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Sets the grid dimensions based on the requirements of the FFT.
!!!Unless you set the grid yourself, this generally should be called
!!!before creating a new FFT object.
!!!IN:
!!!   grid :: grid object to set
!!!   ngr  :: size of the global array
!!!   grdspc :: grid spacing
!!!   narray :: number of grids
!!!   aligned :: whether or not memory will be 16-byte aligned (for FFTW many)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism3d_fft_setgrid(grid,ngr,grdspc,narray,aligned)
    implicit none
    type(rism3d_grid), intent(inout) :: grid
    integer, intent(in) :: ngr(3)
    _REAL_,intent(in) :: grdspc(3)
    integer, intent(in) :: narray
    logical, intent(in) :: aligned
#ifdef MKLdisabled
    integer ::  nkTotal, nkyOff,nky, nrzOff, nrz
#else /*FFTW*/
    integer(C_INTPTR_T) ::  nkTotal, nkyOff,nky, nrzOff, nrz
#endif /*MKLdisabled FFTW*/
    !sites_per_transform : when calling FFTW to do a transform, this
    !                      defines how many arrays are transformed at
    !                      once and has to do with FFTW's 'many'.  This
    !                      only matters for MPI where the memory layout
    !                      has performance issues
    integer ::sites_per_transform,ngk(3)
    !determine array size either from FFTW (MPI) or define them ourselves (serial)
#ifdef MPI
    !The FFTW command wants the size of the global complex array, ngk,
    !where we use the complex data type
    ngk = ngr
    ngk(1) = ngk(1)/2+1
    if(aligned)then
       !aligned memory
       sites_per_transform = narray
    else
       !unaligned memory
       sites_per_transform = 1
    end if
#ifdef MKLdisabled
#else /*FFTW*/
    nkTotal = int( fftw_mpi_local_size_many_transposed(3,&
         int(ngk(3:1:-1),C_INTPTR_T),&
         int(sites_per_transform,C_INTPTR_T), &
         FFTW_MPI_DEFAULT_BLOCK,FFTW_MPI_DEFAULT_BLOCK,&
         grid%mpicomm,&
         nrz,nrzOff,&
         nky,nkyOff),&
         kind(nktotal))
    if(mod(nkTotal,narray) /=0)&
         call rism_report_error("Illegal memory distribution from FFTW")
#endif /*MKLdisabled FFTW*/
    call rism3d_grid_resize(grid,grdspc,&
         ngr,&                                   !r-space global
         (/(ngr(1)+2),ngr(3),ngr(2)/),&          !k-space global
         (/ngr(1),ngr(2),int(nrz,kind(ngr(1)))/),&       !r-space local
         (/(ngr(1)+2),ngr(3),int(nky,kind(ngr(1)))/),&   !k-space local
         (/0,0,int(nrzOff,kind(ngr(1)))/), &                !r-space offset
         (/0,0,int(nkyOff,kind(ngr(1)))/))                  !k-space offset

    !As returned, nkTotal is the number of local complex elements for
    !all solvent sites combined.  However, to be consistent with the
    !serial code, nkTotal will be the local number of _REAL_ data
    !elements for each solvent site
    nkTotal = nkTotal*2/sites_per_transform
#else
    !Numerical recipes memory layout
    !-no transpose
    !-extra 2 X nr(2) X nr(3) in k-space starts at nr(1)*nr(2)*nr(3)+1
    nrz = ngr(3)
    nrzOff=0
    nky = ngr(2)
    nkyOff=0
    nkTotal = (ngr(1)+2)*ngr(2)*ngr(3)
    call rism3d_grid_resize(grid,grdspc,&
         ngr,&                                   !r-space global
         (/(ngr(1)+2),ngr(2),ngr(3)/),&          !k-space global
         ngr,&       !r-space local
         (/(ngr(1)+2),ngr(2),ngr(3)/),&   !k-space local
         (/0,0,0/), &                !r-space offset
         (/0,0,0/))                  !k-space offset
#endif /*MPI*/
#  ifdef RISM_DEBUG
    write(0,*) "FFTW Plans Made"
    write(0,*) "nz_local", nrz
    write(0,*) "nz_start", nrzOff
    write(0,*) "nky_local", nky
    write(0,*) "nky_start", nkyOff
    write(0,*) "total_work", nkTotal
#  endif /*RISM_DEBUG*/
  end subroutine rism3d_fft_setgrid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! create a new object to handle 3D forward and backward transforms.
!!!IN:
!!!   this :: fft object
!!!   planner :: FFTW planner type: FFTW_ESTIMATED, FFTW_MEASURE, FFTW_PATIENT, 
!!!              FFTW_EXHAUSTIVE
!!!   localtrans :: transpose data locally.  Not used yet.
!!!   aligned ::  Needed for SIMD.
!!!               .true.  - assume that the memory is 16-byte aligned.  
!!!               .false. - assume that memory is not aligned.
!!!   grid :: grid object that provides the dimensions of all arrays
!!!   inr  :: input r-space array (fwd input)
!!!   outk :: output k-space array (fwd output)
!!!   ink  :: input k-space array (bwd input)
!!!   outr :: output r-space array (bwd output)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism3d_fft_new(this,&
       planner, localtrans, aligned,&
       grid,inr,outr)
    use safemem
    implicit none

    type(rism3d_fft),intent(inout) :: this
    integer, intent(in) :: planner
    type(rism3d_grid), target, intent(in) :: grid
    logical, intent(in) :: localtrans, aligned
    !inr and outr are intent(inout) but an Intel 10.1 compiler bug
    !won't let us say so
    _REAL_,pointer :: inr(:,:),outr(:,:)
    !sites_per_transform : when calling FFTW to do a transform, this
    !                      defines how many arrays are transformed at
    !                      once and has to do with FFTW's 'many'.  This
    !                      only matters for MPI where the memory layout
    !                      has performance issues
    !ncmplx :: number of k-space complex numbers in each dimension for the global 
    !      array. This is typically (/nk(1)/2,nk(2),nk(3)/)
    integer ::sites_per_transform, ncmplx(3)
    _REAL_, pointer :: test
#ifdef MKLdisabled
    integer :: status, commit_status
    integer :: stride(4)
#else /*FFTW*/
    complex(kind(1d0)), pointer :: outk(:,:),ink(:,:)
    outk=> r2c_pointer(inr)
    ink=> r2c_pointer(outr)
#endif /*MKLdisabled FFTW*/

#ifdef RISM_DEBUG
    write(6,*) 'FFT_MAKE_PLANS'
    call flush(6)
#endif /*RISM_DEBUG*/
    this%localtrans = localtrans
    this%aligned=0
    if(.not.aligned)&
         this%aligned = FFT_UNALIGNED
    this%grid => grid
    ncmplx = (/this%grid%ngk(1)/2,this%grid%ngk(2),this%grid%ngk(3)/)
    this%narray = ubound(inr,2)
    if(this%aligned == FFT_ALIGNED)then
       !aligned memory
       sites_per_transform = this%narray
#ifdef MPI
       this%wrk_trans=>safemem_realloc(this%wrk_trans,this%narray*product(this%grid%ngr),.false.)
#endif
    else
       !unaligned memory
       sites_per_transform = 1
    end if
#ifdef MPI
    this%planfwd =  fftw_mpi_plan_many_dft_r2c(3,&
         int(this%grid%ngr(3:1:-1),C_INTPTR_T),&
         int(sites_per_transform,C_INTPTR_T), &
         FFTW_MPI_DEFAULT_BLOCK,&
         FFTW_MPI_DEFAULT_BLOCK,&
!!$           int(ubound(inr,1),C_INTPTR_T), &
!!$           int(ubound(outk,1),C_INTPTR_T),&
         inr, outk,&
         this%grid%mpicomm,ior(ior(planner,FFTW_MPI_TRANSPOSED_OUT),this%aligned))
!!$           comm,FFTW_ESTIMATE)
    this%planbwd =  fftw_mpi_plan_many_dft_c2r(3,&
         int(this%grid%ngr(3:1:-1),C_INTPTR_T),&
         int(sites_per_transform,C_INTPTR_T), &
         FFTW_MPI_DEFAULT_BLOCK,&
         FFTW_MPI_DEFAULT_BLOCK,&
         !           int(ubound(ink,1),C_INTPTR_T),&
!!$           int(1,C_INTPTR_T),&
!!$           int(1,C_INTPTR_T), &
         ink, outr,&
         this%grid%mpicomm,ior(ior(planner,FFTW_MPI_TRANSPOSED_IN), this%aligned))      
!!$           comm,FFTW_ESTIMATE)      
#else    
    this%wrk => safemem_realloc(this%wrk,2*product(this%grid%ngr(2:3)),.false.)
#  ifdef MKLdisabled
    call mkl_check(DftiCreateDescriptor(this%planfwd,DFTI_DOUBLE,DFTI_REAL,3,this%grid%ngr))
    call mkl_check(DftiSetValue(this%planfwd,DFTI_INPUT_STRIDES,&
         (/0,1,this%grid%ngk(1), product(this%grid%ngk(1:2))/)))
    call mkl_check(DftiSetValue(this%planfwd,DFTI_OUTPUT_STRIDES,&
         (/0,1,this%grid%ngk(1)/2, product(this%grid%ngk(1:2))/2/)))
    call mkl_check(DftiSetValue(this%planfwd,DFTI_NUMBER_OF_TRANSFORMS,sites_per_transform))
    call mkl_check(DftiSetValue(this%planfwd,DFTI_INPUT_DISTANCE,product(this%grid%ngk)))
    call mkl_check(DftiSetValue(this%planfwd,DFTI_OUTPUT_DISTANCE,product(this%grid%ngk)))
    call mkl_check(DftiSetValue(this%planfwd,DFTI_CONJUGATE_EVEN_STORAGE,DFTI_COMPLEX_COMPLEX))
    call mkl_check(DftiCommitDescriptor(this%planfwd))

    call mkl_check(DftiCopyDescriptor(this%planfwd,this%planbwd))
    call mkl_check(DftiSetValue(this%planbwd,DFTI_INPUT_STRIDES,&
         (/0,1,this%grid%ngk(1)/2, product(this%grid%ngk(1:2))/2/)))
    call mkl_check(DftiSetValue(this%planbwd,DFTI_OUTPUT_STRIDES,&
         (/0,1,this%grid%ngk(1), product(this%grid%ngk(1:2))/)))
    call mkl_check(DftiCommitDescriptor(this%planbwd))

    call mkl_print_settings(this%planfwd)
    call mkl_print_settings(this%planbwd)
#  else /*FFTW*/
    this%planfwd =  fftw_plan_many_dft_r2c(3,this%grid%ngr(3:1:-1),sites_per_transform, &
         inr(:,1),this%grid%ngk(3:1:-1),1,product(this%grid%ngk),&
         outk(:,1),ncmplx(3:1:-1),1,product(ncmplx),&
         ior(planner,this%aligned))

    this%planbwd =  fftw_plan_many_dft_c2r(3,this%grid%ngr(3:1:-1),sites_per_transform, &
         ink(:,1),ncmplx(3:1:-1),1,product(ncmplx),&
         outr(:,1),this%grid%ngk(3:1:-1),1,product(this%grid%ngk),&
         ior(planner,this%aligned))
#  endif /*MKLdisabled FFTW*/
#endif
  end subroutine rism3d_fft_new

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Frees memory and destroys the object
!!!IN:
!!!   this :: the rism3d_fft object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_fft_destroy(this)
      use safemem
      implicit none
      type(rism3d_fft), intent(inout) :: this
#ifdef MKLdisabled
      if(associated(this%planfwd))&
           call mkl_check(DftiFreeDescriptor(this%planfwd))
      if(associated(this%planbwd))&
           call mkl_check(DftiFreeDescriptor(this%planbwd))
#else /*FFTW*/
      if(C_associated(this%planfwd))then
         call fftw_destroy_plan(this%planfwd)
         this%planfwd=C_NULL_PTR
      end if
      if(C_associated(this%planbwd))then
         call fftw_destroy_plan(this%planbwd)
         this%planbwd=C_NULL_PTR
      end if
#endif /*MKLdisabled FFTW*/
      if(associated(this%grid))&
           nullify(this%grid)
      if(safemem_dealloc(this%wrk_trans) /=0)&
           call rism_report_error("Could not deallocate FFT working memory")
      if(safemem_dealloc(this%wrk) /=0 )&
           call rism_report_error("Could not deallocate FFT working memory")
    end subroutine rism3d_fft_destroy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Interface to forward 3-D Fast Fourier Transform of Real Data.
!!! Normalizing is performed at this stage
!!! datar and datak may point to the same memory
!!!IN:
!!!   this :: the rism3d_fft object
!!!   datar :: local array of _REAL_ data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_fft_fwd(this,datar)
      implicit none
      type(rism3d_fft), intent(inout) :: this
      !datar is intent(inout) but an Intel 10.1 compiler bug won't let
      !us say so
      _REAL_,pointer :: datar(:,:)
      call rlft3i (this,datar, 1)
    end subroutine rism3d_fft_fwd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Interface to backward 3-D Fast Fourier Transform of Real Data.
!!! Normalizing is not performed at this stage
!!! datar and datak may point to the same memory
!!!IN:
!!!   this :: the rism3d_fft object
!!!   datar :: local array of _REAL_ data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_fft_bwd(this,datar)
      implicit none
      type(rism3d_fft), intent(inout) :: this
      !datak is intent(inout) but an Intel 10.1 compiler bug won't let
      !us say so
      _REAL_,pointer :: datar(:,:)
      call rlft3i (this,datar, -1)
    end subroutine rism3d_fft_bwd


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!PRIVATE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef MKLdisabled
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!tests the output status of MKL DFT commands and triggers the
!!!appropriate error if something goes wrong
!!!IN:
!!!   status :: MKL DTF status
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mkl_check(status)
      implicit none
      integer, intent(in) :: status
      if (.NOT. DftiErrorClass(Status, DFTI_NO_ERROR)) then
         call rism_report_error("(a,i4)",trim(DftiErrorMessage(status)),status)
      end if
    end subroutine mkl_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Writes MKL DFT settings to the error unit
!!!IN:
!!!   plan :: MKL plan
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mkl_print_settings(plan)
      implicit none
      type(DFTI_DESCRIPTOR), pointer, intent(in) :: plan
      integer :: unit
      integer :: status, commit_status
      integer :: array(4)
      character(len=64) :: whtspc
      character(len=DFTI_VERSION_LENGTH) :: version
      _REAL_ :: scale

      write(whtspc,'(a64)')" "
      unit = rism_report_geteunit()

      status = DftiGetValue(plan,DFTI_PRECISION,array(1))
      if(array(1) == DFTI_DOUBLE)then
         write(unit,'(a24,a)') "PRECISION:"//whtspc,"DFTI_DOUBLE"
      elseif(array(1) == DFTI_SINGLE)then
         write(unit,'(a24,a)') "PRECISION:"//whtspc,"DFTI_SINGLE"
      else
         write(unit,'(a24,a)') "PRECISION:"//whtspc,"UNKNOWN"
      endif

      status = DftiGetValue(plan,DFTI_FORWARD_DOMAIN,array(1))
      if(array(1) == DFTI_REAL)then
         write(unit,'(a24,a)') "FORWARD DOMAIN:"//whtspc,"DFTI_REAL"
      elseif(array(1) == DFTI_COMPLEX)then
         write(unit,'(a24,a)') "FORWARD DOMAIN:"//whtspc,"DFTI_COMPLEX"
      else
         write(unit,'(a24,a)') "FORWARD DOMAIN:"//whtspc,"UNKNOWN"
      endif

      status = DftiGetValue(plan,DFTI_DIMENSION,array(1))
      write(unit,'(a24,4i6)') "DIMENSION:"//whtspc, array(1)

      status = DftiGetValue(plan,DFTI_LENGTHS,array)
      write(unit,'(a24,4i6)') "LENGTHS:"//whtspc, array(1:3)

      status = DftiGetValue(plan,DFTI_PLACEMENT,array(1))
      if(array(1) == DFTI_INPLACE)then
         write(unit,'(a24,a)') "PLACEMENT:"//whtspc,"DFTI_INPLACE"
      elseif(array(1) == DFTI_NOT_INPLACE)then
         write(unit,'(a24,a)') "PLACEMENT:"//whtspc,"DFTI_NOTINPLACE"
      else
         write(unit,'(a24,a)') "PLACEMENT:"//whtspc,"UNKNOWN"
      endif

      status = DftiGetValue(plan,DFTI_FORWARD_SCALE,scale)
      write(unit,'(a24,1p,e16.8)') "FORWARD SCALE:"//whtspc, scale

      status = DftiGetValue(plan,DFTI_BACKWARD_SCALE,scale)
      write(unit,'(a24,1p,e16.8)') "BACKWARD SCALE:"//whtspc, scale

      status = DftiGetValue(plan,DFTI_NUMBER_OF_USER_THREADS,array(1))
      write(unit,'(a24,4i6)') "NUMBER OF THREADS:"//whtspc, array(1)

      status = DftiGetValue(plan,DFTI_INPUT_STRIDES,array)
      write(unit,'(a24,4i6)') "INPUT STRIDE:"//whtspc, array

      status = DftiGetValue(plan,DFTI_OUTPUT_STRIDES,array)
      write(unit,'(a24,4i6)') "OUTPUT STRIDE:"//whtspc, array

      status = DftiGetValue(plan,DFTI_NUMBER_OF_TRANSFORMS,array(1))
      write(unit,'(a24,4i6)') "NUMBER OF TRANSFORMS:"//whtspc, array(1)

      status = DftiGetValue(plan,DFTI_INPUT_DISTANCE,array(1))
      write(unit,'(a24,4i6)') "INPUT DISTANCE:"//whtspc, array(1)

      status = DftiGetValue(plan,DFTI_OUTPUT_DISTANCE,array(1))
      write(unit,'(a24,4i6)') "OUTPUT DISTANCE:"//whtspc, array(1)

      status = DftiGetValue(plan,DFTI_REAL_STORAGE,array(1))
      if(array(1) == DFTI_REAL_REAL)then
         write(unit,'(a24,a)') "REAL STORAGE:"//whtspc,"DFTI_REAL_REAL"
      else
         write(unit,'(a24,a)') "REAL STORAGE:"//whtspc,"UNKNOWN"
      endif

      status = DftiGetValue(plan,DFTI_CONJUGATE_EVEN_STORAGE,array(1))
      if(array(1) == DFTI_COMPLEX_REAL)then
         write(unit,'(a24,a)') "CONJUGATE EVEN STORAGE:"//whtspc,"DFTI_COMPLEX_REAL"
      elseif(array(1) == DFTI_COMPLEX_COMPLEX)then
         write(unit,'(a24,a)') "CONJUGATE EVEN STORAGE:"//whtspc,"DFTI_COMPLEX_COMPLEX"
      else
         write(unit,'(a24,a)') "CONJUGATE EVEN STORAGE:"//whtspc,"UNKNOWN"
      endif

      status = DftiGetValue(plan,DFTI_PACKED_FORMAT,array(1))
      if(array(1) == DFTI_CCE_FORMAT)then
         write(unit,'(a24,a)') "PACKED FORMAT:"//whtspc,"DFTI_CCE_FORMAT"
      elseif(array(1) == DFTI_CCS_FORMAT)then
         write(unit,'(a24,a)') "PACKED FORMAT:"//whtspc,"DFTI_CCS_FORMAT"
      elseif(array(1) == DFTI_PACK_FORMAT)then
         write(unit,'(a24,a)') "PACKED FORMAT:"//whtspc,"DFTI_PACK_FORMAT"
      elseif(array(1) == DFTI_PERM_FORMAT)then
         write(unit,'(a24,a)') "PACKED FORMAT:"//whtspc,"DFTI_PERM_FORMAT"
      else
         write(unit,'(a24,a)') "PACKED FORMAT:"//whtspc,"UNKNOWN"
      endif

!!$      status = DftiGetValue(plan,DFTI_WORKSPACE,array(1))
!!$      if(array(1) == DFTI_ALLOW)then
!!$         write(unit,'(a24,a)') "WORKSPACE:"//whtspc,"DFTI_ALLOW"
!!$      elseif(array(1) == DFTI_AVOID)then
!!$         write(unit,'(a24,a)') "WORKSPACE:"//whtspc,"DFTI_AVOID"
!!$      else
!!$         write(unit,'(a24,a)') "WORKSPACE:"//whtspc,"UNKNOWN"
!!$      endif

      status = DftiGetValue(plan,DFTI_ORDERING,array(1))
      if(array(1) == DFTI_ORDERED)then
         write(unit,'(a24,a)') "ORDERING:"//whtspc,"DFTI_ORDERED"
      elseif(array(1) == DFTI_BACKWARD_SCRAMBLED)then
         write(unit,'(a24,a)') "ORDERING:"//whtspc,"DFTI_BACKWARD_SCRAMBLED"
      else
         write(unit,'(a24,a)') "ORDERING:"//whtspc,"UNKNOWN"
      endif

      status = DftiGetValue(plan,DFTI_COMMIT_STATUS,array(1))
      if(array(1) == DFTI_COMMITTED)then
         write(unit,'(a24,a)') "COMMIT STATUS:"//whtspc,"DFTI_COMMITTED"
      elseif(array(1) == DFTI_UNCOMMITTED)then
         write(unit,'(a24,a)') "COMMIT STATUS:"//whtspc,"DFTI_UNCOMMITTED"
      else
         write(unit,'(a24,a)') "COMMIT STATUS:"//whtspc,"UNKNOWN"
      endif

      status = DftiGetValue(plan,DFTI_VERSION,version)
      write(unit,'(a24,a)') "VERSION:"//whtspc, trim(version)

     write(unit,"(a24,a,i4)") "STATUS:"//whtspc, trim(DftiErrorMessage(status)),status
      
    end subroutine mkl_print_settings

#else /*FFTW*/
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!takes a _REAL_ pointer and returns a complex pointer of the same
!!!size with the leading dimension cut in half.  This is the memory
!!!layout required for FFTW.  This uses Fortan 2003 features.
!!!IN:
!!!   rdata :: real data (ndata,narray)
!!!OUT:
!!!    pointer to complex data (ndata/2,narray)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     function r2c_pointer(rdata) result (cdata)
      _REAL_, pointer, intent(in) :: rdata(:,:)
      complex(kind(1d0)),pointer :: cdata(:,:)
      type(c_ptr) :: cptr
      nullify(cdata)
      cptr = c_loc(rdata(1,1))
      call c_f_pointer(cptr, cdata, [ubound(rdata,1)/2, ubound(rdata,2)])
    end function r2c_pointer
#endif

#ifdef MPI
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Interface to MPI 3-Dimensional Fast Fourier Transform of Real Data.
!!! Normalizing while Forward FFT (KEY=1)
!!! datar and datak may point to the same memory
!!!IN:
!!!   this :: the rism3d_fft object
!!!   datar :: local array of _REAL_ data
!!!   key   :: FWD-> 1 and BWD-> -1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine  rlft3i (this,datar, key)
      implicit none
      type(rism3d_fft), intent(inout) :: this
      !datar is intent(inout) but an Intel 10.1 compiler bug won't let
      !us say so
      _REAL_,pointer :: datar(:,:)
      integer,intent(in) ::  key
      integer :: iarray, ir, i, j, ix, iy, iz
      integer :: nr,nk, dim(3)

#  ifdef RISM3D_DEBUG
      integer :: irank, ierr
      character(len=20) :: label
#  endif
#ifdef MKLdisabled
#else
      complex(kind(1d0)), pointer :: datak(:,:)
      datak => r2c_pointer(datar)
#endif

      !....................... R-space data array size .......................
#  ifdef RISM_DEBUG
      write(0,*) "rlft3i",n_local
      call flush(0)
#  endif /*RISM_DEBUG*/
      nr = ubound(datar,1)
      nk = nr/2

      dim = this%grid%ngr
      dim(1) = dim(1) + 2
      if (key == 1)  then
         datar = datar / product(this%grid%ngr)
#  ifdef RISM3D_DEBUG
      write(label,"(a,i2)") "int", key
      call rism3d_debug_print(datar,.true., label,.false.)
#  endif      
!!$#ifdef MANY
         if(this%aligned == FFT_ALIGNED)then
            !transpose site and spatial indices before and after transform
            do i = 1,this%narray
               do j=1,nr
                  this%wrk_trans(i+(j-1)*this%narray) = datar(j,i)
               end do
            end do
            call dcopy(this%narray*nr,this%wrk_trans,1,datar,1)
!!$         do iarray =1, ubound(datar,2)
!!$            do ir=1, ubound(datar,1)
!!$               write(0,*) "ink",iarray,ir, -1, datar(ir,iarray)
!!$            end do
!!$         end do
#  ifdef MKLdisabled
#  else /*FFTW*/
            call fftw_mpi_execute_dft_r2c(this%planfwd, datar, datak)
#  endif /*MKLdisabled FFTW*/
!!$         do iarray =1, ubound(datar,2)
!!$            do ir=1, ubound(datar,1)
!!$               write(0,*) "OUk",iarray,ir, 1, datar(ir,iarray)
!!$            end do
!!$         end do
            call dcopy(this%narray*nr,datak,1,this%wrk_trans,1)
            do i = 1,this%narray
               do j=1,nk
                  datak(j,i) = cmplx(this%wrk_trans(1+(i-1)*2 + (j-1)*this%narray*2),&
                       this%wrk_trans(2+(i-1)*2 + (j-1)*this%narray*2))
               end do
            end do
         else
            do iarray =1, ubound(datar,2)
#  ifdef MKLdisabled
#  else /*FFTW*/
               call fftw_mpi_execute_dft_r2c(this%planfwd, datar(:,iarray), datak(:,iarray))
#  endif /*MKLdisabled FFTW*/
            end do
         end if
!!$         do iarray =1, ubound(datar,2)
!!$            do ir=1, ubound(datar,1)
!!$               write(0,*) "OUT",iarray,ir, 1, datar(ir,iarray)
!!$            end do
!!$         end do

#  ifdef RISM3D_DEBUG
      write(label,"(a,i2)") "out", key
      call rism3d_debug_print(datar,.false., label,.false.)
#  endif      
      elseif (key == -1)  then
#  ifdef RISM3D_DEBUG
      write(label,"(a,i2)") "int", key
      call rism3d_debug_print(datar,.false., label,.false.)
#  endif      
!!$         do iarray =1, ubound(datar,2)
!!$            write(0,*) "in",ir, -1, sum(datar(:,iarray)), ubound(datar),ubound(datak)
!!$         end do
!!$         do iarray =1, ubound(datar,2)
!!$            do ir=1, ubound(datar,1)
!!$               write(0,*) "int",iarray,ir, -1, datar(ir,iarray)
!!$            end do
!!$         end do
!!$#ifdef MANY
         if(this%aligned==FFT_ALIGNED)then
            do i = 1,this%narray
               do j=1,nk
                  this%wrk_trans(1+(i-1)*2 + (j-1)*this%narray*2) = real(datak(j,i))
                  this%wrk_trans(2+(i-1)*2 + (j-1)*this%narray*2) = aimag(datak(j,i))
               end do
            end do
            call dcopy(this%narray*nr,this%wrk_trans,1,datak,1)
!!$         do iarray =1, ubound(datar,2)
!!$            do ir=1, ubound(datar,1)
!!$               write(0,*) "ink",iarray,ir, -1, datar(ir,iarray)
!!$            end do
!!$         end do
#  ifdef MKLdisabled
#  else /*FFTW*/
            call fftw_mpi_execute_dft_c2r(this%planbwd, datak, datar)
#  endif /*MKLdisabled FFTW*/
!!$         do iarray =1, ubound(datar,2)
!!$            do ir=1, ubound(datar,1)
!!$               write(0,*) "ouk",iarray, ir,-1, datar(ir,iarray)/product(this%grid%ngr)
!!$            end do
!!$         end do
            call dcopy(this%narray*nr,datar,1,this%wrk_trans,1)
            do i = 1,this%narray
               do j=1,nr
                  datar(j,i) = this%wrk_trans(i + (j-1)*this%narray)
               end do
            end do
         else
            do iarray =1, ubound(datar,2)
#  ifdef MKLdisabled
#  else /*FFTW*/
               call fftw_mpi_execute_dft_c2r(this%planbwd, datak(:,iarray), datar(:,iarray))
#  endif /*MKLdisabled FFTW*/
            end do
         end if
#  ifdef RISM3D_DEBUG
      write(label,"(a,i2)") "out", key
      call rism3d_debug_print(datar,.true., label,.false.)
#  endif      
!!$         do iarray =1, ubound(datar,2)
!!$            write(0,*) "out", -1, sum(data(1:this%grid%ngr(1),:,:,iarray))
!!$         end do
      end if
      
    end subroutine rlft3i

#else /*defined(MPI)*/
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Interface to 3-Dimensional Fast Fourier Transform of Real Data.
!!! Normalizing while Forward FFT (KEY=1)
!!! datar and datak may point to the same memory
!!!IN:
!!!   this :: the rism3d_fft object
!!!   datar :: local array of real data
!!!   key   :: FWD-> 1 and BWD-> -1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine  rlft3i (this,datar,key)
      implicit none
      type(rism3d_fft), intent(inout) :: this
      integer,intent(in) ::  key
      !datar is intent(inout) but an Intel 10.1 compiler bug won't let
      !us say so
      _REAL_, pointer :: datar(:,:)
      integer ::  ngr_total, ngk_total, igx, igy, igz, iv, ig1
      integer :: ir, iarray, ix, iy, iz, dim(3)
#  ifdef RISM3D_DEBUG
      integer :: irank, ierr
      character(len=20) :: label
#  endif
#ifdef MKLdisabled
#else
      complex(kind(1d0)), pointer :: datak(:,:)
      datak => r2c_pointer(datar)
#endif
      !....................... R-space data array size .......................
      ngr_total = product(this%grid%ngr)
      ngk_total = (this%grid%ngr(1)+2)*product(this%grid%ngr(2:3))
      dim = this%grid%ngr
      dim(1) = dim(1)+2
      if (key == 1)  then
         do iarray = 1, this%narray
            !complex conjugate
            call dscal(ngk_total/2,-1d0,datar(2,iarray),2)
            !normalize
            call dscal(ngr_total,1d0/dble(ngr_total),datar(1,iarray),1)
            call fromnr(this,datar(:,iarray),this%grid%ngr(1),this%grid%ngr(2),this%grid%ngr(3))
         end do
#  ifdef RISM3D_DEBUG
      write(label,"(a,i2)") "int", key
      call rism3d_debug_print(datar,.true., label,.false.)
#  endif      
#  ifdef MKLdisabled
      if(this%aligned == FFT_UNALIGNED)then
         do iarray =1, ubound(datar,2)
            call mkl_check(DftiComputeForward(this%planfwd, datar(:,iarray)))
         end do
      else
!!$         call mkl_check(DftiComputeForward(this%planfwd, datar))
      end if
#  else /*FFTW*/
      if(this%aligned == FFT_UNALIGNED)then
         do iarray =1, ubound(datar,2)
            call fftw_execute_dft_r2c(this%planfwd, datar(:,iarray), datak(:,iarray))
         end do
      else
         call fftw_execute_dft_r2c(this%planfwd, datar, datak)
      end if
#  endif /*MKLdisabled FFTW*/
#  ifdef RISM3D_DEBUG
      write(label,"(a,i2)") "out", key
      call rism3d_debug_print(datar,.true., label,.false.)
#  endif      
         do iarray = 1, this%narray
            !work with the knowledge that datar and datak are the same data
            call tonr(this,datar(:,iarray),this%grid%ngr(1),this%grid%ngr(2),this%grid%ngr(3))
         end do
      elseif (key == -1)  then
         do iarray = 1, this%narray
            !work with the knowledge that datar and datak are the same data
            call fromnr(this,datar(:,iarray),this%grid%ngr(1),this%grid%ngr(2),this%grid%ngr(3))
         end do
#  ifdef RISM3D_DEBUG
      write(label,"(a,i2)") "int", key
      call rism3d_debug_print(datar,.true., label,.false.)
#  endif      
#  ifdef MKLdisabled
      if(this%aligned == FFT_UNALIGNED)then
         do iarray =1, ubound(datar,2)
            call mkl_check(DftiComputeBackward(this%planbwd, datar(:,iarray)))
         end do
      else
!!$         call mkl_check(DftiComputeBackward(this%planbwd, datar))
      end if
#  else /*FFTW*/
      if(this%aligned == FFT_UNALIGNED)then
         do iarray =1, ubound(datar,2)
            call fftw_execute_dft_c2r(this%planbwd, datak(:,iarray), datar(:,iarray))
         end do
      else
         call fftw_execute_dft_c2r(this%planbwd, datak, datar)
      end if
#  endif /*MKLdisabled FFTW*/
#  ifdef RISM3D_DEBUG
      write(label,"(a,i2)") "out", key
      call rism3d_debug_print(datar,.true., label,.false.)
#  endif      
         do iarray = 1, this%narray
            call tonr(this,datar(:,iarray),this%grid%ngr(1),this%grid%ngr(2),this%grid%ngr(3))
            call dscal(ngk_total/2,-1d0,datar(2,iarray),2)
         end do
      endif
    end subroutine rlft3i
#endif /*defined(MPI)*/

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Convert 3D FFT data from FFTW layout to Numerical Recipes layout.  Both 
!!!layouts require (nx+2)*ny*nz elements.  The difference is how the Nyquist
!!!frequency data in reciprocal space is stored.
!!!
!!!NR Layout:
!!!Data is layed out as a contiguous piece of nx*ny*nz memory.  The Nyquist 
!!!frequency data (2*ny*nz) starts at the end of this block.
!!!
!!!FFTW Layout:
!!!Nyqusit frequency data is stored as two extra elements at the end of the first
!!!dimension.
!!!
!!!The strategy employed is to travel through the array, extracting Nyquist 
!!!frequency data to a temporary array, and shift the data to fill in the gaps
!!!as we go.  Once completed, the temporary array is then copied into data.
!!!This does require allocating temporary memory for the Nyquist 
!!!frequencies but this is less than using RESHAPE which creates a full sized
!!!array and copies data back and forth.  Furthermore, we can put this memory on
!!!the heap and avoid exhausting the stack.
!!!IN: 
!!!   this :: the rism3d_fft object
!!!   data :: 3D FFT array in NR format.  On return, in FFTW format
!!!   n1   :: number of x points
!!!   n2   :: number of y points
!!!   n3   :: number of z points
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine tonr(this, data, n1,n2,n3)
      implicit none
      type(rism3d_fft),intent(inout) :: this
      integer, intent(in) :: n1,n2,n3
      _REAL_, intent(inout) ::  data((n1+2)*n2*n3)
      integer :: i,j,k,i1,j1,k1,ngr,ngk, ispec,ondx, nndx
      ngr = n1*n2*n3
      ngk = (n1+2)*n2*n3

      ispec=1
      do k = 1,n3
         do j=1,n2
            ondx = 1+ (j-1)*(n1+2) + (k-1)*(n1+2)*n2
            nndx = 1+ (j-1)*(n1)   + (k-1)*(n1)*n2
            this%wrk(ispec:ispec+1) = data(ondx+n1:ondx+n1+1)
            data(nndx:nndx+n1-1) = data(ondx:ondx+n1-1) 
            ispec = ispec+2
        end do
      end do
      data(ngr+1:ngk) = this%wrk
!!$      data=RESHAPE((/                                                      &
!!$           (((data(i+(j-1)*(n1+2)+(k-1)*(n1+2)*n2),                    &
!!$           i=1,n1),j=1,n2),k=1,n3),            &
!!$           (((data(n1+i1+(j1-1)*(n1+2)+(k1-1)*(n1+2)*n2),              &
!!$           i1=1,2),j1=1,n2),k1=1,n3)           &
!!$           /), SHAPE=(/ngr+2*n2*n3/))
      data(2:ngr+2*n2*n3:2) = -data(2:ngr+2*n2*n3:2)
    end subroutine tonr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Convert 3D FFT data from Numerical Recipes layout to FFTW layout.  Both 
!!!layouts require (nx+2)*ny*nz elements.  The difference is how the Nyquist
!!!frequency data in reciprocal space is stored.
!!!
!!!NR Layout:
!!!Data is layed out as a contiguous piece of nx*ny*nz memory.  The Nyquist 
!!!frequency data (2*ny*nz) starts at the end of this block.
!!!
!!!FFTW Layout:
!!!Nyqusit frequency data is stored as two extra elements at the end of the first
!!!dimension.
!!!
!!!The strategy employed is to move data from the end of the NR array, n1 elements 
!!!at a time, to the end of the aggregate array, filling in Nyquist frequencies 
!!!as we go.  This does require allocating temporary memory for the Nyquist 
!!!frequencies but this is less than using RESHAPE which creates a full sized
!!!array and copies data back and forth.  Furthermore, we can put this memory on
!!!the heap and avoid exhausting the stack.
!!!IN: 
!!!   this :: the rism3d_fft object
!!!   data :: 3D FFT array in NR format.  On return, in FFTW format
!!!   n1   :: number of x points
!!!   n2   :: number of y points
!!!   n3   :: number of z points
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine fromnr(this,data, n1,n2,n3)
      implicit none

      type(rism3d_fft),intent(inout) :: this
      integer, intent(in) :: n1,n2,n3
      _REAL_,intent(inout) ::  data((n1+2)*n2*n3)
      integer :: i,j,k, ngr,ngk, ispec, ondx,nndx
      ngr = n1*n2*n3
      ngk = (n1+2)*n2*n3
      data(2:ngr+2*n2*n3:2) = -data(2:ngr+2*n2*n3:2)
      this%wrk = data(ngr+1:ngk)
      ispec = size(this%wrk)-1
      do k = 1,n3
         do j=1,n2
            ondx = ngr +1 - (n1 + (j-1)*n1 + (k-1)*n1*n2)
            nndx =  ngk  +1 - (1 + size(this%wrk) - ispec + n1 + (j-1)*n1 + (k-1)*n1*n2)
            data(nndx:nndx+n1-1) = data(ondx:ondx+n1-1) 
            data(nndx+n1:nndx+n1+1) = this%wrk(ispec: ispec+1)
            ispec = ispec-2
         end do
      end do
!!$      data = RESHAPE((/                                                    &
!!$           (((data(i+(j-1)*n1+(k-1)*n1*n2),i=1, n1),                  &
!!$           data(ngr+1+(j-1)*2+(k-1)*2*n2),                         &
!!$           data(ngr+2+(j-1)*2+(k-1)*2*n2),j=1,n2),k=1,n3 )         &
!!$           /),SHAPE=(/ngr+2*n2*n3/))

    end subroutine fromnr

  end module rism3d_fft_c
