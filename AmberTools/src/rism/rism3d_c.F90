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
!3D-RISM solver.  
!This is defines the 3D-RISM type and associated subroutines.  All type elements 
!are public.  In general, read but do not write these variables.  This provides
!an object-orientiented interface without needing a function to access every 
!variable.
!Features of this solver include:
! o Multiple closures: KH, HNC, PSE-n
! o MDIIS accelerated solutions
! o Optional cutoffs
! o Supercell method for long range asymptotics
! o Analytic forces
! o Variable grid size and dynamic memory allocation
! o MPI support
! o units:  energy       [kT]
!           distances    [A]         (Angstroms)
!           site charges [sqrt(kT A)]
!           temperature  [K]
!           density      [#/A^3]
!           mass         [au]
! o To convert [e] to [sqrt(kT A)] *sqrt(COULOMB_CONST_E/KB/temperature)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "../include/dprec.fh"

module rism3d_c
    use rism3d_solu_c
    use rism3d_solv_c
    use rism3d_potential_c
    use rism3d_grid_c
    use rism3d_closure_c
    use rism_report_c
    use rism_timer_c
    use mdiis_c
    use rism3d_fft_c
#ifdef RISM3D_DEBUG
!    use rism3d_debug_c
#endif
    implicit none
#include "def_time.h"

    type rism3d
!       sequence

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Solute/solvent information
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !solu  :: solute object
       type(rism3d_solu) :: solu
       !solv  :: solvent object
       type(rism3d_solv) :: solv

       !pot   :: potential object
       type(rism3d_potential) :: pot
       !grid  :: grid object
       type(rism3d_grid) :: grid
       !closure  :: closure object
       type(rism3d_closure) :: closure
       !closureList :: list of closure names to use in order.  Only
       !               the last closure is used for thermodynamic
       !               output.  This can be used to progressively
       !               increase the order of the closure to aid
       !               convergence.
       character(len=8), pointer :: closureList(:)=> NULL()

       !TIMERS.  Subtimers only account for computation.  We ignore setup etc.
       !timer :: timer for this class.  Activated for all public routines
       !resizeTimer :: time to resize solvent box
       !reorientTimer :: time to reorient solute
       !cuvpropTimer :: time to propagate Cuv solution
       !fftTimer :: specifically times FFT calculation
       !solveTimer :: specifically times rism1d_solve calculation
       !rxrismTimer :: specifically times rxrism calculation
       !r1rismTimer :: specifically times r1rism calculation
       !thermoTimer :: specifically times thermodynamics calculations
       !forceTimer :: specifically times force calculation
       !exchemTimer :: specificall times excess chemical potential calculation
       type(rism_timer) :: timer, resizeTimer, reorientTimer, &
            cuvpropTimer, fftTimer, solveTimer,&
            rxrismTimer, r1rismTimer, thermoTimer, forceTimer, exchemTimer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!private !(should be)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
       !FFTW options
       !fftw_planner :: FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT, FFTW_EXHAUSTIVE
       integer :: fftw_planner=FFT_ESTIMATE
       !fft_aligned :: .true.  - use aligned memory and to enable SIMD; 
       !               .false. - don't use aligned memory
       logical :: fft_aligned=.true.
       !fftw_localtrans :: transpose site number and spatial data locally before and after fft
       logical :: fftw_localtrans=.true.

       !verbose :: 0 - no ouput
       !           1 - memory allocation and steps for convergence
       !           2 - 1 + convergence progress
       integer :: verbose=0

       !This is a bit ugly and there may be a better solution.  We
       !need to keep track of the number of solutions for both charged
       !and un-charged solutes.  When we change between the solutes we
       !set point the nsolutions pointer to the appropriate variable.
       !However, the 'target' attribute is not allow in type
       !definitions so these variables have to be pointers and we have
       !to allocate memory for them.
       !nsolution :: number of times full solutions have been calculated
       !nsolutionChg   :: number of times full solutions with a charged solute 
       !                   have been calculated
       !nsolutionNoChg :: number of times full solutions with an uncharged solute
       !                   have been calculated
       integer, pointer :: nsolution=>NULL(), nsolutionChg=>NULL(), nsolutionNoChg=>NULL()

       !ucenter :: center the solute in the solvation box.  
       !           0 - off
       !           1 - center of mass
       !           2 - center of geometry
       !           3 - center of mass shifted to the nearest grid point
       !           4 - center of geometry shifted to the nearest grid point
       !           For negative numbers the centering translation is only calculated
       !           the for the first solution and used for subsequent calculations.
       !           This allows the solute to drift in the box.
       integer :: ucenter = 1

       !ratucm   :: center of mass (geometry) of the solute [A]
       _REAL_ :: ratucm(3)

       !ncuvsteps :: number of past cuv time steps saves
       integer :: ncuvsteps

       !buffer   :: buffer distance to the edge of the box for the solvent [A]
       !boxfixlen:: fixed box size for 3d-rism.  for PBC calculations, these
       !            should generallly be equal [A]
       _REAL_ :: buffer=12d0,boxfixlen(3)
       !nboxfix(3)    :: number of grid points in each dimension for a fixed box size
       integer ::  nboxfix(3)

       !varbox  :: variable box size
       logical :: varbox = .true.

       !mdiis_method :: MDIIS implementation to use    
       !NVec       :: number of vectors used for MDIIS (consequently, the number of copies of
       !             CUV we need to keep for MDIIS)
       integer :: NVec,mdiis_method
       type(mdiis) :: mdiis_o

       !deloz :: 'step size' for MDIIS
       !mdiis_restart :: restart threshold factor. Ratio of the current residual to the 
       !                 minimum residual in the basis that causes a restart
       _REAL_ :: deloz=0.7d0, mdiis_restart

!!!!!!!!!!!!!!!!!!
!!! MPI Support !!
!!!!!!!!!!!!!!!!!!
       integer :: mpirank=0, mpicomm=0, mpisize=1
!!!!!!!!!!!!!!!!!!!
!!! LARGE ARRAYS !!
!!!!!!!!!!!!!!!!!!!

       !
       ! all arrays are declared as pointers to ensure we can reallocate them as necessary
       !
       
       !xvva       :: solvent chi interpolated for our grid size
       !guv        :: solvent distribution function
       !huv        :: guv - 1
       !cuv        :: solvent direct correlation function and points to the
       !              current active solution in cuvWRK
       !cuvres     :: residual value for cuv calculation and points to the
       !              current active solution in cuvresWRK.
       !cuvWRK     :: Working Cuv memory.  Holds Cuv from previous iterations.
       !cuvresWRK  :: Working Cuvres memory.  Holds Cuvres from previous iterations.
       !oldcuv     :: previous solutions of cuv. Points to oldcuvChg
       !              or oldcuvNoChg depending on the charge state of the
       !calculation
       !oldcuvChg  :: previous solutions for the standard charged system
       !oldcuvNoChg:: previous solutions for the chargeless system.  This is only allocated
       !           :: if _unsetCharges() is called
       _REAL_,pointer :: xvva(:)=>NULL(),&
            oldcuv(:,:,:,:,:)=>NULL(),&
            oldcuvChg(:,:,:,:,:)=>NULL(),&
            oldcuvNoChg(:,:,:,:,:)=>NULL(),&
            cuv(:,:,:,:)=>NULL(), cuvWRK(:,:,:,:,:)=>NULL(),&
            cuvres(:,:)=>NULL(), cuvresWRK(:,:,:)=>NULL()
            
       
       _REAL_,pointer :: guv(:,:)=>NULL(),huv(:,:)=>NULL()

       !fft :: fft object for standard 3D-RISM solution
       type(rism3d_fft) :: fft
    end type rism3d

    interface rism3d_setbox
       module procedure rism3d_setbox_variable, rism3d_setbox_fixed
    end interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Public suroutine and functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    public :: rism3d_new, rism3d_destroy, rism3d_solve, rism3d_force, &
         rism3d_exchem_tot, rism3d_exchem,&
         rism3d_exchemGF_tot, rism3d_exchemGF,&
         rism3d_setbox,&
         rism3d_setclosure, rism3d_setverbosity, rism3d_setcut, rism3d_setmdiis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Private suroutine and functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    private :: resizebox, reallocbox, &
         interpolate_xvva, center_solute, rxrism, r1rism, &
         force_lj, force_coulomb, cuv_propagate, cuv_update
    
  contains
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutines
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Constructor - precalculates the solute solvent terms that are not 
!!!configuration dependent and sets box parameters.
!!!
!!!The solvation box may be fixed size or variable.  For fixed size,
!!!define o_boxlen and o_ng3.  For variable box size, define buffer
!!!and grdspc.  Do not mix these parameters as this will cause the
!!!program to halt.
!!!
!!!If this is an MPI run, supply the MPI communicator.  Only the rank
!!!0 parameters will be used. However, due to the limitations of
!!!pre-Fortran2003, the closure must be the same length on all
!!!processes. The values and number of elements for the closure list
!!!on non-rank 0 processes still do not matter.
!!!
!!!IN:
!!!   this :: new rism3d object
!!!   solu :: 3D-RISM solute object
!!!   solv :: 3D-RISM solvent object
!!!   ucenter :: center the solute in the solvation box.  
!!!   ncuvsteps :: number of past cuv time steps saves
!!!   closure :: list of closures. Closures may be KH, HNC or PSEn
!!!              where n is an integer. Ensure the length attribute is
!!!              the same on all processes.
!!!   cut     :: distance cutoff for potential and force calculations
!!!   mdiis_nvec :: number of MDIIS vectors (previous iterations) to keep
!!!   mdiis_del :: scaling factor (step size) applied to estimated gradient (residual)
!!!   mdiis_method :: which implementation of the algorithm
!!!   o_buffer :: (optional) shortest distance between solute and solvent box boundary
!!!   o_grdspc :: (optional) linear grid spacing for the solvent box in each dimension
!!!   o_boxlen :: (optional) solvent box size in each dimension [A]
!!!   o_ng3    :: (optional) number of grid points in each dimension
!!!   o_mpicomm :: (optional) MPI communicator 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_new(this,solu,solv,ucenter, ncuvsteps,&
         closure,cut,mdiis_nvec, mdiis_del, mdiis_method, mdiis_restart, &
         o_buffer,o_grdspc,o_boxlen, o_ng3, o_mpicomm)
      use rism3d_solu_c
      use rism3d_solv_c
      implicit none
#ifdef MPI
      include 'mpif.h'
#endif /*MPI*/
      type(rism3d), intent(inout) :: this
      type(rism3d_solu), intent(in),target :: solu
      type(rism3d_solv), intent(in),target :: solv
      integer, intent(in) :: ucenter, ncuvsteps
      character(len=*),intent(in) :: closure(:)
      _REAL_, intent(in) :: cut
      integer, intent(in) :: mdiis_nvec,mdiis_method
      _REAL_, intent(in) :: mdiis_del, mdiis_restart
      _REAL_, optional, intent(in) :: o_buffer, o_grdspc(3)
      _REAL_, optional, intent(in) :: o_boxlen(3)
      integer, optional, intent(in) :: o_ng3(3)
      integer, optional, intent(in) :: o_mpicomm
      !temporary copies
      character(len=len(closure)),pointer :: t_closure(:)
      _REAL_ :: t_cut
      integer :: t_mdiis_nvec,t_mdiis_method
      _REAL_ :: t_mdiis_del, t_mdiis_restart
      _REAL_ :: t_buffer, t_grdspc(3)
      _REAL_ :: t_boxlen(3)
      integer :: t_ng3(3)
      integer :: t_mpicomm
      integer :: nclosure
      integer :: err
      
      !MPI set up starts by obtaining rank and size.  Temporary copies
      !of input parameters that do not directly set object variables
      !are made.  These are they broadcast to the rank > 0 processes.
      !Then all processes complete the intitialization proceedure
      !using the temporary copies.  This leaves the input parameters
      !untouched.

      call rism_timer_new(this%timer, "3D-RISM")
      call rism_timer_start(this%timer)
      call rism_timer_new(this%thermoTimer, "Thermodynamics")
      call rism_timer_setParent(this%thermoTimer,this%timer)
      call rism_timer_new(this%forceTimer, "Force")
      call rism_timer_setParent(this%forceTimer,this%thermoTimer)
      call rism_timer_new(this%exchemTimer, "Excess Chemical Potential")
      call rism_timer_setParent(this%exchemTimer,this%thermoTimer)
      call rism_timer_new(this%solveTimer, "Solve 3D-RISM")
      call rism_timer_setParent(this%solveTimer,this%timer)
      call rism_timer_new(this%resizeTimer, "Solvation box resize")
      call rism_timer_setParent(this%resizeTimer,this%solveTimer)
      call rism_timer_new(this%reorientTimer, "Solute reorientation")
      call rism_timer_setParent(this%reorientTimer,this%solveTimer)
      call rism_timer_new(this%cuvpropTimer, "Cuv propagation")
      call rism_timer_setParent(this%cuvpropTimer,this%solveTimer)
      call rism_timer_new(this%rxrismTimer, "RXRISM")
      call rism_timer_setParent(this%rxrismTimer,this%solveTimer)
      call rism_timer_new(this%r1rismTimer, "R1RISM")
      call rism_timer_setParent(this%r1rismTimer,this%rxrismTimer)
      call rism_timer_new(this%fftTimer, "FFT")
      call rism_timer_setParent(this%fftTimer,this%r1rismTimer)

      ! GET RANK AND SIZE
      this%mpicomm = 0
      this%mpisize = 1
      this%mpirank = 0
#ifdef MPI
      if(present(o_mpicomm)) then
         this%mpicomm = o_mpicomm
         if(this%mpicomm == MPI_COMM_NULL)&
              call rism_report_error("RISM3D: received NULL MPI communicator")
         call mpi_comm_rank(this%mpicomm,this%mpirank,err)
         if(err /=0) call rism_report_error&
              ("(a,i8)","RISM3D: could not get MPI rank for communicator ",this%mpicomm)
         call mpi_comm_size(this%mpicomm,this%mpisize,err)
         if(err /=0) call rism_report_error&
              ("(a,i8)","RISM3D: could not get MPI size for communicator ",this%mpicomm)
         call rism_report_mpi(this%mpicomm)
      end if
#endif /*MPI*/
      !MAKE TEMPORARY COPIES
      if(this%mpirank == 0)then
         call rism3d_solu_clone(solu,this%solu)
         call rism3d_solv_clone(solv,this%solv)
         this%ucenter = ucenter
         this%ncuvsteps = ncuvsteps
         nclosure = size(closure)
         t_closure=>safemem_realloc(t_closure,len(t_closure),nclosure)
         t_closure = closure
         t_cut=cut
         t_mdiis_nvec = mdiis_nvec
         t_mdiis_method = mdiis_method
         t_mdiis_del = mdiis_del
         t_mdiis_restart = mdiis_restart
         !check box parameters
         if(present(o_buffer) .and. present(o_grdspc))then
            t_buffer = o_buffer
            t_grdspc = o_grdspc
            if(present(o_boxlen) .or. present(o_ng3))&
                 call rism_report_error("RISM3D: do not set BOXLEN or NG3 for variable box size")
         elseif(present(o_boxlen) .and. present(o_ng3)) then
            t_boxlen = o_boxlen
            t_ng3 = o_ng3
            if(present(o_buffer) .or. present(o_grdspc)) &
                 call rism_report_error("RISM3D: do not set BUFFER or GRDSPC for fixed box size")
         else
            call rism_report_error("RISM3D: not enough parameters for fixed or variable box size")
         end if
      end if
#ifdef MPI
      !BROADCAST PARAMETERS
      !set solu on all processes
      call rism3d_solu_mpi_clone(this%solu,this%mpirank,this%mpicomm)
      !set solv on all processes
      call rism3d_solv_mpi_clone(this%solv,this%mpirank,this%mpicomm)
      !set ucenter on all processes
      call mpi_bcast(this%ucenter, 1, mpi_integer, 0, this%mpicomm,err)
      if(err /=0) call rism_report_error("RISM3D: broadcast UCENTER in constructor failed")
      !set ncuvstpes on all processes
      call mpi_bcast(this%ncuvsteps, 1, mpi_integer, 0, this%mpicomm,err)
      if(err /=0) call rism_report_error("RISM3D: broadcast NCUVSTEPS in constructor failed")
      call mpi_bcast(nclosure,1,mpi_integer,0,this%mpicomm,err)
      if(err /=0) call rism_report_error&
           ("RISM3D interface: could not broadcast PROGRESS")
      if(this%mpirank/=0)&
           t_closure=>safemem_realloc(t_closure,len(t_closure),nclosure)
      call mpi_bcast(t_closure, len(t_closure)*nclosure, mpi_character, 0, this%mpicomm,err)
      if(err /=0) call rism_report_error("RISM3D: broadcast CLOSURE in constructor failed")
      call mpi_bcast(t_cut, 1, mpi_double_precision, 0, this%mpicomm,err)
      if(err /=0) call rism_report_error("RISM3D: broadcast CUT in constructor failed")
      call mpi_bcast(t_mdiis_nvec, 1, mpi_integer, 0, this%mpicomm,err)
      if(err /=0) call rism_report_error("RISM3D: broadcast MDIIS_NVEC in constructor failed")
      call mpi_bcast(t_mdiis_del, 1, mpi_double_precision, 0, this%mpicomm,err)
      if(err /=0) call rism_report_error("RISM3D: broadcast MDIIS_DEL in constructor failed")
      call mpi_bcast(t_mdiis_restart, 1, mpi_double_precision, 0, this%mpicomm,err)
      if(err /=0) call rism_report_error("RISM3D: broadcast MDIIS_RESTART in constructor failed")
      call mpi_bcast(t_mdiis_method, 1, mpi_integer, 0, this%mpicomm,err)
      if(err /=0) call rism_report_error("RISM3D: broadcast MDIIS_METHOD in constructor failed")
      if(present(o_buffer) .and. present(o_grdspc))then
         call mpi_bcast(t_buffer, 1, mpi_double_precision, 0, this%mpicomm,err)
         if(err /=0) call rism_report_error("RISM3D: broadcast BUFFER in constructor failed")
         call mpi_bcast(t_grdspc, 3, mpi_double_precision, 0, this%mpicomm,err)
         if(err /=0) call rism_report_error("RISM3D: broadcast GRDSPC in constructor failed")
      else
         call mpi_bcast(t_boxlen, 3, mpi_double_precision, 0, this%mpicomm,err)
         if(err /=0) call rism_report_error("RISM3D: broadcast BOXLEN in constructor failed")
         call mpi_bcast(t_ng3, 3, mpi_integer, 0, this%mpicomm,err)
         if(err /=0) call rism_report_error("RISM3D: broadcast NG3 in constructor failed")
      end if
#endif /*MPI*/
      !INITIALIZE
      call rism3d_grid_new(this%grid,this%mpicomm)
      call rism_timer_stop(this%timer)
      call rism3d_setmdiis(this,t_mdiis_nvec,t_mdiis_del,t_mdiis_method,t_mdiis_restart)
      call rism_timer_start(this%timer)
      call rism3d_potential_new(this%pot, this%grid, this%solv, this%solu,0d0)
      call rism3d_potential_setTimerParent(this%pot,this%solveTimer)

#ifdef MPI
      call mdiis_new_mpi(this%mdiis_o,this%mdiis_method, &
           this%deloz,0d0,&
           this%MDIIS_restart,&
           this%mpirank, this%mpisize, this%mpicomm)
#else
      call mdiis_new(this%mdiis_o,this%mdiis_method, &
           this%deloz,0d0,&
           this%MDIIS_restart)
#endif /*MPI*/

      call mdiis_setTimerParent(this%mdiis_o,this%r1rismtimer)
      call rism_timer_stop(this%timer)
      call rism3d_setcut(this,t_cut)
      call rism3d_setclosurelist(this,t_closure)
      if(present(o_buffer) .and. present(o_grdspc))then
         call rism3d_setbox(this,t_buffer,t_grdspc)
      else
         call rism3d_setbox(this,t_boxlen,t_ng3)
      end if
      allocate(this%nsolutionChg, this%nsolutionNoChg)
      this%nsolutionChg=0
      this%nsolutionNoChg=0
      this%nsolution=>this%nsolutionChg

      call rism3d_fft_global_init()
      this%fftw_planner=FFT_MEASURE
#if defined(MPI) || defined(MKLdisabled)
      this%fft_aligned=.false.
#else
      this%fft_aligned=.true.
#endif
      this%fftw_localtrans=.true.

#ifdef RISM3D_DEBUG
      call rism3d_debug_new(this%grid,this%solv,this%mpirank,this%mpisize,this%mpicomm)
#endif

    end subroutine rism3d_new

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Set parent for this timer
!!!IN:
!!!   this : rism3d object
!!!   parent : parent timer object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism3d_setTimerParent(this, parent)
    implicit none
    type(rism3d), intent(inout) :: this
    type(rism_timer), intent(inout) :: parent
    call rism_timer_start(this%timer)
    call rism_timer_setParent(this%timer,parent)
    call rism_timer_stop(this%timer)
  end subroutine rism3d_setTimerParent

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Sets the parameters for a variable solvation box.
!!!IN:
!!!  this :: rism3d object
!!!  buffer :: shortest distance between solute and solvent box boundary
!!!  grdspc :: linear grid spacing for the solvent box in each dimension
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_setbox_variable(this,buffer,grdspc)
      implicit none
      type(rism3d),intent(inout) :: this
      _REAL_, intent(in) :: buffer, grdspc(3)
      call rism_timer_start(this%timer)
      this%varbox = .true.
      this%buffer = buffer
      call rism3d_grid_setSpacing(this%grid,grdspc)
      call rism_timer_stop(this%timer)
    end subroutine rism3d_setbox_variable

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Sets the parameters for a variable solvation box.
!!!IN:
!!!  this :: rism3d object
!!!  boxlen :: solvent box size in each dimension in Angstroms
!!!  ng3    :: number of grid points in each dimension
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_setbox_fixed(this,boxlen,ng3)
      implicit none
      type(rism3d),intent(inout) :: this
      _REAL_, intent(in) :: boxlen(3)
      integer, intent(in) :: ng3(3)
      call rism_timer_start(this%timer)
      this%varbox = .false.
      this%boxfixlen = boxlen
      this%nboxfix = ng3
      call rism_timer_stop(this%timer)
    end subroutine rism3d_setbox_fixed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Sets the closure list and sets the current closure to the first one
!!!in the list.  When there is no previous solution to work from, the
!!!solver will use each closure in the list in turn. By choosing the
!!!list to increase in order, it makes it possible to converge
!!!otherwise difficult closures. Only the last closure is used for
!!!thermodynamic output.
!!!IN:
!!!   this :: rism3d object
!!!   closure :: array of closure types (see closure enumeration).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_setclosurelist(this,closure)
      implicit none
      type(rism3d),intent(inout) :: this
      character(len=*), intent(in) :: closure(:)
      call rism_timer_start(this%timer)
      this%closureList => safemem_realloc(this%closureList,len(this%closureList),&
           ubound(closure,1))
      this%closureList = closure
      call rism_timer_stop(this%timer)
      call rism3d_setclosure(this,this%closureList(1))
    end subroutine rism3d_setclosurelist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Sets the closure type
!!!IN:
!!!   this :: rism3d object
!!!   closure :: closure type (see closure enumeration).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_setclosure(this,closure)
      implicit none
      type(rism3d),intent(inout) :: this
      character(len=*), intent(in) :: closure
      call rism_timer_start(this%timer)
      call rism3d_closure_destroy(this%closure)
      call rism3d_closure_new(this%closure,closure,this%pot)
      call rism_timer_stop(this%timer)
    end subroutine rism3d_setclosure

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Sets verbosity of output
!!!IN:
!!!   this :: rism3d object
!!!   verbosity :: 0 - no output
!!!                1 - memory allocation and steps for convergence
!!!                2 - 1 + convergence progress
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_setverbosity(this,verbosity)
      implicit none
      type(rism3d),intent(inout) :: this
      integer, intent(in) :: verbosity
      call rism_timer_start(this%timer)
      this%verbose = verbosity
      call rism_timer_stop(this%timer)
    end subroutine rism3d_setverbosity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Sets the cut off distance for potential and force calculations
!!!IN:
!!!   this :: rism3d object
!!!   cut     :: distance cutoff for potential and force calculations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_setcut(this,cut)
      implicit none
      type(rism3d),intent(inout) :: this
      _REAL_, intent(in) :: cut
      call rism_timer_start(this%timer)
      call rism3d_potential_setCut(this%pot,cut)
      call rism_timer_stop(this%timer)
    end subroutine rism3d_setcut

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Sets MDIIS parameters
!!!IN:
!!!   this :: rism3d object!
!!!   nvec :: number of MDIIS vectors (previous iterations) to keep
!!!   del :: scaling factor (step size) applied to estimated gradient (residual)
!!!   method :: which implementation of the algorithm
!!!   restart :: restart threshold factor. Ratio of the current residual to the 
!!!              minimum residual in the basis that causes a restart
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_setmdiis(this,nvec,del,method,restart)
      implicit none
      type(rism3d),intent(inout) :: this
      integer, intent(in) :: nvec,method
      _REAL_, intent(in) :: del, restart
      call rism_timer_start(this%timer)
      this%NVec = nvec
      this%deloz = del
      this%mdiis_method = method
      this%mdiis_restart = restart
      call rism_timer_stop(this%timer)
    end subroutine rism3d_setmdiis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Sets solute coordinates.
!!!IN:
!!!   this :: rism3d object
!!!   ratu :: coordinates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_setCoord(this,ratu)
      implicit none
      type(rism3d), intent(inout) :: this
      _REAL_, intent(in) :: ratu(:,:)
      call rism_timer_start(this%timer)
      call rism3d_solu_setCoord(this%solu,ratu)
      call rism_timer_stop(this%timer)
    end subroutine rism3d_setCoord

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Sets all solute partial charges to zero, resets MDIIS and wipes out
!!!working memory.
!!!IN:
!!!   this :: rism3d object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_unsetCharges(this)
      implicit none
      type(rism3d), intent(inout) :: this
      integer :: i
      !reset MDIIS.  This makes the working vector index 1
      call mdiis_reset(this%mdiis_o)
      this%cuv=>this%cuvWRK(:,:,:,:,mdiis_getWorkVector(this%mdiis_o))
      this%cuvres=>this%cuvresWRK(:,:,mdiis_getWorkVector(this%mdiis_o))
      !turn off charges
      call rism3d_solu_unsetCharges(this%solu)
      !Use the number of no charge solutions
      this%nsolution=>this%nsolutionNoChg
      !Use no charge previous soluitions
      this%oldcuv=> this%oldcuvNoChg
      !If we have run with no charges before, copy the previous solution
      if(associated(this%oldcuvNoChg))&
           call dcopy(product(ubound(this%oldcuv)),this%oldcuv,1,this%cuv,1)
    end subroutine rism3d_unsetCharges

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Sets all solute partial charges to to their original
!!!values. (Undoes rism3d_unsetCharge().)
!!!IN:
!!!   this :: rism3d object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_resetCharges(this)
      implicit none
      type(rism3d), intent(inout) :: this
      integer :: i
      !reset MDIIS.  This makes the working vector index 1
      call mdiis_reset(this%mdiis_o)
      this%cuv=>this%cuvWRK(:,:,:,:,mdiis_getWorkVector(this%mdiis_o))
      this%cuvres=>this%cuvresWRK(:,:,mdiis_getWorkVector(this%mdiis_o))
      !get back the charges
      call rism3d_solu_resetCharges(this%solu)
      !restore the number of previous solutions
      this%nsolution=>this%nsolutionChg
      !point to previous charged solutions
      this%oldcuv=> this%oldcuvChg
      !If we have run with charges before, copy the previous solution
      if(associated(this%oldcuvNoChg))&
           call dcopy(product(ubound(this%oldcuv)),this%oldcuv,1,this%cuv,1)
    end subroutine rism3d_resetCharges

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the full 3D-RISM solvent distribution.  This is required to 
!!!calculate thermodynamic quantities
!!!IN:
!!!   this :: rism3d object
!!!   ksave :: save itermediate results every ksave interations (0
!!!            means no saves)
!!!   kshow :: print parameter for relaxation steps every kshow
!!!            iteration (0 means no saves)
!!!   maxste :: maximum number of rism relaxation steps
!!!   tol :: convergence tolerance.  There should be one tolerance for
!!!          each closure in the closure list. If there is a mismatch
!!!          then the program dies with an error.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_solve (this,ksave,kshow,maxste,tol)
      implicit none
#if defined(MPI)
      include 'mpif.h'
#endif /*defined(MPI)*/
      type(rism3d), intent(inout) :: this
      integer, intent(in) :: ksave,kshow,maxste
      _REAL_, intent(in) :: tol(:)
      !iclosure :: counter for closures
      integer :: iclosure
      call rism_timer_start(this%solveTimer)

!!!!
!! 0) quick check that the tolerance list is of the correct length
!!!!
      if(ubound(tol,1) /= ubound(this%closureList,1))&
           call rism_report_error("(a,i3,a,i3)",&
           "RISM3D_SOLVE: number of tolerances, ",&
           ubound(tol,1),", is not equal to numer of closures, ",&
           ubound(this%closureList,1))
!!!!
!! 1) Reorient solute along the principal axis and resize the grids if necessary
!!!!
#ifdef RISM_DEBUG
      write(0,*) "(1)"
      write(0,*) "RATU", this%mpirank, this%solu%natom, this%grid%nr, this%grid%nga
      call flush(0)
#endif /*RISM_DEBUG*/

!!!!
!! 1b) get the minimum box size for this frame
!!!!
!      if(this%varbox .or. .not. associated(this%grid%ga))then
      if(this%varbox .or. this%nsolution ==0)then
         call rism_timer_start(this%resizeTimer)
         call timer_start(TIME_RESIZE)
         call resizebox(this)
         call timer_stop(TIME_RESIZE)
         call rism_timer_stop(this%resizeTimer)
      end if
      call timer_start(TIME_REORIENT)
      call rism_timer_start(this%reorientTimer)
      call center_solute(this)
      call rism_timer_stop(this%reorientTimer)
      call timer_stop(TIME_REORIENT)

!!!!
!! 3) Calculate electrostatic and Lennard-Jones potential about the solute
!!!!
#ifdef RISM_DEBUG
      write(0,*) "(3)"
      call flush(0)
#endif /*RISM_DEBUG*/
      call rism3d_potential_calc(this%pot)

!!!!
!! 4) Calculate long range asymptotpics
!!!!
#ifdef RISM_DEBUG
      write(0,*) "(4)"
      call flush(0)
#endif /*RISM_DEBUG*/
      call timer_start(TIME_ASYMP)
      !.....tabulating asymptotic function of Cuv and Huv.....
      call rism3d_potential_asympch(this%pot)
      call timer_stop(TIME_ASYMP)

!!!!
!! 2) Propagate previously saved CUV solutions to create an initial guess for this solution
!!!!
#ifdef RISM_DEBUG
      write(0,*) "(2)", this%grid%nkTotal, this%grid%nga
      call flush(0)
#endif /*RISM_DEBUG*/
      call timer_start(TIME_CUVPROP)
      call rism_timer_start(this%cuvpropTimer)
      call cuv_propagate(this)
      call rism_timer_stop(this%cuvpropTimer)
      call timer_stop(TIME_CUVPROP)

!!!!
!! 5) Calculate 3D-RISM solution
!!!!
#ifdef RISM_DEBUG
      write(0,*) "(5) RXRISM"
      call flush(0)
#endif /*RISM_DEBUG*/
      !................. solving 3D-RISM/HNC(PLHNC) by MDIIS .................
      
     !if the user has provide a list of closure, use it only if this
     !is the first solution (nsolution==0) or solution propagation is
     !turned off (ncuvsteps==0). Otherwise, the current closure will
     !be the last one in the list.
     if(this%nsolution == 0 .or. this%ncuvsteps==0)then
        do iclosure =1, size(this%closureList)
           if(this%verbose >=1)&
                call rism_report_message("|Switching to "//&
                trim(this%closureList(iclosure))//" closure")
           call rism_timer_stop(this%solveTimer)
           call rism3d_setClosure(this,this%closureList(iclosure))
           call rism_timer_start(this%solveTimer)
           call timer_start(TIME_RXRISM)
           call rxrism(this,ksave,kshow,maxste,tol(iclosure))
           call timer_stop(TIME_RXRISM)
           !increment nsolution and ncuvsteps to ensure the previous
           !closure solution is used
           if(iclosure ==1)then
              this%nsolution=this%nsolution+1
              this%ncuvsteps=this%ncuvsteps+1
           end if
        end do
        this%nsolution=this%nsolution-1
        this%ncuvsteps=this%ncuvsteps-1
     else
        call timer_start(TIME_RXRISM)
        call rxrism(this,ksave,kshow,maxste,tol(size(tol)))
        call timer_stop(TIME_RXRISM)
     end if

!!!!
!! 11) Update stored variables
!!!!
#ifdef RISM_DEBUG
      write(0,*) "(11)", ubound(this%oldcuv)
      call flush(0)
#endif /*RISM_DEBUG*/
      call timer_start(TIME_CUVPROP)
      call rism_timer_start(this%cuvpropTimer)
      this%nsolution = this%nsolution+1
      call cuv_update(this)
      call rism_timer_stop(this%cuvpropTimer)
      call timer_stop(TIME_CUVPROP)

#ifdef RISM_DEBUG
      write(0,*) "DONE SOLVE"
      call flush(0)
#endif /*RISM_DEBUG*/
      call rism_timer_stop(this%solveTimer)
    end subroutine rism3d_solve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!THERMODYNAMICS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!calculates the forces on the solute contributed by the solvent according
!!!to 3D-RISM.  In fact, this subroutine calls the appropriate subroutines
!!!to calculate this
!!!IN:
!!!   this :: rism3d object with computed solution
!!!   ff   :: 3D-RISM forces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism3d_force (this,ff)
    implicit none
    type(rism3d):: this
    _REAL_,intent(out) :: ff(3,this%solu%natom)
    !................... Calculation force fields    .......................
    call rism_timer_start(this%forceTimer)
    call timer_start(TIME_FF)
    call rism3d_closure_force(this%closure,ff,this%guv)
    call timer_stop(TIME_FF)
    call rism_timer_stop(this%forceTimer)
  end subroutine rism3d_force

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculate the excess chemical potential of solvation for each solvent species
!!!IN:
!!!   this :: rism3d object with computed solution
!!!   o_lr   :: (optional) (default=.true.) Apply asymptotic long range correction
!!!OUT:
!!!    excess chemical potential of solvation for each solvent species
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism3d_exchem (this,o_lr) result (exchem)
      implicit none
      type(rism3d), intent(inout) :: this
      logical, optional, intent(in) :: o_lr
      logical :: lr
      _REAL_ :: exchem(this%solv%natom)
      call rism_timer_start(this%exchemTimer)

      lr = .true.
      if(present(o_lr)) lr = o_lr
      
      if(lr)then
         exchem = rism3d_closure_aexchem(this%closure,this%huv,this%cuv(:,:,:,:))
      else
         exchem = rism3d_closure_exchem(this%closure,this%huv,this%cuv(:,:,:,:))
      end if
      call rism_timer_stop(this%exchemTimer)
    end function rism3d_exchem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculate the total excess chemical potential of solvation 
!!!IN:
!!!   this :: rism3d object with computed solution
!!!   o_lr   :: (optional) (default=.true.) Apply asymptotic long range correction
!!!OUT:
!!!    total excess chemical potential of solvation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism3d_exchem_tot (this,o_lr) result (exchem)
      implicit none
      type(rism3d), intent(inout) :: this
      logical, optional, intent(in) :: o_lr
      logical :: lr
      _REAL_ :: exchem
      call rism_timer_start(this%exchemTimer)

      lr = .true.
      if(present(o_lr)) lr = o_lr

      call rism_timer_stop(this%exchemTimer)
      exchem = sum(rism3d_exchem(this,o_lr))
    end function rism3d_exchem_tot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculate the excess chemical potential of solvation for each solvent species
!!!with the Gaussian fluctuation correction
!!!IN:
!!!   this :: rism3d object with computed solution
!!!   o_lr   :: (optional) (default=.true.) Apply asymptotic long range correction
!!!OUT:
!!!    Gaussian fluctuation excess chemical potential of solvation for each 
!!!    solvent species
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism3d_exchemGF (this,o_lr) result (exchem)
      implicit none
      type(rism3d), intent(inout) :: this
      logical, optional, intent(in) :: o_lr
      logical :: lr
      _REAL_ :: exchem(this%solv%natom)
      call rism_timer_start(this%exchemTimer)

      lr = .true.
      if(present(o_lr)) lr = o_lr
      if(lr)then
         exchem = rism3d_closure_aexchemGF(this%closure,this%huv,this%cuv(:,:,:,:))
      else
         exchem = rism3d_closure_exchemGF(this%closure,this%huv,this%cuv(:,:,:,:))
      end if
      call rism_timer_stop(this%exchemTimer)
    end function rism3d_exchemGF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculate the total excess chemical potential of solvation with the Gaussian 
!!!fluctuation correction
!!!IN:
!!!   this :: rism3d object with computed solution
!!!   o_lr   :: (optional) (default=.true.) Apply asymptotic long range correction
!!!OUT:
!!!    total Gaussian fluctuation excess chemical potential of solvation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism3d_exchemGF_tot (this,o_lr) result (exchem)
      implicit none
      type(rism3d), intent(inout) :: this
      logical, optional, intent(in) :: o_lr
      logical :: lr
      _REAL_ :: exchem
      call rism_timer_start(this%exchemTimer)

      lr = .true.
      if(present(o_lr)) lr = o_lr
      
      call rism_timer_stop(this%exchemTimer)
      exchem = sum(rism3d_exchemGF(this,o_lr))
    end function rism3d_exchemGF_tot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculate the solvation interaction energy: de = rho sum g*u for
!!!each solvent site.  I.e., the direct intection potential energy of
!!!solute and solvent and not the total solvation energy (see solvEne).
!!!IN:
!!!   this :: rism3d object with computed solution
!!!OUT:
!!!    the contribution of each solvent site to the total solvation interaction 
!!!    energy [kT]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism3d_solvPotEne (this) result (ene)
      implicit none
      type(rism3d), intent(inout) :: this
      _REAL_ :: ene(this%solv%natom)
      call rism_timer_start(this%thermoTimer)
      ene = rism3d_closure_solvPotEne(this%closure,this%guv)
      call rism_timer_stop(this%thermoTimer)
    end function rism3d_solvPotEne

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculate the total solvation interaction energy: de = rho sum g*u for
!!!each solvent site.  I.e., the direct intection potential energy of
!!!solute and solvent and not the total solvation energy (see solvEne).
!!!IN:
!!!   this :: rism3d object with computed solution
!!!OUT:
!!!    the total solvent-solute potential energy [kT]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism3d_solvPotEne_tot (this) result (ene)
      implicit none
      type(rism3d), intent(inout) :: this
      _REAL_ :: ene
      call rism_timer_start(this%thermoTimer)
      ene = sum(rism3d_solvPotEne(this))
      call rism_timer_stop(this%thermoTimer)
    end function rism3d_solvPotEne_tot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculating PMV of solute
!!!IN:
!!!   this :: rism3d object with computed solution
!!!OUT:
!!!   partial molar volume
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function  rism3d_pmv (this) result(pmv)
      implicit none
      type(rism3d), intent(inout) :: this
      _REAL_ :: pmv
      call rism_timer_start(this%thermoTimer)
      pmv = rism3d_closure_pmv(this%closure,this%cuv(:,:,:,:))
      call rism_timer_stop(this%thermoTimer)
    end function rism3d_pmv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculating excess number of each solvent type associated with the solute
!!!IN:
!!!   this :: rism3d object with computed solution
!!!   o_lr   :: (optional) (default=.true.) Apply asymptotic long range correction
!!!OUT:
!!!    excess number of each solvent type associated with the solute
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function  rism3d_exNum (this,o_lr) result(num)
      implicit none
      type(rism3d), intent(inout) :: this
      logical, optional, intent(in) :: o_lr
      logical :: lr
      _REAL_ :: num(this%solv%natom)
      call rism_timer_start(this%thermoTimer)

      lr = .true.
      if(present(o_lr)) lr = o_lr
      
      if(lr)then
         num = rism3d_closure_aexNum(this%closure,this%guv)
      else
         num = rism3d_closure_exNum(this%closure,this%guv)
      end if
      call rism_timer_stop(this%thermoTimer)
    end function rism3d_exNum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!DEALLOCATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! deconstructor - frees all memory
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_destroy(this)
      use safemem
      implicit none
      type(rism3d) :: this
      call rism_timer_destroy(this%fftTimer)
      call rism_timer_destroy(this%r1rismTimer)
      call rism_timer_destroy(this%rxrismTimer)
      call rism_timer_destroy(this%cuvpropTimer)
      call rism_timer_destroy(this%reorientTimer)
      call rism_timer_destroy(this%resizeTimer)
      call rism_timer_destroy(this%solveTimer)
      call rism_timer_destroy(this%exchemTimer)
      call rism_timer_destroy(this%forceTimer)
      call rism_timer_destroy(this%thermoTimer)
      call rism_timer_destroy(this%timer)

      call rism3d_solv_destroy(this%solv)
      call rism3d_solu_destroy(this%solu)
      call rism3d_potential_destroy(this%pot)
      call rism3d_grid_destroy(this%grid)
      call rism3d_closure_destroy(this%closure)
      call mdiis_destroy(this%mdiis_o)

      if(safemem_dealloc(this%xvva)/=0) &
           call rism_report_error("RISM3D: failed to deallocate XVVA")
      if(safemem_dealloc(this%cuvWRK)/=0) &
           call rism_report_error("RISM3D: failed to deallocate CUVWRK")
      if(safemem_dealloc(this%oldcuvChg)/=0) &
           call rism_report_error("RISM3D: failed to deallocate OLDCUVCHG")
      if(safemem_dealloc(this%oldcuvNoChg)/=0) &
           call rism_report_error("RISM3D: failed to deallocate OLDCUVNOCHG")
      if(safemem_dealloc(this%cuvresWRK)/=0) &
           call rism_report_error("RISM3D: failed to deallocate CUVRESWRK")
      if(safemem_dealloc(this%closureList)/=0) &
           call rism_report_error("RISM3D: failed to deallocate CLOSURELIST")
      if(associated(this%nsolutionChg))&
           deallocate(this%nsolutionChg)
      if(associated(this%nsolutionNoChg))&
           deallocate(this%nsolutionNoChg)
      nullify(this%cuv)
      nullify(this%cuvres)
      nullify(this%oldcuv)
      nullify(this%nsolution)

      if(safemem_dealloc(this%guv,o_aligned=.true.)/=0) &
           call rism_report_error("RISM3D: failed to deallocate GUV")
      if(safemem_dealloc(this%huv,o_aligned=.true.)/=0) &
           call rism_report_error("RISM3D: failed to deallocate HUV")
      call rism3d_fft_destroy(this%fft)
      call rism3d_fft_global_finalize()
    end subroutine rism3d_destroy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!using the current orientation of the solute, define the minimum box size
!!!and resize all associated grids
!!!IN:
!!!   this :: rism3d object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine resizebox(this)
      use rism_util, only : isprime, lcm, isfactorable
      implicit none
#if defined(MPI)
      include 'mpif.h'
      integer :: ierr
#endif /*defined(MPI)*/
      type(rism3d), intent(inout) :: this
      integer :: ngr(3)
      integer :: id
      _REAL_ :: boxlen(3)
      integer :: primes(4)=(/2,3,5,7/)
#  ifdef RISM_DEBUG
      write(0,*) "RESIZEBOX"
      call flush(0)
#  endif /*RISM_DEBUG*/

      !to properly distribute the run y and z dimension must have a number of 
      !grid points that is a multiple of mpisize.  Thus, mpisize must be 
      !factorizable by small prime numbers

      if(.not.isfactorable(this%mpisize,primes))then
         call rism_report_error("(a,10i4)","Sorry, 3D-RISM requires that the number "&
              //"of processes be a product of ",primes)
      end if
      !if we have a fixed box size, we simply retain all of the
      !previously calculated box size values
      if(this%varbox) then
         !get minimum box size defined by the buffer
         do id = 1,3
            boxlen(id) = maxval(this%solu%ratu(id,:)) &
                 - minval(this%solu%ratu(id,:)) + 2*this%buffer
         end do

         !round this box size using the prescribed linear density and get the
         !number of grid points required
         ngr = ceiling(boxlen/this%grid%grdspc)
         boxlen = ngr*this%grid%grdspc

         !determine if the number of grid points in each dimension product
         !only has prime factors 2, 3, 5 or 7. If not increment the number of
         !points (in that dimension) until this is true

         !make sure that each dimension is divisable by 2 and that the
         !y- and z-dimensions are divisible by this%mpisize if this%mpisize > 1

         do id = 1,3
            ngr(id) = ngr(id)+mod(ngr(id),2)
         end do
         if(this%mpisize > 2)then
            do id=2,3
               if(mod(ngr(id),this%mpisize) /=0) then
                  ngr(id) = ngr(id)+lcm(this%mpisize,2)&
                       -mod(ngr(id),lcm(this%mpisize,2))
               end if
            end do
         end if

         do id = 1,3
            do while(.not.isfactorable(ngr(id),primes)) 
                  if(this%mpisize > 1 .and. id > 1)then
                     ngr(id) = ngr(id)&
                          +lcm(this%mpisize,2)
                  else
                     ngr(id) = ngr(id)+2
                  end if
            end do
         end do
         boxlen = ngr*this%grid%grdspc

#ifndef MPI
         if(.not. associated(this%cuv) .or. &
              .not. associated(this%oldcuv) .or. &
              ngr(1) /= size(this%cuv,1) .or.&
              size(this%cuv,2) /= ngr(2) .or. size(this%cuv,3) /= ngr(3)) then
            call reallocbox(this,ngr,this%grid%grdspc)
         end if
#else
         if(.not. associated(this%cuv) .or. &
              .not. associated(this%oldcuv) .or. &
              ngr(1) /= size(this%cuv,1) &
              .or. size(this%cuv,2) /= ngr(2) &
              .or. size(this%cuv,3) /= ngr(3)/this%mpisize) then
            call reallocbox(this,ngr,this%grid%grdspc)
         end if
#endif /*defined(MPI)*/
      else !fixed box size
         do id =2,3
            if(this%mpirank == 0 .and. (mod(this%nboxfix(id),this%mpisize) /=0 &
                 .or. mod(this%nboxfix(id),2) /=0) )then
               call rism_report_error("Sorry, MPI 3D-RISM requires that fixed grid sizes be "//&
                "divisible by two and the number of processes in the "//&
                "y and z dimensions")
            end if
         end do
         call reallocbox(this,this%nboxfix,this%boxfixlen/this%nboxfix)
      end if
    end subroutine resizebox

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!using the current box size and resize all associated grids and variables
!!!IN:
!!!   this :: rism3d object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine reallocbox(this,ngr,grdspc)
  use rism3d_fft_c
  use safemem
  implicit none
  type(rism3d) :: this
  integer, intent(in) :: ngr(3)
  _REAL_,intent(in) :: grdspc(3)
  integer :: i, id, irank
  _REAL_ :: memuse
#ifdef RISM_DEBUG
  write(0,*) "reallocbox"
  call flush(0)
#endif /*RISM_DEBUG*/

  !MUST CHECK THIS FOR MPI AND NON_MPI CASE  

  !All grid sizes have been determined so report details about the grid at this point
  if(this%verbose >=1)then
     call rism_report_message("||Setting solvation box to")
     call rism_report_message("(3(a,i10))", "|grid size: ",&
          ngr(1)," X ", ngr(2), " X ", ngr(3))
     call rism_report_message("(3(a,f10.3))", "|box size [A]:  ",&
          ngr(1)*grdspc(1)," X ", ngr(2)*grdspc(2), " X ", ngr(3)*grdspc(3))
     call rism_report_message("(3(a,f10.3))","|Effective buffer [A]:",&
          (ngr(1)*grdspc(1)-(maxval(this%solu%ratu(1,:))&
          - minval(this%solu%ratu(1,:))))/2d0,",  ",&
          (ngr(2)*grdspc(2)-(maxval(this%solu%ratu(2,:))&
          - minval(this%solu%ratu(2,:))))/2d0,",  ",&
          (ngr(3)*grdspc(3)-(maxval(this%solu%ratu(3,:))&
          - minval(this%solu%ratu(3,:))))/2d0)
!!$     memuse = 8d0*(dble(product(ngr))*this%solv%natom*(2*this%NVec+1+this%ncuvsteps)&
!!$          +(dble(product(ngr))+dble(product(ngr(1:2))))*(4+2*this%solv%natom)& !see manual. Does not include FFTW scratch
!!$          +product(ubound(this%solv%xvv))*this%mpisize & !XVV
!!$          +product(ubound(this%solv%fourier_tbl))*this%mpisize & !fourier_table
!!$          +this%solu%natom*7*this%mpisize & !solute
!!$          +this%grid%nga*2*this%mpisize & !XVVA and grid%ga
!!$          +this%grid%nkTotal/2*4*this%mpisize)& !grid%gv and grid%g2
!!$          +4*(this%grid%nkTotal/2*this%mpisize) !grid%indga
!!$#ifdef MPI
!!$     memuse = memuse +8d0*this%grid%nktotal !FFTW scratch
!!$#else
!!$     memuse = memuse +8d0*product(ngr(2:3))*2 !FFTW scratch
!!$#endif /*MPI*/
!!$     call rism_report_message("(a,g12.4,a)","|Projected mimimum distributed memory use: ",&
!!$          memuse/1024d0**3," GB" )
!!$     memuse = memuse + 8d0*(dble(product(ngr))*this%solv%natom*this%ncuvsteps&
!!$          +this%solu%natom*7*this%mpisize)  !solute
!!$
!!$     call rism_report_message("(a,g12.4,a)","|with polar/apolar decomposition:          ",&
!!$          memuse/1024d0**3," GB" )
     call flush(rism_report_getmunit())
  end if

  call rism3d_fft_setgrid(this%grid,ngr,grdspc,this%solv%natom,this%fft_aligned)

  !
  !2) allocation
  !
  !reallocate arrays that require preservation of their contents

  !I THINK WE CAN GET RID OF THE MPI HERE

#if defined(MPI)
  this%cuvWRK => safemem_realloc(this%cuvWRK,this%grid%nr(1),this%grid%nr(2),this%grid%nr(3),&
       this%solv%natom,this%NVec,.true.,.true.)
  if(rism3d_solu_charged(this%solu))then
     this%oldcuvChg => safemem_realloc(this%oldcuvChg,this%grid%nr(1),this%grid%nr(2),&
          this%grid%nr(3),this%solv%natom,this%ncuvsteps,.true.,.true.)
     this%oldcuv => this%oldcuvChg
  else
     this%oldcuvNoChg => safemem_realloc(this%oldcuvNoChg,this%grid%nr(1),this%grid%nr(2),&
          this%grid%nr(3),this%solv%natom,this%ncuvsteps,.true.,.true.)
     this%oldcuv => this%oldcuvNoChg
  end if
#else
  this%cuvWRK => safemem_realloc(this%cuvWRK,this%grid%nr(1),this%grid%nr(2),this%grid%nr(3),&
       this%solv%natom,this%NVec,.true.,.true.)
  if(rism3d_solu_charged(this%solu))then
     this%oldcuvChg => safemem_realloc(this%oldcuvChg,this%grid%nr(1),this%grid%nr(2),&
          this%grid%nr(3),this%solv%natom,this%ncuvsteps,.true.,.true.)
     this%oldcuv => this%oldcuvChg
  else
     this%oldcuvNoChg => safemem_realloc(this%oldcuvNoChg,this%grid%nr(1),this%grid%nr(2),&
          this%grid%nr(3),this%solv%natom,this%ncuvsteps,.true.,.true.)
     this%oldcuv => this%oldcuvNoChg
  end if
!!$  this%oldcuv => safemem_realloc(this%oldcuv,this%grid%nr(1),this%grid%nr(2),&
!!$       this%grid%nr(3),this%solv%natom,this%ncuvsteps,.true.)
#endif /*MPI*/
  !reallocate arrays that do not require preservation of their contents
  this%guv => safemem_realloc(this%guv,this%grid%nkTotal,this%solv%natom,&
       o_preserve=.false.,o_aligned=.true.)
  this%huv => safemem_realloc(this%huv,this%grid%nkTotal,this%solv%natom,&
        o_preserve=.false.,o_aligned=.true.)
  this%cuvresWRK => safemem_realloc(this%cuvresWRK,this%grid%nrTotal,this%solv%natom,&
       this%NVec,.false.)
  this%xvva => safemem_realloc(this%xvva,this%grid%nga* (this%solv%natom) **2,.false.)

  call rism3d_fft_destroy(this%fft)
  call rism3d_fft_new(this%fft,&
       this%fftw_planner, this%fftw_localtrans, this%fft_aligned,&
       this%grid,&
       this%guv,this%huv)

  !updated pointers
  call mdiis_resize(this%mdiis_o,this%cuvWRK, this%cuvresWRK,&
       this%grid%nrTotal*this%solv%natom, this%nvec)
  this%cuv=>this%cuvWRK(:,:,:,:,mdiis_getWorkVector(this%mdiis_o))
  this%cuvres=>this%cuvresWRK(:,:,mdiis_getWorkVector(this%mdiis_o))
  
  !
  !3) the remaining variables are handled by rism_setup_wavevector,interpolate_xvva
  !
  call interpolate_xvva(this,this%solv%xvv,this%xvva)
end subroutine reallocbox

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Prints the maximum amount of memory allocated at any one time so far in the
!!!run
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rism3d_max_memory(this)
  use safemem
  implicit none
#ifdef MPI
  include "mpif.h"
#endif /*MPI*/
  type(rism3d) :: this
  integer*8 :: memstats(10), tmemstats(10)
  integer :: err, irank, outunit
  outunit = rism_report_getMUnit()
  memstats = memStatus()
#ifdef MPI
#  ifdef USE_MPI_IN_PLACE
     if(this%mpirank==0)then
        call MPI_REDUCE(MPI_IN_PLACE, memstats, ubound(memstats,1),MPI_INTEGER8,&
             MPI_SUM,0,this%mpicomm,err)
     else
        call MPI_REDUCE(memstats, memstats, ubound(memstats,1),MPI_INTEGER8,&
             MPI_SUM,0,this%mpicomm,err)
     end if
#  else /*USE_MPI_IN_PLACE*/
     call MPI_REDUCE(memstats, tmemstats, ubound(memstats,1),MPI_INTEGER8,&
          MPI_SUM,0,this%mpicomm,err)
     memstats = tmemstats
#  endif /*USE_MPI_IN_PLACE*/
     if(err/=0) call rism_report_warn("RISM_MAX_MEMORY: MPI_REDUCE failed.")
#endif
  if(this%mpirank==0)then
     write(outunit,'(a)')
     write(outunit,'(a)') "|3D-RISM memory allocation summary"
     write(outunit,'(a)') "|Type          Maximum"
     write(outunit,'(a,f12.5,a)') "|Integer  ",&
          dble(memstats(1))/BYTES_PER_GB," GB"
     write(outunit,'(a,f12.5,a)') "|Real     ",&
          dble(memstats(2))/BYTES_PER_GB," GB"
     write(outunit,'(a,f12.5,a)') "|Logical  ",&
          dble(memstats(3))/BYTES_PER_GB," GB"
     write(outunit,'(a,f12.5,a)') "|Character",&
          dble(memstats(4))/BYTES_PER_GB," GB"
     write(outunit,'(a)') "|------------------------"
     write(outunit,'(a,f12.5,a)') "|Total    ",&
          dble(memstats(5))/BYTES_PER_GB," GB"
  end if
end subroutine rism3d_max_memory

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!interpolate the xvv array for the solvent that has be read in to construct
!!!this%xvva for the specific box size we are using
!!!IN:
!!!  this :: rism3d object
!!!  xvv  :: 1D-RISM Xvv
!!!  xvva :: interpolated result
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine interpolate_xvva (this,xvv,xvva)
      use rism_util, only : poly_interp
      implicit none
      type(rism3d),intent(inout) :: this
      _REAL_, intent(in) :: xvv(:,:,:)
      _REAL_, intent(out) :: xvva(:)
      integer :: iiga, itab, itab1,iv1,iv2
      _REAL_ :: err
      integer :: maxip
      parameter (maxip=5)

#ifdef RISM_DEBUG
      write(0,*) "INTERPOLATE_XVVA", sum(this%solv%fourier_tbl),sum(this%grid%ga), this%grid%nga; call flush(0)
#endif /*RISM_DEBUG*/

      !........................ checking R-grid size .........................
      if (this%grid%ga(this%grid%nga) > this%solv%fourier_tbl(this%solv%nr))  then
         call rism_report_error('(a,1pe16.8,a,1pe16.8)',&
              'DISTVV: bulk solvent Kmax=',this%solv%fourier_tbl(this%solv%nr),&
              'insufficient for 3D-grid Kmax=',this%grid%ga(this%grid%nga))
      elseif (maxip > this%solv%nr)  then
         call rism_report_error('(a,i7,a,i7)', &
              'DISTVV: bulk solvent grid size Nr=',this%solv%nr,&
              'insufficient for interpolation MaxIp=',maxip)
      endif

      !................. interpolating Xvv(K) at 3D-grid |K| .................
      do iiga=1,this%grid%nga
         do itab=1,this%solv%nr-maxip+1
            itab1 = itab
            if (this%solv%fourier_tbl(itab1+maxip/2) > this%grid%ga(iiga))  then
               !itab1 is now lower bracket for interpolation
               exit
            end if
         enddo
         do iv2=1,this%solv%natom
            do iv1=1,this%solv%natom
               call  poly_interp (this%solv%fourier_tbl(itab1:itab1+maxip), &
                    xvv(itab1:itab1+maxip,iv1,iv2),maxip, &
                    this%grid%ga(iiga), &
                    xvva(iiga +(iv1-1)*this%grid%nga + (iv2-1)*this%grid%nga*this%solv%natom), &
                    err)
            enddo
         enddo
      enddo
   end subroutine interpolate_xvva

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! center the solute in the solvent box
!!!IN:
!!!   this :: rism3d object    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine center_solute(this)
      use rism_util, only : calc_cm,translate
      implicit none
      type(rism3d), intent(inout) :: this
      _REAL_ :: weight(this%solu%natom)
      _REAL_ :: gridshift(3)
      integer :: idim
      !if ucenter is /= 0 we want to move the solute to the center of the solvent box
      !However, if ucenter < 0 we only figure out the displacement required the _first_
      !time we see the solute.  Thus, for ucenter <= 0, the solute's CM can move relative
      !to the grid
      if(mod(abs(this%ucenter),2)==1)then
         weight=this%solu%mass
      elseif(mod(abs(this%ucenter),2)==0)then
         weight=1
      end if
      if(this%ucenter > 0 .or. (this%ucenter < 0 .and. this%nsolution == 0))then
         call calc_cm(this%solu%ratu,this%ratucm,weight,this%solu%natom)
      end if
      !move the center to the nearest multiple of the grid spacing
      if(this%ucenter >=3 .and. this%ucenter <=4)then
         gridshift = mod(this%ratucm,this%grid%grdspc)
         do idim=1,3
            !round down for the smallest possible translation (i.e., the closest multiple)
            if(abs(gridshift(idim)) > this%grid%grdspc(idim)/2d0)then
               gridshift(idim) = gridshift(idim)-sign(this%grid%grdspc(idim),gridshift(idim))
            end if
         end do
         this%ratucm = this%ratucm - gridshift
      end if
      if(this%ucenter /= 0) then
         call translate(this%solu%ratu,this%solu%natom,this%grid%boxlen/2d0-1*this%ratucm)
      end if
    end subroutine center_solute

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!subroutines to find the iterative 3D-RISM solution.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Main driver for the 3D-RISM solver
!!!IN:
!!!   this :: rism3d object
!!!   ksave  :: save itermediate results every ksave interations (0 means no saves)
!!!   kshow  :: print parameter for relaxation steps every kshow iteration (0 means no saves)
!!!   maxste :: maximum number of rism relaxation steps
!!!   mdiis_method :: MDIIS implementation to use    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rxrism (this,ksave,kshow,maxste,tol)
      use mdiis_c
      use rism3d_csv
      implicit none
#include "def_time.h" 
#if defined(MPI)
      include 'mpif.h'
#endif /*defined(MPI)*/
      type(rism3d), intent(inout) :: this
      integer, intent(in) :: ksave,kshow,maxste
      _REAL_, intent(in) :: tol
      character(72) ::  cuvsav='rism.csv',guvfile
      integer :: guv_local=77
      !    character*(*)  cuvsav


      integer :: iatv, igx, igy, igz
      logical ::  found, conver=.false.
      integer ::  ig,iv,iv2,istep
      _REAL_ ::  resoz0=0

      !first :: absolute first time in rxrism
      logical, save :: first = .true.

      !irank  :: mpi rank counter
      !stat :: iostat
      integer :: irank,stat, ientry,nentry
      !ierr :: mpi error
      logical :: ierr
      call rism_timer_start(this%rxrismTimer)

      if(this%verbose>=1 .and. size(this%closureList) > 1)&
           call rism_report_message("|Using "//&
           trim(rism3d_closure_type(this%closure))//" closure")

      this%cuvres=0
      !!FIX - delete this...

#ifdef MPI
#else
      if(this%mpirank == 0) then 
         inquire (file=cuvsav,exist=found)
         if (found .and. first .and. ksave/=0)  then
#ifdef RISM_DEBUG
            write(0,*)'reading saved Cuv file:  ',cuvsav
#endif /*RISM_DEBUG*/
            call  rdufma (cuvsav, this%cuv(:,:,:,:),this%grid%nrTotal,&
                 this%solv%natom)
         else
#endif /*MPI*/
            !!ENDFIX
            !..... initial Cuv(R)
            if(this%nsolution == 0 .or. this%ncuvsteps == 0) then
               !       if(.true.) then
               this%cuv=0
               !..... add long-range part,
               !..... because this long-range part is subtracted in next routine..
               if(this%solu%charged)then
                  do iv=1,this%solv%natom
!!$               do igz=1,this%grid%ngr(3)
                     do igz=1,this%grid%nr(3)
                        do igy=1,this%grid%nr(2)
                           do igx=1,this%grid%nr(1)
                              ig = igx + (igy-1)*this%grid%nr(1) + (igz-1)*this%grid%nr(2)*this%grid%nr(1)
                              this%cuv(igx,igy,igz,iv) = this%cuv(igx,igy,igz,iv)&
                                   + this%solv%charge(iv)*this%pot%asymcr(ig)
                           enddo
                        enddo
                     enddo
                  enddo
               end if
            endif
#ifdef MPI
#else
         end if
      end if
#endif /*MPI*/
      !.......................... relaxing UV RISM  ..........................
#ifdef RISM_DEBUG
      write(0,*)'relaxing 3D uv RISM:'
      call flush(0)
#endif
      call mdiis_reset(this%mdiis_o)
      this%cuv=>this%cuvWRK(:,:,:,:,mdiis_getWorkVector(this%mdiis_o))
      this%cuvres=>this%cuvresWRK(:,:,mdiis_getWorkVector(this%mdiis_o))

      do istep=1,maxste
         
         !................... one relaxation step of UV RISM ....................
         call timer_start(TIME_R1RISM)    
         call r1rism(this,resoz0,conver,tol)
         call timer_stop(TIME_R1RISM)    
         
         !............. showing selected and last relaxation steps ..............
         if (kshow /= 0 .and. this%mpirank == 0 .and. this%verbose >= 2)  then
            if (conver .OR. mod(istep,kshow) == 0 .OR. &
                 ksave > 0 .AND. mod(istep,max(ksave,1)) == 0)  then
               call rism_report_message('(a,i5,5x,a,1pg10.3,5x,a,i3)',&
                    ' Step=',istep, 'Resid=',resoz0, 'IS=',getCurrentNVec(this%mdiis_o))
               call rism_report_flush()
            endif
         endif

         !!FIX _ DELETE this
         !.............. saving selected and last relaxation steps ..............
#ifdef MPI
#else
         if (ksave /= 0 .and. first)  then
            if (conver .OR. ksave > 0 .AND. mod(istep,ksave) == 0)  then
               call  wrufma (cuvsav, this%cuv(:,:,:,:),this%grid%nrTotal,this%solv%natom)
            endif
         endif
#endif /*MPI*/
         !!endfix
         !............... exiting relaxation loop on convergence ................
         if (conver)  goto 30
      enddo
      call rism_report_error('(a,i5)','RXRISM: reached limit # of relaxation steps: ',maxste)
30    first = .false.
      if(this%mpirank == 0 .and. this%verbose >= 1) then
         call rism_report_message('(a,i5,a)',"|RXRISM converged in ",istep)!," steps")
      end if
      call rism_timer_stop(this%rxrismTimer)
      return
    end subroutine rxrism

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! One Relaxation Step                           *
!!! for the UV RISM equation with the HNC Closure,             *
!!! Guv(R) = exp( - This%pot%uuv(R) + Tuv(R) - DelHv0 ) + DelHv0        *
!!! Cuv(R) = Guv(R) - 1 - Tvv(R)                       *
!!! Huv(K) = Cuv(K) * (Wvv(K)+Rho*Hvv(K))                 *
!!! TuvRes(R) = Huv(R) - Guv(R) - 1                    *
!!!IN:
!!!  this :: rism3d object
!!!  resoz0 ::
!!!  conver ::
!!!  tol    :: target residual tolerance
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine r1rism(this,resoz0,conver,tol)
      use rism3d_fft_c
      use rism_util, only : checksum
      implicit none
#include "def_time.h" 
#if defined(MPI)
      include 'mpif.h'
#endif /*defined(MPI)*/
      type(rism3d), intent(inout) :: this
      logical, intent(inout) :: conver
    _REAL_, intent(inout) ::  resoz0
    _REAL_, intent(in) :: tol
    integer :: iis
    _REAL_ :: earg, tuv0,tvvr
    integer :: istep

    integer ::  ig1,iga, iv,iv1,iv2, igx, igy, igz, igk
#ifdef FFW_THREADS
    integer :: nthreads,totthreads
    integer,external :: OMP_get_max_threads,OMP_get_num_threads
    logical,external :: OMP_get_dynamic,OMP_get_nested
#endif
    integer :: ierr,irank
    call rism_timer_start(this%r1rismTimer)

#ifdef RISM_DEBUG
    write(0,*)"R1RISM"
    call flush(0)
#endif
    !.....subtract short-range part from Cuv(R)
    !.....short-range part of Cuv(R) is loaded in guv array
#if defined(MPI)
    do iv=1,this%solv%natom
       do igz=1,this%grid%nr(3)
          do igy=1,this%grid%nr(2)
             if(this%solu%charged)then
                do igx=1,this%grid%nr(1)
                   ig1 = igx + (igy-1)*this%grid%nr(1) + (igz-1)*this%grid%nr(2)*this%grid%nr(1)
                   igk = igx + (igy-1)*(this%grid%nr(1)+2) + (igz-1)*this%grid%nr(2)*(this%grid%nr(1)+2)
                   this%guv(igk,iv) = this%cuv(igx,igy,igz,iv) - this%solv%charge(iv)*this%pot%asymcr(ig1)
                enddo
             else
                do igx=1,this%grid%nr(1)
                   igk = igx + (igy-1)*(this%grid%nr(1)+2) + (igz-1)*this%grid%nr(2)*(this%grid%nr(1)+2)
                   this%guv(igk,iv) = this%cuv(igx,igy,igz,iv)
                end do
             end if
             igk = this%grid%nr(1)+1 + (igy-1)*(this%grid%nr(1)+2) + (igz-1)*this%grid%nr(2)*(this%grid%nr(1)+2)
             this%guv(igk:igk +1,iv) =0
          enddo
       enddo
    enddo
#else
    do iv=1,this%solv%natom
       if(this%solu%charged)then
          do igz=1,this%grid%nr(3)
             do igy=1,this%grid%nr(2)
                do igx=1,this%grid%nr(1)
                   ig1 = igx + (igy-1)*this%grid%nr(1) + (igz-1)*this%grid%nr(2)*this%grid%nr(1)
                   this%guv(ig1,iv) = this%cuv(igx,igy,igz,iv) - this%solv%charge(iv)*this%pot%asymcr(ig1)
                enddo
             enddo
          enddo
       else
          do igz=1,this%grid%nr(3)
             do igy=1,this%grid%nr(2)
                do igx=1,this%grid%nr(1)
                   ig1 = igx + (igy-1)*this%grid%nr(1) + (igz-1)*this%grid%nr(2)*this%grid%nr(1)
                   this%guv(ig1,iv) = this%cuv(igx,igy,igz,iv)
                enddo
             enddo
          enddo
       end if
          !zero out extra space
       this%guv(this%grid%nrTotal+1:this%grid%nkTotal,iv) =0d0
    enddo
#endif /*defined(MPI)*/

    !.....short-range part of Cuv(R) FFT>K
    call timer_start(TIME_RISMFFT)    
    call rism_timer_start(this%fftTimer)
!    do iv=1,this%solv%natom
#if defined(MPI)
       call  rism3d_fft_fwd(this%fft,this%guv)
       this%guv(2:this%grid%nkTotal:2,:) =&
            -this%guv(2:this%grid%nkTotal:2,:)
#else
       call  rism3d_fft_fwd(this%fft,this%guv)
#endif /*defined(MPI)*/
    call rism_timer_stop(this%fftTimer)
    call timer_stop(TIME_RISMFFT)    
    !.....add long-range part to Cuv(K) in K-space
    if(this%solu%charged .or. .true.)then
       do iv=1,this%solv%natom
          do ig1=1,this%grid%nkTotal
             this%guv(ig1,iv) = this%guv(ig1,iv) -  this%solv%charge(iv)*this%pot%asymck(ig1)
          enddo
       enddo
    end if

    !.....Huv(k) by RISM
    do iv1=1,this%solv%natom
       do ig1=1,this%grid%nkTotal
          this%huv(ig1,iv1) = 0d0
          iga = this%grid%indga((ig1+1)/2)
          do iv2=1,this%solv%natom
             this%huv(ig1,iv1) = this%huv(ig1,iv1) + &
                  this%guv(ig1,iv2)&
                  *this%xvva(iga + (iv2-1)*this%grid%nga + (iv1-1)*this%grid%nga*this%solv%natom)
          enddo
       enddo
    enddo

    !.....add long-range part of Huv(k) at k=0
    !.....which was estimated by long-range part of Cuv(k) at k=0
    if(this%mpirank==0)then
       do ig1=1,2
          do iv=1,this%solv%natom
             this%huv(ig1,iv) = this%huv(ig1,iv) + this%pot%huvk0(ig1,iv)
          enddo
       enddo
    end if

    !.....subtract long-range part from huv in K-space
    if(this%solv%ionic)then
#if defined(MPI)
       do iv=1,this%solv%natom
          if(this%mpirank==0)then
             do ig1=3,this%grid%nkTotal
                this%huv(ig1,iv) = this%huv(ig1,iv) &
                     + 1d0/this%solv%dielconst * this%solv%charge_sp(iv)*this%pot%asymhk(ig1)
             enddo
          else
             do ig1=1,this%grid%nkTotal
                this%huv(ig1,iv) = this%huv(ig1,iv) &
                     + 1d0/this%solv%dielconst * this%solv%charge_sp(iv)*this%pot%asymhk(ig1)
             enddo
          end if
       enddo
#else
       do iv=1,this%solv%natom
          do ig1=3,this%grid%nkTotal
             this%huv(ig1,iv) = this%huv(ig1,iv) &
                  + 1d0/this%solv%dielconst * this%solv%charge_sp(iv)*this%pot%asymhk(ig1)
          enddo
       enddo
#endif /*defined(MPI)*/
    end if

    !.....short-range part of Huv(K) FFT>R
    call timer_start(TIME_RISMFFT)    
    call rism_timer_start(this%fftTimer)
#if defined(MPI)
       this%huv(2:this%grid%nkTotal:2,:)=&
            -this%huv(2:this%grid%nkTotal:2,:)
       call  rism3d_fft_bwd(this%fft,this%huv)
#else
       call  rism3d_fft_bwd(this%fft,this%huv)
#endif /*defined(MPI)*/
    call rism_timer_stop(this%fftTimer)
    call timer_stop(TIME_RISMFFT)    

    !.....add long-range part to huv in R-space
    if(this%solv%ionic)then
       do iv=1,this%solv%natom
          do igz=1,this%grid%nr(3)
             do igy=1,this%grid%nr(2)
                do igx=1,this%grid%nr(1)
                   ig1 = igx + (igy-1)*this%grid%nr(1) + (igz-1)*this%grid%nr(2)*this%grid%nr(1)
#if defined(MPI)
                   igk = igx + (igy-1)*(this%grid%nr(1)+2) + (igz-1)*this%grid%nr(2)*(this%grid%nr(1)+2)
                   this%huv(igk,iv) = this%huv(igk,iv) &
                        +  this%solv%charge_sp(iv)*this%pot%asymhr(ig1)
#else
                   this%huv(ig1,iv) = this%huv(ig1,iv) &
                        +  this%solv%charge_sp(iv)*this%pot%asymhr(ig1)
#endif /*defined(MPI)*/
                enddo
             enddo
          enddo
       enddo
    end if
    this%cuvres(:,:)=0
#ifdef RISM3D_DEBUG
#  ifdef MPI    
    call rism3d_debug_print(this%huv,.true.,"PRE_CLOSE",.false.)
#  else
    call rism3d_debug_print(this%huv,.true.,"PRE_CLOSE",.true.)
#  endif
#endif    
    call rism3d_closure_guv(this%closure,this%guv,this%huv,this%cuv)
#ifdef RISM3D_DEBUG
#  ifdef MPI    
    call rism3d_debug_print(this%guv,.true.,"POST_CLOSE_GUV",.false.)
    call rism3d_debug_print(this%huv,.true.,"POST_CLOSE_HUV",.false.)
#  else
    call rism3d_debug_print(this%guv,.true.,"POST_CLOSE_GUV",.true.)
    call rism3d_debug_print(this%huv,.true.,"POST_CLOSE_HUV",.true.)
#  endif
#endif    
    do iv=1,this%solv%natom
       do igz=1,this%grid%nr(3)
          do igy=1,this%grid%nr(2)
             do igx=1,this%grid%nr(1)
                ig1 = igx + (igy-1)*this%grid%nr(1) + (igz-1)*this%grid%nr(2)*this%grid%nr(1)
#if defined(MPI)
                igk = igx + (igy-1)*(this%grid%nr(1)+2) + (igz-1)*this%grid%nr(2)*(this%grid%nr(1)+2)
#else
                igk=ig1
#endif /*defined(MPI)*/
                this%cuvres(ig1,iv) = this%guv(igk,iv) - 1d0 - this%huv(igk,iv)
             enddo
          enddo
       enddo
    enddo

    call timer_start(TIME_MDIIS)    
#ifdef RISM3D_DEBUG
    call rism3d_debug_print(this%cuv,"PRE_MDIIS")
    call rism3d_debug_print(this%cuvres,"PRE_MDIIS_res")
#endif    
    call mdiis_advance (this%mdiis_o, resoz0,conver,tol)
    this%cuv=>this%cuvWRK(:,:,:,:,mdiis_getWorkVector(this%mdiis_o))
    this%cuvres=>this%cuvresWRK(:,:,mdiis_getWorkVector(this%mdiis_o))
    call timer_stop(TIME_MDIIS)    
    call rism_timer_stop(this%r1rismTimer)
#ifdef RISM3D_DEBUG
    call rism3d_debug_print(this%cuv,"POST_MDIIS")
    call rism3d_debug_print(this%cuvres,"POST_MDIIS_res")
#endif    
  end subroutine r1rism

!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!subroutines to find the iterative 3D-RISM solution.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PROPAGATE PREVIOUS SOLUTIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculates a new initial guess for CUV based on the final solutions 
!from previous timesteps.  The maximum number of previous time steps to
!used is provided by the user in ncuvsteps.  However, if there are not
!enough previous timesteps only nsolution previous timesteps will be used.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cuv_propagate(this)
    implicit none
    type(rism3d) :: this
    integer :: iv, igr,n
    n=this%grid%nrTotal*this%solv%natom
    if(this%ncuvsteps >= 5 .and. this%nsolution >= 5) then
!!$     this%cuv(:,:,:,:) = 5d0*this%cuv(:,:,:,:) - 10d0*this%oldcuv(:,:,:,:,2)&
!!$          +10d0*this%oldcuv(:,:,:,:,3) - 5d0*this%oldcuv(:,:,:,:,4) + this%oldcuv(:,:,:,:,5)
       call dscal(n,5d0,this%cuv(:,:,:,:),1)
       call daxpy(n,-10d0,this%oldcuv(:,:,:,:,2),1,this%cuv(:,:,:,:),1)
       call daxpy(n,10d0,this%oldcuv(:,:,:,:,3),1,this%cuv(:,:,:,:),1)
       call daxpy(n,-5d0,this%oldcuv(:,:,:,:,4),1,this%cuv(:,:,:,:),1)
       call daxpy(n,1d0,this%oldcuv(:,:,:,:,5),1,this%cuv(:,:,:,:),1)
    elseif(this%ncuvsteps >= 4 .and. this%nsolution >= 4) then
!!$     this%cuv(:,:,:,:) = 4d0*this%cuv(:,:,:,:) - 6d0*this%oldcuv(:,:,:,:,2)&
!!$          +4d0*this%oldcuv(:,:,:,:,3) - this%oldcuv(:,:,:,:,4)
       call dscal(n,4d0,this%cuv(:,:,:,:),1)
       call daxpy(n,-6d0,this%oldcuv(:,:,:,:,2),1,this%cuv(:,:,:,:),1)
       call daxpy(n,4d0,this%oldcuv(:,:,:,:,3),1,this%cuv(:,:,:,:),1)
       call daxpy(n,-1d0,this%oldcuv(:,:,:,:,4),1,this%cuv(:,:,:,:),1)
    elseif(this%ncuvsteps >= 3 .and. this%nsolution >= 3) then
!!$     this%cuv(:,:,:,:) = 3d0*(this%cuv(:,:,:,:) - this%oldcuv(:,:,:,:,2))&
!!$          +this%oldcuv(:,:,:,:,3)
       call dscal(n,3d0,this%cuv(:,:,:,:),1)
       call daxpy(n,-1d0,this%oldcuv(:,:,:,:,2),1,this%cuv(:,:,:,:),1)
       call daxpy(n,1d0,this%oldcuv(:,:,:,:,3),1,this%cuv(:,:,:,:),1)
    elseif(this%ncuvsteps >= 2 .and. this%nsolution >= 2) then
!!$     this%cuv(:,:,:,:) = 2*this%cuv(:,:,:,:) - this%oldcuv(:,:,:,:,2)
       call dscal(n,2d0,this%cuv(:,:,:,:),1)
       call daxpy(n,-1d0,this%oldcuv(:,:,:,:,2),1,this%cuv(:,:,:,:),1)
    elseif(this%ncuvsteps == 0) then
       this%cuv(:,:,:,:) = 0
    end if
  end subroutine cuv_propagate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Updates the values in the this%oldcuv queue.  The oldest value (the
!ncuvstep index) is pushed out, the remainder of the data is shifted 
!and the newest solution is placed in the first index
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cuv_update(this)
    implicit none
    type(rism3d) :: this
    integer :: istep,iv,igr
    !update this%oldcuv
#ifdef RISM_DEBUG
    write(0,*) "CUV_UPDATE"; call flush(0)
#endif /*RISM_DEBUG*/  
    if(this%ncuvsteps == 0) return
    do istep = min(this%ncuvsteps,this%nsolution),2,-1
       call dcopy(this%grid%nrTotal*this%solv%natom,this%oldcuv(:,:,:,:,istep-1),1,&
            this%oldcuv(:,:,:,:,istep),1)
    end do
    call dcopy(this%grid%nrTotal*this%solv%natom,this%cuv(:,:,:,:),1,&
         this%oldcuv(:,:,:,:,1),1)
  end subroutine cuv_update
end module rism3d_c
