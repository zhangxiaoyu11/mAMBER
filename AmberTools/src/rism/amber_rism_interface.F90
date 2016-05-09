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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!General interface between 3D-RISM and SANDER/SFF in the Amber suite.  Except
!!!where noted, all data and subroutines may be called from either SFF or SANDER.
!!!
!!!To make this file interoperable with C while still maintaining the Fortran 95
!!!standard, we must not put subroutines or functions inside of modules. However,
!!!some variables must be available at the global scope either to be accessed
!!!from outside routines (e.g. print statements) or be shared between local
!!!routines (e.g. instances of derived types and parameters for the run).  This 
!!!is accomplished with two modules. AMBER_RISM_INTERFACE is always compiled and
!!!contains RISM specific variables.  SANDER_RISM_INTERFACE is created for
!!!SANDER but not for SFF. In the Fortran world (SANDER), the 
!!!SANDER_RISM_INTERFACE module provides access to subroutines and functions 
!!!with all the benefits of a module. In the C world, there is no module and
!!!functions/subroutines may be directly called.
!!!
!!!Setup and initialization (serial and MPI):
!!!The method of setting up and initializing is somewhat flexible and complicated
!!!in order to maintain the correct output initialization and output sequence
!!!of SANDER.  SANDER reads all of the input files and then prints a summary
!!!from the master node.  RISM must do the same, so RISM_SETPARAM must be called
!!!from the master node in SANDER.  However, it is safe to call it from all
!!!nodes at the same time with the caveat that irism is defined and the same on
!!!all nodes; this done in SFF.  Initializing the calculation
!!!must be done in parallel so RISM_INIT must be call from all processes. 
!!!
!!!Note that it is always safe to call RISM_SETPARAM and RISM_INIT as long as
!!!irism is define on the master node in the relevent data structure.
!!!To summarize:
!!!
!!!In SANDER:
!!!
!!! if(master)then 
!!!   call rism_setparam(mdin,&
!!!        commsander,&
!!!        igb,natom,ntypes,x(L15:L15+natom-1),&
!!!        x(LMASS:LMASS+natom-1),cn1,cn2,&
!!!        ix(i04:i04+ntypes**2-1), ix(i06:i06+natom-1))
!!! endif
!!! call rism_init(commsander)
!!!
!!!In SFF:
!!! rism_setparam_( &rismData, &xvvlen,xvvfile, 
!!!                 &guvlen,guvfile,&huvlen,huvfile,&cuvlen,cuvfile,
!!!                 &uuvlen,uuvfile,&asymplen,asympfile,
!!!                 &quvlen,quvfile,&chgdistlen,chgdistfile,
!!!                 &volfmtlen,volfmt,
!!!                 &comm, 
!!!                 &gb, &(prm->Natom), &(prm->Ntypes), prm->Charges,  
!!!                 prm->Masses,prm->Cn1, prm->Cn2, prm->Iac, prm->Cno); 
!!! rism_init_(&comm);
!!!
!!!MPI and other subroutines and functions:  
!!!All public functions are safe to call from all nodes and usually must be.
!!!Only RISM_THERMO_PRINT and RISM_SETPARAM can be called from just the master.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Module to hold global data and preserve namespace.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module amber_rism_interface
  use rism3d_c
  use rism3d_solv_c
  use rism3d_solu_c
  use fce_c
  use rism_report_c
  use rism_timer_c
  use safemem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Parameter derived type for storing calculation parameters.  This can be use
!!!to transfer parameters from a C program.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type rismprm_t
     sequence
     !solvcut :: cutoff for rism calculations (separate from SANDER non-bond). 
     _REAL_ :: solvcut
     !buffer   :: buffer distance to the edge of the box for the solvent
     _REAL_ :: buffer
     !grdspc   :: grid spacing for all of the grids
     _REAL_ :: grdspc(3)
     !boxlen   :: box size for 3d-rism.  for PBC calculations, these should generallly be equal
     _REAL_ :: solvbox(3)
     !mdiis_del :: 'step size' for MDIIS
     _REAL_ :: mdiis_del
     !mdiis_restart :: restart threshold factor. Ratio of the current residual to the 
     !                 minimum residual in the basis that causes a restart
     _REAL_ :: mdiis_restart
     !fcecut :: FCE cutoff distance
     _REAL_ :: fcecut
     !closureOrder :: for backwards compatibility, we still need to read this
     integer :: closureOrder
     !ng3(3)    :: number of grid points in each dimension
     integer :: ng3(3)
     !irism - use 3d-rism
     integer :: irism
     !asympCorr :: use long range asymptotic corrections for thermodynamic calculations
     logical*4 :: asympCorr
     !mdiis_nvec       :: number of vectors used for MDIIS (consequently, the number of copies of
     !             CUV we need to keep for MDIIS)
     integer :: mdiis_nvec
     !mdiis_method :: MDIIS implementation to use    
     integer :: mdiis_method
     !maxstep :: maximum number of rism relaxation steps
     integer :: maxstep
     !npropagate :: number of past cuv time steps saves
     integer :: npropagate
     !centering :: center the solute in the solvation box.  
     !           0 - off
     !           1 - center of mass
     !           2 - center of geometry
     !           3 - center of mass shifted to the nearest grid point
     !           4 - center of geometry shifted to the nearest grid point
     !           For negative numbers the centering translation is only calculated
     !           the for the first solution and used for subsequent calculations.
     !           This allows the solute to drift in the box.
     integer :: centering
     !zerofrc ::0 - do nothing 
     !          1 - redistribute forces to get zero total force
     integer :: zerofrc
     !apply_rism_force :: if 0, the 3D-RISM solution is calculated but the resulting forces are not
     integer :: apply_rism_force
     !polarDecomp :: do a polar/apolar decomposition of the chemical potential
     integer :: polarDecomp
     !rismnrespa :: size of rism multiple timestep
     integer :: rismnrespa
     !fce_stride :: FCE MTS stride length
     integer :: fcestride
     !fcenbasis :: number of FCE basis vectors
     integer :: fcenbasis
     !fcecrd :: FCE coordinate basis type
     integer ::  fcecrd
     !saveprogress  :: save itermediate results every saveprogress interations (0 means no saves)
     integer :: saveprogress
     integer :: ntwrism
     !verbose :: 0 - no ouput
     !                1 - memory allocation and steps for convergence
     !                2 - 1 + convergence progress
     integer :: verbose
     !progress :: print parameter for relaxation steps every progress
     !            iteration (0 means no saves)
     integer :: progress

     !write_thermo :: calculate and print out thermodynamics.  This is
     !                primarily used by sander but also serves as
     !                padding for alignment for NAB
     integer :: write_thermo

     ! BPR: Make sure the number of INTEGERs is not odd,
     ! to shut up a compiler whinge about misalignment
     ! Note: this should be commented if ever there are grounds
     ! for a new real INTEGER
     integer :: padding     
  end type rismprm_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Derived type to store thermodynamic results.  To facilitate MPI
!!!communication, many values are stored in a common array.  Pointers,
!!!however, are used to access the data.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type rismthermo_t

     !Thermodynamic values that use mpi_buffer for memory
     _REAL_, pointer :: exchem(:)=>NULL(), exchemGF(:)=>NULL(), exchemUC=>NULL(),&
          potUV(:)=>NULL(), solvEne(:)=>NULL(), apol_exchem(:)=>NULL(), &
          apol_exchemGF(:)=>NULL(), apol_exchemUC=>NULL(), apol_potUV(:)=>NULL(),&
          apol_solvEne(:)=>NULL(), PMV=>NULL(), exNum(:)=>NULL()
     !Thermodynamic values that have their own memory and are not passed by MPI
     _REAL_, pointer :: pol_exchem(:)=>NULL(), pol_exchemGF(:)=>NULL(), &
          pol_potUV(:)=>NULL(), pol_solvEne(:)=>NULL()
     !mpi_buffer : all the results are stored in a single array so the
     !         necessary reductions are done in a single communication.
     !         This is also used for the serial calculation for
     !         simplicity.  Order and size: exchem(solv%natom),
     !         exchemGF(solv%natom) exchemUC(1), potUV(solv%natom),
     !         solvEne(solv%natom), apol_exchem(solv%natom),
     !         apol_exchemGF(solv%natom), apol_exchemUC(1),
     !         potUV(solv%natom), solvEne(solv%natom),
     !         excessNum(solv%natom), PMV(1)
     _REAL_, pointer :: mpi_buffer(:)
#if !defined( USE_MPI_IN_PLACE) && defined( MPI )
     _REAL_, pointer :: tmpi_buffer(:)
#endif
  end type rismthermo_t

  !Possible RISM calculation types for MD
  !RISM_NONE :: no RISM calculation
  !RISM_FULL :: full RISM solution
  !RISM_INTERP :: interpolation
  integer, parameter :: RISM_NONE=0, RISM_FULL=1, RISM_INTERP=2

  type(rismprm_t), save :: rismprm
  type(rismthermo_t), save :: rismthermo
  type(rism3d),save :: rism_3d
  type(fce),save :: fce_o
  type(rism3d_solv),save :: solv
  type(rism3d_solu),save :: solu
  type(rism_timer), save :: timer
  type(rism_timer), save :: timer_write
  type(rism_timer), save :: timer_init
  integer :: pa_orient, rmsd_orient
  _REAL_:: ratucm(3)

  integer :: outunit

  !closurelist :: list of closures to use in order.  Only the last closure
  !           is used for thermodynamic output.  This can be used to
  !           progressively increase the order of the closure to aid
  !           convergence.  Closure types: KH, HNC or PSEn where n is
  !           an integer. This is initialized to a length of 10 in
  !           defaults() to allow it to be used with the namelist in
  !           sander
  character(len=8),pointer :: closurelist(:)=>NULL()
  !tolerance :: residual tolerance for the solution of each closure in
  !             the list. On input this can be of length one, two or 
  !             size(closurelist). If length one, use this value for the
  !             final closure and the default for all others (see 
  !             sanity_check()). If length two, use the last value for
  !             the last closure and the first value for all other 
  !             closures. Otherwise, match each value to each closure.
  _REAL_,pointer :: tolerancelist(:)=>NULL()
  !nclosuredefault :: default number of closures to start with in the
  !                   list. Namelists don't support pointers in many
  !                   ways so we need to initialize the closure and
  !                   tolerance lists to some reasonable size
  integer, parameter :: nclosuredefault=10

  !I/O file names:
  !xvvfile     : (input) site-site solvent susceptibility from RISM1D (.xvv)
  !guvfile     : (output) pair distribution function.  Volumetric file.
  !huvfile     : (output) total correlation function.  Volumetric file.
  !cuvfile     : (output) direct correlation function.  Volumetric file.
  !uuvfile     : (output) solute-solvent potential.  Volumetric file.
  !asympfile   : (output) long range asymptotics function.  Volumetric file.
  !quvfile     : (output) charge density distribution.  Volumetric file.
  !chgDistFile : (output) charge distribution. Volumetric file.
  !volfmt      : either 'dx' or 'xyzv'
  character(len=256) :: xvvfile='', guvfile='', huvfile='', cuvfile='',&
       uuvfile='', asympfile='', quvFile='', chgDistFile='', volfmt='dx'

  integer :: mpirank=0,mpisize=1,mpicomm=0

  !working memory for rism_force() so it is not reallocated every time
  !ff :: forces
  !ratu_fce :: coordinates for FCE 
  _REAL_,pointer :: ff(:,:)
#if defined(RISM_CRDINTERP)
  _REAL_,pointer :: ratu_fce(:,:)
#endif /*RISM_CRDINTER*/

  private :: rism_mpi_bcast
end module amber_rism_interface


#ifdef SANDER
module sander_rism_interface
  use amber_rism_interface
  implicit none
  contains
#endif
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Sets all input parameters for 3D-RISM.  This _must_ be called by the head 
!!!node and may be called by all nodes.  If all nodes call this subroutine they
!!!must all agree on the value of rismprm%irism (SANDER) or userData%irism(SFF).
!!!
!!!SANDER prerequisites:
!!!   -Names of 3D-RISM specific I/O files (Xvv, Guv, etc.) are specified on the 
!!!    command line and are read by mdfil.F90 into the variables in 
!!!    AMBER_RISM_INTERFACE.
!!!   -igb should be set to 6 (vacuum electrostatics)
!!!SANDER IN:
!!!   mdin: name of the mdin file that SANDER name lists are read from
!!!
!!!SFF prerequisites:
!!!   -All user options (including file names) are read from mm_options.l.  
!!!    Non-string options are read into a C struct equivalent to rismprm_t.
!!!    String options are supplied as char* array, integer length pairs.
!!!   -gb should be set to 0 (vacuum electrostatics)
!!!SFF IN:
!!!   userdata:    rismprm_t C struct equivalent with use options
!!!   ntol:        number of tolerances in the array
!!!   tol:         array of tolerances read in from the user
!!!   closurelen:  length of the closurechar strings
!!!   nclosure:    number of closures read in
!!!   closurechar: an array of nclosure strings of closurelen characters
!!!   xvvlen:      length of xvvchar array
!!!   xvvchar:     character array for Xvv input file name
!!!   guvlen:      length of guvchar array
!!!   guvchar:     character array for Guv output file name
!!!   huvlen:      length of huvchar array
!!!   huvchar:     character array for Huv output file name
!!!   cuvlen:      length of cuvchar array
!!!   cuvchar:     character array for Cuv output file name
!!!   uuvlen:      length of uuvchar array
!!!   uuvchar:     character array for Uuv output file name
!!!   asymplen:    length of asympchar array
!!!   asympchar:   character array for asymptotics output file name
!!!   quvlen:      length of quvchar array
!!!   quvchar:     character array for Quv output file name
!!!   chgDistlen:  length of chgDistchar array
!!!   chgDistchar: character array for charge distribution output file name
!!!   volfmtlen:   length of volfmtchar array
!!!   volfmtchar:  character array for the format type for volumetric data
!!!
!!!IN:
!!!   comm : MPI communicator
!!!   natom: number of solvent atoms
!!!   ntypes: number of atom solvent types
!!!   charge: solvent atom partial charges in Amber units
!!!   mass: solvent atom masses [AU]
!!!   cn1: cn1 Lennard-Jones parameters for each atom type pair
!!!   cn2: cn2 Lennard-Jones parameters for each atom type pair
!!!   iac: ATOM_TYPE_INDEX
!!!   ico: NONBONDED_PARM_INDEX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rism_setparam(&
#ifdef SANDER
     mdin,&
#else /*.not.SANDER*/
     userData,ntol,tol,&
     closurelen,nclosure,closurechar,xvvlen,xvvchar,&
     guvlen,guvchar,huvlen,huvchar,cuvlen,cuvchar,&
     uuvlen,uuvchar,asymplen,asympchar,quvlen,quvchar,chgDistlen,chgDistchar,&
     volfmtlen,volfmtchar, &
#endif /*SANDER*/
     comm,&
     natom, ntypes, charge, mass, cn1, cn2, iac, ico)
  use amber_rism_interface

  use constants, only : KB
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif /*MPI*/
  integer, intent(in) ::natom, ntypes,iac(ntypes**2),ico(natom)
  _REAL_, intent(in) :: charge(natom),mass(natom),cn1(ntypes*(ntypes+1)/2),&
       cn2(ntypes*(ntypes+1)/2)
#ifdef SANDER
  !  character(*), intent(in) :: xvvfile,mdin
  character(*), intent(in) :: mdin
  integer :: mdin_unit=55
#else /*SANDER*/
  type(rismprm_t), intent(in) :: userData;
  integer, intent(in) :: ntol
  _REAL_, intent(in)  :: tol(ntol)
  integer, intent(in) :: closurelen, nclosure, xvvlen, guvlen, &
       huvlen, cuvlen, uuvlen,&
       asymplen, quvlen, chgDistLen, volfmtlen
  integer(kind=1), intent(in) :: closurechar(closurelen*nclosure),&
       xvvchar(xvvlen+1),guvchar(guvlen+1),&
       huvchar(huvlen+1),cuvchar(cuvlen+1),uuvchar(uuvlen+1),&
       asympchar(asymplen+1),quvchar(quvlen+1),chgDistchar(chgDistlen+1),&
       volfmtchar(volfmtlen+1)
  !  character(xvvlen):: xvvfile
#endif /*SANDER*/
  integer, intent(in) :: comm
  character(len=16) :: whtspc
  integer :: i,stat, err
  !iclosure :: counter for closures
  !iclosurechar :: current index in the closurechar array from sff
  integer :: iclosure, iclosurechar
  logical :: op
  integer ::  un

  !in case this is not the first time init has been called (i.e.,
  !multiple runs) destroy timers and re-create them.
  call rism_timer_destroy(timer_write)
  call rism_timer_destroy(timer_init)
  call rism_timer_destroy(timer)
  call rism_timer_new(timer,"3D-RISM Total")
  call rism_timer_new(timer_init,"3D-RISM initialization",timer)
  call rism_timer_start(timer_init)
  call rism_timer_new(timer_write,"3D-RISM Output",timer)

  write(whtspc,'(a16)')" "

#ifdef MPI
  mpicomm = comm
  if(mpicomm == MPI_COMM_NULL)&
       call rism_report_error("RISM3D interface: received NULL MPI communicator")
  call mpi_comm_rank(mpicomm,mpirank,err)
  if(err /=0) call rism_report_error&
       ("(a,i8)","RISM3D interface: could not get MPI rank for communicator ",mpicomm)
  call mpi_comm_size(mpicomm,mpisize,err)
  if(err /=0) call rism_report_error&
       ("(a,i8)","RISM3D interface: could not get MPI size for communicator ",mpicomm)
#endif /*MPI*/  

  !if this is not a RISM run, we're done
#ifdef SANDER
 if(rismprm%irism == 0) then
    call rism_timer_stop(timer_init)
    return
 end if
#else /*SANDER*/
 if(userData%irism == 0) then
    call rism_timer_stop(timer_init)
    return
 end if
#endif /*SANDER*/

  !rank 0 only
 if(mpirank /= 0) then
    call rism_timer_stop(timer_init)
    return
 end if


 outunit = rism_report_getMUnit()
 call defaults()
 
#ifdef SANDER
 inquire(file=mdin,opened=op, number=un)
 if(op) mdin_unit=un
 open(unit=mdin_unit,file=mdin,status='OLD', form='FORMATTED',iostat=stat)
 if(stat/=0)then
    call rism_report_error('(a,i4)', "opening "//trim(mdin)//"failed. IOSTAT=",stat)
 end if
 call read_namelist(mdin_unit)
 if(.not.op) close(unit=mdin_unit)
#else /*.not.SANDER*/
 if(ntol /=0)then
    tolerancelist => safemem_realloc(tolerancelist,ntol)
    tolerancelist = tol(1:ntol)
 end if
 call update_param(userData)
 iclosurechar=1
 closurelist => safemem_realloc(closurelist,len(closurelist),max(nclosure,1))
 do iclosure=1, nclosure
    !the default in SFF is to pass an empty string, since this would
    !overwrite out defaults, ignore it
    if(closurechar((iclosure-1)*closurelen+1)==0)exit
    call cstr2fstr(closurelist(iclosure),closurechar((iclosure-1)*closurelen+1),closurelen)
 end do
 call cstr2fstr(xvvfile,xvvchar,xvvlen)
 call cstr2fstr(guvfile,guvchar,guvlen)
 call cstr2fstr(huvfile,huvchar,huvlen)
 call cstr2fstr(cuvfile,cuvchar,cuvlen)
 call cstr2fstr(uuvfile,uuvchar,uuvlen)
 call cstr2fstr(asympfile,asympchar,asymplen)
 call cstr2fstr(quvfile,quvchar,quvlen)
 call cstr2fstr(chgDistfile,chgDistchar,chgDistlen)
 if(volfmtlen > 0)&
      call cstr2fstr(volfmt,volfmtchar,volfmtlen)
#endif /*SANDER*/
  call sanity_check() 
#ifdef SANDER
  if(rismprm%irism >=1)then
     write(outunit,'(a)') "3D-RISM:"
     if(rismprm%irism < 1) then
        write(outunit,'(5x,a,i10)') 'irism   =',rismprm%irism
     elseif(rismprm%irism == 1) then
        write(outunit,'(5x,3(a10,"=",100a10))') &
             'closure'//whtspc, closurelist
        write(outunit,'(5x,3(a10,"="f10.5))')&
             'solvcut'//whtspc, rismprm%solvcut,&
             ', buffer'//whtspc, rismprm%buffer
        write(outunit,'(5x,a10,"=",3(f10.5,1x))')&
             'grd_spc'//whtspc, rismprm%grdspc
        write(outunit,'(5x,a10,"=",3(i10,1x))')&
             'ng3'//whtspc, rismprm%ng3
        write(outunit,'(5x,a10,"=",3(f10.5,1x))')&
             'solvbox'//whtspc, rismprm%solvbox
        write(outunit,'(5x,a10,"=",1p,100e10.2)')  &
             'tolerance'//whtspc, tolerancelist
        write(outunit,'(5x,a10,"=",f10.5,a10,"=",i10)')  &
             'mdiis_del'//whtspc, rismprm%mdiis_del,&
             ', mdiis_nvec'//whtspc,rismprm%mdiis_nvec
        write(outunit,'(5x,a10,"=",i10,a10,"=",1p,e10.2)')&
             'mdiis_method'//whtspc, rismprm%mdiis_method, &
             ', mdiis_restart'//whtspc,rismprm%mdiis_restart
        write(outunit,'(5x,3(a10,"=",i10))') &
             'maxstep'//whtspc, rismprm%maxstep, &
             ', npropagate'//whtspc,rismprm%npropagate
        write(outunit,'(5x,3(a10,"=",i10))')&
             'centering'//whtspc, rismprm%centering, &
             ', zerofrc'//whtspc, rismprm%zerofrc
        write(outunit,'(5x,a10,"=",i10,a10,"=",l10)') &
             'apply_rism_force'//whtspc, rismprm%apply_rism_force, &
             ' asympcorr'//whtspc, rismprm%asympcorr
        write(outunit,'(5x,2(a10,"=",i10),a10,"=",f10.5)')&
             'rismnrespa'//whtspc, rismprm%rismnrespa, &
             ', fcestride'//whtspc, rismprm%fcestride, &
             ', fcecut'//whtspc,rismprm%fcecut
        write(outunit,'(5x,3(a10,"=",i10))')&
             'fcenbasis'//whtspc, rismprm%fcenbasis, &
             ', fcecrd'//whtspc, rismprm%fcecrd
        write(outunit,'(5x,2(a10,"=",i10),a10,"=  ",a8)')&
             'polardecomp'//whtspc, rismprm%polardecomp, &
             ' write_thermo'//whtspc, rismprm%write_thermo,&
             ' volfmt'//whtspc, volfmt
        write(outunit,'(5x,3(a10,"=",i10))')&
             'saveprogress'//whtspc, rismprm%saveprogress, &
             ', ntwrism'//whtspc, rismprm%ntwrism,&
             ', verbose'//whtspc, rismprm%verbose
        write(outunit,'(5x,3(a10,"=",i10))') &
             'progress'//whtspc, rismprm%progress
     endif
     call flush(outunit)
  end if
#endif /*SANDER*/

  !initialize 3D-RISM solute and solvent
  call rism3d_solv_new(solv,xvvfile)
 
  call rism3d_solu_new_sander(solu,natom,ntypes,iac,&
      ico,charge,cn1,cn2,mass,solv%temperature)

  call rism_timer_stop(timer_init)
end subroutine rism_setparam

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Performs all of the initialization required for 3D-RISM calculations. All
!!!parameters should have already been set with RISM_SETPARAM. For MPI 
!!!calculations this _must_ be called by all processes. Only the head node irism
!!!value is used. 
!!!IN:
!!!  comm : MPI communicator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rism_init(comm)
  use amber_rism_interface
  use safemem
#ifdef RISM_CRDINTERP
  use fce_c
#endif /*RISM_CRDINTERP*/
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif /*MPI*/  
  integer, intent(in) :: comm

  integer :: err

  !ensure that timers have been created incase RISM_SETPARAM was only
  !called by the master node
  if(trim(timer%name) .ne. "3D-RISM Total")&
       call rism_timer_new(timer,"3D-RISM Total")
  if(trim(timer_init%name) .ne. "3D-RISM initialization")&
       call rism_timer_new(timer_init,"3D-RISM initialization",timer)
  call rism_timer_start(timer_init)
  if(trim(timer_write%name) .ne. "3D-RISM Output")&
       call rism_timer_new(timer_write,"3D-RISM Output",timer)

#ifdef SANDER
  call rism_report_setMUnit(6)
  call rism_report_setWUnit(6)
  call rism_report_setEUnit(6)
#endif /*SANDER*/

#ifdef MPI
  mpicomm = comm
  if(mpicomm == MPI_COMM_NULL)&
       call rism_report_error("RISM3D interface: received NULL MPI communicator")
  call mpi_comm_rank(mpicomm,mpirank,err)
  if(err /=0) call rism_report_error&
       ("(a,i8)","RISM3D interface: could not get MPI rank for communicator ",mpicomm)
  call mpi_comm_size(mpicomm,mpisize,err)
  if(err /=0) call rism_report_error&
       ("(a,i8)","RISM3D interface: could not get MPI size for communicator ",mpicomm)
  call rism_mpi_bcast(mpirank,mpisize,mpicomm)
#endif /*MPI*/

  !STOP HERE IF THIS IS NOT A RISM RUN
  if(rismprm%irism == 0) then
     call rism_timer_stop(timer_init)
     return
  end if

#ifdef RISM_CRDINTERP
  call fce_new(fce_o, solu%natom, rismprm%fcenbasis, rismprm%fcecrd, rismprm%fcecut,mpicomm)
#endif /*RISM_CRDINTERP*/

  !3D-RISM may have already been initialized. In the absence of a
  !subroutine to set all of these parameters individually, we destroy
  !the original instance and re-initialize. Since this is a rare
  !event, it should not add to the expense of the calculation in any
  !practical way
  call rism3d_destroy(rism_3d)
  if(rismprm%buffer >0)then
     call rism3d_new(rism_3d,solu,solv,rismprm%centering, rismprm%npropagate,&
          closurelist,rismprm%solvcut,&
          rismprm%mdiis_nvec,rismprm%mdiis_del,rismprm%mdiis_method, rismprm%mdiis_restart,&
          o_buffer=rismprm%buffer, o_grdspc=rismprm%grdspc, o_mpicomm=mpicomm)
  else
     call rism3d_new(rism_3d,solu,solv,rismprm%centering, rismprm%npropagate,&
          closurelist,rismprm%solvcut,&
          rismprm%mdiis_nvec,rismprm%mdiis_del,rismprm%mdiis_method,rismprm%mdiis_restart,&
          o_boxlen=rismprm%solvbox, o_ng3=rismprm%ng3, o_mpicomm=mpicomm)
  end if
  call rism3d_setverbosity(rism_3d,rismprm%verbose)
  call rism3d_setTimerParent(rism_3d,timer)

  !allocate working memory
  ff => safemem_realloc(ff, 3, rism_3d%solu%natom)

#if defined(RISM_CRDINTERP)
  if(rismprm%fcestride>0)&
       ratu_fce => safemem_realloc(ratu_fce, 3, rism_3d%solu%natom)
#endif /* RISM_CRDINTERP */

  !Free up a bit of memory
  call rism3d_solv_destroy(solv)
  call rism3d_solu_destroy(solu)
  
  call flush(outunit)
  call rism_timer_stop(timer_init)
end subroutine rism_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the RISM calculation type for this step.  I.e., no
!!!calculation (RISM_NONE), full RISM solution (RISM_FULL),
!!!interpolation (RISM_INTERP).
!!!IN:
!!!    irespa     :: respa iteration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function rism_calc_type(irespa) result(calc_type)
  use amber_rism_interface

  implicit none
  integer,intent(in) :: irespa
  integer :: calc_type
  if(mod(irespa,rismprm%rismnrespa) /= 0) then !test for RESPA
     calc_type = RISM_NONE
  elseif(rismprm%fcestride >0 .and. fce_o%nsample >= fce_o%nbasis &
       .and. mod(irespa,rismprm%rismnrespa*rismprm%fcestride) /= 0)then
     calc_type = RISM_INTERP
  else
     calc_type = RISM_FULL
  end if
end function rism_calc_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Driver routine for 3D-RISM.
!!!Calculates the 3D-RISM solution, energy and forces.  Determines if an interpolation 
!!!step can be used.  If so, the interpolation code is called, otherwise a full 
!!!3D-RISM calculation is performed.  Energies are only valid for full 3D-RISM 
!!!calculations
!!!IN:
!!!    ratu_md    :: atom positions for solute.
!!!    frc        :: pre-3D-RISM force.  The 3D-RISM forces are added to this
!!!    epol       :: polarization energy, excess chemical potential
!!!    irespa     :: respa iteration
!!!MODIFIED:
!!!    frc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rism_force(ratu_md,frc, epol, irespa)
  use amber_rism_interface
  use constants, only : KB
  use rism3d_c, only : rism3d_solve
  use rism_util, only: corr_drift, alignorient, translate, calc_cm, calc_rotation
  implicit none
#include "def_time.h" 

#ifdef MPI
  include 'mpif.h'
#endif /*MPI*/

  integer,intent(in) :: irespa
  _REAL_,intent(in) :: ratu_md(3,rism_3d%solu%natom)
  _REAL_,intent(inout) :: frc(3,rism_3d%solu%natom)
  _REAL_,intent(out) :: epol

  ! solvation energy and entropy !
  _REAL_ :: epol_e, epol_ts

  !iclosure :: counter for closures
  integer :: iclosure
  integer :: i, iatu
  integer :: err
  _REAL_ :: mpi_temp, pmv
#if !defined(SANDER)
  integer, external :: rism_calc_type
#endif


  call rism_timer_start(timer)
  epol = 0
  ff=0
  !
  !Test for interpolation, minimum samples and if this is an interpolation step
  !
!  if(mod(irespa,rismprm%rismnrespa) /= 0) then !test for RESPA
  if(rism_calc_type(irespa) == RISM_NONE)then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!No forces this steps. DO NOTHING!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined(RISM_CRDINTERP)
!  elseif(rismprm%fcestride >0 .and. fce_o%nsample >= fce_o%nbasis &
!       .and. mod(irespa,rismprm%rismnrespa*rismprm%fcestride) /= 0)then
  elseif(rism_calc_type(irespa) == RISM_INTERP)then
!!!!!!!!!!!!!!!!
!!!FCE forces!!!
!!!!!!!!!!!!!!!!
     if(rismprm%verbose>=2) call rism_report_message("|LINEAR PROJECTION PREDICT!!!")

     call timer_start(TIME_REORIENT)
     ratu_fce = ratu_md
     call orient(rism_3d%solu,ratu_fce,rism_3d%nsolution)
     call timer_stop(TIME_REORIENT)
     !linproj predict
     call timer_start(TIME_CRDINTERP)
     if(rismprm%apply_rism_force==1)&
          call fce_force(fce_o,ff,ratu_fce)
     call timer_stop(TIME_CRDINTERP)

!!$     call timer_start(TIME_REORIENT)
!!$     call unorient(rism_3d%solu,ratu_md)
!!$     call timer_stop(TIME_REORIENT)

     if(rismprm%zerofrc==1)then
#ifdef MPI
        call corr_drift(ff,rism_3d%solu%mass,rism_3d%solu%natom,&
             mpirank,mpisize,mpicomm)
#else
        call corr_drift(ff,rism_3d%solu%mass,rism_3d%solu%natom)
#endif /*MPI*/
     end if
     if(rismprm%fcestride >1)then
        if(rismprm%verbose>=2) call rism_report_message("|IMPULSE FORCE - INTERP!!!")
        ff=rismprm%rismnrespa*ff
     end if

#endif /*RISM_CRDINTERP*/
  else
     if(rismprm%verbose>=2) call rism_report_message("|FULL RISM!!!")
!!!!!!!!!!!!!!!!!!!!!!!!
!!!Full RISM SOLUTION!!!
!!!!!!!!!!!!!!!!!!!!!!!!
     call rism3d_setCoord(rism_3d,ratu_md)
     call rism3d_solve(rism_3d,rismprm%saveprogress,rismprm%progress,&
          rismprm%maxstep,tolerancelist)
#ifdef SANDER
     !ugly, ugly hack.  SANDER runs force on the first frame twice.
     !This messes up the solution propagation.  Here we set the
     !solution counter back one to ignore one of the duplicate
     !solutions
     if(irespa == 1)then
        rism_3d%nsolution = 1
     end if
#else /*SANDER*/
#endif /*SANDER*/
     if(rismprm%apply_rism_force==1)then
        call rism3d_force(rism_3d,ff)
        if(rismprm%zerofrc==1)then
#ifdef MPI
           call corr_drift(ff,rism_3d%solu%mass,rism_3d%solu%natom,&
                mpirank,mpisize,mpicomm)
#else
           call corr_drift(ff,rism_3d%solu%mass,rism_3d%solu%natom)
#endif /*MPI*/
        end if
        !convert to [kcal/mol/A]
        ff = ff*KB*rism_3d%solv%temperature
     end if
     !get the excess chemical potential
     call timer_start(TIME_EXCHEM)
     epol = rism3d_exchem_tot(rism_3d,rismprm%asympCorr)*KB*rism_3d%solv%temperature

     call timer_stop(TIME_EXCHEM)

#ifdef RISM_CRDINTERP
     if(rismprm%fcestride >0 .and. rismprm%apply_rism_force==1)then
        call timer_start(TIME_REORIENT)
        ratu_fce = ratu_md
        call orient(rism_3d%solu,ratu_fce,rism_3d%nsolution)
        call timer_stop(TIME_REORIENT)
        call timer_start(TIME_SAVECRDINTERP)
        call fce_update(fce_o,ff,ratu_fce)
        call timer_stop(TIME_SAVECRDINTERP)
!!$        call timer_start(TIME_REORIENT)
!!$        call unorient(rism_3d%solu,ratu_fce)
!!$        call timer_stop(TIME_REORIENT)
     endif
#endif /*RISM_CRDINTERP*/

     !         if(rismnrespa >1 .and. (.not. interpcrd>0 .or. .not. nsample >= fcenbasis))then
     if(rismprm%rismnrespa >1)then
        if(rismprm%verbose>=2) call rism_report_message("|IMPULSE FORCE!!!")
        ff=rismprm%rismnrespa*ff
     end if

  end if
  if(rismprm%apply_rism_force==1) frc = frc+ff
  call flush(outunit)
  call rism_timer_stop(timer)
end subroutine rism_force

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates and stores thermodynamics for 3D-RISM, optionally
!!!printing distribution files.  Decompostion is performed based off
!!!of the global parameters 'polarDecomp' and 'entropicDecomp'.  If
!!!polarDecomp==.true. repeat the most recent calculation with solute
!!!charges turned off.  Report the usual quanities but also decompose
!!!the chemical potential into polar and non-polar terms. May be
!!!combined with entropicDecomp.  If entropicDecomp ==
!!!.true. calculate the temperature derivative of the most recent
!!!calculation.  Report the energetic and entropic components of the
!!!excess chemical potential.  May be combined with polarDecomp.
!!!
!!!Since performing polar decomposition destroys solvent distributions
!!!for both standard solutions and entropic decompositions,
!!!distributions must be output here.  This is also necessary since
!!!the temperature derivative solution is only calculated here.
!!!
!!!IN:
!!!  writedist : .true. to write distributions
!!!  step : time step number
!!!SIDEEFECTS:
!!!  The first time through, memory is allocated.  This must be destroyed later.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rism_solvdist_thermo_calc(writedist,step)
  use amber_rism_interface
  use constants, only : kb,COULOMB_CONST_E
  implicit none
#ifdef MPI
  include "mpif.h"
#endif /*MPI*/
  logical*4, intent(in) :: writedist
  integer,intent(in) :: step
  integer :: err
  call rism_timer_start(timer_write)
       !setup memory space
     rismthermo%mpi_buffer => safemem_realloc(rismthermo%mpi_buffer,9*rism_3d%solv%natom+3)
     !initialize for case that some values are not calculated
     rismthermo%mpi_buffer = huge(1d0)
#if !defined( USE_MPI_IN_PLACE) && defined( MPI )
     rismthermo%tmpi_buffer => safemem_realloc(rismthermo%tmpi_buffer,9*rism_3d%solv%natom+3)
#endif  
     
     rismthermo%exchem => rismthermo%mpi_buffer(1:rism_3d%solv%natom) 
     rismthermo%exchemGF => rismthermo%mpi_buffer(rism_3d%solv%natom+1:2*rism_3d%solv%natom) 
     rismthermo%potUV => rismthermo%mpi_buffer(2*rism_3d%solv%natom+2:3*rism_3d%solv%natom+1) 
     rismthermo%solvEne => rismthermo%mpi_buffer(3*rism_3d%solv%natom+2:4*rism_3d%solv%natom+1) 
     rismthermo%apol_exchem => &
          rismthermo%mpi_buffer(4*rism_3d%solv%natom+2:5*rism_3d%solv%natom+1) 
     rismthermo%apol_exchemGF => &
          rismthermo%mpi_buffer(5*rism_3d%solv%natom+2:6*rism_3d%solv%natom+1) 
     rismthermo%apol_potUV => &
          rismthermo%mpi_buffer(6*rism_3d%solv%natom+3:7*rism_3d%solv%natom+2) 
     rismthermo%apol_solvEne => &
          rismthermo%mpi_buffer(7*rism_3d%solv%natom+3:8*rism_3d%solv%natom+2) 
     rismthermo%exNum => &
          rismthermo%mpi_buffer(8*rism_3d%solv%natom+3:9*rism_3d%solv%natom+2) 
     rismthermo%pmv => rismthermo%mpi_buffer(9*rism_3d%solv%natom+3)

     !calculate thermodynamics
     call rism_timer_stop(timer_write)
     rismthermo%exchem = &
          rism3d_exchem(rism_3d,rismprm%asympCorr)*KB*rism_3d%solv%temperature
     rismthermo%exchemGF = &
          rism3d_exchemGF(rism_3d,rismprm%asympCorr)*KB*rism_3d%solv%temperature
     rismthermo%potUV = rism3d_solvPotEne(rism_3d)*KB*rism_3d%solv%temperature
     rismthermo%pmv= rism3d_pmv(rism_3d)
     rismthermo%exNum = rism3d_exNum(rism_3d,rismprm%asympCorr)

     if(writedist .and. step>=0)&
          call rism_writeSolvDistF(rism_3d,step)
     call rism_timer_start(timer_write)

     !By doing the the second RISM calculation here output from the
     !RISM routines does not break up the thermodynamic output below
     if(rismprm%polarDecomp==1)then
        !memory
         rismthermo%pol_exchem=>safemem_realloc(rismthermo%pol_exchem,rism_3d%solv%natom)
         rismthermo%pol_exchemGF=>safemem_realloc(rismthermo%pol_exchemGF,rism_3d%solv%natom)
         rismthermo%pol_potUV=>safemem_realloc(rismthermo%pol_potUV,rism_3d%solv%natom)

        !redo calculation with charges off
        call rism_timer_stop(timer_write)
        call rism3d_unsetCharges(rism_3d)
        call rism3d_solve(rism_3d,rismprm%saveprogress,rismprm%progress,&
          rismprm%maxstep,tolerancelist)

        !calculate thermodynamics
        rismthermo%apol_exchem = &
             rism3d_exchem(rism_3d,rismprm%asympCorr)*KB*rism_3d%solv%temperature
        rismthermo%apol_exchemGF = &
             rism3d_exchemGF(rism_3d,rismprm%asympCorr)*KB*rism_3d%solv%temperature
        rismthermo%apol_potUV = rism3d_solvPotEne(rism_3d)*KB*rism_3d%solv%temperature

        call rism3d_resetCharges(rism_3d)
        call rism_timer_start(timer_write)
     end if

     !parallel communication
#ifdef MPI
#  ifdef USE_MPI_IN_PLACE
     if(mpirank==0)then
        call MPI_REDUCE(MPI_IN_PLACE, rismthermo%mpi_buffer, &
             ubound(rismthermo%mpi_buffer,1),MPI_DOUBLE_PRECISION,&
             MPI_SUM,0,mpicomm,err)
     else
        call MPI_REDUCE(rismthermo%mpi_buffer, rismthermo%mpi_buffer, &
             ubound(rismthermo%mpi_buffer,1),MPI_DOUBLE_PRECISION,&
             MPI_SUM,0,mpicomm,err)
     end if
#  else /*USE_MPI_IN_PLACE*/
     call MPI_REDUCE(rismthermo%mpi_buffer, rismthermo%tmpi_buffer, &
          ubound(rismthermo%mpi_buffer,1),MPI_DOUBLE_PRECISION,&
          MPI_SUM,0,mpicomm,err)
     rismthermo%mpi_buffer = rismthermo%tmpi_buffer
#  endif /*USE_MPI_IN_PLACE*/
     if(err/=0) call rism_report_error("RISM_THERMO: MPI_REDUCE failed.")
#endif
  call rism_timer_stop(timer_write)
end subroutine rism_solvdist_thermo_calc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!The format for NAB programs will consist of one line for each
!!!catgory of data: energy, volume and excess.  Each is preceeded with
!!!"rism_[cat]:" where "cat" is a abreviation of the
!!!category.
!!!IN:
!!!   description : if .true., output a description of the table that will
!!!                  be printed but do not print out any values
!!!   soluPot     : total and component potential energy of the solute.
!!!                 Expected order is: total, LJ, elec, bond, angle, dihedral,
!!!                 H-bond, LJ-14, elec-14, restraints, 3D-RISM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rism_thermo_print(description,soluPot)
  use amber_rism_interface
  use constants, only : kb,COULOMB_CONST_E
  implicit none
#ifdef MPI
  include "mpif.h"
#endif /*MPI*/
  logical*4, intent(in) :: description
#ifdef SANDER
  _REAL_, intent(in) :: soluPot(25)
#else
  _REAL_, intent(in) :: soluPot(11)
#endif
  !catFmt    : Format for category (calculation type, e.g. excess chemical potential)
  !catbarFmt : Format for category (calculation type, e.g. excess chemical potential) with comment bar and units
  !titleFmt  : Format for column headings
  !valFmt    : Format for floating point values
  character(len=32) :: catFmt="(a15)", catbarFMT="('|',a15,' ',a10)", titleFmt="(a17)", valFmt='(1p,1x,e16.8e3)' 
  !whtspc : long string of whitespace that can be used to effect a
  !         left-justified string.  Otherwise strings are right-justified.
  !         Simply concatenate this to the end of the string you wish
  !         left-justified.
  character(len=64) :: whtspc
  integer :: iv, il, iu, err

  call rism_timer_start(timer_write)

  write(whtspc,'(a64)')" "

  if(description)then
     if(mpirank==0)then

        !thermodynamics key

        !free energy-type properties
        write(outunit,'(a)') "|3D-RISM thermodynamic data key"
        call flush(outunit)
        write(outunit,catbarFmt,advance='no') 'solute_epot'//whtspc, '[kcal/mol]'//whtspc
        write(outunit,'(12'//titleFmt//')') 'Total', 'LJ', 'Coulomb', 'Bond', &
             'Angle', 'Dihedral', 'H-Bond', 'LJ-14', 'Coulomb-14', 'Restraints', '3D-RISM'

       write(outunit,catbarFmt,advance='no') 'rism_exchem'//whtspc, '[kcal/mol]'//whtspc
        write(outunit,titleFmt,advance='no') 'Total'
        do iv=1, rism_3d%solv%natom
           write(outunit,titleFmt,advance='no') 'ExChem_'//trim(rism_3d%solv%atomname(iv))
        end do
        write(outunit,'(a)')

        write(outunit,catbarFmt,advance='no') 'rism_exchGF'//whtspc, '[kcal/mol]'//whtspc
        write(outunit,titleFmt,advance='no') 'Total'
        do iv=1, rism_3d%solv%natom
           write(outunit,titleFmt,advance='no') 'ExChem_GF_'//trim(rism_3d%solv%atomname(iv))
        end do
        write(outunit,'(a)')

        write(outunit,catbarFmt,advance='no') 'rism_potUV'//whtspc, '[kcal/mol]'//whtspc
        write(outunit,titleFmt,advance='no') 'Total'
        do iv=1, rism_3d%solv%natom
           write(outunit,titleFmt,advance='no') 'UV_potential_'//trim(rism_3d%solv%atomname(iv))
        end do
        write(outunit,'(a)')


        !polar decomposition
        if(rismprm%polarDecomp==1)then
           write(outunit,catbarFmt,advance='no') 'rism_polar'//whtspc, '[kcal/mol]'//whtspc
           write(outunit,titleFmt,advance='no') 'Total'
           do iv=1, rism_3d%solv%natom
              write(outunit,titleFmt,advance='no') 'polar_'//trim(rism_3d%solv%atomname(iv))
           end do
           write(outunit,'(a)')
           
           write(outunit,catbarFmt,advance='no') 'rism_apolar'//whtspc, '[kcal/mol]'//whtspc
           write(outunit,titleFmt,advance='no') 'Total'
           do iv=1, rism_3d%solv%natom
              write(outunit,titleFmt,advance='no') 'apolar_'//trim(rism_3d%solv%atomname(iv))
           end do
           write(outunit,'(a)')
           
           write(outunit,catbarFmt,advance='no') 'rism_polGF'//whtspc, '[kcal/mol]'//whtspc
           write(outunit,titleFmt,advance='no') 'Total'
           do iv=1, rism_3d%solv%natom
              write(outunit,titleFmt,advance='no') 'polarGF_'//trim(rism_3d%solv%atomname(iv))
           end do
           write(outunit,'(a)')
           
           write(outunit,catbarFmt,advance='no') 'rism_apolGF'//whtspc, '[kcal/mol]'//whtspc
           write(outunit,titleFmt,advance='no') 'Total'
           do iv=1, rism_3d%solv%natom
              write(outunit,titleFmt,advance='no') 'apolarGF_'//trim(rism_3d%solv%atomname(iv))
           end do
           write(outunit,'(a)')

           !solute-solvent potential energy 
           write(outunit,catbarFmt,advance='no') 'rism_pol_potuv'//whtspc, '[kcal/mol]'//whtspc
           write(outunit,titleFmt,advance='no') 'Total'
           do iv=1, rism_3d%solv%natom
              write(outunit,titleFmt,advance='no') 'polar_UV_pot_'//trim(rism_3d%solv%atomname(iv))
           end do
           write(outunit,'(a)')
           
           write(outunit,catbarFmt,advance='no') 'rism_apol_potuv'//whtspc, '[kcal/mol]'//whtspc
           write(outunit,titleFmt,advance='no') 'Total'
           do iv=1, rism_3d%solv%natom
              write(outunit,titleFmt,advance='no') 'apolar_UV_pot_'//trim(rism_3d%solv%atomname(iv))
           end do
           write(outunit,'(a)')
        end if

        !Thermodynamic properties not related to free energy
        write(outunit,catbarFmt,advance='no') 'rism_volume'//whtspc, '[A^ 3]'//whtspc
        write(outunit,titleFmt,advance='no') 'PMV'
        write(outunit,'(a)')

        write(outunit,catbarFmt,advance='no') 'rism_exNumb'//whtspc, '[#]'//whtspc
        write(outunit,titleFmt,advance='no') ""
        do iv=1, rism_3d%solv%natom
           write(outunit,titleFmt,advance='no') 'ExNum_'//trim(rism_3d%solv%atomname(iv))
        end do
        write(outunit,'(a)')

        write(outunit,catbarFmt,advance='no') 'rism_exChrg'//whtspc, '[e]'//whtspc
        write(outunit,titleFmt,advance='no') 'Total'
        do iv=1, rism_3d%solv%natom
           write(outunit,titleFmt,advance='no') 'ExChg_'//trim(rism_3d%solv%atomname(iv))
        end do
        write(outunit,'(a)')

     endif
  else
     !DATA: free energy-based properties
     if(mpirank==0 .and. associated(rismthermo%mpi_buffer))then
        write(outunit,catFmt,advance='no') 'solute_epot'//whtspc
#ifdef SANDER
        !cast the SANDER array into the order of the NAB array
        write(outunit,'(12'//trim(valFmt)//')') (/soluPot(1),soluPot(2),&
             soluPot(3),soluPot(5),soluPot(6), soluPot(7), soluPot(13), &
             soluPot(8),soluPot(9),soluPot(10),soluPot(24)/)
#else
        write(outunit,'(12'//trim(valFmt)//')') soluPot
#endif
        write(outunit,catFmt,advance='no') 'rism_exchem'//whtspc
        write(outunit,valFmt,advance='no') sum(rismthermo%exchem)
        do iv=1, rism_3d%solv%natom
           write(outunit,valFmt,advance='no') rismthermo%exchem(iv)
        end do
        write(outunit,'(a)')

        write(outunit,catFmt,advance='no') 'rism_exchGF'//whtspc
        write(outunit,valFmt,advance='no') sum(rismthermo%exchemGF)
        do iv=1, rism_3d%solv%natom
           write(outunit,valFmt,advance='no') rismthermo%exchemGF(iv)
        end do
        write(outunit,'(a)')

        write(outunit,catFmt,advance='no') 'rism_potUV'//whtspc
        write(outunit,valFmt,advance='no') sum(rismthermo%potUV)
        do iv=1, rism_3d%solv%natom
           write(outunit,valFmt,advance='no') rismthermo%potUV(iv)
        end do
        write(outunit,'(a)')

        !polar decomposition
        if(rismprm%polarDecomp==1)then
           !these values require MPI reductions before they can be calculated
           rismthermo%pol_exchem = rismthermo%exchem - rismthermo%apol_exchem
           rismthermo%pol_exchemGF = rismthermo%exchemGF - rismthermo%apol_exchemGF
           rismthermo%pol_potUV = rismthermo%potUV - rismthermo%apol_potUV
           write(outunit,catFmt,advance='no') 'rism_polar'//whtspc
           write(outunit,valFmt,advance='no') sum(rismthermo%pol_exchem)
           do iv=1, rism_3d%solv%natom
              write(outunit,valFmt,advance='no') rismthermo%pol_exchem(iv)
           end do
           write(outunit,'(a)')
           
           write(outunit,catFmt,advance='no') 'rism_apolar'//whtspc
           write(outunit,valFmt,advance='no') sum(rismthermo%apol_exchem)
           do iv=1, rism_3d%solv%natom
              write(outunit,valFmt,advance='no') rismthermo%apol_exchem(iv)
           end do
           write(outunit,'(a)')
           
           write(outunit,catFmt,advance='no') 'rism_polGF'//whtspc
           write(outunit,valFmt,advance='no') sum(rismthermo%pol_exchemGF)
           do iv=1, rism_3d%solv%natom
              write(outunit,valFmt,advance='no') rismthermo%pol_exchemGF(iv)
           end do
           write(outunit,'(a)')
           
           write(outunit,catFmt,advance='no') 'rism_apolGF'//whtspc
           write(outunit,valFmt,advance='no') sum(rismthermo%apol_exchemGF)
           do iv=1, rism_3d%solv%natom
              write(outunit,valFmt,advance='no') rismthermo%apol_exchemGF(iv)
           end do
           write(outunit,'(a)')

           write(outunit,catFmt,advance='no') 'rism_pol_potUV'//whtspc
           write(outunit,valFmt,advance='no') sum(rismthermo%pol_potuv)
           do iv=1, rism_3d%solv%natom
              write(outunit,valFmt,advance='no') rismthermo%pol_potuv(iv)
           end do
           write(outunit,'(a)')

           write(outunit,catFmt,advance='no') 'rism_apol_potUV'//whtspc
           write(outunit,valFmt,advance='no') sum(rismthermo%apol_potuv)
           do iv=1, rism_3d%solv%natom
              write(outunit,valFmt,advance='no') rismthermo%apol_potuv(iv)
           end do
           write(outunit,'(a)')

        end if
        !non-free energy-based properties
        write(outunit,catFmt,advance='no') 'rism_volume'//whtspc
        write(outunit,valFmt,advance='no') rismthermo%pmv
        write(outunit,'(a)')

        write(outunit,catFmt,advance='no') 'rism_exNumb'//whtspc
        write(outunit,titleFmt,advance='no') ""
        do iv=1, rism_3d%solv%natom
           write(outunit,valFmt,advance='no') rismthermo%exNum(iv)
        end do
        write(outunit,'(a)')

        write(outunit,catFmt,advance='no') 'rism_exChrg'//whtspc
        write(outunit,valFmt,advance='no') sum(rismthermo%exNum*rism_3d%solv%charge)&
             *sqrt((KB *rism_3d%solv%temperature)/COULOMB_CONST_E)
        do iv=1, rism_3d%solv%natom
           write(outunit,valFmt,advance='no') rismthermo%exNum(iv)*rism_3d%solv%charge(iv)&
                *sqrt((KB *rism_3d%solv%temperature)/COULOMB_CONST_E)
        end do
        write(outunit,'(a)')

     end if
  end if
  call flush(outunit)
  call rism_timer_stop(timer_write)
end subroutine rism_thermo_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Prints out the heirarchical timer summary
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rism_printTimer()
  use amber_rism_interface
  implicit none
  call rism_timer_start(timer)
  call rism_timer_summary(timer,'|',outunit,mpicomm)
  call rism_timer_stop(timer)
end subroutine rism_printTimer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Finalizes all of the 3D-RISM objects and frees memory.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rism_finalize()
  use amber_rism_interface
  use fce_c, only : fce_destroy
  use safemem
  implicit none
  integer :: err
  integer*8 :: memstats(10)
  call rism_timer_destroy(timer)
  if(rismprm%irism==1)then
     call rism3d_destroy(rism_3d)
     call fce_destroy(fce_o)

     if(safemem_dealloc(ff)/=0)&
          call rism_report_error("Deallocation in Amber-RISM interface failed")
#if defined(RISM_CRDINTERP)     
     if(safemem_dealloc(ratu_fce)/=0)&
          call rism_report_error("Deallocation in Amber-RISM interface failed")
#endif /*RISM_CRDINTER*/

     if(safemem_dealloc(closurelist) /= 0)&
          call rism_report_error("Deallocation in Amber-RISM interface failed")

     !Thermodynamics memory
     if(safemem_dealloc(rismthermo%mpi_buffer) /=0 )&
          call rism_report_error("Dealloc failed in rism_thermo")
#if !defined( USE_MPI_IN_PLACE) && defined( MPI )
     if(safemem_dealloc(rismthermo%tmpi_buffer) /=0 )&
          call rism_report_error("Dealloc failed in rism_thermo")
#endif
     if(safemem_dealloc(rismthermo%pol_exchem) /= 0)&
          call rism_report_error("Dealloc failed in rism_thermo")
     if(safemem_dealloc(rismthermo%pol_exchemGF) /= 0)&
          call rism_report_error("Dealloc failed in rism_thermo")
     if(safemem_dealloc(rismthermo%pol_exchemGF) /= 0)&
          call rism_report_error("Dealloc failed in rism_thermo")
     if(safemem_dealloc(rismthermo%pol_solvEne) /= 0)&
          call rism_report_error("Dealloc failed in rism_thermo")
     if(safemem_dealloc(rismthermo%pol_potUV) /= 0)&
          call rism_report_error("Dealloc failed in rism_thermo")

!!$     memstats = memStatus()
!!$     write(0,'(a)') "3D-RISM memory allocation summary"
!!$     write(0,'(a)') "Type         Current         Maximum"
!!$     write(0,'(a,i12,a,f12.5,a)') "Integer  ",memstats(1)," B ",&
!!$          dble(memstats(6))/BYTES_PER_GB," GB"
!!$     write(0,'(a,i12,a,f12.5,a)') "Real     ",memstats(2)," B ",&
!!$          dble(memstats(7))/BYTES_PER_GB," GB"
!!$     write(0,'(a,i12,a,f12.5,a)') "Logical  ",memstats(3)," B ",&
!!$          dble(memstats(8))/BYTES_PER_GB," GB"
!!$     write(0,'(a,i12,a,f12.5,a)') "Character",memstats(4)," B ",&
!!$          dble(memstats(9))/BYTES_PER_GB," GB"
!!$     write(0,'(a)') "---------------------------------------"
!!$     write(0,'(a,i12,a,f12.5,a)') "Total    ",memstats(5)," B ",&
!!$          dble(memstats(10))/BYTES_PER_GB," GB"
  end if
end subroutine rism_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Prints the maximum amount of memory allocated at any one time so far in the
!!!run
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rism_max_memory()
  use amber_rism_interface
  use safemem
  implicit none
#ifdef MPI
  include "mpif.h"
#endif /*MPI*/
  integer*8 :: memstats(10), tmemstats(10)
  integer :: err, irank
  outunit = rism_report_getMUnit()
  memstats = memStatus()
#ifdef MPI
#  ifdef USE_MPI_IN_PLACE
     if(mpirank==0)then
        call MPI_REDUCE(MPI_IN_PLACE, memstats, ubound(memstats,1),MPI_INTEGER8,&
             MPI_SUM,0,mpicomm,err)
     else
        call MPI_REDUCE(memstats, memstats, ubound(memstats,1),MPI_INTEGER8,&
             MPI_SUM,0,mpicomm,err)
     end if
#  else /*USE_MPI_IN_PLACE*/
     call MPI_REDUCE(memstats, tmemstats, ubound(memstats,1),MPI_INTEGER8,&
          MPI_SUM,0,mpicomm,err)
     memstats = tmemstats
#  endif /*USE_MPI_IN_PLACE*/
     if(err/=0) call rism_report_warn("RISM_MAX_MEMORY: MPI_REDUCE failed.")
#endif
  if(mpirank==0)then
     write(outunit,'(a)')
     write(outunit,'(a)') "|3D-RISM memory allocation summary"
     write(outunit,'(a)') "|Type          Maximum"
     write(outunit,'(a,f12.5,a)') "|Integer  ",&
          dble(memstats(6))/BYTES_PER_GB," GB"
     write(outunit,'(a,f12.5,a)') "|Real     ",&
          dble(memstats(7))/BYTES_PER_GB," GB"
     write(outunit,'(a,f12.5,a)') "|Logical  ",&
          dble(memstats(8))/BYTES_PER_GB," GB"
     write(outunit,'(a,f12.5,a)') "|Character",&
          dble(memstats(9))/BYTES_PER_GB," GB"
     write(outunit,'(a)') "|------------------------"
     write(outunit,'(a,f12.5,a)') "|Total    ",&
          dble(memstats(10))/BYTES_PER_GB," GB"
  end if
end subroutine rism_max_memory

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! I/O: performs RISM related I/O for files that only deal with RISM data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Provides access to writeSolvDist for non-Fortran code
!!!IN:
!!!   step :: step number used as a suffix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rism_writeSolvDistC(step)
  use amber_rism_interface
  implicit none
  integer, intent(in) :: step
!  write(6,*) "writeSolvDistC stub"
  call rism_writeSolvDistF(rism_3d,step)
end subroutine rism_writeSolvDistC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Outputs current solvent distributions to file.  Each distribution is written
!!!in a separate file with the step number before the suffix.  File names are taken 
!!!from guvfile, huvfile and cuvfile.
!!!IN:
!!!   this :: 3D-RISM object
!!!   step :: step number used as a suffix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rism_writeSolvDistF(this,step)
  use constants, only : COULOMB_CONST_E, KB
  use amber_rism_interface
  use rism3d_opendx
  use rism3d_xyzv
  use safemem
  implicit none
#if defined(MPI)
  include 'mpif.h'
#endif
  type(rism3d),intent(in) :: this
  integer, intent(in) :: step
  integer :: iv, ivv, nsolv, igx,igy,igz,ios, i,j,k
  character(len=16) :: cstep
  character(len=64) :: suffix
  _REAL_, pointer::work(:,:,:)=>NULL()
#ifdef MPI
  integer :: err
#endif /*MPI*/  
  call rism_timer_start(timer_write)

#ifdef RISM_DEBUG
  write(6,*) "printrism"
#endif
#ifdef MPI
  if(len_trim(guvfile) /= 0 .or. len_trim(huvfile) /= 0)then
     work => safemem_realloc(work, this%grid%ngr(1), this%grid%ngr(2), this%grid%nr(3))
  end if
#endif /*MPI*/
  !................... outputting Guv and Cuv profiles ...................
  do iv=1,this%solv%natom
     write(cstep,'(i16)') step
     cstep = adjustl(cstep)
     suffix = '.'//trim(rism_3d%solv%atomname(iv))//'.'//trim(cstep)
     if(volfmt.eq.'dx')then
        suffix = trim(suffix)//'.dx'
     else
        suffix = trim(suffix)//'.xyzv'
     end if

     !GUV
     if (len_trim(guvfile) /= 0)  then
#  if defined(MPI)
        do k=1,this%grid%nr(3)
           do j=1,this%grid%ngr(2)
              do i=1,this%grid%ngr(1)
                 work(i,j,k) = this%guv(i+(j-1)*(this%grid%ngr(1)+2)+(k-1)*(this%grid%ngr(1)+2)*this%grid%ngr(2),iv)
              end do
           end do
        end do
        if(volfmt .eq. 'dx')then
           call writeDX(trim(guvfile)//suffix,&
                work,&
                this%grid%boxlen,this%grid%nr,this%grid%ngr(3),this%ratucm,&
                mpirank,mpisize,mpicomm)
        else
           call rism3d_xyzv_write(trim(guvfile)//suffix,work,this%grid%grdspc,&
                this%grid%nr,(this%ratucm-this%grid%boxlen/2d0),&
                mpirank,mpisize,mpicomm)
        end if
#  else
        if(volfmt .eq. 'dx')then
           call writeDX(trim(guvfile)//suffix,this%guv(:,iv),this%grid%boxlen,&
                this%grid%nr,this%grid%ngr(3),this%ratucm)
        else
           call rism3d_xyzv_write(trim(guvfile)//suffix,this%guv(:,iv),this%grid%grdspc,&
                this%grid%nr,(this%ratucm-this%grid%boxlen/2d0))
        end if
#  endif /*defined(MPI)*/
     endif

     !HUV
     if (len_trim(huvfile) /= 0)  then
#  if defined(MPI)
        do k=1,this%grid%nr(3)
           do j=1,this%grid%ngr(2)
              do i=1,this%grid%ngr(1)
                 work(i,j,k) =this%huv(i+(j-1)*(this%grid%ngr(1)+2)+(k-1)*(this%grid%ngr(1)+2)*this%grid%ngr(2),iv)
              end do
           end do
        end do
        if(volfmt .eq. 'dx')then
           call writeDX(trim(huvfile)//suffix,&
                work,&
                this%grid%boxlen,this%grid%nr,this%grid%ngr(3),this%ratucm,&
                mpirank,mpisize,mpicomm)
        else
           call rism3d_xyzv_write(trim(huvfile)//suffix,work,this%grid%grdspc,&
                this%grid%nr,(this%ratucm-this%grid%boxlen/2d0),&
                mpirank,mpisize,mpicomm)
        end if
#  else
        if(volfmt .eq. 'dx')then
           call writeDX(trim(huvfile)//suffix,this%huv(:,iv),this%grid%boxlen,&
                this%grid%nr,this%grid%ngr(3),this%ratucm)
        else
           call rism3d_xyzv_write(trim(huvfile)//suffix,this%huv(:,iv),this%grid%grdspc,&
                this%grid%nr,(this%ratucm-this%grid%boxlen/2d0))
        end if
#  endif /*defined(MPI)*/
     endif

     !CUV
     if (len_trim(cuvfile) /= 0)  then
        if(volfmt .eq. 'dx')then
           call writeDX(trim(cuvfile)//suffix,this%cuv(:,:,:,iv),&
                this%grid%boxlen,this%grid%nr,this%grid%ngr(3),this%ratucm,&
                mpirank,mpisize,mpicomm)
        else
           call rism3d_xyzv_write(trim(cuvfile)//suffix,this%cuv(:,:,:,iv),this%grid%grdspc,&
                this%grid%nr,(this%ratucm-this%grid%boxlen/2d0),&
                mpirank,mpisize,mpicomm)
        end if           
     endif

     !UUV
     if (len_trim(uuvfile) /= 0)  then
        if(volfmt .eq. 'dx')then
           call writeDX(trim(uuvfile)//suffix,this%pot%uuv(:,:,:,iv),&
                this%grid%boxlen,this%grid%nr,this%grid%ngr(3),this%ratucm,&
                mpirank,mpisize,mpicomm)
        else
           call rism3d_xyzv_write(trim(uuvfile)//suffix,this%pot%uuv(:,:,:,iv),this%grid%grdspc,&
                this%grid%nr,(this%ratucm-this%grid%boxlen/2d0),&
                mpirank,mpisize,mpicomm)
        end if           
     endif
  end do
  !asymptotics files.  There are four but it is a simple loop
  if (len_trim(asympfile) /= 0 .and. this%solu%charged)  then
     write(cstep,'(i16)') step
     cstep = adjustl(cstep)
     suffix = '.'//trim(cstep)
     if(volfmt .eq. 'dx')then
        suffix = trim(suffix)//'.dx'
        !h(r)
        if(this%solv%ionic)then
           call writeDX(trim(asympfile)//"hr"//suffix,this%pot%asymhr,this%grid%boxlen,this%grid%nr,this%grid%ngr(3),this%ratucm,&
                mpirank,mpisize,mpicomm)
        else
           call rism_report_warn("Not writing long-range asymptotics of h(r); not used for non-ionic solvent.")
        end if
        !c(r)
        call writeDX(trim(asympfile)//"cr"//suffix,this%pot%asymcr,this%grid%boxlen,this%grid%nr,this%grid%ngr(3),this%ratucm,&
             mpirank,mpisize,mpicomm)
     else
        suffix = trim(suffix)//'.xyzv'
        !h(r)
        if(this%solv%ionic)then
             call rism3d_xyzv_write(trim(asympfile)//"hr"//suffix,this%pot%asymhr,this%grid%grdspc,&
             this%grid%nr,(this%ratucm-this%grid%boxlen/2d0),&
             mpirank,mpisize,mpicomm)
          else
             call rism_report_warn("Not writing long-range asymptotics of h(r); not used for non-ionic solvent.")
          end if
        !c(r)
        call rism3d_xyzv_write(trim(asympfile)//"cr"//suffix,this%pot%asymcr,this%grid%grdspc,&
             this%grid%nr,(this%ratucm-this%grid%boxlen/2d0),&
             mpirank,mpisize,mpicomm)
     end if
     !h(k)
     !in (k) space the array is transposed (distributed over the y axis) so writeDX doesn't
     !know how to handle it
!     call writeDX(unit,this%pot%asymhk,this%grid%boxlen,uboundthis%grid%nr,this%grid%ngr(3),this%ratucm,&
!          mpirank,mpisize,mpicomm)
!     call writeDX(trim(asympfile)//"hk"//suffix,this%pot%asymhk,&
!          this%grid%boxlen,(/this%grid%ngr(1)+2,this%ny_local, this,this%,this%grid%ngr(2),this%ratucm)
     !c(k)
     !in (k) space the array is transposed (distributed over the y axis) so writeDX doesn't
     !know how to handle it
!     call writeDX(unit,this%pot%asymck,this%grid%boxlen,uboundthis%grid%nr,this%grid%ngr(3),this%ratucm,&
!          mpirank,mpisize,mpicomm)
!     call writeDX(trim(asympfile)//"ck"//suffix,this%pot%asymck,&
!          this%grid%boxlen,ubound(this%pot%asymck),this%grid%ngr(2),this%ratucm)
  elseif(.not.this%solu%charged)then
     call rism_report_warn("Not writing long-range asymptotics; not used for uncharged solute.")
  endif
  !................... outputting charge density distribution ...................
  if (len_trim(quvfile) /= 0 .or. len_trim(chgDistFile) /= 0)  then
     write(cstep,'(i16)') step
     cstep = adjustl(cstep)
     suffix = '.'//trim(cstep)
     if(volfmt.eq.'dx')then
        suffix = trim(suffix)//'.dx'
     else
        suffix = trim(suffix)//'.xyzv'
     end if
     work => safemem_realloc(work, this%grid%ngr(1), this%grid%ngr(2), this%grid%nr(3))
     work=0
     !sum the contributions from each solvent type at each grid point and convert units to [e/A^3]
     do iv = 1, this%solv%natom
#  if defined(MPI)
        do k=1,this%grid%nr(3)
           do j=1,this%grid%ngr(2)
              do i=1,this%grid%ngr(1)
                 work(i,j,k) = work(i,j,k) &
                      + this%guv(i+(j-1)*(this%grid%ngr(1)+2)+(k-1)*(this%grid%ngr(1)+2)*this%grid%ngr(2),iv) &
                      *sqrt((KB *this%solv%temperature)/COULOMB_CONST_E)&
                      *this%solv%charge(iv)*this%solv%rho(iv)
              end do
           end do
        end do
#  else
        call DAXPY(this%grid%nrtotal,sqrt((KB *this%solv%temperature)/COULOMB_CONST_E)&
             *this%solv%charge(iv)*this%solv%rho(iv),this%guv(1,iv),1,work,1)
#  endif /*defined(MPI)*/
     end do
     if(len_trim(quvfile) /= 0)then
        if(volfmt .eq. 'dx')then
           call writeDX(trim(quvfile)//suffix,work,this%grid%boxlen,&
                this%grid%nr,this%grid%ngr(3),this%ratucm,&
                mpirank,mpisize,mpicomm)
        else
           call rism3d_xyzv_write(trim(quvfile)//suffix,work,this%grid%grdspc,&
                this%grid%nr,(this%ratucm-this%grid%boxlen/2d0),&
                mpirank,mpisize,mpicomm)
        end if
     end if
     if(len_trim(chgDistFile) /= 0)  then
        call DSCAL(this%grid%nrtotal,this%grid%voxel,work,1)
        if(volfmt .eq. 'dx')then
           call writeDX(trim(chgDistfile)//suffix,work,this%grid%boxlen,&
                this%grid%nr,this%grid%ngr(3),this%ratucm,&
                mpirank,mpisize,mpicomm)
        else
           call rism3d_xyzv_write(trim(chgDistfile)//suffix,work,this%grid%grdspc,&
                this%grid%nr,(this%ratucm-this%grid%boxlen/2d0),&
                mpirank,mpisize,mpicomm)
        end if
     end if
  endif

  if(safemem_dealloc(work)/=0) call rism_report_error("RISM_WRITESOLVEDIST: failed to deallocate WORK")
  call rism_timer_stop(timer_write)
end subroutine rism_writeSolvDistF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!INTERPOLATION RESTART FILE I/O
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if 0
#if defined(RISM_CRDINTERP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Writes the interpolation restart file
!!IN:
!!   ratu    :: the position of each solute atom for each step 
!!   frc     :: the force of each solute atom for each step
!!   this%solu%natom    :: number of solute atoms
!!   nstep   :: number of steps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FIX - modify this to use the NetCDF format
subroutine rismrestrtwrit(this,ratu,frc,nstep)
  use amber_rism_interface
  use rism_util, only : freeUnit
  implicit none
#include "files.h"
  type(rism3d) :: this
  _REAL_,intent(in) :: ratu(3,this%solu%natom,nstep),frc(3,this%solu%natom,nstep)
  integer,intent(in) :: nstep
  integer :: istep
  integer :: unit
  !this is performed by the master node only:
  if(mpirank == 0) then
#ifdef RISM_DEBUG
     write(6,*) 'entering rismrestrtwrit'
     call flush(6)
#endif /*RISM_DEBUG*/
     unit = freeUnit()
     if (len_trim(rismcrdfil) /= 0)  then
        open(unit=unit,file=rismcrdfil,status='new', form='FORMATTED',iostat=stat)
        write(unit,'(i8)') nstep
        do istep = 1,nstep
           call corpac(ratu(1:3,1:this%solu%natom,istep),1,this%solu%natom*3,unit,.true.)
        end do
        close(unit)
     end if
     if (len_trim(rismfrcfil) /= 0)  then
        open(unit=unit,file=rismfrcfil,status='new', form='FORMATTED',iostat=stat)
        write(unit,'(i8)') nstep
        do istep = 1,nstep
           call corpac(frc(1:3,1:this%solu%natom,istep),1,this%solu%natom*3,unit,.true.)
        end do
        close(unit)
     end if
#ifdef RISM_DEBUG
     write(6,*) 'done rismrestrtwrit'
     call flush(6)
#endif /*RISM_DEBUG*/
  end if
end subroutine rismrestrtwrit
!!endfix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Reads in the interpolation restart file
!!IN:
!!   ratu    :: the position of each atom for each step read in
!!   frc     :: the force of each atom for each step read in
!!   this%solu%natom    :: number of solute atoms
!!   nstep   :: maximum number of steps to read
!!   nsample :: will hold the total number of steps read in
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FIX - modify this to use the NetCDF format
subroutine rismrestrtread(this,ratu,frc,nstep,nsample)
  use amber_rism_interface
  use rism_util, only : freeUnit
  implicit none
#include "files.h"
#if defined(MPI)
  include 'mpif.h'
#endif /*defined(MPI)*/
  type(rism3d) :: this
  _REAL_,intent(out) :: ratu(3,this%solu%natom,nstep),frc(3,this%solu%natom,nstep)
  integer,intent(in) :: nstep
  integer,intent(out) :: nsample
  integer :: istep,iatu,csteps,fsteps,err
  integer :: unit
  !this is performed by the master node only:
  if(mpirank == 0) then
#ifdef RISM_DEBUG
     write(6,*) "Reading interpolation restart files..."
     write(6,*) 'entering rismrestrtread'
     call flush(6)
#endif /*RISM_DEBUG*/
     unit=freeeUnit()
     if (len_trim(rismcrdrstfil) /= 0)  then
        open(unit=unit,file=rismcrdfil,status='old', form='FORMATTED',iostat=stat)
        read(unit,'(i12)') csteps
        !discard the extra samples at the beginning
#ifdef RISM_DEBUG
        write(6,*) csteps, nstep, csteps-nstep
        call flush(6)
#endif /*RISM_DEBUG*/
        do istep = 1,csteps-nstep
           !            read(unit,'(10f8.3)') (ratu(1:3,iatu,1),iatu=1,this%solu%natom)
           read(unit,'(10f21.16)') (ratu(1:3,iatu,1),iatu=1,this%solu%natom)
        end do
        do istep = 1,min(csteps,nstep)
           !            read(unit,'(10f8.3)') (ratu(1:3,iatu,istep),iatu=1,this%solu%natom)
           read(unit,'(10f21.16)') (ratu(1:3,iatu,istep),iatu=1,this%solu%natom)
        end do
        close(unit)
     endif
     if (len_trim(rismfrcrstfil) /= 0)  then
        open(unit=unit,file=rismfrcfil,status='old', form='FORMATTED',iostat=stat)
        read(unit,'(i12)') fsteps
        !discard the extra samples at the beginning
        do istep = 1,fsteps-nstep
           !            read(unit,'(10f8.3)') (frc(1:3,iatu,1),iatu=1,this%solu%natom)
           read(unit,'(10f21.16)') (frc(1:3,iatu,1),iatu=1,this%solu%natom)
        end do
        do istep = 1,min(fsteps,nstep)
           !            read(unit,'(10f8.3)') (frc(1:3,iatu,istep),iatu=1,this%solu%natom)
           read(unit,'(10f21.16)') (frc(1:3,iatu,istep),iatu=1,this%solu%natom)
        end do
        close(unit)
        if(csteps /= fsteps)then
           call rism_report_error('RISM interpolation restart files have different numbers of entries.')
        end if
        nsample = min(csteps,nstep)
     end if
#ifdef RISM_DEBUG
     write(6,*) nsample," samples read"
     write(6,*) 'done rismrestrtread'
     call flush(6)
#endif /*RISM_DEBUG*/
  end if
  !
  !broadcast results to the other processes: nsamples, ratu and frc
  !
!!$#if defined(MPI)
!!$      CALL MPI_BCAST(nsample,1,MPI_INTEGER,0,mpicomm,err)
!!$      CALL MPI_BCAST(ratu,3*this%solu%natom*nstep,MPI_DOUBLE_PRECISION,0,mpicomm,err)
!!$      CALL MPI_BCAST(frc,3*this%solu%natom*nstep,MPI_DOUBLE_PRECISION,0,mpicomm,err)
!!$#endif /*defined(MPI)*/

end subroutine rismrestrtread
#endif /*RISM_CRDINTERP*/
#endif  /* 0 */

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!PRIVATE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Broadcasts initialization information about the system from the master to all
!!!the other processes.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef MPI
subroutine rism_mpi_bcast(rank,size,comm)
  use amber_rism_interface

  implicit none
  include 'mpif.h'
  integer, intent(in) :: rank, size, comm
  integer :: err
  integer :: nclosure
  !Private subroutine so there should be no timer

  mpirank = rank
  mpisize = size
  mpicomm = comm

  !could be done by passing creating and passing an MPI derived type
  !but this is done once so it is not worth the effort
  call mpi_bcast(rismprm%irism,1,mpi_integer,0,mpicomm,err)
  if (rismprm%irism == 1) then
     !broadcast the entire rismprm object.  Not everything actually needs
     !to be transferred but there are not that many exceptions.  This is
     !could be done by passing creating and passing an MPI derived type.
     !However, this is only being done once so it is not really worth it.
     call mpi_bcast(rismprm%solvcut,1,mpi_double_precision,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast SOLVCUT")
     call mpi_bcast(rismprm%buffer,1,mpi_double_precision,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast BUFFER")
     call mpi_bcast(rismprm%grdspc,3,mpi_double_precision,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast GRDSPC")
     call mpi_bcast(rismprm%solvbox,3,mpi_double_precision,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast SOLVBOX")
     call mpi_bcast(rismprm%mdiis_del,1,mpi_double_precision,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast MDIIS_DEL")
     call mpi_bcast(rismprm%mdiis_restart,1,mpi_double_precision,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast MDIIS_RESTART")
     call mpi_bcast(rismprm%fcecut,1,mpi_double_precision,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast FCECUT")
     call mpi_bcast(rismprm%ng3,3,mpi_integer,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast NG3")
     call mpi_bcast(rismprm%polarDecomp,1,mpi_integer,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast POLARDECOMP")
     call mpi_bcast(rismprm%asympCorr,1,mpi_integer,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast ASYMPCORR")
     call mpi_bcast(rismprm%maxstep,1,mpi_integer,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast MAXSTEP")
     call mpi_bcast(rismprm%npropagate,1,mpi_integer,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast NPROPAGATE")
     call mpi_bcast(rismprm%centering,1,mpi_integer,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast CENTERING")
     call mpi_bcast(rismprm%zerofrc,1,mpi_integer,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast ZEROFRC")
     call mpi_bcast(rismprm%apply_rism_force,1,mpi_integer,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast APPLY_RISM_FORCE")
     call mpi_bcast(rismprm%rismnrespa,1,mpi_integer,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast RISMNRESPA")
     call mpi_bcast(rismprm%fcestride,1,mpi_integer,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast FCESTRIDE")
     call mpi_bcast(rismprm%fcenbasis,1,mpi_integer,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast FCENBASIS")
     call mpi_bcast(rismprm%fcecrd,1,mpi_integer,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast FCECRD")
     call mpi_bcast(rismprm%saveprogress,1,mpi_integer,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast SAVEPROGRESS")
     call mpi_bcast(rismprm%ntwrism,1,mpi_integer,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast NTWRISM")
     call mpi_bcast(rismprm%verbose,1,mpi_integer,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast VERBOSE")
     call mpi_bcast(rismprm%progress,1,mpi_integer,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast PROGRESS")

     if(mpirank==0)&
          nclosure=ubound(closurelist,1)
     call mpi_bcast(nclosure,1,mpi_integer,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast PROGRESS")
     if(mpirank/=0)then
        closurelist=>safemem_realloc(closurelist,len(closurelist),nclosure)
        tolerancelist=>safemem_realloc(tolerancelist,nclosure)
     end if
     call mpi_bcast(closurelist,len(closurelist)*nclosure,mpi_character,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast CLOSURE")
     call mpi_bcast(tolerancelist,nclosure,mpi_double_precision,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast TOLERANCE")

#ifdef SANDER
     call mpi_bcast(rismprm%write_thermo,1,mpi_integer,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast WRITE_THERMO")
#endif

     !these are not being used currently.
     call mpi_bcast(pa_orient,1,mpi_integer,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast PA_ORIENT")
     call mpi_bcast(rmsd_orient,1,mpi_integer,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast RMSD_ORIENT")

     !I/O
     !special output files that all nodes write to
     call mpi_bcast(guvfile,len(guvfile),mpi_character,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast GUVFILE")
     call mpi_bcast(huvfile,len(huvfile),mpi_character,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast HUVFILE")
     call mpi_bcast(cuvfile,len(cuvfile),mpi_character,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast CUVFILE")
     call mpi_bcast(uuvfile,len(uuvfile),mpi_character,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast UUVFILE")
     call mpi_bcast(asympfile,len(asympfile),mpi_character,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast ASYMPFILE")
     call mpi_bcast(quvfile,len(quvfile),mpi_character,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast QUVFILE")
     call mpi_bcast(chgdistfile,len(chgdistfile),mpi_character,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast CHGDISTFILE")
     call mpi_bcast(volfmt,len(volfmt),mpi_character,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("RISM3D interface: could not broadcast VOLFMT")

  end if
end subroutine rism_mpi_bcast
#endif /*MPI*/

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Sets default values for 3D-RISM paramters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine defaults()
  use amber_rism_interface
  implicit none

  closurelist => safemem_realloc(closurelist,len(closurelist),nclosuredefault)
  closurelist(2:)               = ''
  closurelist(1)                = 'KH'
  rismprm%closureOrder      = 1
  rismprm%asympCorr         = .true.
  rismprm%polarDecomp       = 0

  !solvation box
  rismprm%solvcut         = -1
  rismprm%buffer          = 14d0
  rismprm%grdspc          = 0.5d0
  rismprm%ng3             = -1
  rismprm%solvbox         = -1d0

  !convergence
  tolerancelist => safemem_realloc(tolerancelist,1)
  tolerancelist                 = 1d-5
  rismprm%mdiis_del         = 0.7d0
  rismprm%mdiis_nvec        = 5
  rismprm%mdiis_method      = 2
  rismprm%mdiis_restart     = 10d0
  rismprm%maxstep           = 10000
  rismprm%npropagate        = 5

  !imin = 1 (minimization)
  rismprm%centering          = 1
  rismprm%zerofrc          = 1

  !imin = 5 (trajectory analysis)
  rismprm%apply_rism_force = 1
  pa_orient        = 0
  rmsd_orient      = 0

  !imin = 0 (MD)
  rismprm%rismnrespa       = 1
#ifdef RISM_CRDINTERP
  rismprm%fcestride        = 0
  rismprm%fcecut           = 9999d0
  rismprm%fcenbasis        = 10
  rismprm%fcecrd           = 0
#endif

  !output
  rismprm%saveprogress     = 0
  rismprm%ntwrism          = -1
  rismprm%verbose          = 0
  rismprm%progress         = 1
  volfmt                   = 'dx'

#ifdef SANDER
  rismprm%write_thermo=1
#endif
end subroutine defaults

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Transfers RISM settings from a rismprm_t type (user) to the rismprm parameter set
!!!for the 3D-RISM calculation.  For each possible setting the new value is used
!!!IFF the value > -9999
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine update_param(user)
  use amber_rism_interface
  implicit none
  type(rismprm_t), intent(in) :: user
  rismprm%irism = user%irism
  !asympCorr is a logical, so we can't perform the test
  rismprm%asympCorr = user%asympCorr
  if(user%closureOrder > -9999) rismprm%closureOrder = user%closureOrder
  if(user%polarDecomp > -9999) rismprm%polarDecomp = user%polarDecomp
  
  if(user%buffer > -9999) rismprm%buffer = user%buffer
  if(user%solvcut > 0) then
     rismprm%solvcut = user%solvcut
  else
     rismprm%solvcut = rismprm%buffer     
  end if
  if(user%grdspc(1) > -9999) rismprm%grdspc = user%grdspc
  if(user%ng3(1) > -9999) rismprm%ng3 = user%ng3
  if(user%solvbox(1) > -9999) rismprm%solvbox = user%solvbox
  if(user%mdiis_del > -9999) rismprm%mdiis_del = user%mdiis_del
  if(user%mdiis_nvec > -9999) rismprm%mdiis_nvec = user%mdiis_nvec
  if(user%mdiis_method > -9999) rismprm%mdiis_method = user%mdiis_method
  if(user%mdiis_restart > -9999) rismprm%mdiis_restart = user%mdiis_restart
  if(user%maxstep > -9999) rismprm%maxstep = user%maxstep
  if(user%npropagate > -9999) rismprm%npropagate = user%npropagate
  if(user%centering > -9999) rismprm%centering = user%centering
  if(user%zerofrc > -9999) rismprm%zerofrc = user%zerofrc
  if(user%apply_rism_force > -9999) rismprm%apply_rism_force = user%apply_rism_force
  if(user%rismnrespa > -9999) rismprm%rismnrespa = user%rismnrespa
  if(user%fcestride > -9999) rismprm%fcestride = user%fcestride
  if(user%fcecut > -9999) rismprm%fcecut = user%fcecut
  if(user%fcenbasis > -9999) rismprm%fcenbasis = user%fcenbasis
  if(user%fcecrd > -9999) rismprm%fcecrd = user%fcecrd
  if(user%saveprogress > -9999) rismprm%saveprogress = user%saveprogress
  if(user%ntwrism > -9999) rismprm%ntwrism = user%ntwrism
  if(user%verbose > -9999) rismprm%verbose = user%verbose
  if(user%progress > -9999) rismprm%progress = user%progress
end subroutine update_param

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Reads the RISM namelist from the input file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_namelist(mdin_unit)
  use amber_rism_interface
  implicit none
  !mdin_unit :: unit number for mdin file
  character(len=8) :: closure(nclosuredefault)
  _REAL_ :: tolerance(nclosuredefault)
  integer :: closureOrder
  integer, intent(in) :: mdin_unit
  logical :: asympCorr
  integer :: polarDecomp
  _REAL_ :: solvcut
  _REAL_ :: buffer
  _REAL_ :: grdspc(3)
  integer ::  ng3(3)
  _REAL_ :: solvbox(3)
  _REAL_ :: mdiis_del
  integer :: mdiis_nvec
  integer :: mdiis_method
  _REAL_ :: mdiis_restart
  integer :: maxstep
  integer :: npropagate
  integer :: centering
  integer :: zerofrc
  integer :: apply_rism_force
  integer :: rismnrespa
  integer :: fcestride
  _REAL_ :: fcecut
  integer ::  fcenbasis
  integer ::  fcecrd
  integer :: saveprogress
  integer :: ntwrism
  integer :: verbose
  integer :: progress
#ifdef SANDER
  integer :: write_thermo
#endif
  namelist /rism/ &
       !closure
       closure, closureOrder,&
       !thermodynamics
       asympCorr, &
       !solvation box
       buffer, grdspc, solvcut, &
       ng3, solvbox,&
       !convergence
       tolerance, mdiis_del, mdiis_nvec, mdiis_method, mdiis_restart, maxstep, npropagate,&
       !minimization
       centering,zerofrc,&
       !imin=5
       apply_rism_force, pa_orient, rmsd_orient,&
       !md
       rismnrespa,&
#ifdef RISM_CRDINTERP
       fcestride, fcecut, fcenbasis, fcecrd, &
#endif /*RISM_CRDINTERP*/
  !output
#ifdef SANDER
       write_thermo, &
#endif
       saveprogress, ntwrism, verbose, progress, volfmt
  call flush(0)
  
  closure = closurelist
  tolerance = tolerancelist
  closureOrder = rismprm%closureOrder
  asympCorr = rismprm%asympCorr
  polarDecomp = rismprm%polarDecomp
  solvcut = rismprm%solvcut
  buffer = rismprm%buffer
  grdspc= rismprm%grdspc
  ng3 = rismprm%ng3
  solvbox = rismprm%solvbox
  mdiis_del = rismprm%mdiis_del
  mdiis_nvec = rismprm%mdiis_nvec
  mdiis_method = rismprm%mdiis_method
  mdiis_restart = rismprm%mdiis_restart
  maxstep = rismprm%maxstep
  npropagate = rismprm%npropagate
  centering = rismprm%centering
  zerofrc = rismprm%zerofrc
  apply_rism_force = rismprm%apply_rism_force
  rismnrespa = rismprm%rismnrespa
  fcestride = rismprm%fcestride
  fcecut = rismprm%fcecut
  fcenbasis = rismprm%fcenbasis
  fcecrd = rismprm%fcecrd
  saveprogress = rismprm%saveprogress
  ntwrism = rismprm%ntwrism
  verbose= rismprm%verbose
  progress = rismprm%progress
#ifdef SANDER
  write_thermo = rismprm%write_thermo
#endif

  !resize tolerance to the size of closure
  tolerancelist => safemem_realloc(tolerancelist,size(closurelist))
  tolerancelist(2:)=HUGE(1d0)

  rewind(mdin_unit)
  read(mdin_unit,nml=rism)
  
  closurelist = closure
  tolerancelist = tolerance
  rismprm%closureOrder = closureOrder
  rismprm%asympCorr=asympCorr
  rismprm%polarDecomp = polarDecomp
  !solvation box
  rismprm%buffer=buffer
  rismprm%grdspc=grdspc
  rismprm%solvcut=solvcut
  rismprm%ng3=ng3
  rismprm%solvbox=solvbox
  !convergence
  rismprm%mdiis_del=mdiis_del
  rismprm%mdiis_nvec=mdiis_nvec
  rismprm%mdiis_method=mdiis_method
  rismprm%mdiis_restart=mdiis_restart
  rismprm%maxstep=maxstep
  rismprm%npropagate=npropagate
  !minimization
  rismprm%centering=centering
  rismprm%zerofrc=zerofrc
  !imin=5
  rismprm%apply_rism_force=apply_rism_force
  pa_orient=pa_orient
  rmsd_orient=rmsd_orient
  !md
  rismprm%rismnrespa=rismnrespa
#ifdef RISM_CRDINTERP
  rismprm%fcestride=fcestride
  rismprm%fcecut=fcecut
  rismprm%fcenbasis=fcenbasis
  rismprm%fcecrd=fcecrd
#endif /*RISM_CRDINTERP*/
  !output
  rismprm%saveprogress=saveprogress
  rismprm%ntwrism=ntwrism
  rismprm%verbose=verbose
  rismprm%progress=progress
#ifdef SANDER
  rismprm%write_thermo = write_thermo
#endif

  !set the RISM cutoff if not set by the user
  if(rismprm%solvcut < 0) then
     rismprm%solvcut = rismprm%buffer
  end if

end subroutine read_namelist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Checks  user input to ensure that it is not completely crazy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine sanity_check()
  use amber_rism_interface
  use rism_util, only : caseup
  use array_util, only : array_index
  implicit none
  character(len=32) :: fmt
  !iclosure :: counter for closures
  integer :: iclosure
!!$  if(rismprm%asympCorr > 1 .or. rismprm%asympCorr < 0) then
!!$     call rism_report_error('(a,i2)','RISM asymptotic correction must be 0 or 1:  ',rismprm%asympCorr)
!!$  endif

  !ensure that the cutoff is set to a reasonable value.  This can
  !happen if buffer is < 0
  if(rismprm%solvcut <0) &
       call rism_report_error('solvcut must >= 0.')

  !ensure that solvbox and ng3 have been set if buffer < 0
  if(rismprm%buffer <0)then
     if(minval(rismprm%ng3) < 0) &
          call rism_report_error('if buffer < 0 ng3 (grid size) must be set.')
     if(minval(rismprm%solvbox) < 0) &
          call rism_report_error('if buffer < 0 solvbox must be set.')
  end if

  !ensure that an apropriate file format has been chosen for volumetric output
  if(.not.(volfmt .eq. "dx" .or. volfmt .eq. "xyzv"))then
     call rism_report_error("Only 'dx' and 'xyzv' volumetric data formats are supported")
  end if
  
  !make sure that a valid centering method is used and doesn't conflict with other options
  if(rismprm%centering >4 .or. rismprm%centering < -4)&
       call rism_report_error("Centering must be between -4 and 4")
  if((rismprm%centering >2 .or. rismprm%centering < -2) .and. rismprm%fcestride > 0)&
       call rism_report_error("CENTERING must be between -2 and 2 when FCESTRIDE > 0")

  !resize closure list to the appropriate size
  if(len_trim(closurelist(size(closurelist))) == 0)&
     closurelist => safemem_realloc(closurelist,len(closurelist),array_index(closurelist,'')-1)

  !resize and setup the tolerance list
  !1)get rid of extraneous values
  if(array_index(tolerancelist,HUGE(1d0))>0)&
       tolerancelist=>safemem_realloc(tolerancelist,&
       array_index(tolerancelist,HUGE(1d0))-1)
  !2)if there is only one closure, use only the last tolerance
  if(size(closurelist) == 1)then
     tolerancelist(1) = tolerancelist(size(tolerancelist))
     tolerancelist=>safemem_realloc(tolerancelist,1)
     
  !3)if there is one tolerance, the default for intermediate closures is 1
  elseif(size(tolerancelist) == 1)then
     tolerancelist=>safemem_realloc(tolerancelist,size(closurelist))
     tolerancelist(size(tolerancelist)) = tolerancelist(1)
     tolerancelist(1:size(tolerancelist)-1) = 1d0

  !4)if there are two tolerances, the first is for intermediate closures
  elseif(size(tolerancelist) == 2)then
     tolerancelist=>safemem_realloc(tolerancelist,size(closurelist))
     tolerancelist(size(tolerancelist)) = tolerancelist(2)
     tolerancelist(1:size(tolerancelist)-1) = tolerancelist(1)

  !5)otherwise, there should be the same number of tolerances and closures
  elseif(size(tolerancelist) /= size(closurelist))then
     call rism_report_error("number of tolerances must be 1, 2 or the number of closures.")
  end if

  !if a closure number is given, map it to a name
  do iclosure=1,ubound(closurelist,1)
     if(trim(closurelist(iclosure)) .eq. "0")then
        closurelist(iclosure)="HNC"
     elseif(trim(closurelist(iclosure)) .eq. "1") then
        closurelist(iclosure)="KH"
     elseif(trim(closurelist(iclosure)) .eq. "2") then
        closurelist(iclosure)="PSEN"
     end if
     !if the old method of indicating the PSE order has been use (closureOrder) then
     !reformat to the new method
     call caseup(closurelist(iclosure))
     if(trim(closurelist(iclosure)).eq."PSEN" .or. trim(closurelist(iclosure)).eq."PSE")then
        !check if 'closureOrder' is used with a list of closures
        if(iclosure > 1)&
             call rism_report_error("'closureOrder' is depricated and not compatible "&
             //"with closure lists")
        write(fmt,'(a,i4,a)') '(a,i',int(log10(dble(rismprm%closureOrder)))+1,')'
        write(closurelist,fmt) "PSE", rismprm%closureOrder
     end if
  end do
end subroutine sanity_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!modifies the position of the solute using the centering method 
!!!requested.  (This used to also rotate the solute in the box but this
!!!has been disabled.)
!!!If centering is /= 0 we want to move the solute to the center of the solvent box
!!!However, if centering < 0 we only figure out the displacement required the _first_
!!!time we see the solute.  Thus, for centering <= 0, the solute's CM can move relative
!!!to the grid
!!!IN:
!!!   solu :: solute object
!!!   ratu :: the x,y,z position of each solute atom.  This is modified 
!!!           to place the solute in the center of the box
!!!   nsolution :: number of times a RISM solution has been calculated.
!!!                Methods < 0 only are used it nsolution ==0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine orient(solute,ratu, nsolution)
  use amber_rism_interface
  use rism_util, only : calc_cm,translate
  use constants, only : PI
  implicit none
  type(rism3d_solu) :: solute
  _REAL_,intent(inout) :: ratu(3,solute%natom)
  integer, intent(in) :: nsolution

  _REAL_ :: cm(3), weight(solute%natom)
  !counters
  integer :: iu, id
  !................ setting solute site positions in box .................
#ifdef RISM_DEBUG
  write(6,*)"RISM REORIENT"
  call flush(6)
#endif /*RISM_DEBUG*/


  if(abs(rismprm%centering)==1)then
     weight=solute%mass
  elseif(abs(rismprm%centering)==2)then
     weight=1
  end if
  if(rismprm%centering > 0 .or. (rismprm%centering < 0 .and. nsolution == 0))then
     call calc_cm(ratu,ratucm,weight,solute%natom)
  end if
  if(rismprm%centering /= 0) then
     call translate(ratu,solute%natom,-1*ratucm)
  end if
end subroutine orient

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!return the solute to its original position
!!!IN:
!!!   solute :: solute object
!!!   ratu  :: the x,y,z position of each solute atom.  This is modified 
!!!            to return the solute to its original MD position and orientation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine unorient(solute,ratu)
  use amber_rism_interface
  use rism_util, only: calc_cm, translate
  implicit none
  type(rism3d_solu) :: solute
  _REAL_,intent(inout) :: ratu(3,solute%natom)

  _REAL_ :: cm(3),weight(solute%natom)
  !counters
  integer :: iu, id

  !we don't know if additional translations have been performed on the
  !solute.  1) The first step is to calculate the current COM and
  !translate this to the origin.  2) Then the system and forces are
  !rotated according to qback.  3) Finally, using ratucm, we translate
  !the system back to its original MD COM

#ifdef RISM_DEBUG
  write(6,*)"RISM UNORIENT", ratucm
  write(6,*) qback
  call flush(6)
#endif /*RISM_DEBUG*/
  if(abs(rismprm%centering)==1)then
     weight=solute%mass
     call calc_cm(ratu,cm,weight,solute%natom)
     call translate(ratu,solute%natom,-1*cm)
     call calc_cm(ratu,cm,weight,solute%natom)
  elseif(abs(rismprm%centering)==2)then
     weight=1
     call calc_cm(ratu,cm,weight,solute%natom)
     call translate(ratu,solute%natom,-1*cm)
     call calc_cm(ratu,cm,weight,solute%natom)
  end if
  if(rismprm%centering /=0)then
     call translate(ratu,solute%natom,ratucm)
  end if

end subroutine unorient

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Writes the contents of a C string (array of chars) to a Fortran string.  Will
!!!not write past the end of the Fortran string
!!!IN:
!!!   fstr : Fortran string to write to.
!!!   cstr : C string to write from.
!!!   nchar  : Number of chars in cstr (not including null).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cstr2fstr(fstr,cstr,nchar)
  implicit none
  character(len=*), intent(out) :: fstr
  integer, intent(in) :: nchar
  integer(kind=1), intent(in) :: cstr(nchar+1)
  integer :: i
  fstr = ""
  do i = 1, min(nchar,len(fstr))
     if(cstr(i)==0)exit
     fstr = trim(fstr)//char(cstr(i),1)
  end do
end subroutine cstr2fstr


#ifdef SANDER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the least common multiple of integers a and b
!!!IN:
!!!  a : integer
!!!  b : integer
!!!OUT:
!!! least common multiple
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function mylcm(a,b)
  use rism_util, only: lcm
  implicit none
  integer :: mylcm,a,b
  mylcm = lcm(a,b)
end function mylcm
#else /*SANDER*/
!!!stubs for SANDER timers for non-SANDER executables
subroutine timer_start( label )
integer label
end subroutine timer_start

subroutine timer_stop( label )
integer label
end subroutine timer_stop
#endif /*SANDER*/

#ifdef SANDER
end module sander_rism_interface
#endif
