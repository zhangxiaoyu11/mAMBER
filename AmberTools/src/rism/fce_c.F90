!<compile=optimized>

!The 3D-RISM-KH software found here is copyright (c) 2010-2011 by Andriy Kovalenko,
!Tyler Luchko and David A. Case.
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
!!!Object for force-coordinate-extrapolation (FCE) multiple time step (MTS). 
!!!
!!!N total previous frames are held in memory for both position and
!!!solvation force.  When we are extrapolating the force for a given
!!!time step we first express the current position as a linear
!!!combination of the previous steps:
!!!
!!!  R^(N+1) ~ \sum^N_i a_i R^i
!!!
!!!the a_i that best satisfy this equation are then used to calculate
!!!an approximation of the solvation:
!!!
!!!  F^(N+1) ~ \sum^N_i a_i F^i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module fce_c
  implicit none
  type fce
     sequence
     !enumerated values for coordinate type to use as the basis for extrapolation
     integer :: CRDPOS=0, CRDDIST=1, CRDXYZ=2

     !nbasis  :: number of force evaluations before we start using interpolation
     !crd     :: 0-position, 1-distance, 2-xyz distance
     !nsample :: number of samples collected
     !natom   :: number of atoms
     integer :: nbasis=-1, crd=-1, nsample=0, natom=-1

     !mpirank - MPI rank
     !mpicomm - MPI mpicomm
     !mpisize - number of MPI processes
     !atom0 - first atom in MPI decomposition
     !atomF - final atom in MPI decomposition
     integer :: mpirank=0, mpicomm=0, mpisize=1
     integer :: atom0, atomf
     
     !cut :: distance cutoff used for creating the basis set
     _REAL_ :: cut

     !cutlist :: for each atoms, list of atoms to include in force extrapolation (natom+1,natom).
     !           the first row element indicates how many atoms follow in this row.
     integer,pointer :: cutlist(:,:) => NULL()
     !refit :: TBA
     logical,pointer :: refit(:) => NULL()
     
     !force :: previous forces (dimensions:natom:nbasis)
     !coord :: previous coordinates (dimensions:natom:nbasis)
     _REAL_,pointer :: force(:,:,:) => NULL(), coord(:,:,:) => NULL()
     !! coeff :: the a_i coefficients in the above equation
     _REAL_, pointer :: coeff(:,:) => NULL()

  end type fce

  interface fce_new
     module procedure fce_new_std
  end interface
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!public subroutines and functions  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
public fce_new, fce_destroy, fce_update, fce_force

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!private subroutines and functions  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
private nlist

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!public subroutines and functions  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Constructor. If the this is MPI, then only the rank 0 process parameters will
!!!be used.
!!!IN:
!!!   this :: FCE object
!!!   natom :: number of atoms
!!!   nbasis :: number of basis vectors
!!!   crd :: coordinate orientation method. 0-position, 1-distance, 2-xyz distance
!!!   cut :: cut off distance for dependence
!!!   o_mpicomm :: (optional) MPI communicator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fce_new_std(this, natom, nbasis, crd, cut, o_mpicomm)
    use safemem
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif /*MPI*/
    type(fce),intent(inout) :: this
    integer, intent(in) :: natom, nbasis, crd
    _REAL_, intent(in) :: cut
    integer, optional, intent(in) :: o_mpicomm
    integer :: err
#ifdef MPI
    this%mpicomm = 0
    if(present(o_mpicomm)) this%mpicomm = o_mpicomm
    if(this%mpicomm == MPI_COMM_NULL)&
       call rism_report_error("FCE: received NULL MPI communicator")
    call mpi_comm_rank(this%mpicomm,this%mpirank,err)
    if(err /=0) call rism_report_error&
         ("(a,i8)","FCE: could not get MPI rank for communicator ",this%mpicomm)
    call mpi_comm_size(this%mpicomm,this%mpisize,err)
    if(err /=0) call rism_report_error&
         ("(a,i8)","FCE: could not get MPI size for communicator ",this%mpicomm)
#endif /*MPI*/
    if(this%mpirank==0)then
       this%natom = natom
       this%nbasis = nbasis
       this%crd = crd
       this%cut = cut**2
       
       this%cutlist => safemem_realloc(this%cutlist,natom,natom,.false.)
       this%force => safemem_realloc(this%force,3,natom,nbasis,.false.)
       this%coord => safemem_realloc(this%coord,3,natom,nbasis,.false.)
       this%coeff => safemem_realloc(this%coeff,3,natom,.false.)
       
    end if

#ifdef MPI
    !first distribute the pieces of information needed to allocate memory
    call mpi_bcast(this%natom,1,mpi_integer,0,this%mpicomm,err)
    if(err /=0) call rism_report_error&
         ("FCE: could not broadcast NATOM")
    call mpi_bcast(this%nbasis,1,mpi_integer,0,this%mpicomm,err)
    if(err /=0) call rism_report_error&
         ("FCE: could not broadcast NBASIS")
    call mpi_bcast(this%crd,1,mpi_integer,0,this%mpicomm,err)
    if(err /=0) call rism_report_error&
         ("FCE: could not broadcast CRD")
    call mpi_bcast(this%cut,1,mpi_double_precision,0,this%mpicomm,err)
    if(err /=0) call rism_report_error&
         ("FCE: could not broadcast CUT")

    !non-master processes should now allocate memory
    if(this%mpirank /= 0) then
       this%cutlist => safemem_realloc(this%cutlist,this%natom,this%natom,.false.)
       this%force => safemem_realloc(this%force,3,this%natom,this%nbasis,.false.)
       this%coord => safemem_realloc(this%coord,3,this%natom,this%nbasis,.false.)
       this%coeff => safemem_realloc(this%coeff,3,this%natom,.false.)
    end if

    !Arrays contain no data yet so there is nothing to transfer
        
#endif /*MPI*/  
    
    !set local atom range for this process
    this%atom0 = this%natom/this%mpisize*this%mpirank+1
    this%atomF = min(int(this%natom*dble(this%mpirank+1)/dble(this%mpisize)),this%natom)

  end subroutine fce_new_std

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Destroyer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fce_destroy(this)
    use safemem
    implicit none
    type(fce),intent(inout) :: this
    integer :: err
    this%natom = -1
    this%nbasis = -1
    this%crd = -1
    
    err = safemem_dealloc(this%cutlist)
    err = safemem_dealloc(this%force)
    err = safemem_dealloc(this%coord)
    err = safemem_dealloc(this%coeff)
  end subroutine fce_destroy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!update the force and coordinate basis vectors
!!!IN:
!!!   this  : FCE object
!!!   force : forces on atoms.  In the case of MPI, it is assumed that the forces
!!!           are distributed across the processes and must be reduced internally
!!!   coord : coordinates of the atoms, assume to the the same for all processes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fce_update(this, force, coord)
    use safemem
    implicit none
#if defined(MPI)
      include 'mpif.h'
#endif /*defined(MPI)*/
    type(fce), intent(inout) :: this
    _REAL_, intent(in) :: force(3,this%natom), coord(3,this%natom)

    !workspace for the case of MPI 1.1
#ifdef MPI    
    integer :: err
#  ifndef USE_MPI_IN_PLACE
    _REAL_, pointer :: tforce(:,:)=>NULL()
#  endif /*ifdef USE_MPI_IN_PLACE*/
#endif /*MPI*/

#ifdef RISM_DEBUG
    write(6,*) "FCE_UPDATE"; call flush(6)
    write(6,*) "FCE", this%atom0, this%atomF, this%mpisize, this%mpirank, this%mpicomm, this%cut
#endif /*RISM_DEBUG*/

    !
    !For now we will shift the data through the storage arrays as more data is
    !added.  Later, this should be modified to just update a pointer to the most
    !recent entry
    !
    if(this%nsample > 0)then
       this%coord(:,:,2:min(this%nsample+1,this%nbasis)) = this%coord(:,:,1:min(this%nsample, this%nbasis-1))
       this%force(:,:,2:min(this%nsample+1,this%nbasis)) = this%force(:,:,1:min(this%nsample, this%nbasis-1))
    end if
    this%coord(:,:,1) = coord
    this%force(:,:,1) = force

    !
    !reduce the atoms locally
    !
#ifdef MPI    
#  ifdef USE_MPI_IN_PLACE
    call mpi_allreduce(MPI_IN_PLACE,this%force(:,:,1),3*this%natom,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
    if(err /=0) call rism_report_error&
         ("FCE: could not reduce force")
#  else /*ifdef USE_MPI_IN_PLACE*/
    tforce => safemem_realloc(tforce,3,this%natom,.false.)
    call mpi_allreduce(this%force(:,:,1),tforce,3*this%natom,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)    
    if(err /=0) call rism_report_error&
         ("FCE: could not reduce force")
    this%force(:,:,1) = tforce
    if(safemem_dealloc(tforce) /=0) call rism_report_error("FCE_UPDATE: deallocate TFORCE failed")
#  endif /*ifdef USE_MPI_IN_PLACE*/
#endif /*MPI*/


    this%nsample = min(this%nsample+1,this%nbasis)
  end subroutine fce_update

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Using previous forces and coordinates, predicts solvation forces
!!! for the current set of coordinates.  Like  linprojpredict, except 
!!! we are finding the forces on each atom individually and rotate the
!!! solute for each to optimize the prediction.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fce_force(this,force,coord)
    use safemem
    implicit none
    type(fce), intent(inout) :: this
    _REAL_, intent(out) :: force(3,this%natom)
    _REAL_, intent(in) :: coord(3,this%natom)

    integer :: iatm,jatm,id,err
    
    !LAPACK variables
    integer :: M,N,NRHS, LDA,LDB, LWORK,RANK=0, INFO=0
    integer, pointer :: JPVT(:)=>NULL()
    !not really sure what value this should be...
    _REAL_ :: RCOND=0d0
    _REAL_, pointer :: A(:,:)=>NULL(),B(:,:)=>NULL(),work(:)=>NULL()

    _REAL_,external :: ddot

!!$    write(6,*) "FCE_FORCE", this%atom0, this%atomF, this%cut

    force=0
    N = this%nbasis
    JPVT=>safemem_realloc(JPVT,N,.false.)
    JPVT=0
    !the cutoff and non-cutoff situations are handled differently
    if(this%cut>0) then
       !update cutoff list if requested
       call nlist(this)

       !iterate over each atom and calculate the extrapolated force
       do iatm = this%atom0, this%atomF
!          call orient(this)

          M = 3*(this%cutlist(1,iatm)+1)
          LDA = M
          LDB = max(M,N)
          NRHS = 1

          A=>safemem_realloc(A,M,N,.false.)
          B=>safemem_realloc(B,max(N,M),1,.false.)

          !the target atom is never in the cutlist so put it in the first index
          B(1:3,1) = coord(:,iatm)
          A(1:3,:) = this%coord(:,iatm,:)

          !now add the atoms from within the cutoff
          do jatm=2,this%cutlist(1,iatm)+1
             B(jatm*3-2: jatm*3,1) = coord(:,this%cutlist(jatm,iatm))
             A(jatm*3-2: jatm*3,:) = this%coord(:,this%cutlist(jatm,iatm),:)
          end do

          !first determine the amount of memory required
          LWORK=-1
          WORK => safemem_realloc(work,1,.false.)
          call dgelsy(M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, WORK, LWORK, INFO )
          
          !perform the calculation
          LWORK=WORK(1)
          WORK => safemem_realloc(WORK,LWORK,.false.)
          call dgelsy(M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, WORK, LWORK, INFO )

          !the force in each dimension is now a dot product
          do id=1,3
             force(id,iatm) = ddot(N,this%force(id,iatm,:),1,B(1:N,1),1)
          end do
!!$          write(6,*) "FCE FORCE",iatm, force(:,iatm)
 !         call unorient(this)
       end do
    else !cutoff
       
    end if !cutoff
    err = safemem_dealloc(A)
    err = safemem_dealloc(B)
    err = safemem_dealloc(WORK)
    err = safemem_dealloc(JPVT)
  end subroutine fce_force

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!private subroutines and functions  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Very simple method to generate a cutoff list.  The resulting matrix is 2D is 
!!!lists the number of neighbour atoms for each atom (first row) and then the 
!!!atom numbers themselves.
!!!E.g.
!!!   1 2 1
!!!   2 1 2
!!!   0 3 0
!!!   0 0 0
!!!For this three atom system, atoms 1 and 2, 3 and 2 are neighbours.  I.e. 
!!!atoms 1 and 3 each have one neighbour and 2 has two.  The number of 
!!!neighbours is list in the first row.  The subsequent rows list the ids of the
!!!neighbours.
!!!IN:
!!!   this : FCE object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine nlist(this)
  implicit none
  type(fce),intent(inout) :: this
  integer :: id,iatm1,iatm2
#ifdef RISM_DEBUG
  write(6,*) "NLIST"; call flush(6)
#endif /*RISM_DEBUG*/

  this%cutlist=0
  do iatm1 = 1,this%natom
     do iatm2 = iatm1+1,this%natom
        if( sum((this%coord(1:3,iatm1,1) - this%coord(1:3,iatm2,1))**2) < this%cut)then
           this%cutlist(1,iatm1) = this%cutlist(1,iatm1) +1
           this%cutlist(1,iatm2) = this%cutlist(1,iatm2) +1
           this%cutlist(this%cutlist(1,iatm1)+1,iatm1) = iatm2
           this%cutlist(this%cutlist(1,iatm2)+1,iatm2) = iatm1
        end if
     end do
  end do
end subroutine nlist

end module fce_c
