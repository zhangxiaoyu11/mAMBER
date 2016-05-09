!+ Specification and control of Amber's parallel implementation.

! The parallel programs use MPI in a replicated data paradigm.
! The usual algorithms to distribute coordinates and forces,
! which are activated when the preprocessor name noBTREE is undefined,
! employ a recursive doubling scheme based on ideas from the book
! by Fox, et al., and Robert van de Geijn's papers.

#ifndef noBTREE
#  undef  MPI_MAX_PROCESSORS
#  define MPI_MAX_PROCESSORS 128
integer logtwo(MPI_MAX_PROCESSORS)
data logtwo/ 0,1,0,2,0,0,0,3,0,0,0,0,0,0,0,4, &
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5, &
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6, &
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7 /
#else
#  ifndef MPI_MAX_PROCESSORS
#    define MPI_MAX_PROCESSORS 512
#  endif
#endif
integer commworld, commsander, commmaster
integer worldrank, sanderrank, masterrank
integer worldsize, sandersize, mastersize
integer groupmaster, worldmaster
logical ng_sequential
integer numtasks,mytaskid
integer iparpt  ! the atom partition among the processors
! processor i owns atoms from iparpt(i) + 1 to iparpt(i+1)
integer iparpt3,rcvcnt,rcvcnt3
logical mpi_orig
integer ierr,notdone
dimension iparpt(0:MPI_MAX_PROCESSORS)
dimension iparpt3(0:MPI_MAX_PROCESSORS)
dimension rcvcnt(0:MPI_MAX_PROCESSORS)
dimension rcvcnt3(0:MPI_MAX_PROCESSORS)

common/parallel/numtasks,mytaskid,notdone, &
      iparpt,iparpt3,rcvcnt,rcvcnt3,mpi_orig
common/parallel_multi/commworld, commsander, commmaster, &
      worldrank, sanderrank, masterrank, &
      worldsize, sandersize, mastersize, &
      groupmaster, worldmaster, ng_sequential

