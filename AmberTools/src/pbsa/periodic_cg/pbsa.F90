#include "copyright.h"
#include "dprec.h"
#include "assert.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Main program
!  Setup MPI and file handling. 
!  Call pbsa to perform calculations.

program pbsadrv

   implicit none

#  include "files.h"

#ifdef MPI
#  include "mpif.h"
#  include "parallel.h"

   integer ierror, nodeid, masterid, i
   character ext*3


   call mpi_init(ierror)


   !        SET UP THE WORLD COMMUNICATORS AND ASSOCIATED VARIABLES

   !   call MPI_COMM_DUP(MPI_COMM_WORLD, CommWorld, Ierror)
   commworld = mpi_comm_world
   call mpi_comm_rank( commworld, worldrank,  ierror)
   call mpi_comm_size( commworld, worldsize,  ierror)
   call mpi_barrier(commworld, ierror)


   !        CALL MDFIL TO SET THE NAMES OF THE INPUT AND OUTPUT FILES AND
   !        PROCESS OTHER COMMAND LINE ARGUMENTS.

   mytaskid = worldrank
   numtasks = worldsize
   if (worldrank == 0) then
      ng_sequential = .true.
      numgroup = 1
      ! mdfil may modify these depending on the command line options.
      call mdfil
   end if

   !        BROADCAST and VALIDATE THE NUMBER OF GROUPS

   call mpi_bcast(numgroup, 1, mpi_integer, 0, commworld, ierror)

   if (numgroup > worldsize) then
      if (worldrank == 0) then
         write(6,*) 'Error: specified more groups (', numgroup, &
               ') than the number of processors (', numtasks, ') !'
      end if
      call mexit(6,1)
   end if

   !        ...todo: handle -gpes "file" which will alter this requirement

   if ( mod(worldsize, numgroup ) /= 0 ) then
      if (worldrank == 0) then
         write(6,*) 'Error: the number of processors ', &
               'is not a multiple of the number of groups!'
      end if
      call mexit(6,1)
   end if


   !        BROADCAST THE GROUP FILE INFORMATION

   call mpi_bcast(groups, len(groups), mpi_character, 0, &
         commworld, ierror)


   !       PROCESSOR ALLOCATION:  (todo -gpes "files")

   call mpi_bcast(ng_sequential, 1, mpi_logical, 0, commworld, ierror)

   if ( ng_sequential ) then
      nodeid = worldrank / (worldsize / numgroup)
   else
      nodeid = mod(worldrank, numgroup)
   end if

   call mpi_barrier(commworld, ierror)


   !        CREATE A COMMUNICATOR FOR EACH GROUP OF -ng NumGroup PROCESSORS

   commsander = mpi_comm_world
   sandersize = worldsize
   sanderrank = worldrank
   if (numgroup > 1) then
      commsander = mpi_comm_null
      call mpi_comm_split(commworld, nodeid, worldrank, &
            commsander, ierror)
      if (commsander == mpi_comm_null) then
         if (worldrank == 0) then
            write(6,'(a,i5,a,i5)') 'Error: NULL Communicator on PE ', &
                  worldrank, ' from group ', nodeid
         end if
         call mexit(6,1)
      end if
      call mpi_comm_size(commsander, sandersize, ierror)
      call mpi_comm_rank(commsander, sanderrank, ierror)
   end if

   !        DEFINE A COMMUNICATOR (CommMaster) THAT ONLY TALKS BETWEEN THE LOCAL
   !        "MASTER" IN EACH GROUP.  THIS IS EQUIVALENT TO A SanderRank .eq. 0

   masterid = 0
   masterrank = MPI_UNDEFINED
   mastersize = 0
   if (numgroup > 1) then
      commmaster = mpi_comm_null
      if(sanderrank /= 0) then
         masterid = MPI_UNDEFINED
      end if

      call mpi_comm_split(commworld, masterid, worldrank, &
            commmaster, ierror)
      ! will this be emitted when using the default MPI error handler ?
      if (ierror /= MPI_SUCCESS) then
         write(6,*) 'Error: MPI_COMM_SPLIT error ', ierror, &
               ' on PE ', worldrank
      end if

      if(commmaster /= mpi_comm_null) then
         call mpi_comm_size(commmaster, mastersize, ierror)
         call mpi_comm_rank(commmaster, masterrank, ierror)
      end if
   end if


   !        SETUP mytaskid AND numtasks FOR EACH GROUP

   mytaskid = sanderrank
   numtasks = sandersize


   !        DETERMINE AND COMMUNICATE THE FILE INFORMATION TO THE MASTERS

   if (numgroup > 1) then
      if (groups == ' ') then
         if (commmaster /= mpi_comm_null) then

            ! Original/default method: Use the specified command line 
            ! names with 000, 001, ... appended.  This is used if a 
            ! groupfile is NOT specified.

            call mpi_bcast(mdin,   len(mdin),   mpi_character, 0, &
                  commmaster, ierror)
            call mpi_bcast(mdout,  len(mdout),  mpi_character, 0, &
                  commmaster, ierror)
            call mpi_bcast(inpcrd, len(inpcrd), mpi_character, 0, &
                  commmaster, ierror)
            call mpi_bcast(parm,   len(parm),   mpi_character, 0, &
                  commmaster, ierror)
            call mpi_bcast(restrt, len(restrt), mpi_character, 0, &
                  commmaster, ierror)
            call mpi_bcast(refc,   len(refc),   mpi_character, 0, &
                  commmaster, ierror)
            call mpi_bcast(mdvel,  len(mdvel),  mpi_character, 0, &
                  commmaster, ierror)
            call mpi_bcast(mden,   len(mden),   mpi_character, 0, &
                  commmaster, ierror)
            call mpi_bcast(mdcrd,  len(mdcrd),  mpi_character, 0, &
                  commmaster, ierror)
            call mpi_bcast(mdinfo, len(mdinfo), mpi_character, 0, &
                  commmaster, ierror)
            call mpi_bcast(vecs,   len(vecs),   mpi_character, 0, &
                  commmaster, ierror)
            call mpi_bcast(freqe,  len(freqe),  mpi_character, 0, &
                  commmaster, ierror)
            call mpi_bcast(rstdip, len(rstdip), mpi_character, 0, &
                  commmaster, ierror)
            call mpi_bcast(inpdip, len(inpdip), mpi_character, 0, &
                  commmaster, ierror)
            call mpi_bcast(mddip,  len(mddip),  mpi_character, 0, &
                  commmaster, ierror)
            call mpi_bcast(radii,  len(radii),  mpi_character, 0, &
                  commmaster, ierror)
            call mpi_bcast(owrite, 1, mpi_character, 0, &
                  commmaster, ierror)
         end if  ! (commmaster /= mpi_comm_null)

         call mpi_barrier(commworld, ierror)

         write(ext, '(i3.3)') nodeid

         write(inpcrd(index(inpcrd, " "):), '(A3)') ext
         write(  parm(index(parm,   " "):), '(A3)') ext
         write(  mdin(index(mdin,   " "):), '(A3)') ext

         write( mdout(index(mdout,  " "):), '(A3)') ext
         write(restrt(index(restrt, " "):), '(A3)') ext
         write( mdcrd(index(mdcrd,  " "):), '(A3)') ext
         write(  mden(index(mden,  " "):), '(A3)') ext
         write(mdinfo(index(mdinfo, " "):), '(A3)') ext

      else  ! (groups == ' ')

         !     Use the -groupfile command line option to specify the input file
         !     names for each processor.  This is done assuming *each* master has
         !     access to the groupfile.


         if (sanderrank == 0) then

            call amopen(5,groups,'O','F','R')
            i = 0

            do while (i < numgroup )

               read(5,'(a512)') groupbuffer
               if ( groupbuffer(1:1) /= '#' .and. &
                     groupbuffer(1:1) /= '/' ) then

                  if (i == masterrank) then
                     call mdfil
                  end if
                  i = i + 1
               end if
            end do

         end if

         call mpi_barrier(commworld, ierror)

      end if  ! (groups == ' ')


      !        PRINT SUMMARY OF MULTISANDER RUN

      if (worldrank == 0) then
         write(6,*) ''
         write(6,*) ' RUNNING MULTISANDER VERSION OF SANDER AMBER8'
         write(6,*) '    Total processors = ', worldsize
         write(6,*) '    Number of groups = ', numgroup
         if ( .not. ng_sequential) then
            write(6,*) '    Allocation of processors is non-sequential.'
         end if
         write(6,*) ' '
      end if

   end if  ! (numgroup > 1)

#else /* MPI */

      numgroup = 1
      call mdfil()

#endif /* MPI */

   call pbsa()

#ifdef MPI

   !       Clean up and exit

   if(numgroup > 1 .and. commmaster /= mpi_comm_null) then
      call mpi_comm_free(commmaster, ierror)
      commmaster = mpi_comm_null
   end if
   if(numgroup > 1 .and. commsander /= mpi_comm_null) then
      call mpi_comm_free(commsander, ierror)
      commsander = mpi_comm_null
   end if
#endif

   call mexit(6,0)

end program pbsadrv
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ The Molecular Dynamics Module of the AMBER
subroutine pbsa()

   use genborn
   use decomp, only : allocate_int_decomp, allocate_real_decomp, &
                      deallocate_int_decomp, deallocate_real_decomp
   use fastwt

   implicit none

   logical belly, erstop
   integer ier,ifind,jn,ncalls
   character(len=4) itest
   logical ok
#  include "files.h"
#  include "memory.h"
#  include "box.h"
#  include "md.h"
#  include "parms.h"
#  include "extra.h"
#ifdef MPI
   !     =========================== AMBER/MPI ===========================
#  include "parallel.h"
#  ifdef MPI_DOUBLE_PRECISION
#    undef MPI_DOUBLE_PRECISION
#  endif
#  include "mpif.h"
#  ifdef MPI_BUFFER_SIZE
   integer*4 mpibuf(mpi_buffer_size)
#  endif
   _REAL_ ener(30),vir(4)
   !     ========================= END AMBER/MPI =========================
#endif
#  include "nonper.h"
#  include "def_time.h"

   logical do_list_update
   data do_list_update / .true. /

   _REAL_ ene(51)
   integer native,nr3,nr

   ! nmrcal vars
   _REAL_ f,enmr,devdis,devang,devtor,ag,bg,cg
   integer numphi,nttyp,nhb

   ! runmin/trajene var
   _REAL_ carrms

   ! dipole momemt stuff
   integer ngrp

   character(len=8) initial_date, setup_end_date, final_date
   character(len=10) initial_time, setup_end_time, final_time
   _REAL_ time0, time1

   _REAL_,  dimension(:), allocatable :: x, r_stack
   integer, dimension(:), allocatable :: ix, ipairs, i_stack
   character(len=4), dimension(:), allocatable :: ih
   

   !     ---- HERE BEGIN THE EXECUTABLE STATEMENTS ----

   ! Initialize the cpu timer. Needed for machines where returned cpu times
   ! are relative.

   call date_and_time( initial_date, initial_time )
   call wallclock( time0 )
   call init_timers()

#ifdef MPI
   !     =========================== AMBER/MPI ===========================

   !     Parallel initialization (setup is now down outside of sander)

   !     Make PE 0 the master
   master = mytaskid == 0

#  ifndef noBTREE
   !  BTREE is selected by default if noBTREE is not specified
   !     The number of processes is required to be a power of two.
   !
   if ( master .and. numtasks > 1 ) then
      if ( numtasks > MPI_MAX_PROCESSORS .or. &
            logtwo(numtasks) <= 0 ) then  ! assume short-circut logical or
         write(0,*) 'The number of processors must be a power of 2', &
               ' and no greater than ', MPI_MAX_PROCESSORS, &
               ', but is ', numtasks
         call mexit(6,1)
      end if
   end if
#  endif /* BTREE */
#  ifdef MPI_BUFFER_SIZE
   call mpi_buffer_attach(mpibuf, mpi_buffer_size*4, ierr)
#  endif

   !     ========================= END AMBER/MPI =========================

#else   /* not MPI follows */

   !     in the single-threaded version, the one process is master
   master = .true.
#endif  /* MPI */

   erstop = .false.

   !     --- generic packing scheme ---

   nwdvar = 1
   native = 32
#ifdef ISTAR2

   !     --- Int*2 packing scheme ---

   nwdvar = 2
#endif  /*ISTAR2*/
   numpk = nwdvar
   nbit = native/numpk

   !     ----- Only the master node (only node when single-process)
   !           performs the initial setup and reading/writing -----

   call timer_start(TIME_TOTAL)
   if (master) then

      !        --- initialize stack ---

      call rstack_setup()
      call istack_setup()

      !        ---- first, initial reads to determine memry sizes:

      call mdread1()
      call amopen(8,parm,'O','F','R')
      call rdparm1(8)


      !        --- now, we can allocate memory:

      call locmem

      !     --- dynamic memory allocation:

      allocate( x(lastr), ix(lasti), ipairs(lastpr), ih(lasth), stat = ier )
      REQUIRE( ier == 0 )

      lastrst = 1
      if( igb < 10 ) call allocate_gb( natom )
      if( icfe > 0 )  lastrst = lastrst + 3*natom
      allocate( r_stack(1:lastrst), stat = ier )
      REQUIRE( ier == 0 )

      lastist = 1
      allocate( i_stack(1:lastist), stat = ier )
      REQUIRE( ier == 0 )

      if( idecomp > 0 ) then
         call allocate_int_decomp(idecomp, natom, nres)
      endif

      write(6,'(/,a,5x,a)') '|','Memory Use     Allocated'
      write(6,'(a,5x,a,i14)') '|', 'Real      ', lastr
      write(6,'(a,5x,a,i14)') '|', 'Hollerith ', lasth
      write(6,'(a,5x,a,i14)') '|', 'Integer   ', lasti
      write(6,'(a,5x,a,i14)') '|', 'Max Pairs ', lastpr
      write(6,'(a,5x,a,i14)') '|', 'Max Rstack', lastrst
      write(6,'(a,5x,a,i14)') '|', 'Max Istack', lastist
      write(6,'(a,5x,a,i14,a)') '|', '  Total   ', &
           (8*(lastr+lastrst) + 4*(lasth+lasti+lastpr+lastist))/1024, ' kbytes'


      !        --- finish reading the prmtop file and other user input:
      call rdparm2(x,ix,ih,ipairs,8,i_stack)
      call mdread2(x,ix,ih,ipairs,r_stack,i_stack)

      !        --- alloc memory for decomp module that needs info from mdread2
      if( idecomp == 1 .or. idecomp == 2 ) then
         call allocate_real_decomp(nres)
      else if( idecomp == 3 .or. idecomp == 4 ) then
         call allocate_real_decomp(npdec*npdec)
      end if

      !        ----- EVALUATE SOME CONSTANTS FROM MDREAD SETTINGS -----

      nr = nrp
      nr3 = 3*nr
      belly = ibelly > 0

      !        --- seed the random number generator ---

      call amrset(ig)

      if (nbit < 32 .and. nr > 32767) then
         write(6, *) '  Too many atoms for 16 bit pairlist -'
         write(6, *) '    Recompile without ISTAR2'
         call mexit(6, 1)
      end if

      if (ntp > 0.and.iabs(ntb) /= 2) then
         write(6,*) 'Input of NTP/NTB inconsistent'
         call mexit(6, 1)
      end if

      !        ----- READ COORDINATES AND VELOCITIES -----

      call timer_start(TIME_RDCRD)
      call getcor(nr,x(lcrd),x(lvel),x(lforce),ntx,box,irest,t)

      !        ----- SET THE INITIAL VELOCITIES -----

      if (ntx <= 3) call setvel(nr,x(lvel),x(lwinv),tempi,init,0,0.0)

      if (belly) call bellyf(natom,ix(ibellygp),x(lvel))
      call timer_stop(TIME_RDCRD)

      !        --- Call FASTWAT, which will tag those bonds which are part
      !            of 3-point water molecules. Constraints will be effected
      !            for these waters using a fast analytic routine -- dap.

      call timer_start(TIME_FASTWT )

      call fastwat(ih(m04),nres,ix(i02),ih(m02), &
            nbonh,nbona,ix(iibh),ix(ijbh),ibelly,ix(ibellygp), &
            iwtnm,iowtnm,ihwtnm,jfastw,ix(iifstwt), &
            ix(iifstwr),ibgwat,ienwat,iorwat, &
            6,natom)
      call timer_stop(TIME_FASTWT)

      call getwds(ih(m04)   ,nres      ,ix(i02)   ,ih(m02)   , &
            nbonh     ,nbona     ,0         ,ix(iibh)  ,ix(ijbh)  , &
            iwtnm     ,iowtnm    ,ihwtnm    ,jfastw    ,ix(iicbh) , &
            req       ,x(lwinv)  ,rbtarg    ,ibelly  ,ix(ibellygp), &
            6)

      !        --- OPEN THE DATA DUMPING FILES AND POSITION IT DEPENDING
      !            ON THE TYPE OF RUN -----

      call open_dump_files
      if (master) call amflsh(6)

      !        --- end of master process setup ---
   end if  ! (master)


#ifdef MPI
   !     =========================== AMBER/MPI ===========================

   !     NOTE: in the current AMBER/MPI implementation, two means of
   !     running in parallel within sander are supported. The value
   !     of mpi_orig determines which approach is used.
   !     This is turned on when minimization (imin .ne. 0) is requested,
   !     and is otherwise off.

   !     When running the mpi_orig case, a variable notdone is now
   !     set by the master and determines when to exit the force()
   !     loop.  When the master has finished calling force, the
   !     master changes notdone to 0 and broadcasts the data one more
   !     time to signal end of the loop.  force() is modified so that
   !     in the mpi_orig case, an initial broadcast is done to receive
   !     the value from the master to decide whether to do the work or
   !     simply exit.

   !     ...set up initial data and send all needed data to other nodes,
   !     now that the master has it

   nr = nrp
   nr3 = 3*nr
   belly = ibelly > 0

   !     First, broadcast parameters in memory.h, so that all processors
   !     will know how much memory to allocate:

   call mpi_bcast(natom,BC_MEMORY,mpi_integer,0,commsander,ierr)
   call mpi_barrier(commsander,ierr)

   !     ---allocate memory on the non-master nodes:

   if( .not.master ) then
      allocate( x(1:lastr), stat = ier )
      REQUIRE( ier == 0 )

      allocate( ix(1:lasti), stat = ier )
      REQUIRE( ier == 0 )

      allocate( ipairs(1:lastpr), stat = ier )
      REQUIRE( ier == 0 )

      allocate( ih(1:lasth), stat = ier )
      REQUIRE( ier == 0 )

      call allocate_gb( natom )

      allocate( r_stack(1:lastrst), stat = ier )
      REQUIRE( ier == 0 )

      allocate( i_stack(1:lastist), stat = ier )
      REQUIRE( ier == 0 )
   end if  ! ( .not.master )

   call startup_groups(ierr)
   call startup(x,ix,ih)

   call amrset(ig) ! different initialization for each node

   !    ---------------- Old parallel for minimization ----------------------

   if (imin /= 0 .or. plevel == 0) then
      mpi_orig = .true.
      notdone = 1
   else
      mpi_orig = .false.
   end if
   if (mpi_orig .and. .not.master) then

      !        ...all nodes only do the force calculations (JV)

      do while( notdone == 1 )
         call force(x,ix,ih,ipairs,x(lcrd),x(lforce),ener,vir, &
               r_stack,i_stack, x(l96), x(l97), x(l98), do_list_update)
      end do

      goto 999  ! deallocate and return
   end if
   !    ----------------------------------------------------------------------

   if (master) write(6, '(a,i4,a,/)') &
         '|  Running AMBER/MPI version on ',numtasks, ' nodes'
   if (master .and. numgroup > 1) write(6, '(a,i4,a,i4,a,i4,a)') &
         '|  MULTISANDER: ', numgroup, 'groups. ', &
         numtasks, ' processors out of ', worldsize, ' total.'
   if(master)call amflsh(6)

   !     ========================= END AMBER/MPI =========================
#endif /* MPI */

   call date_and_time( setup_end_date, setup_end_time )


   ! ----------------------------------------------------------------------
   ! Now do the dynamics or minimization.
   ! ----------------------------------------------------------------------

   if( master ) write(6,'(/80(1H-)/''   4.  RESULTS'',/80(1H-)/)')

   ! Input flag imin determines the type of calculation: MD, minimization, ...
   select case ( imin )
   case ( 0 )
      !        --- Dynamics:

      call timer_start(TIME_RUNMD)
      call runmd(x,ix,ih,ipairs, &
            x(lcrd),x(lwinv),x(lmass),x(lforce), &
            x(lvel),x(lvel2),x(l45),x(lcrdr), &
            x(l50),x(l95),ix(i70),x(l75),erstop,r_stack,i_stack)
      call timer_stop(TIME_RUNMD)

      if (master) call amflsh(6)

      if (erstop) then
         ! This error condition stems from subroutine shake;
         ! furthermore, it seems that erstop can never be true since shake
         ! can never return with its third last argument, niter, equal to 0.
         ! SRB, Sep 24, 2003
         if (master) then
            write(6, *) 'FATAL ERROR'
         end if
         call mexit(6,1)
      end if

   case ( 1 )
      !        --- Minimization:

      ! Input flag ntmin determines the method of minimization
      select case ( ntmin )
      case ( 0, 1, 2 )
         call runmin(x,ix,ih,ipairs,x(lcrd),x(lforce),x(lvel), &
               ix(iibh),ix(ijbh),x(l50),x(lwinv),ix(ibellygp), &
               x(l95),ene,r_stack,i_stack, carrms)
      case default
         ! invalid ntmin
         ! ntmin input validation occurs in mdread.f
         ASSERT( .false. )
      end select

      if (master) call minrit(x(lcrd))  ! Write the restart file

   case ( 5 )
      !       ---carlos modified for reading trajectories (trajene option)

      write (6,*) "POST-PROCESSING OF TRAJECTORY ENERGIES"

      !       ---read trajectories and calculate energies for each frame

      call trajene(x,ix,ih,ipairs, ene,ok,r_stack,i_stack)

      if (.not.ok) then
         write (6,*) 'error in trajene()'
         call mexit(6,1)
      end if

   case default
      ! invalid imin
      ! imin input validation should be transferred to mdread.f
      write(6,'(/2x,a,i3,a)') 'Error: Invalid IMIN (',imin,').'
      ASSERT( .false. )
   end select

   !     -- calc time spent running vs setup

   call timer_stop(TIME_TOTAL)
   call wallclock( time1 )
   call date_and_time( final_date, final_time )
#ifdef NO_DETAILED_TIMINGS
#else
   call profile_time( time1 - time0 )
#endif

#ifdef MPI
   !     =========================== AMBER/MPI ===========================

   !     Set and broadcast notdone in mpi_orig case to inform
   !     other nodes that we are finished calling force(). (tec3)

   if ( mpi_orig ) then
      notdone = 0
      call mpi_bcast(notdone,1,mpi_integer,0, commsander,ierr)
   end if

   !     ========================= END AMBER/MPI =========================
#endif

   if( master ) then
      call close_dump_files

      !     --- write out final times

      write(6,'(12(a))') '|           Job began  at ', initial_time(1:2), &
           ':', initial_time(3:4), ':', initial_time(5:10), '  on ',&
           initial_date(5:6), '/', initial_date(7:8), '/', initial_date(1:4)
      write(6,'(12(a))') '|           Setup done at ', setup_end_time(1:2),  &
           ':', setup_end_time(3:4), ':', setup_end_time(5:10), '  on ', &
           setup_end_date(5:6), '/',setup_end_date(7:8),'/',setup_end_date(1:4)
      write(6,'(12(a))') '|           Run   done at ', final_time(1:2),  &
           ':', final_time(3:4), ':', final_time(5:10), '  on ', &
           final_date(5:6), '/', final_date(7:8), '/', final_date(1:4)
      call nwallclock( ncalls )
      write(6, '(''|'',5x,''wallclock() was called'',I8,'' times'')') ncalls
   end if

999 continue  !     --- dynamic memory deallocation:

   if (master) then  
      if( idecomp > 0 ) then
         call deallocate_real_decomp()
         call deallocate_int_decomp(idecomp)
      endif
   endif
   deallocate( i_stack, stat = ier )
   REQUIRE( ier == 0 )
   deallocate( r_stack, stat = ier )
   REQUIRE( ier == 0 )
   if( igb < 10 ) then
      call deallocate_gb( )
   else
      call pb_free( )
   end if
   deallocate( ih, stat = ier )
   REQUIRE( ier == 0 )
   deallocate( ipairs, stat = ier )
   REQUIRE( ier == 0 )
   deallocate( ix, stat = ier )
   REQUIRE( ier == 0 )
   deallocate( x, stat = ier )
   REQUIRE( ier == 0 )

   return

end subroutine pbsa

!*********************************************************************
!               SUBROUTINE TRAJENE
!*********************************************************************
! carlos add trajene routine for processing trajectory energies

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine trajene here]
subroutine trajene(x,ix,ih,ipairs,ene,ok, r_stack,i_stack)

   implicit none
   integer ipairs(*),i_stack(*)
   _REAL_ r_stack(*),carrms

   _REAL_ x(*),ene(*)
   integer ix(*),ih(*)
   logical ok
   integer member,j

#  include "memory.h"
#  include "extra.h"
#  include "box.h"

   member=0

   !     loop over trajectory file, exiting only on error or end of file

   do while ( .true. )

      !       --- read next coordinate set from trajectory

      read(12,110,end=1000,err=1010) (x(j),j=lcrd,lcrd+natom*3-1)

      if (ifbox > 0) read(12,110,end=1000,err=1010)

      !       --- uncomment this to force box read if prmtop doesn't have it
      !           but traj does:
      !        READ(12,110,END=1000,ERR=1010)
      110 format(10f8.3)

      member=member+1

      write (6,'(a,i6)') 'minimizing coord set #',member

      call runmin(x,ix,ih,ipairs,x(lcrd),x(lforce),x(lvel), &
            ix(iibh),ix(ijbh),x(l50),x(lwinv),ix(ibellygp), &
            x(l95),ene,r_stack,i_stack,carrms)

      write (6,364) ene(1),carrms
      364 format ('minimization completed, ENE=',1x,e12.6, &
            1x,'RMS=',1x,e12.6)

      !       ---loop for next coordinate set

   end do

   !     ---end of trajectory file

   1000 write (6,*) "Trajectory file ended"
   ok=.true.
   return

   1010 write (6,*) "Error in trajectory file"
   ok=.false.
   return

end subroutine trajene
!-----------------------------------------------------------------
#if defined (MPI)
#error MPI option is not supported
#endif
