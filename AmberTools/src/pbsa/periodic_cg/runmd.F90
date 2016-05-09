! <compile=optimized>
#include "copyright.h"
#include "dprec.h"
#include "assert.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ main driver routine for molecular dynamics
subroutine runmd(xx,ix,ih,ipairs,x,winv,amass,f, &
      v,vold,xr,xc,conp,skip,nsp,tma,erstop,r_stack,i_stack)
   

   use decomp, only : jgroup, index
   use fastwt
   implicit none
   integer   ipairs(*)
   _REAL_ xx(*)
   integer   ix(*)
   character(len=4) ih(*)
   _REAL_ r_stack(*)
   integer i_stack(*)
#ifdef MPI
   !     =========================== AMBER/MPI ===========================
   !     NOTE: this routine contains MPI functionality to update the
   !     positions and velocities for all the atoms on a given node.
   !     After the positions are updated, communication (where necessary)
   !     is performed to update the coordinates and velocities on each
   !     processor.  Also note that wrappers to all the I/O routines are
   !     present to prevent all nodes except the master from executing
   !     I/O.  In addition, extra timing variables and calls are defined
   !     to allow profiling of this routine.
   
#  include "parallel.h"
#ifdef MPI_DOUBLE_PRECISION
#undef MPI_DOUBLE_PRECISION
#endif
#  include "mpif.h"
   _REAL_ ekcmtt(3),ekmpi(3),ekmpit(3)
   !     ========================= END AMBER/MPI =========================
#endif
   
#  include "files.h"
#  include "md.h"
#  include "box.h"
#  include "memory.h"
#  include "extra.h"
#  include "def_time.h"
#  include "constants.h"
#  include "pb_md.h"
#  include "random.h"   
#  include "parms.h"
 

   character(len=6)fnam

   ! Constant pH
   integer :: icpselres, icpselstat ! randomly selected residue and state
   type (rand_gen_state) :: cnstph_rand_gen

   logical resetvelo
   integer nshak
   _REAL_ ekgs,eold3,eold4,etot_save,ekpbs
   
   logical do_list_update
   logical skip(*),belly,lout,loutfm,erstop,vlim,onstep
   _REAL_ x(*),winv(*),amass(*),f(*),v(*),vold(*), &
         xr(*),xc(*),conp(*),vol
   _REAL_ enert(51),enert2(51),ener(51),vir(4),ekcmt(4)
   _REAL_ enert_old(51),enert2_old(51),ecopy(51),edvdl(51), &
         enert_tmp(51),enert2_tmp(51),etot_start,edvdl_r(51)
   _REAL_ pres(4),rmu(3),fac(3),vircopy(3),clfac
   _REAL_ tma(*)

   _REAL_ tspan,atempdrop,fln,scaltp
   _REAL_ vel,vel2,vcmx,vcmy,vcmz,vmax,vx,vy,vz
   _REAL_ winf,aamass,rterm,ekmh,ekph,ekpht,wfac,rsd,ekav
   _REAL_ fit,fiti,fit2
   _REAL_ gammai,c_implic,c_explic,c_ave,sdfac,ekins0
   _REAL_ dtx,dtxinv,dt5,factt,ekin0,ekinp0,dtcp,dttp
   _REAL_ rndf,rndfs,rndfp,ibelsv,onet,boltz2,pconv,tempsu
   _REAL_ xcm(3),acm(3),ocm(3),vcm(3),ekcm,ekrot

   integer nsp(*)
   integer idumar(4)
   integer l_temp
   integer i,j,im,m,i3,nitp,nits
   integer nstep,nrep,nrek,nren,iend,istart3,iend3
   integer nrx,nr,nr3,ntcmt,izero,istart
   logical ixdump,ivarch,itdump
   
   equivalence (scaltp,ener(5)),(vol,ener(10))
   equivalence (pres(1),ener(11)),(ekcmt(1),ener(15))
   equivalence (vir(1),ener(19))
   integer nvalid
   _REAL_ eke,eket
   _REAL_ extent

   _REAL_ xcen,ycen,zcen,extents(3,2),centertest
   _REAL_, allocatable, dimension(:) :: frcti
   integer ier
   
   _REAL_ small
   data small/1.0d-7/
   data nren/51/
   
   !--- VARIABLES FOR DIPOLE PRINTING ---
   integer prndipngrp
   integer prndipfind
   character(len=4) prndiptest
   !--- END VARIABLES FOR DIPOLE PRINTING ---

   !  Runmd operates in kcal/mol units for energy, amu for masses,
   !     and angstoms for distances.  To convert the input time parameters
   !     from picoseconds to internal units, multiply by 20.455
   !     (which is 10.0*sqrt(4.184)).
   
   !==========================================================================
   
   !     ----- INITIALIZE SOME VARIABLES -----
   
   if( master ) call amopen(7,mdinfo,'U','F','W')
   vlim = vlimit > small
   ntcmt = 0
   izero = 0
   belly = ibelly > 0
   lout = .true.
   loutfm = ioutfm <= 0
   nr = nrp
   nr3 = 3*nr
   ekmh = 0.d0
#ifdef MPI
   if ( mpi_orig ) then
      istart = 1
      iend = natom
   else
      istart = iparpt(mytaskid) + 1
      iend = iparpt(mytaskid+1)
   end if
#else
   istart = 1
   iend = natom
#endif
   istart3 = 3*istart -2
   iend3 = 3*iend

   if( icfe /= 0 ) then
      allocate( frcti( nr3 ), stat = ier )
      REQUIRE( ier == 0 )
   end if
      
   ! If NTWPRT.NE.0, only print the atoms up to this value
   nrx  = nr3
   if (ntwprt > 0) nrx = ntwprt*3
   
   ! Cleanup the velocity if belly run
   if(belly) call bellyf(nr,ix(ibellygp),v)
   
   !=======================================================================
   
   ! Determine system degrees of freedom (for T scaling, reporting)
   
   ! Call DEGCNT to get the actual number of degrees of freedom for the
   ! solute and solvent. This call returns the correct numbers for belly
   ! simulations and simulations with separate solute/solvent scaling -- dap
   
   ! "IDUMAR" is dummy array. Used since this routine also used w/ GIBBS.
   
   call degcnt(ibelly,nr,ix(ibellygp),nsolut,nbonh,nbona,0, &
         ix(iibh),ix(ijbh),ix(iiba),ix(ijba),idumar, &
         idumar,ntc,idumar,0,0,0, &
         idumar,ibelsv,rndfp,rndfs)
   
   ! RNDFP = # degrees of freedom for solute
   ! RNDFS = # degrees of freedom for solvent
   ! RNDF = total number of degrees of freedom.
   
   !    modify RNDFP to reflect NDFMIN (set in mdread)
   
   rndfp = rndfp - ndfmin
   rndf = rndfp+rndfs
   
   ! End of degrees of freedom stuff
   
   !=======================================================================
   
   onet = 1.d0/3.d0
   
   boltz2 = 8.31441d-3 * 0.5d0
   pconv = 1.6604345d+04  ! factor to convert the pressure kcal/mole to bar
   
   !     ---convert to kcal/mol units
   
   boltz2 = boltz2/4.184d0   ! k-sub-B/2
   dtx = dt*20.455d+00
   dtxinv = 1.0d0 / dtx
   dt5 = dtx * 0.5d0
   pconv = pconv*4.184d0
   
   ! FAC() are #deg freedom * kboltz / 2
   ! multiply by T to get expected kinetic energy
   ! FAC(1) is for total system
   
   fac(1) = boltz2*rndf
   fac(2) = boltz2*rndfp
   if(rndfp < 0.1d0) fac(2) = 1.d-6

   fac(3) = boltz2*rndfs
   if(rndfs < 0.1d0) fac(3) = 1.d-6
   factt = rndf/(rndf+ndfmin)
   
   ! these are "desired" kinetic energies based on
   ! # degrees freedom and target temperature
   ! they will be used for calculating the velocity scaling factor
   
   ekin0  = fac(1)*temp0
   ekinp0 = fac(2)*temp0
   
   ekins0 = fac(3)*temp0

   !     LN setup:
   
   gammai = gamma_ln/20.455d0
   c_implic = 1.d0/(1.d0+gammai*dt5)
   c_explic = 1.d0 - gammai*dt5
   c_ave    = 1.d0+gammai*dt5
   sdfac = sqrt( 4.d0*gammai*boltz2*temp0/dtx )

   !     Constant pH setup 
   ! (separate stream of random numbers so choices stay sync'ed between MPI nodes)
   !
   if (icnstph /= 0) then
      call amrset_gen(cnstph_rand_gen, ig)
   end if
   
   if (ntt == 1) dttp = dt/tautp
   
   nrek = 4
   nrep = 15
   
   nvalid = 0
   nstep = 0
   fit = 0.d0
   fiti = 0.d0
   fit2 = 0.d0

   do i = 1,nren
      ener(i) = 0.0d0
      enert(i) = 0.0d0
      enert2(i) = 0.0d0
      enert_old(i) = 0.d0
      enert2_old(i) = 0.d0
      edvdl(i) = 0.d0
      edvdl_r(i) = 0.d0
   end do
   
   ener(5) = 1.d0
   ener(6) = 1.d0
   do m = 1,3
      ener(m+6) = box(m)
   end do
   
   nitp = 0
   nits = 0
   
   !=======================================================================
   !     ----- MAKE A FIRST DYNAMICS STEP -----
   !=======================================================================
   !  init = 3:  general startup if not continuing a previous run
   
   if( init == 3 ) then
      
      ! ----- CALCULATE THE FORCE -----

      npbstep = nstep

      !   ---   set irespa to get full energies calculated on step "0":
      irespa = 0
      !       TIME_force is started and stopped inside force
      call force(xx,ix,ih,ipairs,x,f,ener(23),vir, &
            r_stack,i_stack, xx(l96), xx(l97), xx(l98), do_list_update)
      if (icnstph /= 0 ) call cnstphwrite(ix(icpresst),0,ix(icptrsct))
      
      ! This FORCE call does not count as a "step". CALL NMRDCP to decrement
      ! local NMR step counter
      nstep = nstep - 1

      irespa = 1
      
      ! Reset quantities depending on TEMP0 and TAUTP (which may have been
      ! changed by MODWT during FORCE call).
      ! Recalculate target kinetic energies.
      
      ekinp0 = fac(2) * temp0
      ekins0 = fac(3) * temp0
      ekin0 = fac(1) * temp0

      if (ntt == 1) dttp = dt / tautp
      
      ntnb = 0
      i3 = 0
      tempsu = 0.0d0
      
      do j = 1,nrp
         winf = winv(j) * dt5
         aamass = amass(j)
         do m = 1,3
            i3 = i3+1
            rterm = v(i3)*v(i3) * aamass
            tempsu = tempsu + rterm
            v(i3) = v(i3) - f(i3) * winf
            if (vlim) v(i3) = sign(min(abs(v(i3)),vlimit),v(i3))
         end do
      end do
      
      do im=1,0
         v(nr3+im) = v(nr3+im) - f(nr3+im) * dt5 / 100.0
         tempsu = tempsu + 100.0 * v(nr3+im)*v(nr3+im)
      end do
      
      ! store kinetic energies: ENER() are actual KE while FAC() are target values
      ! ENER(2): total system, target is ekin0
      ! ENER(3): non-solvent, target is ekinp0
      ! ENER(4): solvent, target is ekins0
      
      ener(3) = tempsu * 0.5d0
      ener(2) = ener(3)
      ener(1) = ener(2)+ener(23)

      if(ntt == 1) then
         ekmh = max(ener(3),fac(1)*10.d0)
      end if
      
      ! endif for init=3 check
      
   end if  ! ( init == 3 )
   
   !-------------------------------------------------------------------------
   ! init = 4:  continuation of a previous trajectory
   !            this code also done for init=3
   !
   ! Note: if the last printed energy from the previous trajectory was
   !       at time "t", then the restrt file has velocities at time
   !       t + 0.5dt, and coordinates at time t + dt
   !-------------------------------------------------------------------------
   
   ekmh = 0.0d0
   i3 = 0
   do j = 1,nrp
      aamass = amass(j)
      do m = 1,3
         i3 = i3+1
         rterm = v(i3)*v(i3) * aamass
         ekmh = ekmh + rterm
      end do
   end do
   
   do im=1,0
      ekmh = ekmh + 100.0*v(nr3+im)*v(nr3+im)
   end do
   ekmh = ekmh * 0.5d0
   
   do i=1,nr3+0
      vold(i) = v(i)
   end do
   
   if (init == 4) then
   else
      
      !-------------------------------------------------------------------
      !           PRINT THE INITIAL ENERGIES AND TEMPERATURES
      !-------------------------------------------------------------------
      
      if (nstep <= 0 .and. master) then
         rewind(7)
         call prntmd(nstep,nitp,nits,t,ener,fac,7,.false.)
         call amflsh(7)
      end if
      if (nstlim == 0) return
      init = 4
   end if
   
   !=======================================================================
   !     ----- MAIN LOOP FOR PERFORMING THE DYNAMICS STEP -----
   !           (at this point, the coordinates are a half-step "ahead"
   !           of the velocities; the variable EKMH holds the kinetic
   !           energy at these "-1/2" velocities, which are stored in
   !           the array VOLD.)
   !=======================================================================
   
   260 continue
   onstep = mod(irespa,nrespa) == 0
   
   !---------------------------------------------------------------
   !  ---Step 1a: do some setup for pressure calculations:
   !---------------------------------------------------------------
   
   ! Constant pH setup
   if ((icnstph /= 0) .and. (mod(irespa,ntcnstph) == 0)) then
      call cnstphbeginstep(ix(icpstinf),ix(icpresst), &
            ix(icptrsct), &
            xx(lcpene),xx(lcpcrg),xx(l190),icpselres,icpselstat, cnstph_rand_gen)
   end if

   !--------------------------------------------------------------
   !  ---Step 1b: Get the forces for the current coordinates:
   !--------------------------------------------------------------
   
   npbstep = nstep

   !     TIME_force is started and stopped inside force
   call force(xx,ix,ih,ipairs,x,f,ener(23),vir, &
         r_stack,i_stack,xx(l96),xx(l97),xx(l98),do_list_update)

   ! Constant pH transition evaluation
   if ((icnstph /= 0) .and. (mod(irespa,ntcnstph) == 0)) then
      call cnstphendstep(ix(icpstinf),ix(icpresst),ix(icpptcnt), &
            ix(icptrsct), xx(lcpene),xx(lcpcrg),xx(l190),xx(l15),ener(39), &
            icpselres,icpselstat,cnstph_rand_gen)
      call cnstphwrite(ix(icpresst),icpselres,ix(icptrsct))
   end if

   ! Reset quantities depending on TEMP0 and TAUTP (which may have been
   ! changed by MODWT during FORCE call).
   
   ekinp0 = fac(2)*temp0
   ekin0 = fac(1)*temp0
   if (ntt == 1) dttp = dt/tautp
   
   !  Pressure coupling:
   
   !----------------------------------------------------------------
   !  ---Step 1c: do randomization of velocities, if needed:
   !----------------------------------------------------------------
   ! ---Assign new random velocities every Vrand steps, if ntt=2

   resetvelo=.false.
   if (vrand /= 0 .and. ntt == 2) then
      if (mod((nstep+1),vrand) == 0) resetvelo=.true.
   end if

   if (resetvelo) then
      if (master) then
         write (6,'(a,i8)') 'Setting new random velocities at step ', &
               nstep + 1
         call setvel(natom,v,winv,temp0*factt,init,0,100.0)
         if (ibelly > 0) call bellyf(nr,ix(ibellygp),v)
      end if
# ifdef MPI
      call mpi_bcast(v, 3*natom, MPI_DOUBLE_PRECISION, &
            0, commsander, ierr)
# endif
      
   end if  ! (resetvelo)
   
   call timer_start(TIME_VERLET)
   
   !-----------------------------------------------------
   !  ---Step 2: Do the velocity update:
   !-----------------------------------------------------
   
   if( gammai == 0.d0 ) then
      
      !       ---Newtonian dynamics:
      
      i3 = 3*(istart-1)
      do j=istart,iend
         wfac = winv(j) * dtx
         v(i3+1) = v(i3+1) + f(i3+1)*wfac
         v(i3+2) = v(i3+2) + f(i3+2)*wfac
         v(i3+3) = v(i3+3) + f(i3+3)*wfac
         i3 = i3+3
      end do
      
   else  !  gamma_ln .ne. 0, which also implies ntt=3 (see mdread.f)
      
      !       ---simple model for Langevin dynamics, basically taken from
      !          Loncharich, Brooks and Pastor, Biopolymers 32:523-535 (1992),
      !          Eq. 11.  (Note that the first term on the rhs of Eq. 11b
      !          should not be there.)
      
      ! In order to generate the same sequence of pseudorandom numbers that you
      ! would using a single processor you have to go through the atoms 
      ! in order.  The unused results are thrown away

      i3 = 3*(istart-1)
      do j=1,natom
         if( j<istart .or. j>iend ) then
            call gauss( 0.d0, 1.d0, fln )
            call gauss( 0.d0, 1.d0, fln )
            call gauss( 0.d0, 1.d0, fln )
            cycle
         end if
         wfac = winv(j) * dtx
         aamass = amass(j)
         rsd = sdfac*sqrt(aamass)
         call gauss( 0.d0, rsd, fln )
         v(i3+1) = (v(i3+1)*c_explic + (f(i3+1)+fln)*wfac) * c_implic
         call gauss( 0.d0, rsd, fln )
         v(i3+2) = (v(i3+2)*c_explic + (f(i3+2)+fln)*wfac) * c_implic
         call gauss( 0.d0, rsd, fln )
         v(i3+3) = (v(i3+3)*c_explic + (f(i3+3)+fln)*wfac) * c_implic
         i3 = i3+3
      end do
   end if  ! ( gammai == 0.d0 )
   
   !     --- consider vlimit
   
   if (vlim) then
      vmax = 0.0d0
      do i=istart3,iend3
         vmax = max(vmax,abs(v(i)))
         v(i) = sign(min(abs(v(i)),vlimit),v(i))
      end do

      ! Only violations on the master node are actually reported
      ! to avoid both MPI communication and non-master writes.
      if (vmax > vlimit) then
         if (master) then
            write(6,*) 'vlimit exceeded for step ',nstep,'; vmax = ',vmax
         end if
      end if
   end if
   
   do im=1,0
      v(nr3+im) = (v(nr3+im) + f(nr3+im)*dtx/100.0)
   end do
   
   !-------------------------------------------------------------------
   !   Step 3: update the positions, putting the "old" positions into F:
   !-------------------------------------------------------------------
   
   do i3 = istart3,iend3
      f(i3) = x(i3)
      x(i3) = x(i3)+v(i3)*dtx
   end do
   do i = 1,0
      f(nr3+i) = x(nr3+i)
      x(nr3+i) = x(nr3+i)+v(nr3+i)*dtx
   end do
   call timer_stop(TIME_VERLET)
   
   if (ntc /= 1) then

      !-------------------------------------------------------------------
      !   Step 4a: if shake is being used, update the new positions to fix
      !      the bond lengths.
      !-------------------------------------------------------------------
   
      call timer_start(TIME_SHAKE)
      call shake(nrp,nbonh,nbona,0,ix(iibh),ix(ijbh),ix(ibellygp), &
            winv,conp,skip,f,x,nitp,belly,ix(iifstwt))
      if(nitp == 0) then
         erstop = .true.
         goto 480
      end if
      call quick3(f,x,ix(iifstwr),natom,nres,ix(i02))

      !-----------------------------------------------------------------
      !   Step 4b:   Now fix the velocities and calculate KE
      !-----------------------------------------------------------------
      
      !     ---re-estimate the velocities from differences in positions:
      v(istart3:iend3) = (x(istart3:iend3)-f(istart3:iend3)) * dtxinv
      call timer_stop(TIME_SHAKE)
   end if
   call timer_start(TIME_VERLET)
   
   if( ntt == 1 .or. onstep ) then
      
      !-----------------------------------------------------------------
      !   Step 4c: get the KE, either for printing or for Berendsen:
      !-----------------------------------------------------------------
      
      eke = 0.d0
      ekph = 0.d0
      ekpbs = 0.d0
      i3 = 3*(istart-1)
      do j=istart,iend
         aamass = amass(j)
         do m = 1,3
            i3 = i3+1
            if( gammai == 0.d0 ) then

               eke = eke + aamass*0.25d0*(v(i3)+vold(i3))**2

               ! try pseudo KE from Eq. 4.7b of Pastor, Brooks & Szabo,
               ! Mol. Phys. 65, 1409-1419 (1988):
               ekpbs = ekpbs + aamass*v(i3)*vold(i3)

               ekph = ekph + aamass*v(i3)**2
            else
               eke = eke + aamass*0.25d0*c_ave*(v(i3)+vold(i3))**2
            endif
         end do
      end do
      
#ifdef MPI
      
      !  ---   sum up the partial kinetic energies:
      
      if ( .not. mpi_orig .and. numtasks > 1 ) then
         ekmpi(1) = eke
         ekmpi(2) = ekph
         ekmpi(3) = ekpbs
         call mpi_allreduce(ekmpi,ekmpit,3, &
               MPI_DOUBLE_PRECISION,mpi_sum,commsander,ierr)
         eke = ekmpit(1)
         ekph = ekmpit(2)
         ekpbs = ekmpit(3)
      end if
#endif
      
      !         --- all processors handle the "extra" variables:
      
      do im=1,0
         eke = eke + 100.0*0.25d0*(v(nr3+im)+vold(nr3+im))**2
         ekpbs = ekpbs + 100.0*v(nr3+im)*vold(nr3+im)
         ekph = ekph + 100.0*v(nr3+im)**2
      end do
      
      eke = eke * 0.5d0
      ekph = ekph * 0.5d0
      ekpbs = ekpbs * 0.5d0
      
      if( ntt == 1 ) then
         
         !           --- following is from T.E. Cheatham, III and B.R. Brooks,
         !               Theor. Chem. Acc. 99:279, 1998.
         
         scaltp = sqrt(1.d0 + 2.d0*dttp*(ekin0-eke)/(ekmh+ekph))

         do i3 = istart3,iend3
            v(i3) = v(i3)*scaltp
         end do
         do im=1,0
            v(nr3+im) = v(nr3+im)*scaltp
         end do
      end if  ! (ntt == 1 )
      
   end if  ! ( ntt == 1 .or. onstep; end of step 4c )
   
   !-----------------------------------------------------------------
   !  --- put current velocities into VOLD
   !-----------------------------------------------------------------
   
   vold(istart3:iend3) = v(istart3:iend3)
   do im=1,0
      vold(nr3+im) = v(nr3+im)
   end do
   
   !-----------------------------------------------------------------
   !   Step 5:  several tasks related to dumping of trajectory information
   !-----------------------------------------------------------------

   !  --- determine if restart writing is imminent and
   !      requires xdist of v and dipole information in parallel runs:
   
   ivarch = .false.
   ixdump = .false.
   itdump = .false.
   if( ntwv > 0 ) then
      if( mod( nstep+1, ntwv ) == 0 ) ivarch = .true.
   end if
   if( ntwr /= 0 ) then
      if( mod( nstep+1, ntwr ) == 0 ) then
         ivarch = .true.
         ixdump = .true.
      end if
   end if
   if( nstep+1 >= nstlim ) then
      ivarch = .true.
      ixdump = .true.
   end if
   if( mod(nstep+1,nscm) == 0 ) ivarch =.true.
   
   !  --- determine if trajectory writing is imminent:
   
   if(ntwx > 0)then
      itdump= mod(nstep+1,ntwx) == 0
   else
      itdump=.false.
   end if

#ifdef MPI

   !-----------------------------------------------------------------
   !  --- now distribute the coordinates, and if necessary, dipoles and vel:
   !-----------------------------------------------------------------

   call timer_barrier( commsander )
   call timer_stop_start(TIME_VERLET,TIME_DISTCRD)

   if ( .not. mpi_orig .and. numtasks > 1 ) call xdist(x)

      call timer_stop(TIME_DISTCRD)
#endif  /* MPI */
   
   !-------------------------------------------------------------------
   !   Step 6: zero COM velocity if requested; used for preventing
   !   ewald "block of ice flying thru space" phenomenon, or accumulation
   !   of rotational momentum in vacuum simulations
   !-------------------------------------------------------------------
   
   if (mod(nstep+1,nscm) == 0) then
      if (mod(nstep,nsnb) == 0) ntnb = 1
      
      if( ifbox == 0 ) then
        
         !    ---Non-periodic simulation: remove both translation and rotation.
         !       Back the coords up 1/2 step, so that the correspond to the
         !       velocities; temporarily store in the F() array:
         
         f(1:nr3) = x(1:nr3) - v(1:nr3)*dt5
         
         !     --- now compute the com motion, remove it, and recompute (just
         !         to check that it is really gone.....)
         
         call cenmas(natom,f,v,tmass,tmassinv,amass,ekcm, &
               xcm,vcm,acm,ekrot,ocm,4)
         call stopcm(natom,f,v,xcm,vcm,ocm)
         call cenmas(natom,f,v,tmass,tmassinv,amass,ekcm, &
               xcm,vcm,acm,ekrot,ocm,4)
         
      else
         
         !    ---Periodic simulation: just remove the translational velocity:
         
         vcmx = 0.d0
         vcmy = 0.d0
         vcmz = 0.d0

         j = 1
         do i = 1, 3*natom,3
            aamass = amass(j)
            vcmx = vcmx + aamass * v(i)
            vcmy = vcmy + aamass * v(i+1)
            vcmz = vcmz + aamass * v(i+2)
            j = j + 1
         end do
         
         vcmx = vcmx * tmassinv
         vcmy = vcmy * tmassinv
         vcmz = vcmz * tmassinv
         
         vel2 = vcmx*vcmx + vcmy*vcmy + vcmz*vcmz
         atempdrop = 0.5d0 * tmass * vel2 / fac(1)
         vel = sqrt(vel2)
         if ( master ) &
               write (6,'(a,f15.6,f9.2,a)') &
               'check COM velocity, temp: ',vel,atempdrop, &
               '(Removed)'

         do i = 1, 3*natom, 3
            v(i)   = v(i)   - vcmx
            v(i+1) = v(i+1) - vcmy
            v(i+2) = v(i+2) - vcmz
         end do
      end if  ! ( ifbox == 0 )
   end if  ! (mod(nstep+1,nscm) == 0)
   
   !  Also zero out the non-moving velocities if a belly is active:
   if (belly) call bellyf(nr,ix(ibellygp),v)
   
   !-------------------------------------------------------------------
   !  Step 7: scale coordinates if constant pressure run:
   !-------------------------------------------------------------------
   
   ener(4) = ekpbs + ener(23)  ! Pastor, Brooks, Szabo conserved quantity
                               ! for harmonic oscillator: Eq. 4.7b of Mol.
                               ! Phys. 65:1409-1419, 1988
   ener(3) = eke
   ener(2) = ener(3)
   if (ntt == 1 .and. onstep) then
      ekmh = max(ekph,fac(1)*10.d0)
   end if
   
   !     ---if velocities were reset, the KE is not accurate; fudge it
   !        here to keep the same total energy as on the previous step.
   !        Note that this only affects printout and averages for Etot
   !        and KE -- it has no effect on the trajectory, or on any averages
   !        of potential energy terms.
   
   if( resetvelo ) ener(2) = etot_save - ener(23)
   
   !     --- total energy is sum of KE + PE:
   
   ener(1) = ener(2)+ener(23)
   etot_save = ener(1)
   
   !-------------------------------------------------------------------
   !  Step 8:  update the step counter and the integration time:
   !-------------------------------------------------------------------
   
   nstep = nstep+1
   t = t+dt
   
   !     ---full energies are only calculated every nrespa steps
   !     nvalid is the number of steps where all energies are calculated
   
   if ( onstep )then
      nvalid = nvalid + 1
      enert(1:nren) = enert(1:nren)+ener(1:nren)
      enert2(1:nren) = enert2(1:nren) + ener(1:nren)**2
      if( nvalid == 1 ) etot_start = ener(1)
   end if
     
   ntnb = 0
   if (mod(nstep,nsnb) == 0) ntnb = 1
   lout = mod(nstep,ntpr) == 0 .and. onstep
   irespa = irespa + 1
     
   ! reset pb-related flags
#ifdef MPI
   if(mytaskid == 0)then
#endif   
      if ( igb >= 10 ) then
         if ( mod(nstep+1,npbgrid) == 0 .and. nstep+1 /= nstlim ) pbgrid = .true.
         if ( mod(nstep+1,ntpr) == 0 .or. nstep+1 == nstlim ) pbprint = .true.
         if ( mod(nstep+1,nsnbr) == 0 .and. nstep+1 /= nstlim ) ntnbr = 1
         if ( mod(nstep+1,nsnba) == 0 .and. nstep+1 /= nstlim ) ntnba = 1
      end if
#ifdef MPI
   endif
#endif   

   !-------------------------------------------------------------------
   !  Step 9:  output from this step if required:
   !-------------------------------------------------------------------
   
   !     ...only the master needs to do the output
   
   if (master) then
      
      !        -- restrt:
      
      if (ixdump) then
         if( iwrap == 0 ) then
            nr = nrp
            call mdwrit(nstep,nrp,nr,nres,ntxo,ntr,ntb, &
                 x,v,xx(lcrdr),box,ih(m04),ih(m02),ix(i02),t)
         else
            
            ! --- use temp. array to hold coords. so that the master's values
            !     are always identical to those on all other nodes:
            
            call get_stack(l_temp,3*natom)
            do m=1,3*natom
               r_stack(l_temp + m - 1) = x(m)
            end do
            
            nr = nrp
            call mdwrit(nstep,nrp,nr,nres,ntxo,ntr,ntb, &
                  r_stack(l_temp),v,xx(lcrdr),box,ih(m04),ih(m02),ix(i02),t)
            call free_stack(l_temp)
         end if  ! ( iwrap == 0 )

         if (icnstph /= 0) &
            call cnstphwriterestart(ix(icpstinf),ix(icpresst),ix(icpptcnt), &
                  ix(icptrsct), xx(lcpene),xx(lcpcrg))

      end if  ! (ixdump)
      
      !     -- Coordinate archive:
      
      if (itdump) then
         if( iwrap == 0 ) then
            call corpac(x,1,nrx,12,loutfm)
         else
            call get_stack(l_temp,3*natom)
            do m=1,3*natom
               r_stack(l_temp + m - 1) = x(m)
            end do
            call corpac(r_stack(l_temp),1,nrx,12,loutfm)
            call free_stack(l_temp)
         end if
      end if  ! (itdump)
      
      !     Velocity archive:
      
      if (ntwv > 0) then
         if (mod(nstep,ntwv) == 0) call corpac(v,1,nrx,13,loutfm)
      end if
      
      !     Energy archive:
      
      if (ntwe > 0) then
         if (mod(nstep,ntwe) == 0.and.onstep) &
               call mdeng(15,nstep,t,ener,fac,ntp)
      end if
      
      !     General printed output:
      
      if (lout) then
         rewind(7)
         call prntmd(nstep,nitp,nits,t,ener,fac,7,.false.)
         call amflsh(7)
      end if
      
      !     output running average
      
      if ( ntave > 0 )then
         if ( mod(nstep,ntave) == 0 .and. onstep )then
            write(6,542)
            tspan = ntave/nrespa
            do m = 1,nren
               enert_tmp(m) = enert(m)-enert_old(m)
               enert2_tmp(m) = enert2(m)-enert2_old(m)
               enert_old(m) = enert(m)
               enert2_old(m) = enert2(m)
               enert_tmp(m) = enert_tmp(m)/tspan
               enert2_tmp(m) = enert2_tmp(m)/tspan - &
                     enert_tmp(m)*enert_tmp(m)
               if ( enert2_tmp(m) < 0.d0 )enert2_tmp(m) = 0.d0
               enert2_tmp(m) = sqrt(enert2_tmp(m))
            end do
            write(6,540) ntave/nrespa
            call prntmd(nstep,izero,izero,t,enert_tmp,fac,0,.false.)
            write(6,550)
            call prntmd(nstep,izero,izero,t,enert2_tmp,fac,0,.true.)
            if( icfe > 0 ) then
               write(6,541) ntave/nrespa
               edvdl_r(1:51) = edvdl_r(1:51)/tspan
               call prntmd(nstep,izero,izero,t,edvdl_r,fac,0,.false.)
               edvdl_r(1:51) = 0.d0
            end if
            write(6,542)
         end if
      end if  ! ( ntave > 0 )
      
      !     --- end masters output ---
      
   end if  ! (master)
   
   !=======================================================================
   
   !  ---major cycle back to new step unless we have reached our limit:
   
   call timer_stop(TIME_VERLET)
   if (nstep < nstlim) goto 260
   480 continue
   
   !=======================================================================
   
   !     ----- PRINT AVERAGES -----
   
   !=======================================================================
   
   
   if (master) then
      tspan = nvalid
      if (nvalid > 0) then
         do m = 1,nren
            enert(m) = enert(m)/tspan
            enert2(m) = enert2(m)/tspan - enert(m)*enert(m)
            if(enert2(m) < 0.d0) enert2(m) = 0.d0
            enert2(m) =  sqrt(enert2(m))
            edvdl(m) = edvdl(m)/tspan
         end do
         
         write(6,540) nvalid
         call prntmd(nstep,izero,izero,t,enert,fac,0,.false.)
         write(6,550)
         call prntmd(nstep,izero,izero,t,enert2,fac,0,.true.)
         
         if( icfe > 0 ) then
            write(6,541) nvalid
            call prntmd(nstep,izero,izero,t,edvdl,fac,0,.false.)
         end if
         
         !          print Born radii statistics
         
         if ((rbornstat == 1).and.(igb /= 0)) then
            write(6,580) nstep
            write(6,590)
            do m = 1,natom
               xx(l188-1+m) = xx(l188-1+m)/tspan
               xx(l189-1+m) = xx(l189-1+m)/tspan - &
                     xx(l188-1+m)*xx(l188-1+m)
               xx(l189-1+m) = sqrt(xx(l189-1+m))
               write(6,600) m, xx(l186-1+m), xx(l187-1+m), &
                     xx(l188-1+m), xx(l189-1+m)
            end do
         end if
         
         do m = 2,nrek
            enert(m) = enert(m)/fac(m-1)
            enert2(m) = enert2(m)/fac(m-1)
         end do
         temp = enert(2)
      end if  ! (nvalid > 0)
      
   end if  ! (master)
   
   540 format(/5x,' A V E R A G E S   O V E R ',i7,' S T E P S',/)
   541 format(/5x,' DV/DL, AVERAGES OVER ',i7,' STEPS',/)
   542 format('|',79('='))
   550 format(/5x,' R M S  F L U C T U A T I O N S',/)
   580 format('STATISTICS OF EFFECTIVE BORN RADII OVER ',i7,' STEPS')
   590 format('ATOMNUM     MAX RAD     MIN RAD     AVE RAD     FLUCT')
   600 format(i4,2x,4f12.4)
   return
end subroutine runmd 
