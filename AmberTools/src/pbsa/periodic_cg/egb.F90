! <compile=optimized>
#include "copyright.h"
#include "assert.h"
#include "dprec.h"

module genborn

_REAL_, private, dimension(:), allocatable :: r2x,rjx,vectmp1,vectmp2, &
                             vectmp3,vectmp4,vectmp5,sumdeijda,psi
integer, private, dimension(:), allocatable :: jj,k_vals,j_vals
logical, private, dimension(:), allocatable :: skipv

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ allocates scratch space for egb()
subroutine allocate_gb( natom )

   implicit none
   integer, intent(in) :: natom
   integer ier                     

   allocate( r2x(natom), rjx(natom), vectmp1(natom), vectmp2(natom),   &
             vectmp3(natom), vectmp4(natom), vectmp5(natom),           &
             sumdeijda(natom), psi(natom), jj(natom),  &
             skipv(0:natom),k_vals(natom),j_vals(natom), stat = ier )
   REQUIRE( ier == 0 )
   return

end subroutine allocate_gb


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ deallocates scratch space for egb()
subroutine deallocate_gb( )

   implicit none
   integer ier

   ! assume that if r2x is allocated then all are allocated
   if ( allocated( r2x ) ) then
      deallocate( j_vals, k_vals, skipv, jj, psi, sumdeijda, vectmp5, &
            vectmp4, vectmp3, vectmp2, vectmp1, rjx, r2x, stat = ier )
      REQUIRE( ier == 0 )
   else
      ASSERT( .false. )  ! cannot deallocate un-allocated array
   end if
   return

end subroutine deallocate_gb


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ handles generalized Born functionality, plus reg. nonbon, plus surface area
subroutine egb(x,xx,ix,f,rborn,fs,reff,charge,iac,ico,numex, &
      natex,dcharge,cut,ntypes,natom,natbel, &
      epol,eelt,evdw,esurf,dvdl,vdwrad,ineighbor,p1,p2,p3,p4, &
      rbmax,rbmin,rbave,rbfluct)

   
   !--------------------------------------------------------------------------
   
   !     Compute nonbonded interactions with a generalized Born model,
   !     getting the "effective" Born radii via the approximate pairwise method
   !     Use Eqs 9-11 of Hawkins, Cramer, Truhlar, J. Phys. Chem. 100:19824
   !     (1996).  Aside from the scaling of the radii, this is the same
   !     approach developed in Schaefer and Froemmel, JMB 216:1045 (1990).
   
   !     The input coordinates are in the "x" array, and the forces in "f"
   !     get updated; energy components are returned in "epol", "eelt" and
   !     "evdw".
   
   !     Input parameters for the generalized Born model are "rborn(i)", the
   !     intrinsic dielectric radius of atom "i", and "fs(i)", which is
   !     set (in routine mdread) to (rborn(i) - offset)*si.
   
   !     Input parameters for the "gas-phase" electrostatic energies are
   !     the charges, in the "charge()" array.
   
   !     Input parameters for the van der Waals terms are "cn1()" and "cn2()",
   !     containing LJ 12-6 parameters, and "asol" and "bsol" containing
   !     LJ 12-10 parameters.  (The latter are not used in 1994 are more
   !     forcefields.)  The "iac" and "ico" arrays are used to point into
   !     these matrices of coefficients.
   
   !     The "numex" and "natex" arrays are used to find "excluded" pairs of
   !     atoms, for which gas-phase electrostatics and LJ terms are skipped;
   !     note that GB terms are computed for all pairs of atoms.
   
   !     If gbsa=1, then an approximate surface-area dependent term is
   !     computed, with the resulting energy placed into "esurf".  The
   !     algorithm is from J. Weiser, P.S. Shenkin, and W.C. Still,
   !     "Approximate atomic sufraces from linear combinations of pariwise
   !     overlaps (LCPO)", J. Computat. Chem. 20:217 (1999).
   
   !     The code also supports a multiple-time-step facility:
   
   !       pairs closer than sqrt(cut_inner) are evaluated every nrespai steps;
   !         "   between sqrt(cut_inner) and sqrt(cut) are evaluated
   !                        every nrespa steps
   !         "   beyond sqrt(cut) are ignored
   
   !       the forces arising from the derivatives of the GB terms with respect
   !          to the effective Born radii are evaulated every nrespa steps
   
   !       the surface-area dependent term is evaluated every nrespa steps
   
   !       the effective radii are only updated every nrespai steps
   
   !     (Be careful with the above: what seems to work is dt=0.001,
   !     nrespai=2, nrespa=4; anything beyond this seems dangerous.)
   
   !     Written 1999-2000, primarily by D.A. Case, with help from C. Brooks,
   !       T. Simonson, R. Sinkovits  and V. Tsui.  The LCPO implementation
   !       was written by V. Tsui.
   
   !     Vectorization and optimization 1999-2000, primarily by C. P. Sosa,
   !       T. Hewitt, and D. A. Case.  Work presented at CUG Fall of 2000.
   !--------------------------------------------------------------------------
   
   use icosasurf, only : icosa_init, icosa_sphere_approx
   use decomp, only: decsasa, decpair
   implicit none
   _REAL_ fasti,slowi
   
#ifdef MPI
#  include "parallel.h"
#  ifdef MPI_DOUBLE_PRECISION
#    undef MPI_DOUBLE_PRECISION
#  endif
#  include "mpif.h"
#  ifdef CRAY_PVP
#    define MPI_DOUBLE_PRECISION MPI_REAL8
#  endif
#endif
#  include "md.h"
#  include "parms.h"
#  include "def_time.h"
#  include "constants.h"
   
   logical onstep,onstepi,oncpstep
   
   _REAL_ x,xx,f,rborn,fs,reff,charge,cut, &
         epol,eelt,evdw,esurf,vdwrad,p1,p2,p3,p4, &
         dcharge,rbmax,rbmin,rbave,rbfluct
   _REAL_ totsasa,cutxyz,extdieli,intdieli, &
         lastxj,lastyj,lastzj,xi,yi,zi,ri,ri1i,sumi,xij,yij,zij, &
         dij1i,r2,dij,sj,sj2,theta,uij,frespa,si,sumaij,sumajk, &
         sumaijajk,sumdaijddijdxi,sumdaijddijdyi,sumdaijddijdzi, &
         sumdaijddijdxiajk,sumdaijddijdyiajk,sumdaijddijdziajk, &
         xj,yj,zj,rij,tmpaij,aij,daijddij,daijddijdxj, daijddijdyj, &
         daijddijdzj,sumajk2,sumdajkddjkdxj,sumdajkddjkdyj, &
         sumdajkddjkdzj,p3p4aij,xk,yk,zk,rjk2,djk1i,rjk,vdw2dif, &
         tmpajk,ajk,dajkddjk,dajkddjkdxj,dajkddjkdyj,dajkddjkdzj, &
         daidxj,daidyj,daidzj,ai,daidxi,daidyi,daidzi,qi,dumx, &
         dumy,dumz,de,rj,temp1,fgbi,rinv,r2inv,qiqj,rb2,fgbk,expmkf, &
         dl,e,temp4,temp5,temp6,eel,r6inv,f6,f12,r10inv,f10, &
         dedx,dedy,dedz,qi2h,temp7,dij2i,datmp,dij3i, &
         qid2h,si2,rj1i,dvdl,reff_i,thi,thi2,log_r3, &
         beta,gamma,gamfac,fgbiorig
   
   integer count,count2,icount,nearest,ineighbor(*),max_count
   integer ix,iac,ico,numex,natex,ntypes,natom,natbel
   integer i,j,k,kk1,kk2,maxi,num_j_vals,jjj,count2_fin,num_k_vals, &
           iexcl,iaci,jexcl,jexcl_last,jjv,ic,kk
   integer icase, j1, j2, j3
   _REAL_ f_x,f_y,f_z,f_xi,f_yi,f_zi
   _REAL_ dumbo, tmpsd

   ! variables needed for icosa surface area calculation
   integer ineighborpt

   ! variables needed for smooth integration cutoff in Reff:
   _REAL_ rgbmax2, rgbmax1i, rgbmax2i,tmpcs,rgbmaxpsmax2
   !     _REAL_  datmp2
   
   dimension x(*),xx(*),ix(*),f(*),rborn(*),charge(*),iac(*), &
         ico(*),numex(*),natex(*),fs(*),reff(*), &
         vdwrad(*),p1(*),p2(*),p3(*),p4(*), &
         dcharge(*),rbmax(*),rbmin(*),rbave(*),rbfluct(*)
   
   !   FGB taylor coefficients follow
   !   from A to H :
   !   1/3 , 2/5 , 3/7 , 4/9 , 5/11
   !   4/3 , 12/5 , 24/7 , 40/9 , 60/11

   _REAL_  ta
   _REAL_  tb
   _REAL_  tc
   _REAL_  td
   _REAL_  tdd
   _REAL_  te
   _REAL_  tf
   _REAL_  tg
   _REAL_  th
   _REAL_  thh
   parameter ( ta   = third )
   parameter ( tb   = two / five )
   parameter ( tc   = three / seven )
   parameter ( td   = four / nine )
   parameter ( tdd  = five / eleven )
   parameter ( te   = four / three )
   parameter ( tf   = twelve / five )
   parameter ( tg   = three * eight / seven )
   parameter ( th   = five * eight / nine )
   parameter ( thh  = three * four * five / eleven )

   epol = zero
   eelt = zero
   evdw = zero
   esurf = zero
   totsasa = zero
   
   onstep = mod(irespa,nrespa) == 0
   onstepi = mod(irespa,nrespai) == 0
   oncpstep = icnstph /= 0 .and. mod(irespa, ntcnstph) == 0
   
   if( .not.onstepi ) return
   
   cutxyz = sqrt( cut )
   extdieli = one/extdiel
   intdieli = one/intdiel
   maxi = natom
   if(natbel > 0) maxi = natbel

   ! Smooth "cut-off" in calculating GB effective radii.  
   ! Implementd by Andreas Svrcek-Seiler and Alexey Onufriev. 
   ! The integration over solute is performed up to rgbmax and includes 
   ! parts of spheres; that is an atom is not just "in" or "out", as 
   ! with standard non-bonded cut.  As a result, calclated effective 
   ! radii are less than rgbmax. This saves time, and there is no 
   ! discontinuity in dReff/drij.

   ! Only the case rgbmax > 5*max(sij) = 5*fsmax ~ 9A is handled; this is 
   ! enforced in mdread().  Smaller values would not make much physical
   ! sense anyway.
      
   rgbmax2 = rgbmax*rgbmax
   rgbmax1i = 1.d0/rgbmax
   rgbmax2i = rgbmax1i*rgbmax1i
   rgbmaxpsmax2 = (rgbmax+fsmax)**2

   !---------------------------------------------------------------------------
   !      Step 1: loop over pairs of atoms to compute the effective Born radii.
   !---------------------------------------------------------------------------
   
   if( irespa > 1 .and. mod(irespa,nrespai) /= 0 ) goto 90
   call timer_start(TIME_GBRAD1)
   
   reff(1:natom) = zero
   
#ifdef MPI
   do i=mytaskid+1,natom,numtasks
#else
   do i=1,natom
#endif
      xi = x(3*i-2)
      yi = x(3*i-1)
      zi = x(3*i)
      
      reff_i = reff(i)
      ri = rborn(i)-offset
      ri1i = one/ri
      si = fs(i)
      si2 = si*si
      
      !  Here, reff_i will sum the contributions to the inverse effective
      !  radius from all of the atoms surrounding atom "i"; later the
      !  inverse of its own intrinsic radius will be added in
      
      icount = 0
      do j=i+1,natom
         
         xij = xi - x(3*j-2)
         yij = yi - x(3*j-1)
         zij = zi - x(3*j  )
         r2 = xij*xij + yij*yij + zij*zij
         if( r2 > rgbmaxpsmax2 ) cycle
         
         icount = icount + 1
         jj(icount) = j
         r2x(icount) = r2
         
      end do
      
      call vdinvsqrt( icount, r2x, vectmp1 )
      
      kk1 = 0
      kk2 = 0
      !dir$ ivdep
      do k = 1, icount
         
         
         j = jj(k)
         r2 = r2x(k)
         sj =  fs(j)

!        don't fill the remaining vectmp arrays if atoms don't see each other:
         dij1i = vectmp1(k)
         dij = r2*dij1i
         if (dij > rgbmax+si .and. dij > rgbmax+sj) cycle
         rj = rborn(j) - offset

         if( dij <= 4.d0*sj) then
            kk1 = kk1 + 1
            vectmp2(kk1) = dij + sj
            if( dij > ri+sj ) then
               vectmp4(kk1) = dij - sj
            else if ( dij > abs(ri-sj) ) then
               vectmp4(kk1) = ri
            else if ( ri < sj ) then
               vectmp4(kk1) = sj - dij
            else
               vectmp4(kk1) = one
            end if
         end if
         
         if( dij <= 4.d0*si) then
            kk2 = kk2 + 1
            vectmp3(kk2) = dij + si
            if( dij > rj+si) then
               vectmp5(kk2) = dij - si
            else if ( dij > abs(rj-si) ) then
               vectmp5(kk2) = rj
            else if ( rj < si ) then
               vectmp5(kk2) = si - dij
            else
               vectmp5(kk2) = one
            end if
         end if
         
      end do  !  k = 1, icount
      
      call vdinv( kk1, vectmp2, vectmp2 )
      call vdinv( kk2, vectmp3, vectmp3 )
      vectmp4(1:kk1) = vectmp2(1:kk1)*vectmp4(1:kk1)
      vectmp5(1:kk2) = vectmp3(1:kk2)*vectmp5(1:kk2)
      call vdln( kk1, vectmp4, vectmp4 )
      call vdln( kk2, vectmp5, vectmp5 )
      
      kk1 = 0
      kk2 = 0
      do k = 1, icount
         
         j = jj(k)
         r2 = r2x(k)
         
         rj = rborn(j) - offset
         rj1i = one/rj
         sj =  fs(j)
         
         sj2 = sj * sj
         
         xij = xi - x(3*j-2)
         yij = yi - x(3*j-1)
         zij = zi - x(3*j  )
         
         dij1i = vectmp1(k)
         dij = r2*dij1i
         
         
         if (dij <= rgbmax + sj) then
         
            if ((dij > rgbmax - sj)) then
               uij = 1./(dij -sj)
               reff_i = reff_i - 0.125d0 * dij1i * (1.d0 + 2.d0 * dij *uij + &
                     rgbmax2i * (r2 - 4.d0 * rgbmax * dij - sj2) + &
                     2.d0 * log((dij-sj)*rgbmax1i))

             else if( dij > 4.d0*sj ) then
            
               dij2i = dij1i*dij1i
               tmpsd = sj2*dij2i
               dumbo = ta+tmpsd* (tb+tmpsd* (tc+tmpsd* (td+tmpsd* tdd)))

               reff_i = reff_i - tmpsd*sj*dij2i*dumbo
            
            !     ---following are from the Appendix of Schaefer and Froemmel,
            !        J. Mol. Biol. 216:1045-1066, 1990, divided by (4*Pi):

            else if( dij > ri+sj ) then
            
               kk1 = kk1 + 1
               reff_i = reff_i - half*( sj/(r2-sj2) + half*dij1i*vectmp4(kk1) )
            
            !-----------------------------------------------------------------
            
            else if ( dij > abs(ri-sj) ) then
            
               kk1 = kk1 + 1
               theta = half*ri1i*dij1i*(r2 + ri*ri -sj2)
               reff_i = reff_i - fourth*( ri1i*(two-theta) &
                  - vectmp2(kk1) + dij1i*vectmp4(kk1) )
            
            !-----------------------------------------------------------------
            
            else if ( ri < sj ) then
               kk1 = kk1 + 1
               reff_i = reff_i - half*( sj/(r2-sj2) + two*ri1i &
                  + half*dij1i*vectmp4(kk1) )
            
            !-----------------------------------------------------------------
            
            else
               kk1 = kk1 + 1
            end if  ! ( dij > 4.d0*sj )
         
         end if

         !   --- Now the same thing, but swap i and j:
         
         if (dij > rgbmax +si) cycle
             
         if (dij > rgbmax - si) then
            uij = 1./(dij -si)
            reff(j) = reff(j) - 0.125d0 * dij1i * (1.d0 + 2.d0 * dij *uij + &
                     rgbmax2i * (r2 - 4.d0 * rgbmax * dij - si2) + &
                     2.d0 * log((dij-si)*rgbmax1i))

         else if( dij > 4.d0*si ) then
            
            dij2i = dij1i*dij1i
            tmpsd = si2*dij2i
            dumbo = ta+tmpsd* (tb+tmpsd* (tc+tmpsd* (td+tmpsd* tdd)))
            reff(j) = reff(j) - tmpsd*si*dij2i*dumbo
            
         else if( dij > rj+si ) then
            
            kk2 = kk2 + 1
            reff(j) = reff(j) - half*( si/(r2-si2) + &
                  half*dij1i*vectmp5(kk2) )
            
            !-----------------------------------------------------------------
            
         else if ( dij > abs(rj-si) ) then
            
            kk2 = kk2 + 1
            theta = half*rj1i*dij1i*(r2 + rj*rj -si2)
            reff(j) = reff(j) - fourth*( rj1i*(two-theta) &
                  - vectmp3(kk2) + dij1i*vectmp5(kk2) )
            
            !-----------------------------------------------------------------
            
         else if ( rj < si ) then
            
            kk2 = kk2 + 1
            reff(j) = reff(j) - half*( si/(r2-si2) + two*rj1i &
                  + half*dij1i*vectmp5(kk2) )
            
            !-----------------------------------------------------------------
            
         else
            kk2 = kk2 + 1
         end if  ! ( dij > 4.d0*si )
      end do                    !  k = 1, icount
      
      ! we are ending the do-i-loop, reassign the scalar to the original array:
      
      reff(i) = reff_i
      
   end do  !  i = 1,natom
   
#ifdef MPI
   call timer_stop_start(TIME_GBRAD1,TIME_GBRADDIST)
   
   !       collect the (inverse) effective radii from other nodes:
   
   if( numtasks > 1 ) then
      call mpi_allreduce( reff, vectmp1, natom, &
            MPI_DOUBLE_PRECISION,mpi_sum,commsander,ierr)
      reff(1:natom) = vectmp1(1:natom)
   end if
   call timer_stop_start(TIME_GBRADDIST,TIME_GBRAD1)
#endif
   
   if( igb == 2 .or. igb == 5) then
      
      !       apply the new Onufriev "gbalpha, gbbeta, gbgamma" correction:
      
      do i=1,natom
         ri = rborn(i)-offset
         ri1i = one/ri
         psi(i) = -ri*reff(i)
         reff(i) = ri1i - tanh( (gbalpha + gbgamma*psi(i)*psi(i) &
               - gbbeta*psi(i) )*psi(i) )/rborn(i)
         reff(i) = one/reff(i)
      end do
      
   else
      
      !       "standard" GB, including the "diagonal" term here:
      
      do i=1,natom
         ri = rborn(i)-offset
         ri1i = one/ri
         reff(i) = one/(reff(i) + ri1i)
      end do
   end if
   
   if ( rbornstat == 1 ) then
      do i=1,natom
         rbave(i) = rbave(i) + reff(i)
         rbfluct(i) = rbfluct(i) + reff(i)*reff(i)
         if ( rbmax(i) <= reff(i) ) rbmax(i) = reff(i)
         if ( rbmin(i) >= reff(i) ) rbmin(i) = reff(i)
      end do
   end if
   call timer_stop(TIME_GBRAD1)
   
   90 continue
   
   !--------------------------------------------------------------------------
   
   !     Step 2: loop over all pairs of atoms, computing the gas-phase
   !             electrostatic energies, the LJ terms, and the off-diagonal
   !             GB terms.  Also accumulate the derivatives of these off-
   !             diagonal terms with respect to the inverse effective radii,
   !             sumdeijda(k) will hold  sum over i,j>i ( deij/dak ),  where
   !             "ak" is the inverse of the effective radius for atom "k".
   
   !             Update the forces with the negative derivatives of the
   !             gas-phase terms, plus the derivatives of the explicit
   !             distance dependence in Fgb, i.e. the derivatives of the
   !             GB energy terms assuming that the effective radii are
   !             constant.
   
   !--------------------------------------------------------------------------
   
   sumdeijda(1:natom) = zero
   call timer_start(TIME_GBFRC)
   
   !       Note: this code assumes that the belly atoms are the first natbel
   !             atoms...this is checked in mdread.
   
#ifdef MPI
   iexcl = 1
   do i=1,mytaskid
      iexcl = iexcl + numex(i)
   end do
   do i=mytaskid+1,maxi,numtasks
#else
   iexcl = 1
   do i=1,maxi
#endif
      xi = x(3*i-2)
      yi = x(3*i-1)
      zi = x(3*i  )
      qi = charge(i)
      ri = reff(i)
      iaci = ntypes * (iac(i) - 1)
      jexcl = iexcl
      jexcl_last = iexcl + numex(i) -1
      dumx = zero
      dumy = zero
      dumz = zero
      
      !         -- check the exclusion list for eel and vdw:
      
      do k=i+1,natom
         skipv(k) = .false.
      end do
      do jjv=jexcl,jexcl_last
         skipv(natex(jjv))=.true.
      end do
      
      icount = 0
      do j=i+1,natom
         
         xij = xi - x(3*j-2)
         yij = yi - x(3*j-1)
         zij = zi - x(3*j  )
         r2 = xij*xij + yij*yij + zij*zij
         if( r2 > cut ) cycle
         if( .not.onstep .and. r2 > cut_inner ) cycle
         
         icount = icount + 1
         jj(icount) = j
         r2x(icount) = r2
         rjx(icount) = reff(j)
         
      end do
      
      vectmp1(1:icount) = -r2x(1:icount)/(four*ri*rjx(1:icount))
      call vdexp( icount, vectmp1, vectmp1 )
      vectmp3(1:icount) = r2x(1:icount) + rjx(1:icount)*ri*vectmp1(1:icount)
      call vdinvsqrt( icount, vectmp3, vectmp2 )
      call vdinvsqrt( icount, r2x, vectmp3 )
      
      if( kappa /= zero ) then
         vectmp4(1:icount) = -kappa/vectmp2(1:icount)
         call vdexp( icount, vectmp4, vectmp4 )
      end if
      
      !dir$ ivdep
      do k=1,icount
         j = jj(k)
         de = zero
         r2 = r2x(k)
         rj = rjx(k)
         temp1 = vectmp1(k)
         fgbi = vectmp2(k)
         rinv = vectmp3(k)
         
         xij = xi - x(3*j-2)
         yij = yi - x(3*j-1)
         zij = zi - x(3*j  )
         
         frespa = nrespai
         if( onstep .and. r2 > cut_inner ) frespa = nrespa
         
         r2inv = rinv*rinv
         
         qiqj = qi * charge(j)
         rb2 = ri*rj
         if( kappa == zero ) then
            fgbk = zero
            expmkf = extdieli
         else
            fgbk = -kappa/fgbi
            expmkf = vectmp4(k)*extdieli
         end if
         dl = intdieli - expmkf
         
         e = -qiqj*dl*fgbi
         epol = epol + e
         if(idecomp == 1 .or. idecomp == 2) then
            call decpair(1,i,j,e)
         else if(idecomp == 3 .or. idecomp == 4) then
            call decpair(-1,i,j,e)
         end if
#ifdef CHARGE_MIXING
         if( onstep .and. icfe /= 0 ) then
            dvdl = dvdl - dl*fgbi*(charge(i)*dcharge(j) &
                  + charge(j)*dcharge(i) )
         end if
#endif
         if (oncpstep) then
            dvdl = dvdl - (dl*fgbi*dcharge(i) &
                  *dcharge(j)+e)
         end if
         
         temp4 = fgbi*fgbi*fgbi
         
         !   [here, and in the gas-phase part, "de" contains -(1/r)(dE/dr)]
         
         temp6 = -qiqj*temp4*(dl + fgbk*expmkf)
         de = temp6*(one - fourth*temp1)
         temp5 = half*temp1*temp6*(rb2 + fourth*r2)
         
         sumdeijda(i) = sumdeijda(i) + ri*temp5
         sumdeijda(j) = sumdeijda(j) + rj*temp5
         
         !    -- skip exclusions for remaining terms:
         
         if( .not. skipv(j) ) then
            
            !   -- gas-phase Coulomb energy:
            
            eel = intdieli*qi*charge(j)*rinv
            eelt = eelt + eel
            if(idecomp == 1 .or. idecomp == 2) then
               call decpair(2,i,j,eel)
            else if(idecomp == 3 .or. idecomp == 4) then
               call decpair(-2,i,j,eel)
            end if
            de = de + eel*r2inv
#ifdef CHARGE_MIXING
            if( onstep .and. icfe /= 0 ) then
               dvdl = dvdl + intdieli*rinv*(charge(i)*dcharge(j) &
                     + charge(j)*dcharge(i) )
            end if
#endif
            if (oncpstep) then
               dvdl = dvdl + (intdieli*rinv*dcharge(i) &
                     *dcharge(j) - eel)

            end if
            
            !    -- van der Waals energy:
            
            ic = ico( iaci + iac(j) )
            if( ic > 0 ) then
               !                                    6-12 potential:
               r6inv = r2inv*r2inv*r2inv
               f6 = cn2(ic)*r6inv
               f12 = cn1(ic)*(r6inv*r6inv)
               evdw = evdw + (f12 - f6)
               if(idecomp == 1 .or. idecomp == 2) then
                  call decpair(3,i,j,f12-f6)
               else if(idecomp == 3 .or. idecomp == 4) then
                  call decpair(-3,i,j,f12-f6)
               end if
               de = de + (twelve*f12 - six*f6)*r2inv
               
#ifdef HAS_10_12
               
               !    ---The following could be commented out if the Cornell
               !       et al. force field was always used, since then all hbond
               !       terms are zero.
               
            else
               !                                    10-12 potential:
               r10inv = r2inv*r2inv*r2inv*r2inv*r2inv
               f10 = bsol(-ic)*r10inv
               f12 = asol(-ic)*r10inv*r2inv
               !     --- put old hbond term into epol, where it "used to" be:
               evdw = evdw + f12 -f10
               if(idecomp == 3 .or. idecomp == 2) then
                  call decpair(1,i,j,f12-f10)
               else if(idecomp == 3 .or. idecomp == 4) then
                  call decpair(-1,i,j,f12-f10)
               end if
               de = de + (twelve*f12 - ten*f10)*r2inv
#endif
               
            end if  ! ( ic > 0 )
         end if  ! ( .not. skipv(j) )
         
         !    -- derivatives:
         
         de = de*frespa
         dedx = de * xij
         dedy = de * yij
         dedz = de * zij
         dumx = dumx + dedx
         dumy = dumy + dedy
         dumz = dumz + dedz
         f(3*j-2) = f(3*j-2) - dedx
         f(3*j-1) = f(3*j-1) - dedy
         f(3*j  ) = f(3*j  ) - dedz
      end do
      
      f(3*i-2) = f(3*i-2) + dumx
      f(3*i-1) = f(3*i-1) + dumy
      f(3*i  ) = f(3*i  ) + dumz
#ifdef MPI
      do k=1,numtasks
         iexcl = iexcl + numex(i+k-1)
      end do
#else
      iexcl = iexcl + numex(i)
#endif
   end do  !  i=1,maxi
   call timer_stop(TIME_GBFRC)
   
   !--------------------------------------------------------------------------
   
   !    Step 3:  Finally, do the reduction over the sumdeijda terms:, adding
   !             into the forces those terms that involve derivatives of
   !             the GB terms (including the diagonal or "self" terms) with
   !             respect to the effective radii.  This is done by computing
   !             the vector dai/dxj, and using the chain rule with the
   !             previously-computed sumdeijda vector.
   
   !             Also, compute a surface-area dependent term if igbsa=1
   
   !             Do these terms only at "nrespa" multiple-time step intervals;
   !             (when igb=2 or 5, one may need to do this at every step)
   
   !--------------------------------------------------------------------------
   
   if( onstep ) then
      count=0
      frespa = nrespa
#ifdef MPI
      
      !       -- first, collect all the sumdeijda terms:
      
      call timer_start(TIME_GBRADDIST)
      if( numtasks > 1 ) then
         call mpi_allreduce(sumdeijda,vectmp1,natom, &
               MPI_DOUBLE_PRECISION,mpi_sum,commsander,ierr)
         sumdeijda(1:natom) = vectmp1(1:natom)
      end if
      call timer_stop(TIME_GBRADDIST)
#endif
      call timer_start(TIME_GBRAD2)
      
      !       -- diagonal epol term, plus off-diag derivs wrt alpha == reff^-1:
      
#ifdef MPI
      do i=mytaskid+1,maxi,numtasks
#else
      do i=1,maxi
#endif
         
         f_xi = zero
         f_yi = zero
         f_zi = zero
         qi = charge(i)
         expmkf = exp( -kappa * reff(i) )*extdieli
         dl = intdieli - expmkf
         qi2h = half*qi*qi
         qid2h = qi2h*dl
         epol = epol - qid2h/reff(i)
         if(idecomp == 1 .or. idecomp == 2) then
            call decpair(1,i,i,-qid2h/reff(i))
         else if(idecomp == 3 .or. idecomp == 4) then
            call decpair(-1,i,i,-qid2h/reff(i))
         end if
#ifdef CHARGE_MIXING
         if( onstep .and. icfe /= 0 ) then
            dvdl = dvdl - dl*charge(i)*dcharge(i)/reff(i)
         end if
#endif
         if(oncpstep) then
            dvdl = dvdl - (half*dl*dcharge(i)* &
                  dcharge(i)-qid2h)/reff(i)
         end if

         temp7 = -sumdeijda(i) + qid2h - kappa*qi2h*expmkf*reff(i)
         
         xi = x(3*i-2)
         yi = x(3*i-1)
         zi = x(3*i)
         ri = rborn(i)-offset
         ri1i = one/ri
         iaci = ntypes * (iac(i) - 1)
         
         if( igb == 2 .or. igb == 5 ) then
            
            !         --- new onufriev: we have to later scale values by a
            !             alpha,beta,gamma -dependent factor:
            
            ri = rborn(i) - offset
            thi = tanh( (gbalpha + gbgamma*psi(i)*psi(i) &
                  - gbbeta*psi(i) )*psi(i) )
            thi2 = (gbalpha + three*gbgamma*psi(i)*psi(i) &
                  - two*gbbeta*psi(i) ) &
                  *(one - thi*thi)*ri/rborn(i)
         end if
         
         icount = 0
         do j=1,natom
            if( i == j ) cycle
            
            xij = xi - x(3*j-2)
            yij = yi - x(3*j-1)
            zij = zi - x(3*j  )
            r2 = xij*xij + yij*yij + zij*zij
            if( r2 > rgbmaxpsmax2 ) cycle
            ! pairlist contains only atoms within rgbmax + safety margin
            
            icount = icount + 1
            jj(icount) = j
            r2x(icount) = r2
            
         end do
         call vdinvsqrt( icount, r2x, vectmp1 )
         
         kk1 = 0
         do k=1,icount
            j = jj(k)
            r2 = r2x(k)
            sj =  fs(j)

            dij1i = vectmp1(k)
            dij = r2*dij1i
            sj2 = sj * sj
            
            if( dij > 4.d0*sj ) goto 260
            kk1 = kk1 + 1
            vectmp3(kk1) = dij + sj
            if( dij > ri+sj ) then
               vectmp2(kk1) = r2 - sj2
               vectmp4(kk1) = dij - sj
            else if ( dij > abs(ri-sj) ) then
               vectmp2(kk1) = dij + sj
               vectmp4(kk1) = ri
            else if ( ri < sj ) then
               vectmp2(kk1) = r2 - sj2
               vectmp4(kk1) = sj - dij
            else
                vectmp2(kk1) = one
               vectmp4(kk1) = one
            end if
            260 continue
         end do
         
         call vdinv( kk1, vectmp2, vectmp2 )
         call vdinv( kk1, vectmp3, vectmp3 )
         vectmp4(1:kk1) = vectmp4(1:kk1)*vectmp3(1:kk1)
         call vdln( kk1, vectmp4, vectmp4 )
         
         kk1 = 0
         do k=1,icount
            j = jj(k)
            j3 = 3*j
            r2 = r2x(k)
            xij = xi - x(j3-2)
            yij = yi - x(j3-1)
            zij = zi - x(j3  )
            
            dij1i = vectmp1(k)
            dij = r2*dij1i
            sj =  fs(j)
            if (dij > rgbmax +sj) cycle
            sj2 = sj * sj
            
            !           datmp will hold (1/r)(dai/dr):
            
            dij2i = dij1i*dij1i
            dij3i = dij2i*dij1i

            if (dij > rgbmax - sj ) then 

               temp1 = 1./(dij-sj)
               datmp = 0.125d0 * dij3i * ((r2 + sj2) * &
                       (temp1*temp1 - rgbmax2i) - 2.d0 * log(rgbmax*temp1))
            
            else if( dij > 4.d0*sj ) then
               
               tmpsd = sj2*dij2i
               dumbo = te+tmpsd* (tf+tmpsd* (tg+tmpsd* (th+tmpsd* thh)))
               datmp = tmpsd*sj*dij2i*dij2i*dumbo
               
               !     ---check accuracy of above Taylor series:
               !      kk1 = kk1 + 1
               !      datmp2 = vectmp2(kk1)*sj*(-Half*dij2i + vectmp2(kk1)) + &
               !              Fourth*dij3i*vectmp4(kk1)
               !      if( abs( datmp/datmp2 - 1.d0 ) .gt. 0.00001) &
               !             write(6,*) i,j, datmp, datmp2
               
            else if( dij > ri+sj ) then
               
               kk1 = kk1 + 1
               datmp = vectmp2(kk1)*sj*(-half*dij2i + vectmp2(kk1)) + &
                     fourth*dij3i*vectmp4(kk1)
               
            else if ( dij > abs(ri-sj) ) then
               kk1 = kk1 + 1
               datmp = -fourth*(-half*(r2 - ri*ri + sj2)*dij3i*ri1i*ri1i &
                     + dij1i*vectmp2(kk1)*(vectmp2(kk1) - dij1i) &
                     - dij3i*vectmp4(kk1) )
               
            else if ( ri < sj ) then
               kk1 = kk1 + 1
               datmp = -half*(sj*dij2i*vectmp2(kk1) &
                     - two*sj*vectmp2(kk1)*vectmp2(kk1) &
                     - half*dij3i*vectmp4(kk1) )
               
            else
               kk1 = kk1 + 1
               datmp = zero
            end if  ! ( dij > 4.d0*sj )
            
            datmp = -datmp*frespa*temp7
            if( igb == 2 .or. igb == 5 ) datmp = datmp*thi2
            
            f_x = xij*datmp
            f_y = yij*datmp
            f_z = zij*datmp
            f(j3-2) = f(j3-2) + f_x
            f(j3-1) = f(j3-1) + f_y
            f(j3  ) = f(j3  ) + f_z
            f_xi = f_xi - f_x
            f_yi = f_yi - f_y
            f_zi = f_zi - f_z
            
         end do  !  k=1,icount
         
         f(3*i-2) = f(3*i-2) + f_xi
         f(3*i-1) = f(3*i-1) + f_yi
         f(3*i  ) = f(3*i  ) + f_zi
         
         !       --- Define neighbor list ineighbor for calc of LCPO areas ---
         
         if ( gbsa > 0 ) then
            do k=1,icount
               j = jj(k)
               dij = sqrt(r2x(k))
               if ((vdwrad(i) + vdwrad(j)) > dij ) then
                  if ( (gbsa == 2) .or. &  !  Consider all atoms for icosa, only non-H's for LCPO
                       (vdwrad(i) > 2.5).and.(vdwrad(j) > 2.5) ) then
                     count=count+1
                     ineighbor(count) = j
                  end if
               end if
            end do
         end if
         
         if ( gbsa > 0 ) then
            count=count+1
            ineighbor(count)=0
         end if
         
         if ( gbsa == 2 ) then

            !       --- Calc surface area with icosasurf algo;
            !           does not provide forces ==> only for single point calculations ---

            if(i == 1) then
               ineighborpt = 1
               call icosa_init(2, 3, 0.0d0)
            end if
            totsasa = totsasa + icosa_sphere_approx(i,x, &
                                                    vdwrad,ineighborpt,ineighbor, &
                                                    idecomp,xx,ix)
            
         end if  !  ( gbsa == 2 )

      end do   ! end loop over atom i
      
      call timer_stop(TIME_GBRAD2)

      if ( gbsa == 1 ) then
         
         !       --- calculate surface area by computing LCPO (Still) over
         !           all atoms ---
         
#        include "gbsa.h"
         
      end if  !  ( gbsa == 1 )
      esurf = surften*totsasa
      
   end if  !  i=mytaskid+1,maxi,numtasks
   return
end subroutine egb 

end module genborn
