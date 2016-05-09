! <compile=optimized>
#include "copyright.h"
#include "is_copyright.h"
#include "dprec.h"
#include "is_def.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ partition of atoms into internal and external portions
subroutine pb_atmpart( verbose,pbprint,natom,ibgwat,ienwat,inatm,outatm,ipres,outflag,xctr,yctr,zctr,rdiel,sepbuf,x )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Authors:
   ! Lijiang Yang, Luo Research Group, UC-Irvine
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

   ! Passed variables

   logical verbose, pbprint
   integer natom, ibgwat, ienwat, inatm, outatm, ipres(*)
   integer outflag(natom)
   _REAL_ xctr, yctr, zctr, rdiel, sepbuf
   _REAL_ x(3,natom)

   ! Local variables

   integer ires, num
   _REAL_ sepr, atmr, xtmp, ytmp, ztmp
 
   sepr = (rdiel-sepbuf)**2
   outflag = 0
 
   ! internal portion always include protein atoms
   ! external portion only include water atoms
 
   inatm = ipres(ibgwat-1)
   outatm = 0

   do ires = ibgwat, ienwat
      num = ipres(ires)
      xtmp = x(1,num); ytmp = x(2,num); ztmp = x(3,num)
      atmr = (xtmp-xctr)**2 + (ytmp-yctr)**2 + (ztmp-zctr)**2
      if ( atmr > sepr ) then
         outatm = outatm + 3
         outflag(num) = 1
         outflag(num+1) = 1
         outflag(num+2) = 1
      else
         inatm = inatm + 3
      end if
   end do

   if ( verbose .and. pbprint ) then
      write(6,'(a,2i6,a,f6.3)') ' Atoms are partitioned into two regions', inatm, outatm, ' with a buffer of', sepbuf
   end if

end subroutine pb_atmpart
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Residue-based nblist for PBMD
subroutine pb_reslist( verbose,pbprint,maxnbr,natom,nres,ibgwat,ienwat,ntypes,ipres,iac,ico,natex,nshrt,iar1pb,iprlong,cutres,x)
    
   implicit none
    
   ! Passed variables
   
   logical verbose, pbprint
   integer maxnbr, natom, nres, ibgwat, ienwat, ntypes, ipres(*), iac(*), ico(*), natex(*), nshrt(0:natom)
   integer iar1pb(6,0:natom),iprlong(*)
   _REAL_ cutres
   _REAL_ x(3,*)
     
   ! Local variables
     
   integer i, j, ii, jj, ijp, itotal, jtotal, ijtotal, lastii, lastjj1, lastjj2, ip1, ip2, jp1, jp2
   integer lrp, iprpt, lpair
   integer iwe(natom), iwh(natom), iwa(natom), irp(MAXNEI)
   _REAL_ xi, yi, zi, dx, dy, dz, d2, cutlng
    
   iprpt = 0
   lpair = 0
   
   ! long cutoff to search close res/res pairs

   cutlng = (sqrt(cutres) + 10.d0)**2
        
   do i = 1, natom
      iar1pb(5,i) = 0
      iar1pb(6,i) = 0
      iwe(i) = 0
   end do
        
   ! loop over nonsolvent residues ii
   if (ibgwat /= 0) then      
      lastii = ibgwat - 1   ! no. res of solute/ions
      lastjj1 = ibgwat - 1   ! no. res of solute/ions
      lastjj2 = ienwat
   else
      lastii = nres
      lastjj1 = nres
      lastjj2 = -999
   endif
   do ii = 1, lastii
      ip1 = ipres(ii)
      ip2 = ipres(ii+1)-1
      itotal = ip2 - ip1 + 1
      lrp = 0
           
      ! always include ii == jj
        
      lrp = lrp+1
      irp(lrp) = ii
        
      ! loop over nonsolvent residues jj
        
      do jj = ii+1, lastjj1
           
         jp1 = ipres(jj)
         jp2 = ipres(jj+1)-1
         jtotal = jp2 - jp1 + 1
           
         ! see if the approximated backbone centers of the 2 residues fall
         ! within the long cutoff
           
         dx = (x(1,ip1)+x(1,ip2))/TWO - (x(1,jp1)+x(1,jp2))/TWO
         dy = (x(2,ip1)+x(2,ip2))/TWO - (x(2,jp1)+x(2,jp2))/TWO
         dz = (x(3,ip1)+x(3,ip2))/TWO - (x(3,jp1)+x(3,jp2))/TWO
         d2 = dx*dx + dy*dy + dz*dz
         if ( d2 .gt. cutlng ) cycle
          
         ! see if the 2 residues fall within the final residue cutoff
          
         ijtotal = itotal*jtotal - 1
         do ijp = 0, ijtotal
            i = ijp/jtotal + ipres(ii)
            j = mod(ijp, jtotal) + ipres(jj)
            dx = x(1,i) - x(1,j)
            dy = x(2,i) - x(2,j)
            dz = x(3,i) - x(3,j)
            d2 = dx*dx + dy*dy + dz*dz
            if ( d2 > cutres ) cycle
            lrp = lrp + 1
            irp(lrp) = jj
            exit
         enddo
           
      end do  !  jj = ii+1, lastjj1
        
      ! loop over solvent residues jj
        
      do jj = ibgwat, lastjj2
         xi = x(1,ipres(jj))
         yi = x(2,ipres(jj))
         zi = x(3,ipres(jj))
          
         ! see if the approximated backbone centers of the 2 residues fall
         ! within the long cutoff
           
         dx = xi - (x(1,ip1)+x(1,ip2))/TWO
         dy = yi - (x(2,ip1)+x(2,ip2))/TWO
         dz = zi - (x(3,ip1)+x(3,ip2))/TWO
         d2 = dx*dx + dy*dy + dz*dz
         if ( d2 > cutlng ) cycle
          
         ! see if the 2 residues fall within the final residue cutoff
          
         do i = ip1, ip2
            dx = xi - x(1,i)
            dy = yi - x(2,i)
            dz = zi - x(3,i)
            d2 = dx*dx + dy*dy + dz*dz
            if ( d2 > cutres ) cycle
            lrp = lrp + 1
            irp(lrp) = jj
            exit
         end do
       
      end do  ! jj = ibgwat, lastjj2
      
      if (lrp > MAXNEI) then
         write(6, *) 'PB bomb in pb_reslist(): MAXNEI too short' 
         call mexit(6, 1)
      end if
      
      ! pack into nblist, taking care of exclusions
       
      call packlist(ii,ip1,ip2,lrp,lpair,iprpt,ntypes,irp,ipres,iac,ico,natex,nshrt,iwa,iwh,iwe,iprlong)
       
   end do  ! do ii = 1, lastii
    
   ! loop over solvent residues ii, use first atom only
   
   if (ibgwat /= 0) then      
      lastii = nres - 1
   else
      lastii = -999
   endif
   lastjj1 = nres
   do ii = ibgwat, lastii
      ip1 = ipres(ii)
      ip2 = ipres(ii+1)-1
      i = ipres(ii)
      xi = x(1,i)
      yi = x(2,i)
      zi = x(3,i)
      lrp = 0 
                
      ! loop over solvent residues jj, use first atom only
                
      do jj = ii+1, lastjj1
         j = ipres(jj)
         dx = xi - x(1,j)
         dy = yi - x(2,j)
         dz = zi - x(3,j)
         d2 = dx*dx + dy*dy + dz*dz
         if ( d2 <= cutres ) then
            lrp = lrp + 1
            irp(lrp) = jj
         endif
      end do
       
      if (lrp > MAXNEI) then
         write(6, *) 'PB bomb in pb_reslist(): MAXNEI too short' 
         call mexit(6, 1)
      end if

      ! pack into nblist, taking care of exclusions
       
      call packlist(ii,ip1,ip2,lrp,lpair,iprpt,ntypes,irp,ipres,iac,ico,natex,nshrt,iwa,iwh,iwe,iprlong)
       
   end do
   if ( verbose .and. pbprint ) write(6,'(2x,a,i9)') 'NB-update: residue-based nb list', lpair

contains       
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Pack nblist by taking care of excluded atoms
subroutine packlist( ii,ip1,ip2,lrp,lpair,iprpt,ntypes,irp,ipres,iac,ico,natex,nshrt,iwa,iwh,iwe,iprlong )
      
   implicit none
    
   ! Passed variables
   
   integer ii, ip1, ip2, lrp, lpair, iprpt, ntypes
   integer irp(*),ipres(*),iac(*),ico(*),natex(*),nshrt(0:natom),iwa(*),iwh(*),iwe(*),iprlong(*)
      
   ! Local variables 
   integer jj, i, j, k, kp, jp1, jp2, jrp, lpr, lhb, iaci
    
   do i = ip1, ip2
      iaci = ntypes*(iac(i)-1)
      lpr = 0
      lhb = 0
      
      do kp = nshrt(i-1)+1, nshrt(i)
         k = max0(1,natex(kp))
         iwe(k) = i
      end do
        
      do jrp = 1, lrp
         jj = irp(jrp)
         if (ii == jj) then
            jp1 = i+1
         else
            jp1 = ipres(jj)
         end if
         jp2 = ipres(jj+1)-1
           
         do j = jp1, jp2
            if (iwe(j) == i) cycle
            
            if( ico(iaci+iac(j)) > 0) then
               lpr = lpr+1 ! non-bonded pair
               iwa(lpr) = j
            else
               lhb = lhb+1 ! h-bonded pair
               iwh(lhb) = j
            end if
         end do  ! j = jp1, jp2
      end do  ! jrp = 1, lrp
          
      lpair = lpair + lpr+lhb
      iar1pb(5,i) = lpr
      iar1pb(6,i) = lhb
      if ( lpair <= maxnbr ) then
         do j = 1, lpr
            iprpt = iprpt + 1
            iprlong(iprpt) = iwa(j)
         end do
         do j = 1, lhb
            iprpt = iprpt + 1
            iprlong(iprpt) = iwh(j)
         end do
      else
         write (6, *) 'PB bomb in pb_reslist(): maxnbr too small'
         call mexit(6, 1)
      end if
   end do  ! i = ip1, ip2

end subroutine packlist

end subroutine pb_reslist
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Atom-based nblist for PBMD
subroutine pb_atmlist( verbose,pbprint,maxnba,natom,ntypes,iac,ico,natex,nshrt,nex,iex,iar1pb,iprlong,iprshrt,&
                       cutnb,cutsa,cutfd,cn1,cn2,cn1pb,cn2pb,cn3pb,cg,x )
    
   implicit none
    
   ! Passed variables
    
   integer maxnba,natom,ntypes,natex(*),nshrt(0:natom),nex(*),iex(32,*),iac(*),ico(*),iar1pb(6,0:natom)
   integer iprlong(*),iprshrt(*)
   _REAL_ x(3,*), cn1(*), cn2(*), cg(*)
   _REAL_ cn1pb(*), cn2pb(*), cn3pb(*)
   _REAL_ cutnb, cutsa, cutfd
    
   ! Local variables
    
   logical verbose, pbprint 
   integer iaci, ic
   integer ilast, iatm, jatm
   integer i, j, jp, lpr, npr, lpair
   integer cntr, eclose, pclose, sclose, nclose
   integer tmpex(MAXNEI), tmppb(MAXNEI), tmpsa(MAXNEI), tmpnb(MAXNEI)
   _REAL_ xi, yi, zi, dx, dy, dz, d2
   _REAL_ cgi
   
   ! set up zeroth atom for limits
    
   iar1pb(1, 0) = 0
   iar1pb(2, 0) = 0
   iar1pb(3, 0) = 0
   iar1pb(4, 0) = 0
   iar1pb(1, natom) = 0
   iar1pb(2, natom) = 0
   iar1pb(3, natom) = 0
   iar1pb(4, natom) = 0
   nex(1:natom) = 0
    
   ! auxiliary pointer
      
   lpair = 0
    
   ! this is the global index of atom-based pair
    
   cntr = 0
   ilast = natom - 1
   do i = 1, ilast
      iaci = ntypes*(iac(i)-1)
      cgi  = -cg(i)
        
      ! save excluded pair into tmpex
        
      eclose = 0
      do jp = nshrt(i-1) + 1, nshrt(i)
         j = natex(jp)
         if (j == 0)  cycle
         eclose = eclose + 1
         tmpex(eclose) = j
      end do
       
      ! save nonboneded pairs into tmppb and tmpnb
       
      xi = x(1, i)
      yi = x(2, i)
      zi = x(3, i)
      lpr = iar1pb(5, i)
      npr = lpr + iar1pb(6, i)
      pclose = 0
      sclose = 0
      nclose = 0
      do jp = 1, lpr
         j = iprlong(jp + lpair)
         dx = xi - x(1, j)
         dy = yi - x(2, j)
         dz = zi - x(3, j)
         d2 = dx**2 + dy**2 + dz**2
         if (d2 <= cutfd) then
            pclose = pclose + 1
            tmppb(pclose) = j
         else if (d2 <= cutsa) then
            sclose = sclose + 1
            tmpsa(sclose) = j
         else if (d2 <= cutnb) then
            nclose = nclose + 1
            tmpnb(nclose) = j
         end if
      end do
      
      ! save h-bonding pairs into tmppb, these are not needed in nb calculation
      
      do jp = lpr+1, npr
         j = iprlong(jp + lpair)
         dx = xi - x(1, j)
         dy = yi - x(2, j)
         dz = zi - x(3, j)
         d2 = dx**2 + dy**2 + dz**2
         if (d2 <= cutfd) then
            pclose = pclose + 1
            tmppb(pclose) = j
         end if
      end do
      
      ! now pack them into the new atom-based nblist
      
      if (eclose > MAXNEI .or. pclose > MAXNEI .or. sclose > MAXNEI .or. nclose > MAXNEI) then
         write(6, *) 'PB bomb in pb_atmlist(): MAXNEI too short'
         call mexit(6, 1)
      end if
      if (eclose + pclose + sclose + nclose + cntr > maxnba) then
         write(6, '(4x,a,i8)') 'PB Bomb in pb_atmlist(): maxnba too short'
         call mexit(6, 1)
      end if
      do j = 1, eclose
         cntr = cntr + 1
         iprshrt(cntr) = tmpex(j)
         iatm = i
         jatm = tmpex(j)
         nex(iatm) = nex(iatm) + 1
         iex(nex(iatm),iatm) = jatm
         nex(jatm) = nex(jatm) + 1
         iex(nex(jatm),jatm) = iatm
      end do
      iar1pb(1, i) = cntr ! for exclusion and FDCOUL calculation
      do j = 1, pclose
         cntr = cntr + 1
         iprshrt(cntr) = tmppb(j)
         ic = ico(iaci+iac(tmppb(j)))
         if ( ic > 0 ) then ! normal nonbonded pairs
            cn1pb(cntr) = cn1(ic)
            cn2pb(cntr) = cn2(ic)
         else               ! h-bonded pairs
            cn1pb(cntr) = ZERO
            cn2pb(cntr) = ZERO
         end if
         cn3pb(cntr) = cgi*cg(tmppb(j))
      end do
      iar1pb(2, i) = cntr ! for FDCOUL calculation
      do j = 1, sclose
         cntr = cntr + 1
         iprshrt(cntr) = tmpsa(j)
         ic = ico(iaci+iac(tmpsa(j)))
         cn1pb(cntr) = cn1(ic)
         cn2pb(cntr) = cn2(ic)
         cn3pb(cntr) = cgi*cg(tmpsa(j))
      end do
      iar1pb(3, i) = cntr ! for SA calculation
      do j = 1, nclose
         cntr = cntr + 1
         iprshrt(cntr) = tmpnb(j)
         ic = ico(iaci+iac(tmpnb(j)))
         cn1pb(cntr) = cn1(ic)
         cn2pb(cntr) = cn2(ic)
         cn3pb(cntr) = cgi*cg(tmpnb(j))
      end do
      iar1pb(4, i) = cntr ! for VDW calculation

      lpair = lpair + npr
   end do  !  i = 1, ilast
   if ( verbose .and. pbprint ) write(6,'(2x,a,i9)') 'NB-update: atom-based nb list', cntr
   
end subroutine pb_atmlist 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Set up FDPB grid for a de novo call
subroutine pb_setgrd( verbose, prnt, initial, ifcap, natom, xcap, ycap, zcap, cutcap, x )
    
   use poisson_boltzmann
   use solvent_accessibility, only : radi
   implicit none
    
   ! Passed variables
    
   logical verbose, prnt, initial
   integer ifcap, natom, totsavxmymzm
   _REAL_ xcap, ycap, zcap, cutcap
   _REAL_ x(3, *)
    
   ! Local variables
    
   integer l, alloc_err(32)
    
   ! get center and dimension information of the current molecule
    
   call setgrd( verbose, prnt, initial, ifcap, natom, xcap, ycap, zcap, cutcap, x )
    
   ! if allocate working arrays for fdpb
    
   alloc_err(1:32) = 0
   if ( .not. initial ) then
      deallocate(   sphi, stat = alloc_err(1 ) )
      deallocate(    sad, stat = alloc_err(2 ) )
      deallocate(    sbv, stat = alloc_err(3 ) )
      deallocate(    spv, stat = alloc_err(4 ) )
      deallocate(    phi, stat = alloc_err(5 ) )
      deallocate(   epsx, stat = alloc_err(6 ) )
      deallocate(   epsy, stat = alloc_err(7 ) )
      deallocate(   epsz, stat = alloc_err(8 ) )
      deallocate( iepsav, stat = alloc_err(9 ) )
      deallocate(  insas, stat = alloc_err(10) )
      deallocate( atmsas, stat = alloc_err(11) )
      !deallocate( fedgex, stat = alloc_err(12) )
      !deallocate( fedgey, stat = alloc_err(13) )
      !deallocate( fedgez, stat = alloc_err(14) )
      !deallocate( fatomx, stat = alloc_err(15) )
      !deallocate( fatomy, stat = alloc_err(16) )
      !deallocate( fatomz, stat = alloc_err(17) )
      deallocate(     ad, stat = alloc_err(18) )
      deallocate(     rd, stat = alloc_err(19) )
      deallocate(    am1, stat = alloc_err(20) )
      deallocate(    am2, stat = alloc_err(21) )
      deallocate(    am3, stat = alloc_err(22) )
! WJ
      deallocate(    am4, stat = alloc_err(28) )
      deallocate(    am5, stat = alloc_err(29) )
      deallocate(    am6, stat = alloc_err(30) )
!
      deallocate(     bv, stat = alloc_err(23) )
      deallocate(     zv, stat = alloc_err(24) )
      deallocate(     pv, stat = alloc_err(25) )
      deallocate(     tv, stat = alloc_err(26) )
      deallocate(     xs, stat = alloc_err(27) )
      if ( alloc_err( 1)+alloc_err( 2)+alloc_err( 3)+alloc_err( 4)+alloc_err( 5)+&
           alloc_err( 6)+alloc_err( 7)+alloc_err( 8)+alloc_err( 9)+alloc_err(10)+&
           alloc_err(11)+alloc_err(12)+alloc_err(13)+alloc_err(14)+alloc_err(15)+&
           alloc_err(16)+alloc_err(17)+alloc_err(18)+alloc_err(19)+alloc_err(20)+&
           alloc_err(21)+alloc_err(22)+alloc_err(23)+alloc_err(24)+alloc_err(25)+&
           alloc_err(26)+alloc_err(27)+alloc_err(28)+alloc_err(29)+alloc_err(30) /= 0 ) then
         write(6, *) 'PB bomb in pb_setgrd(): Deallocation aborted', alloc_err(1:30)
         call mexit(6, 1)
      end if
   end if
   allocate(   sphi(  1:savxmymzm(1     )), stat = alloc_err(1 ) )
   allocate(    sad(  1:savxmymzm(1     )), stat = alloc_err(2 ) )
   allocate(    sbv(  1:savxmymzm(nfocus)), stat = alloc_err(3 ) )
   allocate(    spv(  1:savxmymzm(nfocus)), stat = alloc_err(4 ) )
   allocate(    phi(  1:savxmymzm(nfocus)), stat = alloc_err(5 ) )
   allocate(   epsx(  1:savxmymzm(nfocus)), stat = alloc_err(6 ) )
   allocate(   epsy(  1:savxmymzm(nfocus)), stat = alloc_err(7 ) )
   allocate(   epsz(  1:savxmymzm(nfocus)), stat = alloc_err(8 ) )
   allocate( iepsav(4,1:savxmymzm(nfocus)), stat = alloc_err(9 ) )
   allocate(  insas(  1:savxmymzm(nfocus)), stat = alloc_err(10) )
   allocate( atmsas(  1:savxmymzm(nfocus)), stat = alloc_err(11) )
   !allocate( fedgex(2,1:savxmymzm(nfocus)), stat = alloc_err(12) )
   !allocate( fedgey(2,1:savxmymzm(nfocus)), stat = alloc_err(13) )
   !allocate( fedgez(2,1:savxmymzm(nfocus)), stat = alloc_err(14) )
   !allocate( fatomx(2,1:savxmymzm(nfocus)), stat = alloc_err(15) )
   !allocate( fatomy(2,1:savxmymzm(nfocus)), stat = alloc_err(16) )
   !allocate( fatomz(2,1:savxmymzm(nfocus)), stat = alloc_err(17) )
   allocate( ad (-savxmym(nfocus):savxmymzm(nfocus)+savxmym(nfocus)         ), stat = alloc_err(18) )
   allocate( rd (-savxmym(nfocus):savxmymzm(nfocus)+savxmym(nfocus)         ), stat = alloc_err(19) )
   allocate( am1(-savxmym(nfocus):savxmymzm(nfocus)+savxmym(nfocus)         ), stat = alloc_err(20) )
   allocate( am2(-savxmym(nfocus):savxmymzm(nfocus)+savxmym(nfocus)         ), stat = alloc_err(21) )
   allocate( am3(-savxmym(nfocus):savxmymzm(nfocus)+savxmym(nfocus)         ), stat = alloc_err(22) )
! WJ
   allocate( am4(1:savxmymzm(nfocus)         ), stat = alloc_err(28) )
   allocate( am5(1:savxmymzm(nfocus)         ), stat = alloc_err(29) )
   allocate( am6(1:savxmymzm(nfocus)         ), stat = alloc_err(30) )
!
   allocate( bv (-savxmym(nfocus):savxmymzm(nfocus)+savxmym(nfocus)         ), stat = alloc_err(23) )
   allocate( zv (-savxmym(nfocus):savxmymzm(nfocus)+savxmym(nfocus)         ), stat = alloc_err(24) )
   allocate( pv (-savxmym(nfocus):savxmymzm(nfocus)+savxmym(nfocus)         ), stat = alloc_err(25) )
   allocate( tv (-savxmym(nfocus):savxmymzm(nfocus)+savxmym(nfocus)         ), stat = alloc_err(26) )
   totsavxmymzm = 0
   do l = 1, nfocus
! WJ
!     totsavxmymzm = totsavxmymzm + savxmymzm(l)
      totsavxmymzm = totsavxmymzm + 2*savxmymzm(l)
!
   end do
   allocate( xs (               1:totsavxmymzm                              ), stat = alloc_err(27) )
   if ( alloc_err( 1)+alloc_err( 2)+alloc_err( 3)+alloc_err( 4)+alloc_err( 5)+&
        alloc_err( 6)+alloc_err( 7)+alloc_err( 8)+alloc_err( 9)+alloc_err(10)+&
        alloc_err(11)+alloc_err(12)+alloc_err(13)+alloc_err(14)+alloc_err(15)+&
        alloc_err(16)+alloc_err(17)+alloc_err(18)+alloc_err(19)+alloc_err(20)+&
        alloc_err(21)+alloc_err(22)+alloc_err(23)+alloc_err(24)+alloc_err(25)+&
        alloc_err(26)+alloc_err(27)+alloc_err(28)+alloc_err(29)+alloc_err(30) /= 0 ) then
      write(6, *) 'PB bomb in pb_setgrd(): Allocation aborted', alloc_err(1:30)
      call mexit(6, 1)
   end if

   ! initialize saved phi map
    
   xs = ZERO
    
   ! save fine grid limits for checking out-of-grid atoms
   
   gxmin = savgox(nfocus) + 3*savh(nfocus);
   gxmax = savgox(nfocus) + (savxm(nfocus)-2)*savh(nfocus)
   gymin = savgoy(nfocus) + 3*savh(nfocus);
   gymax = savgoy(nfocus) + (savym(nfocus)-2)*savh(nfocus)
   gzmin = savgoz(nfocus) + 3*savh(nfocus);
   gzmax = savgoz(nfocus) + (savzm(nfocus)-2)*savh(nfocus)

contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ set up FDPB grid for a de novo call
subroutine setgrd( verbose, prnt, initial, ifcap, natom, xcap, ycap, zcap, cutcap, x )
    
   implicit none
    
   ! Passed variables
    
   logical verbose, prnt, initial
   integer ifcap, natom 
   _REAL_ xcap, ycap, zcap, cutcap
   _REAL_ x(3,*)
    
   ! Local variables
    
   integer l, iatm
   _REAL_ xlength, ylength, zlength
   _REAL_ xbox, ybox, zbox, htmp
    
   ! set bounding box center for all atoms
    
   if ( verbose .and. prnt ) then
      write(6, *)
      write(6, *)
      write(6, *) '======== Setting up Grid Parameters ========'
      write(6, *) 'Using bounding box for grid setup'
   end if
   if ( ifcap == 0 ) then
      xmin = 9999.0; ymin = 9999.0; zmin = 9999.0
      xmax = -9999.0; ymax = -9999.0; zmax = -9999.0
      do iatm = 1, natom
         if ( x(1,iatm)-radi(iatm) .lt. xmin ) xmin = x(1,iatm)-radi(iatm)
         if ( x(1,iatm)+radi(iatm) .gt. xmax ) xmax = x(1,iatm)+radi(iatm)
         if ( x(2,iatm)-radi(iatm) .lt. ymin ) ymin = x(2,iatm)-radi(iatm)
         if ( x(2,iatm)+radi(iatm) .gt. ymax ) ymax = x(2,iatm)+radi(iatm)
         if ( x(3,iatm)-radi(iatm) .lt. zmin ) zmin = x(3,iatm)-radi(iatm)
         if ( x(3,iatm)+radi(iatm) .gt. zmax ) zmax = x(3,iatm)+radi(iatm)
      enddo
      xbox = (xmax + xmin)/TWO; ybox = (ymax + ymin)/TWO; zbox = (zmax + zmin)/TWO

      ! rounding to nearest h unit for easy restarting

      htmp = savh(nfocus)
      xbox = nint(xbox/htmp)*htmp; ybox = nint(ybox/htmp)*htmp; zbox = nint(zbox/htmp)*htmp
   else
      xmin = xcap - cutcap; xmax = xcap + cutcap
      ymin = ycap - cutcap; ymax = ycap + cutcap
      zmin = zcap - cutcap; zmax = zcap + cutcap
      xbox = xcap; ybox = ycap; zbox = zcap
   end if

   if ( verbose .and. prnt ) then
      write(6, '(1x,a,3f10.3)') 'Bounding Box Center:  ', xbox, ybox, zbox
      write(6, '(1x,a,3f10.3)') 'Xmin, Xmax, Xmax-Xmin:', xmin, xmax, xmax-xmin
      write(6, '(1x,a,3f10.3)') 'Ymin, Ymax, Ymax-Ymin:', ymin, ymax, ymax-ymin
      write(6, '(1x,a,3f10.3)') 'Zmin, Zmax, Zmax-Zmin:', zmin, zmax, zmax-zmin
   end if
   if ( initial ) then
      do l = 1, nfocus
         cxbox(l) = xbox; cybox(l) = ybox; czbox(l) = zbox
         if ( verbose .and. prnt ) write(6, '(a,1x,i5,1x,3f10.3)') &
            '   beginning box center at level ', l, cxbox(l), cybox(l), czbox(l)
      end do
   else
      do l = 1, nfocus
         savxbox(l) = cxbox(l); savybox(l) = cybox(l); savzbox(l) = czbox(l)
         if ( verbose .and. prnt ) write(6, '(a,1x,i5,1x,3f10.3)') &
            '   previous box center at level ', l, savxbox(l), savybox(l), savzbox(l)
      end do
   end if
    
   ! set first grid
    
   ! at least 50% fill for boundary condition
    
   xlength = (xmax-xmin)*fillx; ylength = (ymax-ymin)*filly; zlength = (zmax-zmin)*fillz
! debugging for xiang
!  xlength = 20d0; ylength = 20d0;  zlength = 20d0
!
   if ( outphi .and. phiform == 0 ) then
      xlength = max(xlength, ylength, zlength)
      ylength = xlength; zlength = xlength
   endif
   ! debugging for xiang
!   savxm(1) = int(xlength/savh(1))!+1
!   savym(1) = int(ylength/savh(1))!+1
!   savzm(1) = int(zlength/savh(1))!+1
   savxm(1) = nint(xlength/savh(1))
   savym(1) = nint(ylength/savh(1))
   savzm(1) = nint(zlength/savh(1))
   savxm(1) = 2*nint( REAL(savxm(1))*HALF ) + 1
   savym(1) = 2*nint( REAL(savym(1))*HALF ) + 1
   savzm(1) = 2*nint( REAL(savzm(1))*HALF ) + 1
!
   savxmym(1) = savxm(1)*savym(1)
   savxmymzm(1) = savxmym(1)*savzm(1)
   if ( verbose .and. prnt ) write(6, '(a,i5,1x,3i5)') &
      ' Grid dimension at level ', 1, savxm(1),savym(1),savzm(1)
    
   if ( initial ) then
      savgox(1) = - REAL(savxm(1)+1) * savh(1) * HALF + xbox
      savgoy(1) = - REAL(savym(1)+1) * savh(1) * HALF + ybox
      savgoz(1) = - REAL(savzm(1)+1) * savh(1) * HALF + zbox
   else
      cxbox(1) = savxbox(1) + nint( (xbox - savxbox(1))/savh(1) )*savh(1)
      cybox(1) = savybox(1) + nint( (ybox - savybox(1))/savh(1) )*savh(1)
      czbox(1) = savzbox(1) + nint( (zbox - savzbox(1))/savh(1) )*savh(1)
      if ( verbose .and. prnt ) write(6, '(a,i5,1x,3f10.3)') &
         ' Box center corrected at level ', 1, cxbox(1), cybox(1), czbox(1)
      savgox(1) = - REAL(savxm(1)+1) * savh(1) * HALF + cxbox(1)
      savgoy(1) = - REAL(savym(1)+1) * savh(1) * HALF + cybox(1)
      savgoz(1) = - REAL(savzm(1)+1) * savh(1) * HALF + czbox(1)
   end if
   if ( verbose .and. prnt ) write(6, '(a,i5,1x,3f10.3)') &
      ' Grid origin corrected at level ', 1, savgox(1),savgoy(1),savgoz(1)
    
   ! set remaining grids
   ! focus grids HALF grid scale away from the surface
   ! use bounding box for origin
    
   xlength = xmax-xmin; ylength = ymax-ymin; zlength = zmax-zmin
   if ( outphi .and. phiform == 0 ) then
      xlength = max(xlength, ylength, zlength)
      ylength = xlength; zlength = xlength
   endif
   do l = 2, nfocus
      savxm(l) = nint( xlength/savh(l) ) + nbuffer
      savym(l) = nint( ylength/savh(l) ) + nbuffer
      savzm(l) = nint( zlength/savh(l) ) + nbuffer
      savxm(l) = 2*nint( REAL(savxm(l))*HALF ) + 1
      savym(l) = 2*nint( REAL(savym(l))*HALF ) + 1
      savzm(l) = 2*nint( REAL(savzm(l))*HALF ) + 1
      savxmym(l) = savxm(l)*savym(l)
      savxmymzm(l) = savxmym(l)*savzm(l)
      if ( verbose .and. prnt ) write(6, '(a,i5,1x,3i5)') &
         ' Grid dimension at level ', l, savxm(l), savym(l), savzm(l)
      if ( initial ) then
         savgox(l) = - REAL(savxm(l)+1) * savh(l) * HALF + xbox
         savgoy(l) = - REAL(savym(l)+1) * savh(l) * HALF + ybox
         savgoz(l) = - REAL(savzm(l)+1) * savh(l) * HALF + zbox
      else
         cxbox(l) = savxbox(l) + nint( (xbox - savxbox(l))/savh(l) )*savh(l)
         cybox(l) = savybox(l) + nint( (ybox - savybox(l))/savh(l) )*savh(l)
         czbox(l) = savzbox(l) + nint( (zbox - savzbox(l))/savh(l) )*savh(l)
         if ( verbose .and. prnt ) write(6, '(a,i5,1x,3f10.3)') &
            ' Box center corrected at level ', l, cxbox(l), cybox(l), czbox(l)
         savgox(l) = - REAL(savxm(l)+1) * savh(l) * HALF + cxbox(l)
         savgoy(l) = - REAL(savym(l)+1) * savh(l) * HALF + cybox(l)
         savgoz(l) = - REAL(savzm(l)+1) * savh(l) * HALF + czbox(l)
      end if
      if ( verbose .and. prnt ) write(6, '(a,i5,1x,3f10.3)') &
         ' Grid origin corrected at level ', l, savgox(l), savgoy(l), savgoz(2)
   end do
    
   ! do some insanity checking for electrostatic focussing ...
    
   do l = 1, nfocus-1
      if ( savxmymzm(l) > savxmymzm(nfocus) ) then
         write(6, '(a,i5)') 'PB Bomb in setgrd(): coarse grid too large at level', l
         write(6, '(a,f6.3)') 'reset fscale to a larger number', savh(l)/savh(l+1)
         call mexit(6,1)
      end if
   end do
   do l = 2, nfocus
      if ( savgox(l)+savxm(l)*savh(l) >= savgox(l-1)+savxm(l-1)*savh(l-1) .or.&
           savgoy(l)+savym(l)*savh(l) >= savgoy(l-1)+savym(l-1)*savh(l-1) .or.&
           savgoz(l)+savzm(l)*savh(l) >= savgoz(l-1)+savzm(l-1)*savh(l-1) ) then
         write(6, '(a,i5)') 'PB Bomb in setgrd(): focusing grid too large', l
         write(6, '(a,f6.3)') 'reset fillratio to a larger number', fillratio
         call mexit(6,1)
      end if
   end do
    
   ! if requested offseting grid, do it here
    
   if ( offx + offy + offz /= ZERO ) then
      do l = 1, nfocus
         savgox(l) = savgox(l) + offx
         savgoy(l) = savgoy(l) + offy
         savgoz(l) = savgoz(l) + offz
         if ( verbose .and. prnt ) write(6, '(a,i5,1x,3f10.3)') &
         ' Grid origin at level after offset', l, savgox(l), savgoy(l), savgoz(2)
      end do
   end if
    
end subroutine setgrd

end subroutine pb_setgrd
