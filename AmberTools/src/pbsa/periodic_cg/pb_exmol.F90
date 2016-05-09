! <compile=optimized>
#include "copyright.h"
#include "is_copyright.h"
#include "dprec.h"
#include "is_def.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ molecular volume assignment
module pb_exmol
    
   implicit none
    
contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ driver for molecular volume assignment
!  solvent exluded surface approach
subroutine pb_exmol_ses( pbverbose,natom,sprob,stern,istrng,&
                         level,xm,ym,zm,xmymzm,h,gox,goy,goz,&
                         gcrd,acrd,radi,&
                         narcdot,maxarc,marc,m2narc,fstarc,arcatm,arccrd,savarc,&
                         atmsas,insas,itv,zv )
    
   ! Passed variables
    
   logical pbverbose
   integer natom
   _REAL_ sprob,stern,istrng
   integer level,xm,ym,zm,xmymzm,i,j,k
   _REAL_ h,gox,goy,goz
   _REAL_ gcrd(3,*),acrd(3,*),radi(*)
   integer narcdot,maxarc,marc(*),m2narc(maxarc,*),fstarc(*),arcatm(2,*)
   _REAL_ arccrd(3,*),savarc(3,*)
   integer atmsas(xmymzm),insas(xmymzm), itv(xmymzm)
   _REAL_ zv(xmymzm)
    
   ! Local variables
    
   integer iatm
   _REAL_ xi, yi, zi
   _REAL_ range1, r, rh
    
   rh = ONE/h
    
   ! initialize the grid owner list, which records the closest atom to the grid
   ! initialize the grid volume flag to 0, i.e. in the solvent region
    
   ! part 0, if salt presents compute ion exclusion map
   ! mark grid points within the stern layer as 1
   ! itv() saves the ion exclusion map
    
   if ( istrng /= ZERO .and. stern /= sprob ) then
    
   itv(1:xmymzm) = 0
   atmsas(1:xmymzm) = 0
   zv(1:xmymzm) = 9999.0d0
   do iatm = 1, natom
      range1 = radi(iatm)
      if ( range1 == ZERO ) cycle
      range1 = (range1+stern)*rh
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      call exvwsph( 1, xm, ym, zm, itv, atmsas, zv )
   end do
    
   end if
    
   ! part a, mark grid points within the solvent accessible surface as 1
   ! zv() records the distance of the point to the vdw surface of the closest atom
   ! also, if sas is used for the stern layer save it right after sas calculation.
    
   insas(1:xmymzm) = 0
   atmsas(1:xmymzm) = 0
   zv(1:xmymzm) = 9999.0d0
   do iatm = 1, natom
      range1 = radi(iatm)
      if ( range1 == ZERO ) cycle
      range1 = (range1+sprob)*rh; r = radi(iatm)*rh;
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      call exsasph( 1, xm, ym, zm, insas, atmsas, zv )
   end do
    
   if ( istrng /= ZERO .and. stern == sprob ) itv(1:xmymzm) = insas(1:xmymzm)
    
   ! part b, mark grid points within the van der waals surface as 2
   ! zv() records the distance of the point to the center of the closest atom
    
   zv(1:xmymzm) = 9999.0d0
   do iatm = 1, natom
      range1 = radi(iatm)
      if ( range1 == ZERO ) cycle
      range1 = range1*rh
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      call exvwsph( 2, xm, ym, zm, insas, atmsas, zv )
   end do
    
   ! part c, mark 1 grid points that don't belong to the reentry region as -2
    
   call contact( xm, ym, zm, insas, atmsas )
    
   ! part d, mark 1 grid points within any solvent reentry probes as -1
    
   range1 = sprob*rh
   zv(1:xmymzm) = 9999.0d0
   do iatm = 1, narcdot
      xi = (arccrd(1,iatm) - gox)*rh; yi = (arccrd(2,iatm) - goy)*rh; zi = (arccrd(3,iatm) - goz)*rh
      call exresph( -1, xm, ym, zm, insas, atmsas, zv )
   end do
    
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Mark grid points within atomic sasurf
subroutine exsasph( dielsph,xm,ym,zm,insph,inatm,dst )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Mark grid points within an atom  (dielectric constant dielsph) as dielsph
   ! Modified from UHBD (Comp. Phys. Comm. 91:57-95, 1995) routines excrx() and
   ! exsrfx() by Michael Gilson and Malcom Davis.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
   ! Passed variables
    
   integer dielsph
   integer xm, ym, zm    
   integer insph(xm,ym,zm)
   integer inatm(xm,ym,zm)
   _REAL_ dst(xm,ym,zm)
    
   ! Local variables
    
   integer i, j, k
   integer lowi, lowj, lowk
   integer highi, highj, highk
   _REAL_ range2, range3
   _REAL_ d, d2, r2

   r2 = r*r
   lowk = int(zi - range1) + 1; highk = int(zi + range1)
   do k = lowk, highk

      range2 = sqrt(range1**2-(zi-REAL(k))**2)
      lowj = int(yi - range2) + 1; highj = int(yi + range2)
      do j = lowj, highj
          
         range3 = sqrt(range2**2-(yi-REAL(j))**2)
         if ( range3 > ZERO ) then
             
            lowi = int(xi - range3) + 1; highi = int(xi + range3)
            do i = lowi, highi
               d2 = (i-xi)**2 + (j-yi)**2 + (k-zi)**2

               ! outside the vdw sphere

               if ( d2 > r2 ) then
                  d = sqrt(d2); d = d - r
                  if ( insph(i,j,k) == dielsph ) then
                     if ( d < dst(i,j,k) ) then
                        inatm(i,j,k) = iatm; dst(i,j,k) = d
                     end if
                     cycle
                  else
                     insph(i,j,k) = dielsph;
                     inatm(i,j,k) = iatm; dst(i,j,k) = d
                  end if

               ! within the vdw sphere

               else
                  if ( insph(i,j,k) == dielsph ) cycle
                  insph(i,j,k) = dielsph; inatm(i,j,k) = iatm; dst(i,j,k) = ZERO
               end if
            end do  ! i = lowi, highi
             
         end if  ! ( range3 > ZERO )
          
      end do  ! j = lowj, highj
       
   end do  ! k = lowk, highk
             
end subroutine exsasph
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Mark grid points within atomic vdw surf
subroutine exvwsph( dielsph,xm,ym,zm,insph,inatm,dst )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Mark grid points within an atom (dielectric constant dielsph) as dielpsh.
   ! Modified from UHBD (Comp. Phys. Comm. 91:57-95, 1995) routines excrx() and
   ! exsrfx() by Michael Gilson and Malcom Davis.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
   ! Passed variables
    
   integer dielsph
   integer xm, ym, zm
   integer insph(xm,ym,zm)
   integer inatm(xm,ym,zm)
   _REAL_ dst(xm,ym,zm)
    
   ! Local variables
    
   integer  i, j, k
   integer  lowi, lowj, lowk
   integer  highi, highj, highk
   _REAL_ range2, range3, d2
    
   lowk = int(zi - range1) + 1; highk = int(zi + range1)
   do k = lowk, highk
       
      range2 = sqrt(range1**2-(zi-REAL(k))**2)
      lowj = int(yi - range2) + 1; highj = int(yi + range2)
      do j = lowj, highj
          
         range3 = sqrt(range2**2-(yi-REAL(j))**2)
         if ( range3 > ZERO ) then
             
            lowi = int(xi - range3) + 1; highi = int(xi + range3)
            do i = lowi, highi
               d2 = (i-xi)**2 + (j-yi)**2 + (k-zi)**2

               ! it's marked already, make sure it belongs to the cloest atom

               if ( insph(i,j,k) == dielsph ) then

                  if ( d2 < dst(i,j,k) ) then
                     inatm(i,j,k) = iatm; dst(i,j,k) = d2
                  end if
                  cycle

               ! it's not marked yet

               else
                  insph(i,j,k) = dielsph;
                  inatm(i,j,k) = iatm; dst(i,j,k) = d2
               end if
            end do  ! i = lowi, highi
             
         end if  ! ( range3 > ZERO )
          
      end do  ! j = lowj, highj
       
   end do  ! k = lowk, highk
    
end subroutine exvwsph
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Mark contact grid points between vdw and sas surfaces
subroutine contact( xm,ym,zm,insas,atmsas )

   integer xm, ym, zm
   integer insas(xm,ym,zm), atmsas(xm,ym,zm)
    
   integer i, j, k, ip, buffer, iatm, ii, jj, iarc, inside
   _REAL_ xg(3), xi(3), xj(3), xij(3), dx(3), rxij, rd
   _REAL_ cosgij, cosgji, cosaij, cosaji

   _REAL_, parameter :: small = 0.01d0
    
   buffer = 1
   do k = 1+buffer, zm-buffer; do j = 1+buffer, ym-buffer; do i = 1+buffer, xm-buffer
       
      if ( insas(i,j,k) /= 1 ) cycle
       
      xg(1) = gox + i*h; xg(2) = goy + j*h; xg(3) = goz + k*h
       
      ! this is the atom that marked this grid within sasrf, so it will be the
      ! grid's conatct atom if it is marked so.
       
      iatm = atmsas(i,j,k)

      xi(1) = acrd(1,iatm); xi(2) = acrd(2,iatm); xi(3) = acrd(3,iatm)
       
      ! go through all arcs that this atom generates in solvent_accessibility.circle() 
       
      inside = -2
      do ip = 1, marc(iatm)
         iarc = m2narc(ip,iatm)
          
         ! generated by outer loop, i.e. the atom is iatm in solvent_accessibility.circle()
          
         if ( iarc >= fstarc(iatm) ) then
            jj = arcatm(1,iarc)
            xj(1) = acrd(1,jj)
            xj(2) = acrd(2,jj)
            xj(3) = acrd(3,jj)
            cosaij = savarc(1,iarc); cosaji = savarc(2,iarc)
          
         ! generated by inner loop, i.e. the atom is jatm in solvent_accessibility.circle()
          
         else
            xj = xi ! RL, warning ::::: why xj is assigned twice
            ii = arcatm(2,iarc)
            xj(1) = acrd(1,ii)
            xj(2) = acrd(2,ii)
            xj(3) = acrd(3,ii)
            cosaji = savarc(1,iarc); cosaij = savarc(2,iarc)
         end if
         rxij = savarc(3,iarc)
         xij = rxij*(xj - xi)
          
         dx = xg - xi; rd = ONE/sqrt(dx(1)**2+dx(2)**2+dx(3)**2); dx = rd*dx
         cosgij =  (xij(1)*dx(1)+xij(2)*dx(2)+xij(3)*dx(3))-small
          
         dx = xg - xj; rd = ONE/sqrt(dx(1)**2+dx(2)**2+dx(3)**2); dx = rd*dx
         cosgji = -(xij(1)*dx(1)+xij(2)*dx(2)+xij(3)*dx(3))-small
          
         ! if gij < aij .and. gji < aji, this is a reentry grid
          
         if ( cosgij <= cosaij .or. cosgji <= cosaji ) cycle
         inside = 1
         exit
      end do
       
      insas(i,j,k) = inside

   end do; end do; end do
    
end subroutine contact
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Mark grid points within reentry surf
subroutine exresph( dielsph,xm,ym,zm,insph,inatm,dst )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Mark grid points within a renentry sphere (dielectric constant dielsph)
   ! as dielsph. Modified from UHBD (Comp. Phys. Comm. 91:57-95, 1995) routines
   ! excrx() and exsrfx() by Michael Gilson and Malcom Davis.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Passed variables

   integer dielsph
   integer xm, ym, zm
   integer insph(xm,ym,zm)
   integer inatm(xm,ym,zm)
   _REAL_ dst(xm,ym,zm)

   ! Local variables

   integer  i, j, k
   integer  lowi, lowj, lowk
   integer  highi, highj, highk
   _REAL_ range2, range3, d2

   lowk = int(zi - range1) + 1; highk = int(zi + range1)
   do k = lowk, highk
       
      range2 = sqrt(range1**2-(zi-REAL(k))**2)
      lowj = int(yi - range2) + 1; highj = int(yi + range2)
      do j = lowj, highj
          
         range3 = sqrt(range2**2-(yi-REAL(j))**2)
         if ( range3 > ZERO ) then
             
            lowi = int(xi - range3) + 1; highi = int(xi + range3)
            do i = lowi, highi

               ! for a grid point in the contact region that is solvent accessible 
               ! do nothing

               if ( insph(i,j,k) == -2 .or. insph(i,j,k) == -3 ) cycle

               ! for a grid point that is within a solvent sphere, need to remember its
               ! cloest atom for future projection

               d2 = (i-xi)**2 + (j-yi)**2 + (k-zi)**2
               if ( insph(i,j,k) == dielsph ) then
                  if ( d2 < dst(i,j,k) ) then
                     inatm(i,j,k) = iatm; dst(i,j,k) = d2
                  end if
               else
                  insph(i,j,k) = dielsph
                  inatm(i,j,k) = iatm; dst(i,j,k) = d2
               end if
            end do  ! i = lowi, highi
             
         end if  ! ( range3 > ZERO )
          
      end do  ! j = lowj, highj
       
   end do  ! k = lowk, highk
    
end subroutine exresph

end subroutine pb_exmol_ses
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ driver for molecular volume assignment
!  modified van der Waals approach
subroutine pb_exmol_mvdw( pbverbose,natom,xm,ym,zm,xmymzm,h,&
                          gcrd,radipd,&
                          atmsas,insas,zv )
    
   ! Passed variables
    
   logical pbverbose
   integer natom
   integer xm,ym,zm,xmymzm
   _REAL_ h
   _REAL_ gcrd(3,*),radipd(*)
   integer atmsas(xmymzm),insas(xmymzm)
   _REAL_ zv(xmymzm)
    
   ! Local variables
    
   integer iatm
   _REAL_ xi, yi, zi
   _REAL_ range1, rh
   _REAL_ wf0, wf1, exposure, increase
    
   ! part a: dynamic radii are shared with the spline surface, and are assigned
   ! in sa_driver() 
   !
   ! part b: compute van der Waals volume
    
   ! initialize the grid owner list, i.e. which atom owns the grid
   ! initialize the grid volume flag to -2, i.e. in the solvent region
    
   atmsas(1:xmymzm) = 0
   insas(1:xmymzm) = -2
    
   ! mark volume within vdw srf as 2, outside -2, so there are only contact-type
   ! boundary grid points
    
   rh = ONE/h
   zv(1:xmymzm) = 9999.0d0
   do iatm = 1, natom
      if ( radipd(iatm) == ZERO ) cycle
      range1 = radipd(iatm)*rh; xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      call exmvwsph( 2, xm, ym, zm, insas, atmsas, zv(1) )
   end do
    
   ! note parts c & d from Lu and Luo, JCP, 2004 are removed for clarity of the code

contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Mark grid points within atomic vdw surf
subroutine exmvwsph( dielsph,xm,ym,zm,insph,inatm,dst )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Mark grid points within an atom (dielectric constant dielsph) as dielpsh.
   ! Modified from UHBD (Comp. Phys. Comm. 91:57-95, 1995) routines excrx() and
   ! exsrfx() by Michael Gilson and Malcom Davis.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
   ! Passed variables
    
   integer dielsph
   integer xm, ym, zm    
   integer insph(xm,ym,zm)
   integer inatm(xm,ym,zm)
   _REAL_ dst(xm,ym,zm)
    
   ! Local variables
    
   integer  i, j, k
   integer  lowi, lowj, lowk
   integer  highi, highj, highk
   _REAL_ range2, range3, d2
    
   lowk = int(zi - range1) + 1; highk = int(zi + range1)
   do k = lowk, highk
       
      range2 = sqrt(range1**2-(zi-REAL(k))**2)
      lowj = int(yi - range2) + 1; highj = int(yi + range2)
      do j = lowj, highj
          
         range3 = sqrt(range2**2-(yi-REAL(j))**2)
         if ( range3 > ZERO ) then
             
            lowi = int(xi - range3) + 1; highi = int(xi + range3)
            do i = lowi, highi
               d2 = (i-xi)**2 + (j-yi)**2 + (k-zi)**2
                
               ! it's marked already, make sure it belongs to the cloest atom
                
               if ( insph(i,j,k) == dielsph ) then
                
                  if ( d2 < dst(i,j,k) ) then
                     inatm(i,j,k) = iatm; dst(i,j,k) = d2
                  end if
                  cycle
                
               ! it's not marked yet
                
               else
                  insph(i,j,k) = dielsph;
                  inatm(i,j,k) = iatm; dst(i,j,k) = d2
               end if
            end do  ! i = lowi, highi
             
         end if  ! ( range3 > ZERO )
          
      end do  ! j = lowj, highj
       
   end do  ! k = lowk, highk
    
end subroutine exmvwsph

end subroutine pb_exmol_mvdw
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ driver for molecular volume assignment
!  spline approach
subroutine pb_exmol_spline( pbverbose,natom,sprob,&
                         level,xm,ym,zm,xmymzm,h,gox,goy,goz,&
                         gcrd,radi,radip,radipd,&
                         nn,dash,spcoef,&
                         atmsas,insas,zv,tv )
    
   ! Passed variables
    
   logical pbverbose
   integer natom
   _REAL_ sprob
   integer level,xm,ym,zm,xmymzm,grid_reentry(3,30000),i,j,k,l,count_num1(20),count_num2(20)
   _REAL_ h,gox,goy,goz
   _REAL_ gcrd(3,*),radi(*),radip(*),radipd(*)
   integer nn
   _REAL_ dash(nn), spcoef(nn-1,4)
   integer atmsas(xmymzm),insas(xmymzm)
   _REAL_ zv(xmymzm),tv(xmymzm)
    
   ! Local variables
    
   integer iatm, num_reentry, flag
   _REAL_ xi, yi, zi
   _REAL_ range1, r, d2g, rh

   !if(level == 3) then
   !open(101,file='benchmark.dat')
   !open(109,file='contri-SES.dat')
   !open(110,file='contri-noSES.dat')
   !read(101,*)num_reentry
   !write(*,*)num_reentry
   !do i=1,num_reentry
   !read(101,*)grid_reentry(1,i),grid_reentry(2,i),grid_reentry(3,i)
   !enddo
   !rh = ONE/h
   !count_num1(1:20) = ZERO
   !count_num2(1:20) = ZERO
   !endif
    
   rh = ONE/h
   
   ! initialize the grid owner list, which records atom owns the grid
   ! initialize the grid volume flag to -1, i.e. in the solvent region
    
   atmsas(1:xmymzm) = 0
   insas(1:xmymzm) = -1
   
   ! part a, mark grid points within sa srf as 1

   do iatm = 1, natom
      range1 = radip(iatm)
      if ( range1 == ZERO ) cycle
      range1=range1*rh; xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      call exsasph( 1, xm, ym, zm, insas )
   end do

   ! part b, mark grid points within vdw srf as 2

   do iatm = 1, natom
      range1 = radipd(iatm)
      if ( range1 == ZERO ) cycle
      range1 = range1*rh; xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      call exvwsph( 2, xm, ym, zm, insas )
   end do

   ! part c, accumulate density at grid points that are 1 but not 2
   ! zv saves the final density

   zv(1:xmymzm) = ZERO
   tv(1:xmymzm) = ZERO
   do iatm = 1, natom
      range1 = radip(iatm)
      if ( range1 == ZERO ) cycle
      range1 = (range1+sprob)*rh; r = radi(iatm)*rh; d2g = ONE/( (TWO*sprob)*rh )
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      call density_calc( nn, dash, spcoef, xm, ym, zm, insas, zv(1), tv(1) )
   end do
   if(level == 3)then
   do k = 1, zm
   do j = 1, ym
   do i = 1, xm
      flag = 1+(i-1)+xm*(j-1)+xm*ym*(k-1)
      if ( insas(flag) == 1  ) then
           count_num2(int(tv(flag))) = count_num2(int(tv(flag))) + 1
           do l=1, num_reentry
           if(i==grid_reentry(1,l).and.j==grid_reentry(2,l).and.k==grid_reentry(3,l))then
           count_num1(int(tv(flag))) = count_num1(int(tv(flag))) + 1
           count_num2(int(tv(flag))) = count_num2(int(tv(flag))) - 1
           endif
           enddo
      endif
   end do
   end do
   end do
   write(109,'(20I6)')count_num1(1:20)
   write(110,'(20I6)')count_num2(1:20)
   endif

    do iatm = 1, natom
      range1 = radipd(iatm)
      if ( range1 == ZERO ) cycle
      range1 = range1*rh; xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      call exvwsph( 3, xm, ym, zm, insas )
   end do

   ! part f, mark grid points within vdw srf as 2

   do iatm = 1, natom
      range1 = radi(iatm)
      if ( range1 == ZERO ) cycle
      range1 = range1*rh; xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      call exvwsph( 2, xm, ym, zm, insas )
   end do


   ! part d, translate density to grid flags
   ! of contact v.s. reentry

   call density_convert( xm, ym, zm, insas, zv(1), tv(1) )

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Mark grid points within atomic sasurf
subroutine exsasph( dielsph,xm,ym,zm,insph )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Mark grid points within an atom  (dielectric constant dielsph) of index iatm
   ! as dielpsh. Modified from UHBD (Comp. Phys. Comm. 91:57-95, 1995) routines
   ! excrx() and exsrfx() by Michael Gilson and Malcom Davis.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   ! Passed variables

   integer  dielsph, xm, ym, zm
   integer  insph(xm,ym,zm)
 
   ! Local variables
    
   integer  i, j, k
   integer  lowi, lowj, lowk
   integer  highi, highj, highk
   _REAL_ range2, range3

   lowk = int(zi - range1) + 1; highk = int(zi + range1)
   do k = lowk, highk

      range2 = sqrt(range1**2-(zi-REAL(k))**2)
      lowj = int(yi - range2) + 1; highj = int(yi + range2)
      do j = lowj, highj
          
         range3 = sqrt(range2**2-(yi-REAL(j))**2)
         if ( range3 > ZERO ) then
             
            lowi = int(xi - range3) + 1; highi = int(xi + range3)
            do i = lowi, highi
               insph(i,j,k) = dielsph
            end do  ! i = lowi, highi
             
         end if  ! ( range3 > ZERO )
          
      end do  ! j = lowj, highj
       
   end do  ! k = lowk, highk
             
end subroutine exsasph
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Mark grid points within atomic vdw surf
subroutine exvwsph( dielsph,xm,ym,zm,insph )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Mark grid points within an atom (dielectric constant dielsph) of index iatm
   ! as dielpsh. Modified from UHBD (Comp. Phys. Comm. 91:57-95, 1995) routines
   ! excrx() and exsrfx() by Michael Gilson and Malcom Davis.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
   ! Passed variables
    
   integer  dielsph, xm, ym, zm
   integer  insph(xm,ym,zm)
    
   ! Local variables
    
   integer  i, j, k
   integer  lowi, lowj, lowk
   integer  highi, highj, highk
   _REAL_ range2, range3
    
   lowk = int(zi - range1) + 1; highk = int(zi + range1)
   do k = lowk, highk
       
      range2 = sqrt(range1**2-(zi-REAL(k))**2)
      lowj = int(yi - range2) + 1; highj = int(yi + range2)
      do j = lowj, highj
          
         range3 = sqrt(range2**2-(yi-REAL(j))**2)
         if ( range3 > ZERO ) then
             
            lowi = int(xi - range3) + 1; highi = int(xi + range3)
            do i = lowi, highi
               insph(i,j,k) = dielsph
            end do  ! i = lowi, highi
             
         end if  ! ( range3 > ZERO )
          
      end do  ! j = lowj, highj
       
   end do  ! k = lowk, highk
    
end subroutine exvwsph
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ compute density at grid points within mol sas surf
subroutine density_calc( nn,dash,spcoef,xm,ym,zm,insph,dens,atmctr )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Mark grid points within an atom  (dielectric constant dielsph) of index iatm
   ! as dielpsh. Modified from UHBD (Comp. Phys. Comm. 91:57-95, 1995) routines
   ! excrx() and exsrfx() by Michael Gilson and Malcom Davis.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Passed variables

   integer nn, xm, ym, zm
   _REAL_ dash(nn), spcoef(nn-1,4)
   integer  insph(xm,ym,zm)
   _REAL_ dens(xm,ym,zm)
   _REAL_ atmctr(xm,ym,zm)

   ! Local variables
    
   integer  i, j, k
   integer  lowi, lowj, lowk
   integer  highi, highj, highk
   _REAL_ range2, range3
   _REAL_ d, d2, point
   integer  l 
   
   lowk = int(zi - range1) + 1; highk = int(zi + range1)
   do k = lowk, highk

      range2 = sqrt(range1**2-(zi-REAL(k))**2)
      lowj = int(yi - range2) + 1; highj = int(yi + range2)
      do j = lowj, highj
          
         range3 = sqrt(range2**2-(yi-REAL(j))**2)
         if ( range3 > ZERO ) then
             
            lowi = int(xi - range3) + 1; highi = int(xi + range3)
            do i = lowi, highi

               ! no need to spline if it is outside sas

               if ( insph(i,j,k) == 1) then

                  d2 = (i-xi)**2 + (j-yi)**2 + (k-zi)**2; d = sqrt(d2)
                  d = d - r; point = d*d2g ! in the unit of (2*sprob/h)

                  ! spline on nn intervals in the range of [0, 2*sprob/h]

                  do l = 1, nn-1
                     if ( point > dash(l) .and. point <= dash(l+1) ) then
                        dens(i,j,k) = dens(i,j,k)                   +&
                                      spcoef(l,1)                   +&
                                      spcoef(l,2)*(point-dash(l))   +&
                                      spcoef(l,3)*(point-dash(l))**2+&
                                      spcoef(l,4)*(point-dash(l))**3
                     endif
                  enddo
                   
                  ! accumulate contributing atoms
                   
                  atmctr(i,j,k) = atmctr(i,j,k) + ONE
                   
               end if ! ( insph(i,j,k) == 1 )

            end do  ! i = lowi, highi
             
         end if  ! ( range3 > ZERO )
          
      end do  ! j = lowj, highj
       
   end do  ! k = lowk, highk
             
end subroutine density_calc
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ convert grid points with density < 1 to -2.
subroutine density_convert( xm,ym,zm,insas,dens,atmctr )

   integer xm, ym, zm
   integer insas(xm,ym,zm)
   _REAL_ dens(xm,ym,zm),fit_ratio(20)
   _REAL_ atmctr(xm,ym,zm)

   integer i, j, k
   fit_ratio(1) = 0.0
   fit_ratio(2) = 0.646
   fit_ratio(3) = 0.743
   fit_ratio(4) = 0.85348
   fit_ratio(5) = 0.9556739
   fit_ratio(6) = 1.03747
   fit_ratio(7) = 1.17028
   fit_ratio(8) = 1.422610
   fit_ratio(9) = 1.655190
   fit_ratio(10) = 1.7649
   fit_ratio(11) = 1.77995
   fit_ratio(12) = 1.77491
   fit_ratio(13) = 1.773760
   fit_ratio(14) = 1.81975
   fit_ratio(15) = 1.92714
   fit_ratio(16) = 2.05216
   fit_ratio(17) = 0.0
   fit_ratio(18) = 0.0
   fit_ratio(19) = 0.0
   fit_ratio(20) = 0.0

   ! for InsightII display
   !_REAL_ g(3)
   !
   !open (unit=55, file='gauss.dot')
   !write (55, '("DOTS")')
!   open (unit=101, file='check_dens.dat') 
   do k = 1, zm
   do j = 1, ym
   do i = 1, xm
    if ( insas(i,j,k) == 1 .and. atmctr(i,j,k) /= ZERO ) then
          if ( dens(i,j,k)*fit_ratio(atmctr(i,j,k)) <ONE )then
!           if ( dens(i,j,k) <ONE )then
           insas(i,j,k) = -2
           endif
      endif

   end do
   end do
   end do

   ! for InsightII display
   !close(55)

end subroutine density_convert

end subroutine pb_exmol_spline

end module pb_exmol
