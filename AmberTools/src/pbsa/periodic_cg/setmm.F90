#include "copyright.h"
#include "dprec.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine setvar here]
subroutine setvar(ix,belly)
   
   implicit none
   logical belly
   integer i
   
   !     ----- ROUTINE TO DO THE NECESSARY ACCOMODATIONS FOR PROTEIN
   !           BELLY MINIMISATIONS -----
   
#  include "memory.h"
#  include "box.h"
   
   integer ix(*)
   
   !     --- SETUP THE BELLY GROUP ARRAY AND RESIDUE ARRAY for resnba ---
   
   if(belly) goto 190
   do i = 1,natom
      ix(ibellygp+i-1) = 1
   end do
   do i = 1,nres
      ix(i66+i-1) = 1
   end do
   call setbon(nbonh,nbonh,ix(iibh),ix(ijbh),ix(iicbh),ix(ibellygp))
   call setbon(nbona,nbona,ix(iiba),ix(ijba),ix(iicba),ix(ibellygp))
   return
   
   !     ----- DELETE BONDS WHICH ARE IN THE BELLY ALONE -----
   
   190 continue
   call setbon(nbonh,nbonh,ix(iibh),ix(ijbh),ix(iicbh),ix(ibellygp))
   call setbon(nbona,nbona,ix(iiba),ix(ijba),ix(iicba),ix(ibellygp))
   
   !     ----- MAKE THE BOND ARRAYS SEQUENTIAL FOR SHAKE AND FORCE
   !           ROUTINES -----
   
   do i = 1,nbona
      ix(iibh+nbonh+i-1)  = ix(iiba+i-1)
      ix(ijbh+nbonh+i-1)  = ix(ijba+i-1)
      ix(iicbh+nbonh+i-1) = ix(iicba+i-1)
   end do
   iiba = iibh+nbonh
   ijba = ijbh+nbonh
   iicba = iicbh+nbonh
   
   !     ----- DELETE THE BONDS WHICH ARE IN THE BELLY ALONE -----
   
   call setang(ntheth,ntheth,ix(i24),ix(i26),ix(i28),ix(i30), &
         ix(ibellygp))
   call setang(ntheta,ntheta,ix(i32),ix(i34),ix(i36),ix(i38), &
         ix(ibellygp))
   
   !     ----- MAKE THE ANGLE ARRAYS SEQUENTIAL -----
   
   do i = 1,ntheta
      ix(i24+ntheth+i-1) = ix(i32+i-1)
      ix(i26+ntheth+i-1) = ix(i34+i-1)
      ix(i28+ntheth+i-1) = ix(i36+i-1)
      ix(i30+ntheth+i-1) = ix(i38+i-1)
   end do
   i32 = i24+ntheth
   i34 = i26+ntheth
   i36 = i28+ntheth
   i38 = i30+ntheth
   
   !     ----- DELETE THE DIHEDRALS -----
   
   call setdih(nphih,nphih,ix(i40),ix(i42),ix(i44),ix(i46), &
         ix(i48),ix(ibellygp))
   call setdih(nphia,nphia,ix(i50),ix(i52),ix(i54),ix(i56), &
         ix(i58),ix(ibellygp))
   
   !     ----- MAKE THE DIHEDRALS SEQUENTIAL ------
   
   do i = 1,nphia
      ix(i40+nphih+i-1) = ix(i50+i-1)
      ix(i42+nphih+i-1) = ix(i52+i-1)
      ix(i44+nphih+i-1) = ix(i54+i-1)
      ix(i46+nphih+i-1) = ix(i56+i-1)
      ix(i48+nphih+i-1) = ix(i58+i-1)
   end do
   i50 = i40+nphih
   i52 = i42+nphih
   i54 = i44+nphih
   i56 = i46+nphih
   i58 = i48+nphih
   
   !     ----- FIND THE TOTAL NUMBER OF ACTIVE ATOMS AND RESIDUES -----
   
   call setatm(natom,nres,natbel,ix(i02),ix(ibellygp),ix(i66))
   
   !     ----- ALL DONE RETURN -----
   
   return
end subroutine setvar 
