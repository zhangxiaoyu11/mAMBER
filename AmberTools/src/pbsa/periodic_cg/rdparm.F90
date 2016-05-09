#include "copyright.h"
#include "dprec.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine rdparm1 here]
subroutine rdparm1(nf)
   
   implicit none
   
#  include "md.h"
#  include "memory.h"
#  include "parms.h"
#  include "files.h"
#  include "box.h"
   integer nf
   integer i,nspsol,iok
   integer nhparm,idum,nttyp
   integer mbper,mgper,mdper,mbona,mtheta,mphia ! read but ignored
   integer numextra
   character(len=80) fmt
   character(len=80) fmtin,ifmt,afmt,rfmt,type
   ifmt = '(12I6)'
   afmt = '(20A4)'
   rfmt = '(5E16.8)'
   
   !     ----- READ THE MOLECULAR TOPOLOGY -----
   
   nspsol = 0
   
   !     ----- FORMATTED INPUT -----
   
   fmtin = '(A80)'
   type = 'TITLE'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmtin) title
   
   fmtin = ifmt
   type = 'POINTERS'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) natom,ntypes,nbonh,mbona,ntheth,mtheta,nphih,mphia, &
         nhparm,nparm,nnb,nres,nbona,ntheta,nphia, &
         numbnd,numang,nptra,natyp,nphb,ifpert,nbper,ngper,ndper, &
         mbper,mgper,mdper,ifbox,nmxrs,ifcap,numextra,ncopy
   if (nparm == 1) then
      write(6,*) ' *** THIS VERSION WILL NOT ACCEPT TOPOLOGY FILES'
      write(6,*) '     THAT WERE CREATED BY ADDLES, WITH NPARM=1'
      write(6,*) '     USE A VERSION COMPILED WITH -DLES '
      call mexit(6,1)
   end if
   write(6,8118) natom,ntypes,nbonh,mbona,ntheth,mtheta,nphih,mphia, &
         nhparm,nparm,nnb,nres,nbona,ntheta,nphia,numbnd, &
         numang,nptra,natyp,nphb,ifbox,nmxrs,ifcap,numextra,ncopy
   8118 format(t2, &
         'NATOM  = ',i7,' NTYPES = ',i7,' NBONH = ',i7,' MBONA  = ',i7, &
         /' NTHETH = ',i7,' MTHETA = ',i7,' NPHIH = ',i7,' MPHIA  = ',i7, &
         /' NHPARM = ',i7,' NPARM  = ',i7,' NNB   = ',i7,' NRES   = ',i7, &
         /' NBONA  = ',i7,' NTHETA = ',i7,' NPHIA = ',i7,' NUMBND = ',i7, &
         /' NUMANG = ',i7,' NPTRA  = ',i7,' NATYP = ',i7,' NPHB   = ',i7, &
         /' IFBOX  = ',i7,' NMXRS  = ',i7,' IFCAP = ',i7,' NEXTRA = ',i7 &
         ,/' NCOPY  = ',i7/)
   
   !     --- make sure we do not exceed memory limits in commons ---
   
   nttyp = ntypes*(ntypes+1)/2
   if (numbnd > 5000 .or. numang > 900 .or. nptra > 1200 .or. &
         nphb > 200 .or. natyp > 60 .or. nttyp > 1830) then
      write(6,'(/,5x,a)') 'rdparm: a parameter array overflowed'
      write(6,'(/,5x,a)') '       (e.g. the table of dihedral params)'
      call mexit(6, 1)
   end if
   
   if(nbona /= mbona .or. ntheta /= mtheta .or. nphia /= mphia) then
      write(6,*) 'Sander no longer allows constraints in prmtop'
      write(6,*) '...must have nbona=mbona, ntheta=mtheta, nphi=mphi'
      call mexit(6,1)
   end if
   
   return
end subroutine rdparm1 

!=======================================================================

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine rdparm2 here]
subroutine rdparm2(x,ix,ih,ipairs,nf,i_stack)
   implicit none
   integer nf
   _REAL_ x(*)
   integer ix(*),ipairs(*),i_stack(*)
   character(len=4) ih(*)
   
#  include "md.h"
#  include "memory.h"
#  include "parms.h"
#  include "files.h"
#  include "box.h"
   integer nttyp,ntype,i,iok
   integer i1,i2,i3,i4,j,jj,k,l,n,nn
   integer iptres,nspsol,natsm,idum,ip14
   integer l_ib,l_jb,l_bt
   integer l_it,l_jt,l_kt,l_tt
   integer l_id,l_jd,l_kd,l_ld,l_dt
   integer bp,ibp,jbp,btp
   _REAL_ dumd,oldbeta,duma,dumb,dumc
   character(len=4) dumchar
   integer dumint
   _REAL_ dumfloat

   character(len=80) fmt
   character(len=80) fmtin,ifmt,afmt,rfmt,type
   ifmt = '(12I6)'
   afmt = '(20A4)'
   rfmt = '(5E16.8)'
   nttyp = ntypes*(ntypes+1)/2
   ntype = ntypes*ntypes
   
   !     ----- READ THE SYMBOLS AND THE CHARGES AND THE MASSES -----
   
   fmtin = afmt
   type = 'ATOM_NAME'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ih(m04+i-1),i = 1,natom)
   
   fmtin = rfmt
   type = 'CHARGE'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (x(l15+i-1),i = 1,natom)
   
   fmtin = rfmt
   type = 'MASS'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (x(lwinv+i-1),i = 1,natom)
   
   fmtin = ifmt
   type = 'ATOM_TYPE_INDEX'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i04+i-1),i = 1,natom)
   
   fmtin = ifmt
   type = 'NUMBER_EXCLUDED_ATOMS'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i+i08-1),i = 1,natom)
   
   fmtin = ifmt
   type = 'NONBONDED_PARM_INDEX'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i+i06-1),i = 1,ntype)
   
   fmtin = afmt
   type = 'RESIDUE_LABEL'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ih(i+m02-1),i=1,nres)
   
   fmtin = ifmt
   type = 'RESIDUE_POINTER'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i+i02-1),i=1,nres)
   ix(i02+nres) = natom+1
   
   !     ----- READ THE PARAMETERS -----
   
   fmtin = rfmt
   type = 'BOND_FORCE_CONSTANT'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (rk(i),    i = 1,numbnd)
   
   fmtin = rfmt
   type = 'BOND_EQUIL_VALUE'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (req(i),   i = 1,numbnd)
   
   fmtin = rfmt
   type = 'ANGLE_FORCE_CONSTANT'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (tk(i),    i = 1,numang)
   
   fmtin = rfmt
   type = 'ANGLE_EQUIL_VALUE'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (teq(i),   i = 1,numang)
   
   fmtin = rfmt
   type = 'DIHEDRAL_FORCE_CONSTANT'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (pk(i),    i = 1,nptra)
   
   fmtin = rfmt
   type = 'DIHEDRAL_PERIODICITY'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (pn(i),    i = 1,nptra)
   
   fmtin = rfmt
   type = 'DIHEDRAL_PHASE'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (phase(i), i = 1,nptra)
   
   fmtin = rfmt
   type = 'SOLTY'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (solty(i), i = 1,natyp)
   
   fmtin = rfmt
   type = 'LENNARD_JONES_ACOEF'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (cn1(i),   i = 1,nttyp)
   
   fmtin = rfmt
   type = 'LENNARD_JONES_BCOEF'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (cn2(i),   i = 1,nttyp)
   
   !     ----- READ THE BONDING INFORMATIONS -----
   
   fmtin = ifmt
   type = 'BONDS_INC_HYDROGEN'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i+iibh-1),ix(i+ijbh-1),ix(i+iicbh-1), &
         i = 1,nbonh)
   
   fmtin = ifmt
   type = 'BONDS_WITHOUT_HYDROGEN'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt)(ix(i+iiba-1),ix(i+ijba-1),ix(i+iicba-1),i = 1,nbona)
   
   fmtin = ifmt
   type = 'ANGLES_INC_HYDROGEN'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i+i24-1),ix(i+i26-1),ix(i+i28-1),ix(i+i30-1), &
         i = 1,ntheth)
   
   fmtin = ifmt
   type = 'ANGLES_WITHOUT_HYDROGEN'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i+i32-1),ix(i+i34-1),ix(i+i36-1),ix(i+i38-1), &
         i = 1,ntheta)
   
   fmtin = ifmt
   type = 'DIHEDRALS_INC_HYDROGEN'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i+i40-1),ix(i+i42-1),ix(i+i44-1),ix(i+i46-1), &
         ix(i+i48-1),i = 1,nphih)
   
   fmtin = ifmt
   type = 'DIHEDRALS_WITHOUT_HYDROGEN'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i+i50-1),ix(i+i52-1),ix(i+i54-1),ix(i+i56-1), &
         ix(i+i58-1),i = 1,nphia)
   
   fmtin = ifmt
   type = 'EXCLUDED_ATOMS_LIST'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i+i10-1),i=1,nnb)
   
   !     ----- READ THE H-BOND PARAMETERS -----
   
   fmtin = rfmt
   type = 'HBOND_ACOEF'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (asol(i),i=1,nphb)
   
#ifndef HAS_10_12
   do i=1,nphb
      if( asol(i) /= 0.d0 ) then
         write(6,*) 'Found a non-zero 10-12 coefficient, but source', &
               ' was not compiled with -DHAS_10_12.'
         write(6,*) 'If you are using a pre-1994 force field, you', &
               ' will need to re-compile with this flag.'
         call mexit(6,1)
      end if
   end do
#endif
   
   fmtin = rfmt
   type = 'HBOND_BCOEF'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (bsol(i),i=1,nphb)
   
   fmtin = rfmt
   type = 'HBCUT'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (hbcut(i),i=1,nphb)
   
   !     ----- READ ISYMBL,ITREE,JOIN AND IROTAT ARRAYS -----
   
   fmtin = afmt
   type = 'AMBER_ATOM_TYPE'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ih(i+m06-1),i=1,natom)
   
   fmtin = afmt
   type = 'TREE_CHAIN_CLASSIFICATION'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ih(i+m08-1),i=1,natom)
   
   fmtin = ifmt
   type = 'JOIN_ARRAY'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i+i64-1),i=1,natom)
   
   fmtin = ifmt
   type = 'IROTAT'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i+i66-1),i=1,natom)
   
   !     ----- READ THE BOUNDARY CONDITION STUFF -----
   
   nspm = 1
   ix(i70) = natom
   if (ifbox > 0) then
      
      fmtin = ifmt
      type = 'SOLVENT_POINTERS'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) iptres,nspm,nspsol
      
      fmtin = ifmt
      type = 'ATOMS_PER_MOLECULE'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) (ix(i+i70-1),i=1,nspm)
      
      fmtin = rfmt
      type = 'BOX_DIMENSIONS'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) oldbeta,duma,dumb,dumc
      
      !       ---(above values are still read, for backward compatibility, but
      !           ignored. Box info must come from the coord. file or from
      !           the &ewald namelist of the input file)
      
      if( igb /= 0  .or.  ntb == 0 )then
         box(1)=0.0d0
         box(2)=0.0d0
         box(3)=0.0d0
      end if
      
   end if  ! (ifbox > 0)
   
   !     ----- LOAD THE CAP INFORMATION IF NEEDED -----
   
   if(ifcap > 0) then
      fmtin = '(I6)'
      type = 'CAP_INFO'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) natcap
      
      fmtin = '(4E16.8)'
      type = 'CAP_INFO2'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) cutcap,xcap,ycap,zcap
   else
      natcap = 0
      cutcap = 0.0d0
      xcap = 0.0d0
      ycap = 0.0d0
      zcap = 0.0d0
   end if
   
   if( igb /= 0 .and. iok == -1 ) then
      write(0,*) 'GB calculations now require a new-style prmtop file'
      write(6,*) 'GB calculations now require a new-style prmtop file'
      call mexit(6,1)
   end if
   
   if ( igb /= 0 ) then
      fmtin = rfmt
      type = 'RADII'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) (x(l97+i-1),i=1,natom)
      type = 'SCREEN'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) (x(l96+i-1),i=1,natom)
   end if
   
   ! -------------------------------------------------------------------
   !     --- if icfe is set, read in the pert-types and charges
   ! -------------------------------------------------------------------
   
   if( icfe /= 0 ) then
      write(6,'(/a,f6.3)') &
            'Running a free energy calculation with lambda = ',clambda
      
      if( ifpert == 0 ) then
         write(6,'(a)') 'Must use a pert. prmtop file for free energy'
         call mexit(6,1)
      end if
      
      fmtin = ifmt
      type = 'PERT_BOND_ATOMS'
      call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
      read(nf,fmt) (ix(iiba+nbona+i-1),ix(ijba+nbona+i-1),i=1,nbper)
      
      fmtin = ifmt
      type = 'PERT_BOND_PARAMS'
      call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
      read(nf,fmt) (ix(iicba+nbona+i-1),i=1,2*nbper)
      
      fmtin = ifmt
      type = 'PERT_ANGLE_ATOMS'
      call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
      read(nf,fmt) (ix(i32+ntheta+i-1),ix(i34+ntheta+i-1), &
            ix(i36+ntheta+i-1), i=1,ngper)
      
      fmtin = ifmt
      type = 'PERT_ANGLE_PARAMS'
      call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
      read(nf,fmt) (ix(i38+ntheta+i-1),i=1,2*ngper)
      
      fmtin = ifmt
      type = 'PERT_DIHEDRAL_ATOMS'
      call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
      read(nf,fmt) (ix(i50+nphia+i-1),ix(i52+nphia+i-1), &
            ix(i54+nphia+i-1),ix(i56+nphia+i-1),i=1,ndper)
      
      fmtin = ifmt
      type = 'PERT_DIHEDRAL_PARAMS'
      call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
      read(nf,fmt) (ix(i58+nphia+i-1),i=1,2*ndper)
      
      fmtin = afmt
      type = 'PERT_RESIDUE_NAME'
      call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
      read(nf,fmt) (dumchar,i=1,nres)
      
      fmtin = afmt
      type = 'PERT_ATOM_NAME'
      call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
      read(nf,fmt) (dumchar,i=1,natom)
      
      fmtin = afmt
      type = 'PERT_ATOM_SYMBOL'
      call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
      read(nf,fmt) (dumchar,i=1,natom)
      
      fmtin = rfmt
      type = 'ALMPER'
      call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
      read(nf,fmt) (dumfloat,i=1,natom)
      
      fmtin = ifmt
      type = 'IAPER'
      call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
      read(nf,fmt) (dumint,i=1,natom)
      
      fmtin = ifmt
      type = 'PERT_ATOM_TYPE_INDEX'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) (ix(i84-1+i), i=1,natom)
      
      fmtin = rfmt
      type = 'PERT_CHARGE'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) (x(l190-1+i), i=1,natom)
      
#ifdef CHARGE_MIXING
      if( igb > 0 .and. igb /= 10 ) then
         
         !         TI using dV/dlambda:
         !         L190 will hold the charge increments (i.e. state1 - state0 )
         !         L15  will hold the charges for the current value of lambda
         
         do i=1,natom
            x(l190-1+i) = x(l190-1+i) - x(l15-1+i)
            x(l15-1+i) = x(l15-1+i) + clambda*x(l190-1+i)
         end do
         
      else
         write(6,*) 'sander with CHARGE_MIXING requires igb>0'
         call mexit(6,1)
      end if
#endif

   else
      
      if( ifpert /= 0 ) then
         write(6,'(a)') 'Cannot use a pert. prmtop unless icfe>0'
         call mexit(6,1)
      end if
      
   end if  ! ( icfe /= 0 )
   
   !     ----- CALCULATE INVERSE, TOTAL MASSES -----
   
   !       -- save the masses for removal of mass weighted velocity,
   !          leaving the inverse masses in the legacy, Lwinv area
   
   
   tmass = 0.0d0
   !     -- index over molecules
   j = l75-1
   jj = i70-1
   !     -- index over mass->invmass
   k = lwinv-1
   !     -- index over saved mass
   l = lmass-1
   do n = 1,nspm
      j = j + 1
      jj = jj + 1
      x(j) = 0.0d0
      natsm = ix(jj)
      do nn = 1,natsm
         k = k+1
         l = l+1
         
         !         -- sum molecule
         
         x(j) = x(j) + x(k)
         
         !         -- save mass in "new" Lmass area
         
         x(l) = x(k)
         
         !         -- make inverse in "old" Lwinv area
         
         x(k) = 1.0d0 / x(k)
      end do
      tmass = tmass + x(j)
   end do
   tmassinv = 1.0d0 / tmass
   
   !     ----- SCALE THE CHARGES IF DIELC.GT.1.0E0 -----
   
   if (dielc /= 1.0e0) then
      dumd = sqrt(dielc)
      do i = 1,natom
         x(i+l15-1) = x(i+l15-1)/dumd
      end do
   end if
   
   !     ----- INVERT THE HBCUT ARRAY -----
   
   do i = 1,nphb
      if(hbcut(i) <= 0.001e0) hbcut(i) = 1.0d-10
      hbcut(i) = 1.0e0/hbcut(i)
   end do
   
   !     ----- duplicate dihedral pointers for vector ephi -----
   
   call dihdup(nphih,ix(i40),ix(i42),ix(i44),ix(i46),ix(i48),pn)
   call dihdup(nphia,ix(i50),ix(i52),ix(i54),ix(i56),ix(i58),pn)
   
   !     --- pre-calculate some parameters for vector ephi ---
   
   call dihpar(nptra,pk,pn,phase,gamc,gams,ipn,fmn)
   
#ifdef CHARMM
   
   !    ---read in   1-4 parameters at end of parm file:
   
   read(nf,9128) (cn114(i),   i = 1,nttyp)
   read(nf,9128) (cn214(i),   i = 1,nttyp)
   
   !    --- read in   urey-bradley parameters:
   
   read(nf,9128) (rkub(i),i=1,numang)
   read(nf,9128) (rub(i),i=1,numang)
#endif
   
   return
   9108 format(20a4)
   9128 format(5e16.8)
   9138 format(80A1)
end subroutine rdparm2 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine istuff here]
subroutine istuff(i,j,iarray,k)
   
   ! routine to correctly load a strange shaped 2 dim matrix
   
   dimension iarray(15,*)
   iarray(i,j) = k
   return
end subroutine istuff 
