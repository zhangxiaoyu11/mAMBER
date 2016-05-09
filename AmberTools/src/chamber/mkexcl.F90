! <compile=optimized>
#include "copyright.h"

!--------------------------------------------------------------
! CHAMBER: Exclusion list
!
! Description: Routines for building, 1-2, 1-3 and 1-4
!              exclusion lists.
!
! Authors: Mark Williamson (SDSC, 2009)
!          Michael Crowley (NREL, 2009)
!          Ross Walker (SDSC, 2009)
!
!--------------------------------------------------------------

subroutine make_exclusion_list(nonbond_excl_14,nonbond_excl_pointers_14,maxexcl,done)

   !Builds the nonbond exclusion list for CHAMBER.

   use sort2
   use psfprm, only:nbond,natom,outu,verbose_debug,nnb,jb,ib

   implicit none

   integer,intent(in) :: maxexcl

   integer,intent(out),dimension(*) :: nonbond_excl_14,nonbond_excl_pointers_14

   logical, intent(out) :: done

   integer, allocatable, dimension(:) :: itmp,jtmp


   integer nnb14

   integer,allocatable,dimension(:) :: natbon,ipk,jpk,ipk14,jpk14
   integer, allocatable,dimension(:,:) :: iatbon
   logical lex14

   integer modex,npair,npair4,maxwrk,mxwrk4
   integer nx14,iatom,next14
   integer i2,j2,i,j,k,ij,ija,ik,ika,i2l,j2l,ibt,jbt
   integer :: numexi,lasti,index
   integer, parameter :: iatmxb=10 !Max number of bonds to a single atom
   logical ok
!
   done=.false.

   allocate(iatbon(iatmxb, natom))
   allocate(natbon(natom))
   allocate(ipk(maxexcl))
   allocate(jpk(maxexcl))
   allocate(ipk14(maxexcl))
   allocate(jpk14(maxexcl))

! The modex value controls how much goes on the exclusion list. This is the NUMBER_EXCLUDED_ATOMS and
! EXCLUDED_ATOMS_LIST in the prmtop file. The most useful options are:
! 3 = 1-2 and 1-3 interactions are added to the excluded list.
! 5 = 1-2, 1-3 and 1-4 interactions are added to the excluded list.
!
! The charmm force field does not scale 1-4 interactions and therefore there is no reason to exclude
! them from the standard EEL and VDW calculations and instead have a separate loop over unique dihedrals
! like AMBER needs to do for SCEE and SCNB scaling. Except that for VDW charmm has a different set of
! parameters for the A and B coeffs for 1-4 terms so we have to do this inside the 1-4 loop in amber.

      modex=5
!
      npair=0
      maxwrk=maxexcl
!
      npair4=0
      mxwrk4=maxexcl
!
! compile a list of all the specified interactions
!
      if(modex > 1) then
         do i=1,nbond
            if(ib(i) > jb(i)) then
               i2=jb(i)
               j2=ib(i)
            else
               i2=ib(i)
               j2=jb(i)
            endif
            if(i2 > 0) then
               npair=npair+1
               if(npair > maxexcl)then
                  deallocate(iatbon,natbon,ipk,jpk,ipk14,jpk14)
                  return
               endif
               ipk(npair)=i2
               jpk(npair)=j2
            endif
         enddo
      endif
!
!     next make a list of all the possible 1-3 and 1-4 interactions.
!     this is based on the bond list.  first make a list of all bonds
!     for every atom.
!
      if(modex > 2) then
         do i=1,natom
            natbon(i)=0
         enddo
         do i=1,nbond
            ibt=ib(i)
            jbt=jb(i)
            if(ibt > 0 .and. jbt > 0) then
               natbon(ibt)=natbon(ibt)+1
               if(natbon(ibt) > iatmxb) then
                  write(outu,'(" <make_exclusion_list>: too many bonds for atom",i5)')ibt
                  call mexit(outu,1)
               endif
               iatbon(natbon(ibt),ibt)=i
               natbon(jbt)=natbon(jbt)+1
               if(natbon(jbt) > iatmxb) then
                  write(outu,'(" <make_exclusion_list>: too many bonds for atom",i5)') jbt
                  call mexit(outu,1)
               endif
               iatbon(natbon(jbt),jbt)=-i
            endif
         enddo
!
!     now make the unsorted list of 1-3 interactions by taking bonds
!     and extending in a direction one bond.
!
         do i=1,nbond
            ibt=ib(i)
            jbt=jb(i)
            do j=1,natbon(ibt)
               ij=iatbon(j,ibt)
               if(abs(ij) > i) then
                  if(ij > 0) then
                     ija=jb(ij)
                  else
                     ija=ib(abs(ij))
                  endif
                  if(ija > jbt) then
                     i2=jbt
                     j2=ija
                  else
                     i2=ija
                     j2=jbt
                  endif
                  if (i2 > 0 .and. i2 /= j2) then
                     npair=npair+1
                     if(npair> maxwrk) then
                        !We ran out of space - return and let calling routine reallocate.
                        deallocate(iatbon,natbon,ipk,jpk,ipk14,jpk14)
                        return
                     endif
                     ipk(npair)=i2
                     jpk(npair)=j2
                  endif
               endif
            enddo
            do j=1,natbon(jbt)
               ij=iatbon(j,jbt)
               if(abs(ij) > i) then
                  if(ij > 0) then
                     ija=jb(ij)
                  else
                     ija=ib(abs(ij))
                  endif
                  if(ija > ibt) then
                     i2=ibt
                     j2=ija
                  else
                     i2=ija
                     j2=ibt
                  endif
                  if(i2 > 0 .and. i2 /= j2) then
                     npair=npair+1
                     if(npair > maxwrk) then
                        !We ran out of space - return and let calling routine reallocate.
                        deallocate(iatbon,natbon,ipk,jpk,ipk14,jpk14)
                        return
                     endif
                     ipk(npair)=i2
                     jpk(npair)=j2
                  endif
               endif
            enddo
         enddo
      endif
!
!
!     now make the unsorted list of 1-4 interactions by taking bonds
!     and extending in each direction one bond.
!
      if(modex > 3) then
         do i=1,nbond
            ibt=ib(i)
            jbt=jb(i)
            do j=1,natbon(ibt)
               ij=iatbon(j,ibt)
               if(abs(ij) /= i) then
                  if(ij > 0) then
                     ija=jb(ij)
                  else
                     ija=ib(abs(ij))
                  endif
                  do k=1,natbon(jbt)
                     ik=iatbon(k,jbt)
                     if(abs(ik) /= i) then
                        if(ik > 0) then
                           ika=jb(ik)
                        else
                           ika=ib(abs(ik))
                        endif
                        if(ija > ika) then
                           i2=ika
                           j2=ija
                        else
                           i2=ija
                           j2=ika
                        endif
                        if(i2 > 0 .and. i2 /= j2) then
                           npair=npair+1
                           if(npair > maxwrk) then
                              !We ran out of space - return and let calling routine reallocate.
                              deallocate(iatbon,natbon,ipk,jpk,ipk14,jpk14)
                              return
                           endif
                           ipk(npair)=i2
                           jpk(npair)=j2
                           if(modex == 5) then
                              npair4=npair4+1
                              if(npair4 > mxwrk4) then
                                 !We ran out of space - return and let calling routine reallocate.
                                 deallocate(iatbon,natbon,ipk,jpk,ipk14,jpk14)
                                 return
                              endif
                              ipk14(npair4)=i2
                              jpk14(npair4)=j2
                           endif
                        endif
                     endif
                  enddo
               endif
            enddo
         enddo
      endif
!
!     next sort the list of all the possible 1-4 interactions.
!
      if(modex == 5) then
         npair4=npair4+1
         if(npair4 > mxwrk4) then
            !We ran out of space - return and let calling routine reallocate.
            deallocate(iatbon,natbon,ipk,jpk,ipk14,jpk14)
            return
         endif
         ipk14(npair4)=natom
         jpk14(npair4)=natom
         call sorter(npair4,ipk14,jpk14)
      endif
!
! sort the pair list.
!
      call sorter(npair,ipk,jpk)
!
! process the sorted pair list to make inb. check that there are not
! multiple entries.
!
      nnb14=0
      nx14=0
      i2l=0
      j2l=0
      iatom=1
      next14=1
      lex14=.false.
      do i=1,npair
         i2=ipk(i)
         j2=jpk(i)
         if(modex == 5) then
            lex14=i2 == ipk14(next14) .and. j2 == jpk14(next14)
            if(lex14) next14=next14+1
         endif
!
         ok=.true.
         if(ok) then
           if(i2 == i2l .and. j2 == j2l) then
!             this is a repetition. remove 1-4 sign if applied.
              if(.not.lex14 .and. nonbond_excl_14(nnb14) < 0) then
                 nonbond_excl_14(nnb14)=-nonbond_excl_14(nnb14)
                 nx14=nx14-1
              endif
           else
!             this is a new one.
              i2l=i2
              j2l=j2
              if(i2 > iatom) then
                 do j=iatom,i2-1
                    nonbond_excl_pointers_14(j)=nnb14
                 enddo
                 iatom=i2
              endif
              nnb14=nnb14+1
!
              if(nnb14 > maxexcl) then
                 !We ran out of space - return and let calling routine reallocate.
                 deallocate(iatbon,natbon,ipk,jpk,ipk14,jpk14)
                 return
              endif
!
              if(lex14) then
                 nonbond_excl_14(nnb14)=-j2
                 nx14=nx14+1
              else
                 nonbond_excl_14(nnb14)=j2
              endif
           endif
         endif
      enddo
!
      do j=iatom,natom
         nonbond_excl_pointers_14(j)=nnb14
      enddo

      i=nnb14-nx14
      if(verbose_debug) write(outu,432) modex,i,nx14
 432  format(' <makinb> with modex',i4,' found',i7,' exclusions and',i7, &
          ' interactions(1-4)')
!
      if(verbose_debug) then
         write(outu,433) 'nonbond exclusions pointers'
         write(outu,435) (nonbond_excl_pointers_14(j),j=1,natom)
         write(outu,433) 'nonbond exclusions'
         write(outu,435) (nonbond_excl_14(j),j=1,nnb14)
 433     format(' makinb: ',a)
 435     format(12i8)
      endif
!
      allocate(itmp(natom),jtmp(maxexcl))
      itmp=nonbond_excl_pointers_14(1:natom)
      jtmp=nonbond_excl_14(1:maxexcl)
      lasti=0
      index=0
      do i=1,natom
         do j=lasti+1,nonbond_excl_pointers_14(i)
            index=index+1
            if(index > maxexcl)then
               !We ran out of space - return and let calling routine reallocate.
               deallocate(iatbon,natbon,ipk,jpk,ipk14,jpk14)
               return
            endif
            nonbond_excl_14(index)=abs(jtmp(j))
         enddo
         numexi=nonbond_excl_pointers_14(i)-lasti
         if(numexi == 0)then
            numexi=1
            index=index+1
            nonbond_excl_14(index)=0
         endif
         lasti=nonbond_excl_pointers_14(i)
         nonbond_excl_pointers_14(i)=numexi
      enddo  
      nnb14=index
      nnb=nnb14

      deallocate(itmp,jtmp)

      if(verbose_debug) then
         write(outu,433) 'nonbond exclusions pointers'
         write(outu,435) (nonbond_excl_pointers_14(j),j=1,natom)
         write(outu,433) 'nonbond exclusions'
         write(outu,435) (nonbond_excl_14(j),j=1,nnb14)
      endif

      deallocate(iatbon,natbon,ipk,jpk,ipk14,jpk14)
      done=.true.
      return

end subroutine make_exclusion_list

