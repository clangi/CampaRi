!--------------------------------------------------------------------------!
! LICENSE INFO:                                                            !
!--------------------------------------------------------------------------!
!    This file is part of CAMPARI.                                         !
!                                                                          !
!    Version 2.0                                                           !
!                                                                          !
!    Copyright (C) 2014, The CAMPARI development team (current and former  !
!                        contributors)                                     !
!                        Andreas Vitalis, Adam Steffen, Rohit Pappu, Hoang !
!                        Tran, Albert Mao, Xiaoling Wang, Jose Pulido,     !
!                        Nicholas Lyle, Nicolas Bloechliger                !
!                                                                          !
!    Website: http://sourceforge.net/projects/campari/                     !
!                                                                          !
!    CAMPARI is free software: you can redistribute it and/or modify       !
!    it under the terms of the GNU General Public License as published by  !
!    the Free Software Foundation, either version 3 of the License, or     !
!    (at your option) any later version.                                   !
!                                                                          !
!    CAMPARI is distributed in the hope that it will be useful,            !
!    but WITHOUT ANY WARRANTY; without even the implied warranty of        !
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         !
!    GNU General Public License for more details.                          !
!                                                                          !
!    You should have received a copy of the GNU General Public License     !
!    along with CAMPARI.  If not, see <http://www.gnu.org/licenses/>.      !
!--------------------------------------------------------------------------!
! AUTHORSHIP INFO:                                                         !
!--------------------------------------------------------------------------!
!                                                                          !
! MAIN AUTHOR:   Albert Mao                                                !
! CONTRIBUTIONS: Andreas Vitalis                                           !
!                                                                          !
!--------------------------------------------------------------------------!
#include "macros.i"
      
  ! These subroutines and functions and the type definition from
  ! grandensembles.i are supposed to be grouped together into
  ! one module that encapsulates the concept of an integer set.
  ! Unfortunately, in order to fit within the framework Andreas
  ! designed for modules in CAMPARI, the derived type definition
  ! and the methods that operate on it must be separated because
  ! modules cannot use other modules.
subroutine constructset(set, universe)
    use grandensembles
    implicit none
    type(integer_set), intent(inout) :: set
    integer, intent(in) :: universe
    integer :: i
      
    set%universesize = universe
    set%nummembers = 0
    allocate (set%indexof(universe), set%array(universe))
    do i = 1, universe
      set%indexof(i) = i
      set%array(i) = i
    end do
end subroutine constructset
      
subroutine destroyset(set)
    use grandensembles
    implicit none
    type(integer_set), intent(inout) :: set
    
    set%universesize = -1
    set%nummembers = -1
    if (allocated(set%indexof).eqv..true.) then
      deallocate(set%indexof)
    end if
    if (allocated(set%array).eqv..true.) then
      deallocate (set%array)
    end if
end subroutine destroyset
      
subroutine addtoset(set, element)
    use grandensembles
    implicit none
    type(integer_set), intent(inout) :: set
    integer, intent(in) :: element
    integer :: oldelementindex, tomove
    
    if ((element.lt.1) .or. (element.gt.set%universesize)) return
    if (set%indexof(element).le.set%nummembers) return
    
    ! At this point, the element is valid and is not already in the
    ! set
    set%nummembers = set%nummembers + 1
    oldelementindex = set%indexof(element)
    tomove = set%array(set%nummembers)
    set%array(oldelementindex) = tomove
    set%indexof(tomove) = oldelementindex
    set%array(set%nummembers) = element
    set%indexof(element) = set%nummembers
end subroutine addtoset
      
subroutine removefromset(set, element)
    use grandensembles
    implicit none
    type(integer_set), intent(inout) :: set
    integer, intent(in) :: element
    integer :: oldelementindex, tomove
    
    if ((element.lt.1) .or. (element.gt.set%universesize)) return
    if (set%indexof(element).gt.set%nummembers) return
    
    ! At this point, the element is valid and is already in the set
    oldelementindex = set%indexof(element)
    tomove = set%array(set%nummembers)
    set%array(oldelementindex) = tomove
    set%indexof(tomove) = oldelementindex
    set%array(set%nummembers) = element
    set%indexof(element) = set%nummembers
    set%nummembers = set%nummembers - 1
end subroutine removefromset
      
subroutine clearset(set)
    use grandensembles
    implicit none
    type(integer_set), intent(inout) :: set
    
    set%nummembers = 0
end subroutine clearset
    
  pure integer function universecardinality(set)
    use grandensembles
    implicit none
    type(integer_set), intent(in) :: set
    
    universecardinality = set%universesize
end function universecardinality
      
  pure integer function cardinality(set)
    use grandensembles
    implicit none
    type(integer_set), intent(in) :: set
    
    cardinality = set%nummembers
end function cardinality
      
  pure logical function ismember(set, element)
    use grandensembles
    implicit none
    type(integer_set), intent(in) :: set
    integer, intent(in) :: element
    
    ismember = .false.
    if (element.lt.1 .or. element.gt.set%universesize) return
    ismember = (set%indexof(element).le.set%nummembers)
end function ismember
      
  pure logical function isemptyset(set)
    use grandensembles
    implicit none
    type(integer_set), intent(in) :: set
    
    isemptyset = (set%nummembers.eq.0)
end function isemptyset
    
  pure integer function getmember(set, i)
    use grandensembles
    implicit none
    type(integer_set), intent(in) :: set
    integer, intent(in) :: i
    
    getmember = -1
    if ((i.ge.1) .and. (i.le.set%nummembers)) then
      getmember = set%array(i)
    end if
end function getmember
      
subroutine printset(set)
    use grandensembles
    implicit none
    type(integer_set), intent(in) :: set
    
    write(*,*) set%nummembers,'/',set%universesize
end subroutine printset
  ! End of functions from integerset
      
  integer function randomfluctype()
    use grandensembles
    implicit none
    RTYPE random
    integer cardinality, getmember
    integer i
    
    i = floor(random() * cardinality(fluctypes)) + 1
    randomfluctype = getmember(fluctypes, i)
end function randomfluctype

 integer function random_wt_fluctype()
    use grandensembles
    use interfaces
    use movesets
    implicit none
    RTYPE random,rcheck
    integer i
    
    rcheck = random()
    call binary_search(fluclst%nr,fluclst%wt(:),rcheck,i)
    random_wt_fluctype = fluclst%idx(min(fluclst%nr,i+1))
end function random_wt_fluctype 

  integer function randommoloftype(typ, presence)
    use grandensembles
    implicit none
    integer, intent(in) :: typ
!    logical, intent(in) :: presence
    integer, intent(in) :: presence
    RTYPE random
    integer cardinality, getmember
    integer i

    i = floor(random()*cardinality(typesets(typ, presence))) + 1
    randommoloftype = getmember(typesets(typ, presence), i)
end function randommoloftype
      
subroutine insertmol(imol)
    use grandensembles
    use molecule
    implicit none
    integer, intent(in) :: imol
    integer typ
    
    typ = moltypid(imol)
    call addtoset(ispresent, imol)
    call addtoset(typesets(typ,1), imol)
    call removefromset(typesets(typ,2), imol)
end subroutine insertmol
      
subroutine deletemol(imol)
    use grandensembles
    use molecule
    implicit none
    integer, intent(in) :: imol
    integer typ
    
    typ = moltypid(imol)
    call removefromset(ispresent, imol)
    call addtoset(typesets(typ,2), imol)
    call removefromset(typesets(typ,1), imol)
end subroutine deletemol

subroutine transmute(frommol, tomol)
    use grandensembles
    use molecule
    use atoms
    use cutoffs
    implicit none
    integer, intent(in) :: frommol
    integer, intent(in) :: tomol
    RTYPE tvec(3)
    integer k
    
    call deletemol(frommol)
    call insertmol(tomol)
    ! Creates a translation vector pointing from frommol to tomol
    tvec(:) = com(tomol,:) - com(frommol,:)
    ! Swaps the centers of mass of the molecules
    com(frommol,:) = com(frommol,:) + tvec(:)
    com(tomol,:) = com(tomol,:) - tvec(:)
    do k = atmol(frommol,1), atmol(frommol,2)
      x(k) = x(k) + tvec(1)
      y(k) = y(k) + tvec(2)
      z(k) = z(k) + tvec(3)
    end do
    do k = atmol(tomol,1), atmol(tomol,2)
      x(k) = x(k) - tvec(1)
      y(k) = y(k) - tvec(2)
      z(k) = z(k) - tvec(3)
    end do
    ! Updates residue grid points
    if (use_mcgrid) then
      do k = rsmol(frommol,1), rsmol(frommol,2)
        call updateresgp(k)
      end do
      do k = rsmol(tomol,1), rsmol(tomol,2)
        call updateresgp(k)
      end do
    end if
end subroutine transmute

subroutine particleflucstat()
  use grandensembles
  implicit none
  integer cardinality, getmember
  integer i, a, currentnum
      
  do i = 1, cardinality(fluctypes)
    a = getmember(fluctypes, i)
    currentnum = cardinality(typesets(a,1))
    numberhistogram(a,currentnum) = numberhistogram(a,currentnum)+1
  end do
end subroutine particleflucstat

subroutine prt_particlenumhistogram()
!
  use grandensembles
  use molecule
  use mpistuff
!
  implicit none
  integer freeunit, cardinality, getmember
  integer i,a,c,iu,lastbin,maxnum,nft,ii,jj
#ifdef ENABLE_MPI
  integer tl
  character(3) nod
#endif
  character(100) fn
  logical exists, foundnonzero
 
 98   format('#',i7,9999(i8))
 99   format(1000000(i8))

#ifdef ENABLE_MPI
  if (use_REMC.EQV..true.) then
    tl = 3
    call int2str(myrank,nod,tl)
    fn = 'N_'//nod(1:tl)//'_PARTICLENUMHIST.dat'
  else if (use_MPIAVG.EQV..true.) then
    fn = 'PARTICLENUMHIST.dat'
  end if
#else
  fn = 'PARTICLENUMHIST.dat'
#endif
  call strlims(fn,ii,jj)
  inquire(file=fn(ii:jj), exist=exists)
  if(exists) then
    iu = freeunit()
    open(unit=iu,file=fn(ii:jj),status='old')
    close(unit=iu,status='delete')
  end if
  iu=freeunit()
  open(unit=iu,file=fn(ii:jj),status='new')
      
  nft = cardinality(fluctypes)
  maxnum = 0 
  do i = 1, nft
    a = getmember(fluctypes, i)
    lastbin = nmol
    foundnonzero = .false.
    do while ((lastbin.gt.maxnum) .and. .not. foundnonzero)
      if (numberhistogram(a,lastbin).gt.0) then
        foundnonzero = .true.
      else
        lastbin = lastbin - 1
      end if
    end do
    maxnum = max(lastbin, maxnum)
  end do
  write(iu,98) (getmember(fluctypes,i), i=1,nft)
  ! Note: numberhistogram was allocated so that the second index is 0-based.
  do c = 0, maxnum
    write(iu,99) (numberhistogram(getmember(fluctypes,i),c), i=1,n&
 &ft)
  end do
  close(unit=iu)
end subroutine prt_particlenumhistogram
