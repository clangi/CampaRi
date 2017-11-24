!--------------------------------------------------------------------------!
! LICENSE INFO:                                                            !
!--------------------------------------------------------------------------!
!    This file is part of CAMPARI.                                         !
!                                                                          !
!    Version 3.0                                                           !
!                                                                          !
!    Copyright (C) 2017, The CAMPARI development team (current and former  !
!                        contributors)                                     !
!                        Andreas Vitalis, Adam Steffen, Rohit Pappu, Hoang !
!                        Tran, Albert Mao, Xiaoling Wang, Jose Pulido,     !
!                        Nicholas Lyle, Nicolas Bloechliger, Marco Bacci,  !
!                        Davide Garolini, Jiri Vymetal                     !
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

subroutine particleflucstat(doinst,handle)
  use grandensembles
  use molecule, ONLY: nmol
  implicit none
  logical, INTENT(IN)::  doinst
  integer, INTENT(IN):: handle
  integer cardinality, getmember
  integer i, a, currentnum,insts(fluctypes%nummembers)
  character(10) fmstrng
      
  do i = 1, cardinality(fluctypes)
    a = getmember(fluctypes, i)
    currentnum = cardinality(typesets(a,1))
    insts(i) = currentnum
    numberhistogram(a,currentnum) = numberhistogram(a,currentnum)+1
  end do
  if (doinst.EQV..true.) then
    i = 0
    call int2str(max(1,ceiling(log10(1.0*(nmol+0.5)))),fmstrng,i)
    call strlims(fmstrng,i,a)
    write(handle,'(10000(i'//fmstrng(i:a)//',1x))') insts(:)
  end if
end subroutine particleflucstat
!
!---------------------------------------------------------------------
!
subroutine subvol_numstat(inst,handle)
!
  use molecule
  use paircorr
  use system
  use math
  use grandensembles
  use movesets
  use pdb
  use iounit
  use cutoffs, ONLY: mcel_cutoff 
!
  implicit none
!
  logical, INTENT(IN):: inst
  integer, INTENT(IN):: handle
!
  integer i,j,k,kk,imol,imt,cnts(24*nangrps),nsphs
  RTYPE comi(3),dis2,incr
  RTYPE sysmid(3),sphrad,sphvol,sphr2(3)
  character(10) fmstrng
!
  if (use_frame_weights.EQV..true.) then
    incr = 1.0*framewts(curframe)
  else
    incr = 1.0
  end if

!
! does not work correctly with NPT unless bnd_params stay fixed
  if (bnd_shape.eq.1) then
    nsphs = 8
    sysmid(:) = bnd_params(4:6) + 0.25*bnd_params(1:3)
    sphrad = 1.0*minval(bnd_params(1:3))/4.0
  else if (bnd_shape.eq.2) then
    nsphs = 6
    sysmid(:) = bnd_params(1:3)
    sysmid(1) = sysmid(1) - 0.5875*bnd_params(4)
    sphrad = 0.4125*bnd_params(4)
  else if (bnd_shape.eq.3) then
    if (bnd_params(6).le.0.825*bnd_params(4)) then
      nsphs = 4
      sphrad = bnd_params(6)/2.0
      sysmid(3) = bnd_params(3)
    else if (bnd_params(6).le.2.0*0.825*bnd_params(4)) then
      nsphs = 4
      sphrad = 0.4125*bnd_params(4)
      sysmid(3) = bnd_params(3)
    else if (bnd_params(6).le.6.0*0.825*bnd_params(4)) then
      nsphs = 8
      sphrad = 0.4125*bnd_params(4)
      sysmid(3) = bnd_params(3) - 0.25*bnd_params(6)
    else
      sysmid(1:2) = bnd_params(1:2)
      sphrad = bnd_params(4)
      nsphs = min(8,nint(0.5*bnd_params(6)/sphrad))
      sysmid(3) = bnd_params(3) - 0.5*bnd_params(6) + sphrad
    end if
  else
    write(ilog,*) 'Fatal. Encountered unsupported box shape in subvol_numstat(...). This is an omission bug.'
    call fexit()
  end if
  sphvol = (sphrad**3)*(4./3.)*PI
  sphr2(1) = sphrad
  sphr2(2) = sphrad*0.5
  sphr2(3) = sphrad*0.25
  sphr2(:) = sphr2(:)**2
!
!
 66 format(' For spheres ',i3,' to ',i3,' (columns ',i3,'-',i3,' in PARTICLENUMHIST.dat; ',i3,'-',i3,' in PARTICLENOS.dat), the &
 &ratio of subvolume over total volume is ',g12.5,' the observed particle number fluctuations in an ideal gas are expected to be&
 & a factor of ',g12.5,' too small (finite reservoir), and the ratio of diameter and long-range cutoff is ',g12.5,'.')
  if (allocated(numberhistogram).EQV..false.) then
    allocate(numberhistogram(3*nmoltyp, 0:nmol))
    numberhistogram(:,:) = 0
    write(ilog,'(1x,a,i4,a)') 'Remark. Using ',3*nsphs,' subvolumes in particle number fluctuation analysis.'
    do k=1,3
      write(ilog,66) (k-1)*nsphs+1,k*nsphs,(k-1)*nangrps+1,k*nangrps,(k-1)*nangrps*nsphs+1,k*nangrps*nsphs,&
 &sphvol/(1.0*8**(k-1))/ens%insV,(1.0-sphvol/(1.0*8**(k-1))/ens%insV),2*sphrad/(k*2**(k-1)*mcel_cutoff)
    end do
  end if
!
  if (use_dyn.EQV..true.) then
    if ((dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)) then
      if (in_dyncyc.EQV..false.) then
        do i=1,nsolutes
          imol = solutes(i)
          call update_rigidm(imol)
        end do
      end if
    end if
  end if
!
  cnts(:) = 0
  do j=1,nsphs
    if (bnd_shape.eq.1) then
      if ((j.eq.2).OR.(j.eq.4).OR.(j.eq.6).OR.(j.eq.8)) then
        sysmid(3) = sysmid(3) + 0.5*bnd_params(3)
      else if (j.ne.1) then
        sysmid(3) = sysmid(3) - 0.5*bnd_params(3)
      end if
      if ((j.eq.3).OR.(j.eq.7)) then
        sysmid(2) = sysmid(2) + 0.5*bnd_params(2)
      else if (j.eq.5) then
        sysmid(2) = sysmid(2) - 0.5*bnd_params(2)
      end if
      if (j.eq.5) sysmid(1) = sysmid(1) + 0.5*bnd_params(1)
    else if (bnd_shape.eq.2) then
      sysmid(:) = bnd_params(1:3)
      if (mod(j,2).eq.0) then
        sysmid(ceiling((j-0.1)/2.0)) = bnd_params(ceiling((j-0.1)/2.0)) - 0.5875*bnd_params(4)
      else
        sysmid(ceiling((j-0.1)/2.0)) = bnd_params(ceiling((j-0.1)/2.0)) + 0.5875*bnd_params(4)
      end if
    else if (bnd_shape.eq.3) then
      if (bnd_params(6).le.6.0*0.825*bnd_params(4)) then ! 4 or 8 in a z-plane
        sysmid(1:2) = bnd_params(1:2)
        if (mod(j,2).eq.0) then
          sysmid(ceiling((mod(j-1,4)+0.9)/2.0)) = bnd_params(ceiling((mod(j-1,4)+0.9)/2.0)) - bnd_params(4) + sphrad
        else
          sysmid(ceiling((mod(j-1,4)+0.9)/2.0)) = bnd_params(ceiling((mod(j-1,4)+0.9)/2.0)) + bnd_params(4) - sphrad
        end if
        if (j.eq.5) sysmid(3) = sysmid(3) + bnd_params(6)/2.0
      else
        if (j.gt.1) sysmid(3) = sysmid(3) + bnd_params(6)/(1.0*nsphs)
      end if
    end if
    do imol=1,nmol
!    imol = solutes(i)
      imt = an_grp_mol(imol)
      if (use_dyn.EQV..true.) then
        comi(:) = comm(imol,:)
      else
        comi(:) = com(imol,:)
      end if
      call dis_bound2(comi,sysmid,dis2)
      do k=1,3
        if (dis2.le.sphr2(k)) then
          kk = (j-1)*nangrps+imt + (k-1)*nsphs*nangrps
          cnts(kk) = cnts(kk) + 1
        end if
      end do
    end do
  end do
  k = nsphs*nangrps
  do i=1,nangrps
    do j=1,nsphs
      numberhistogram(i,cnts((j-1)*nangrps+i)) = numberhistogram(i,cnts((j-1)*nangrps+i)) + incr
    end do
    do j=nsphs+1,2*nsphs
      numberhistogram(nangrps+i,cnts((j-1)*nangrps+i)) = numberhistogram(nangrps+i,cnts((j-1)*nangrps+i)) + incr
    end do
    do j=2*nsphs+1,3*nsphs
      numberhistogram(2*nangrps+i,cnts((j-1)*nangrps+i)) = numberhistogram(2*nangrps+i,cnts((j-1)*nangrps+i)) + incr
    end do
  end do
  if (inst.EQV..true.) then
    i = 0
    call int2str(max(1,ceiling(log10(1.0*(nmol+0.5)))),fmstrng,i)
    call strlims(fmstrng,i,j)
    write(handle,'(10000(i'//fmstrng(i:j)//',1x))') cnts(1:(nangrps*3*nsphs))
  end if
!
end subroutine subvol_numstat
!
!------------------------------------------------------------------------------------------------------
!
subroutine prt_particlenumhistogram()
!
  use grandensembles
  use molecule
  use mpistuff
  use system, ONLY: ens
!
  implicit none
  integer freeunit, cardinality, getmember
  integer i,a,c,iu,lastbin,maxnum,nft,ii,jj
#ifdef ENABLE_MPI
  character(re_aux(10)) nod
#endif
  character(100) fn
  logical exists, foundnonzero
 
 98   format('#      N',10000(i8))
 99   format(10001(i8))

#ifdef ENABLE_MPI
  if (use_REMC.EQV..true.) then
    call int2str(myrank,nod,re_aux(10))
    fn = 'N_'//nod(1:re_aux(10))//'_PARTICLENUMHIST.dat'
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
!
  maxnum = 0
  lastbin = nmol
  do i = 1,nmoltyp
    foundnonzero = .false.
    do while (foundnonzero.EQV..false.)
      if (maxval(numberhistogram(:,lastbin)).gt.0) then
        foundnonzero = .true.
      else
        lastbin = lastbin - 1
      end if
    end do
    maxnum = lastbin
  end do
  if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
    nft = cardinality(fluctypes) 
    write(iu,98) (getmember(fluctypes,i), i=1,nft)
!   Note: numberhistogram was allocated so that the second index is 0-based.
    do c = 0, maxnum
      write(iu,99) c,(numberhistogram(getmember(fluctypes,i),c), i=1,nft)
    end do
  else
    write(iu,98) (i, i=1,3*nangrps)
!   Note: numberhistogram was allocated so that the second index is 0-based.
    do c = 0, maxnum
      write(iu,99) c,(numberhistogram(i,c), i=1,3*nangrps)
    end do
  end if
  close(unit=iu)
end subroutine prt_particlenumhistogram
