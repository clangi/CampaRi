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
! MAIN AUTHOR:   Andreas Vitalis                                           !
! CONTRIBUTIONS: Rohit Pappu, Nicholas Lyle                                !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
!-----------------------------------------------------------------------
!
!
!             ##################################
!             #                                #
!             #  A SET OF GENERAL TOOLS        #
!             #  (DISORGANIZED TOOLBOX)        #
!             #                                #
!             ##################################
!
!
!-----------------------------------------------------------------------
!
! standard bubble sort of an integer array of size ldim in ascending order
! (O(N^2) average and WCS)
! WARNING: do NOT use this sort algorithm for anything rate-limiting
!
subroutine sort(ldim,list)
!
  implicit none
!
  integer i,ii,j,dummy,vali,valj,ldim
  integer list(*)
!
  do i=1,ldim
    ii = i
    vali = list(ii)
    do j=i+1,ldim
      valj = list(j)
      if (valj.lt.vali) then
        ii = j
        vali = list(ii)
      end if
    end do
    if (ii.ne.i) then
      dummy = list(i)
      list(i) = list(ii)
      list(ii) = dummy
    end if
  end do
!
end
!
!--------------------------------------------------------------------------
!
! merge sort (O(NlogN) average and WCS) 
!
recursive subroutine merge_isort(ldim,up,list,olist,ilo,ihi,idxmap,olist2)
!
  use iounit
!
  implicit none
!
  integer i,j,ldim,naux1,naux2,kn1,kn2,ilo,ihi,ilo2,ihi2
  integer list(*),olist(*)
  integer, OPTIONAL:: idxmap(*),olist2(*)
  logical up
!
  if (ihi.lt.ilo) then
    write(ilog,*) 'Fatal. Called merge_sort with bad settings for ilo and ihi. This is a bug.'
    call fexit()
  end if
  if (ihi.eq.ilo) then
    olist(ilo) = list(ilo)
    if ((present(idxmap).EQV..true.).AND.(present(olist2).EQV..true.)) then
      olist2(ilo) = idxmap(ilo)
    end if
    return
  end if
!
  ihi2 = (ihi-ilo+1)/2 + ilo - 1
  ilo2 = ihi2 + 1
  naux1 = ihi2 - ilo + 1
  naux2 = ihi - ilo2 + 1
!
  if (present(idxmap).EQV..true.) then
    if (present(olist2).EQV..false.) then
      write(ilog,*) 'Fatal. Called merge-sort variant with optional index array but &
 &without then-mandatory output index array. This is a bug.'
      call fexit()
    end if
    call merge_isort(ldim,up,list,olist,ilo,ihi2,idxmap=idxmap,olist2=olist2)
    if (naux1.gt.1) then
      list(ilo:ihi2) = olist(1:naux1)
      idxmap(ilo:ihi2) = olist2(1:naux1)
    end if
    call merge_isort(ldim,up,list,olist,ilo2,ihi,idxmap=idxmap,olist2=olist2)
    if (naux2.gt.1) then
      list(ilo2:ihi) = olist(1:naux2)
      idxmap(ilo2:ihi) = olist2(1:naux2)
    end if
  else
    if (present(olist2).EQV..true.) then
      write(ilog,*) 'Fatal. Called merge-sort variant with optional second output array but &
 &without then-mandatory input index array. This is a bug.'
      call fexit()
    end if
    call merge_isort(ldim,up,list,olist,ilo,ihi2)
    if (naux1.gt.1) list(ilo:ihi2) = olist(1:naux1)
    call merge_isort(ldim,up,list,olist,ilo2,ihi)
    if (naux2.gt.1) list(ilo2:ihi) = olist(1:naux2)
  end if

!
  kn1 = ilo
  kn2 = ilo2
  i = 0
  if ((up.EQV..true.).AND.(present(idxmap).EQV..true.)) then
    do while ((kn1.le.ihi2).AND.(kn2.le.ihi))
      if (list(kn1).le.list(kn2)) then
        i = i + 1
        olist(i) = list(kn1)
        olist2(i) = idxmap(kn1)
        kn1 = kn1 + 1
      else
        i = i + 1
        olist(i) = list(kn2)
        olist2(i) = idxmap(kn2)
        kn2 = kn2 + 1
      end if
    end do
    if (kn1.le.ihi2) then
      j = i + 1
      i = j + ihi2 - kn1 
      olist(j:i) = list(kn1:ihi2)
      olist2(j:i) = idxmap(kn1:ihi2)
    else
      j = i + 1
      i = j + ihi - kn2
      olist(j:i) = list(kn2:ihi)
      olist2(j:i) = idxmap(kn2:ihi)
    end if
  else if ((up.EQV..true.).AND.(present(idxmap).EQV..false.)) then
    do while ((kn1.le.ihi2).AND.(kn2.le.ihi))
      if (list(kn1).le.list(kn2)) then
        i = i + 1
        olist(i) = list(kn1)
        kn1 = kn1 + 1
      else
        i = i + 1
        olist(i) = list(kn2)
        kn2 = kn2 + 1
      end if
    end do
    if (kn1.le.ihi2) then
      j = i + 1
      i = j + ihi2 - kn1 
      olist(j:i) = list(kn1:ihi2)
    else
      j = i + 1
      i = j + ihi - kn2
      olist(j:i) = list(kn2:ihi)
    end if
  else if ((up.EQV..false.).AND.(present(idxmap).EQV..true.)) then
    do while ((kn1.le.ihi2).AND.(kn2.le.ihi))
      if (list(kn1).ge.list(kn2)) then
        i = i + 1
        olist(i) = list(kn1)
        olist2(i) = idxmap(kn1)
        kn1 = kn1 + 1
      else
        i = i + 1
        olist(i) = list(kn2)
        olist2(i) = idxmap(kn2)
        kn2 = kn2 + 1
      end if
    end do
    if (kn1.le.ihi2) then
      j = i + 1
      i = j + ihi2 - kn1 
      olist(j:i) = list(kn1:ihi2)
      olist2(j:i) = idxmap(kn1:ihi2)
    else
      j = i + 1
      i = j + ihi - kn2
      olist(j:i) = list(kn2:ihi)
      olist2(j:i) = idxmap(kn2:ihi)
    end if
  else
    do while ((kn1.le.ihi2).AND.(kn2.le.ihi))
      if (list(kn1).ge.list(kn2)) then
        i = i + 1
        olist(i) = list(kn1)
        kn1 = kn1 + 1
      else
        i = i + 1
        olist(i) = list(kn2)
        kn2 = kn2 + 1
      end if
    end do
    if (kn1.le.ihi2) then
      j = i + 1
      i = j + ihi2 - kn1
      olist(j:i) = list(kn1:ihi2)
    else
      j = i + 1
      i = j + ihi - kn2
      olist(j:i) = list(kn2:ihi)
    end if
  end if
!
  return
!
end
!
!--------------------------------------------------------------------------
!
! merge sort (O(NlogN) average and WCS) 
!
recursive subroutine merge_fsort(ldim,up,list,olist,ilo,ihi,idxmap,olist2)
!
  use iounit
!
  implicit none
!
  integer i,j,ldim,naux1,naux2,kn1,kn2,ilo,ihi,ilo2,ihi2
  RTYPE list(*),olist(*)
  integer, OPTIONAL:: idxmap(*),olist2(*)
  logical up
!
  if (ihi.lt.ilo) then
    write(ilog,*) 'Fatal. Called merge_sort with bad settings for ilo and ihi. This is a bug.'
    call fexit()
  end if
  if (ihi.eq.ilo) then
    olist(ilo) = list(ilo)
    if ((present(idxmap).EQV..true.).AND.(present(olist2).EQV..true.)) then
      olist2(ilo) = idxmap(ilo)
    end if
    return
  end if
!
  ihi2 = (ihi-ilo+1)/2 + ilo - 1
  ilo2 = ihi2 + 1
  naux1 = ihi2 - ilo + 1
  naux2 = ihi - ilo2 + 1
!
  if (present(idxmap).EQV..true.) then
    if (present(olist2).EQV..false.) then
      write(ilog,*) 'Fatal. Called merge-sort variant with optional index array but &
 &without then-mandatory output index array. This is a bug.'
      call fexit()
    end if
    call merge_fsort(ldim,up,list,olist,ilo,ihi2,idxmap=idxmap,olist2=olist2)
    if (naux1.gt.1) then
      list(ilo:ihi2) = olist(1:naux1)
      idxmap(ilo:ihi2) = olist2(1:naux1)
    end if
    call merge_fsort(ldim,up,list,olist,ilo2,ihi,idxmap=idxmap,olist2=olist2)
    if (naux2.gt.1) then
      list(ilo2:ihi) = olist(1:naux2)
      idxmap(ilo2:ihi) = olist2(1:naux2)
    end if
  else
    if (present(olist2).EQV..true.) then
      write(ilog,*) 'Fatal. Called merge-sort variant with optional second output array but &
 &without then-mandatory input index array. This is a bug.'
      call fexit()
    end if
    call merge_fsort(ldim,up,list,olist,ilo,ihi2)
    if (naux1.gt.1) list(ilo:ihi2) = olist(1:naux1)
    call merge_fsort(ldim,up,list,olist,ilo2,ihi)
    if (naux2.gt.1) list(ilo2:ihi) = olist(1:naux2)
  end if

!
  kn1 = ilo
  kn2 = ilo2
  i = 0
  if ((up.EQV..true.).AND.(present(idxmap).EQV..true.)) then
    do while ((kn1.le.ihi2).AND.(kn2.le.ihi))
      if (list(kn1).le.list(kn2)) then
        i = i + 1
        olist(i) = list(kn1)
        olist2(i) = idxmap(kn1)
        kn1 = kn1 + 1
      else
        i = i + 1
        olist(i) = list(kn2)
        olist2(i) = idxmap(kn2)
        kn2 = kn2 + 1
      end if
    end do
    if (kn1.le.ihi2) then
      j = i + 1
      i = j + ihi2 - kn1 
      olist(j:i) = list(kn1:ihi2)
      olist2(j:i) = idxmap(kn1:ihi2)
    else
      j = i + 1
      i = j + ihi - kn2
      olist(j:i) = list(kn2:ihi)
      olist2(j:i) = idxmap(kn2:ihi)
    end if
  else if ((up.EQV..true.).AND.(present(idxmap).EQV..false.)) then
    do while ((kn1.le.ihi2).AND.(kn2.le.ihi))
      if (list(kn1).le.list(kn2)) then
        i = i + 1
        olist(i) = list(kn1)
        kn1 = kn1 + 1
      else
        i = i + 1
        olist(i) = list(kn2)
        kn2 = kn2 + 1
      end if
    end do
    if (kn1.le.ihi2) then
      j = i + 1
      i = j + ihi2 - kn1 
      olist(j:i) = list(kn1:ihi2)
    else
      j = i + 1
      i = j + ihi - kn2
      olist(j:i) = list(kn2:ihi)
    end if
  else if ((up.EQV..false.).AND.(present(idxmap).EQV..true.)) then
    do while ((kn1.le.ihi2).AND.(kn2.le.ihi))
      if (list(kn1).ge.list(kn2)) then
        i = i + 1
        olist(i) = list(kn1)
        olist2(i) = idxmap(kn1)
        kn1 = kn1 + 1
      else
        i = i + 1
        olist(i) = list(kn2)
        olist2(i) = idxmap(kn2)
        kn2 = kn2 + 1
      end if
    end do
    if (kn1.le.ihi2) then
      j = i + 1
      i = j + ihi2 - kn1 
      olist(j:i) = list(kn1:ihi2)
      olist2(j:i) = idxmap(kn1:ihi2)
    else
      j = i + 1
      i = j + ihi - kn2
      olist(j:i) = list(kn2:ihi)
      olist2(j:i) = idxmap(kn2:ihi)
    end if
  else
    do while ((kn1.le.ihi2).AND.(kn2.le.ihi))
      if (list(kn1).ge.list(kn2)) then
        i = i + 1
        olist(i) = list(kn1)
        kn1 = kn1 + 1
      else
        i = i + 1
        olist(i) = list(kn2)
        kn2 = kn2 + 1
      end if
    end do
    if (kn1.le.ihi2) then
      j = i + 1
      i = j + ihi2 - kn1
      olist(j:i) = list(kn1:ihi2)
    else
      j = i + 1
      i = j + ihi - kn2
      olist(j:i) = list(kn2:ihi)
    end if
  end if
!
  return
!
end
!
!--------------------------------------------------------------------------
!
! merge sort (O(NlogN) average and WCS) 
!
recursive subroutine merge_ffsort(ldim,up,list,olist,ilo,ihi,idxmap,olist2)
!
  use iounit
!
  implicit none
!
  integer i,j,ldim,naux1,naux2,kn1,kn2,ilo,ihi,ilo2,ihi2
  real(KIND=4) list(*),olist(*)
  integer, OPTIONAL:: idxmap(*),olist2(*)
  logical up
!
  if (ihi.lt.ilo) then
    write(ilog,*) 'Fatal. Called merge_sort with bad settings for ilo and ihi. This is a bug.'
    call fexit()
  end if
  if (ihi.eq.ilo) then
    olist(ilo) = list(ilo)
    if ((present(idxmap).EQV..true.).AND.(present(olist2).EQV..true.)) then
      olist2(ilo) = idxmap(ilo)
    end if
    return
  end if
!
  ihi2 = (ihi-ilo+1)/2 + ilo - 1
  ilo2 = ihi2 + 1
  naux1 = ihi2 - ilo + 1
  naux2 = ihi - ilo2 + 1
!
  if (present(idxmap).EQV..true.) then
    if (present(olist2).EQV..false.) then
      write(ilog,*) 'Fatal. Called merge-sort variant with optional index array but &
 &without then-mandatory output index array. This is a bug.'
      call fexit()
    end if
    call merge_ffsort(ldim,up,list,olist,ilo,ihi2,idxmap=idxmap,olist2=olist2)
    if (naux1.gt.1) then
      list(ilo:ihi2) = olist(1:naux1)
      idxmap(ilo:ihi2) = olist2(1:naux1)
    end if
    call merge_ffsort(ldim,up,list,olist,ilo2,ihi,idxmap=idxmap,olist2=olist2)
    if (naux2.gt.1) then
      list(ilo2:ihi) = olist(1:naux2)
      idxmap(ilo2:ihi) = olist2(1:naux2)
    end if
  else
    if (present(olist2).EQV..true.) then
      write(ilog,*) 'Fatal. Called merge-sort variant with optional second output array but &
 &without then-mandatory input index array. This is a bug.'
      call fexit()
    end if
    call merge_ffsort(ldim,up,list,olist,ilo,ihi2)
    if (naux1.gt.1) list(ilo:ihi2) = olist(1:naux1)
    call merge_ffsort(ldim,up,list,olist,ilo2,ihi)
    if (naux2.gt.1) list(ilo2:ihi) = olist(1:naux2)
  end if

!
  kn1 = ilo
  kn2 = ilo2
  i = 0
  if ((up.EQV..true.).AND.(present(idxmap).EQV..true.)) then
    do while ((kn1.le.ihi2).AND.(kn2.le.ihi))
      if (list(kn1).le.list(kn2)) then
        i = i + 1
        olist(i) = list(kn1)
        olist2(i) = idxmap(kn1)
        kn1 = kn1 + 1
      else
        i = i + 1
        olist(i) = list(kn2)
        olist2(i) = idxmap(kn2)
        kn2 = kn2 + 1
      end if
    end do
    if (kn1.le.ihi2) then
      j = i + 1
      i = j + ihi2 - kn1 
      olist(j:i) = list(kn1:ihi2)
      olist2(j:i) = idxmap(kn1:ihi2)
    else
      j = i + 1
      i = j + ihi - kn2
      olist(j:i) = list(kn2:ihi)
      olist2(j:i) = idxmap(kn2:ihi)
    end if
  else if ((up.EQV..true.).AND.(present(idxmap).EQV..false.)) then
    do while ((kn1.le.ihi2).AND.(kn2.le.ihi))
      if (list(kn1).le.list(kn2)) then
        i = i + 1
        olist(i) = list(kn1)
        kn1 = kn1 + 1
      else
        i = i + 1
        olist(i) = list(kn2)
        kn2 = kn2 + 1
      end if
    end do
    if (kn1.le.ihi2) then
      j = i + 1
      i = j + ihi2 - kn1 
      olist(j:i) = list(kn1:ihi2)
    else
      j = i + 1
      i = j + ihi - kn2
      olist(j:i) = list(kn2:ihi)
    end if
  else if ((up.EQV..false.).AND.(present(idxmap).EQV..true.)) then
    do while ((kn1.le.ihi2).AND.(kn2.le.ihi))
      if (list(kn1).ge.list(kn2)) then
        i = i + 1
        olist(i) = list(kn1)
        olist2(i) = idxmap(kn1)
        kn1 = kn1 + 1
      else
        i = i + 1
        olist(i) = list(kn2)
        olist2(i) = idxmap(kn2)
        kn2 = kn2 + 1
      end if
    end do
    if (kn1.le.ihi2) then
      j = i + 1
      i = j + ihi2 - kn1 
      olist(j:i) = list(kn1:ihi2)
      olist2(j:i) = idxmap(kn1:ihi2)
    else
      j = i + 1
      i = j + ihi - kn2
      olist(j:i) = list(kn2:ihi)
      olist2(j:i) = idxmap(kn2:ihi)
    end if
  else
    do while ((kn1.le.ihi2).AND.(kn2.le.ihi))
      if (list(kn1).ge.list(kn2)) then
        i = i + 1
        olist(i) = list(kn1)
        kn1 = kn1 + 1
      else
        i = i + 1
        olist(i) = list(kn2)
        kn2 = kn2 + 1
      end if
    end do
    if (kn1.le.ihi2) then
      j = i + 1
      i = j + ihi2 - kn1
      olist(j:i) = list(kn1:ihi2)
    else
      j = i + 1
      i = j + ihi - kn2
      olist(j:i) = list(kn2:ihi)
    end if
  end if
!
  return
!
end
!
!-----------------------------------------------------------------------
!
! sort according to an integer with up to three optional float-list 
! and two additional integer list arguments to be sorted by the same
! index
!
subroutine isort(ldim,up,ilist,flist,ilist2,ilist3,flist2,flist3)
!
  implicit none
!
  integer i,ii,j,ldim,dummi,ilist(*),vali,valj
  RTYPE dummy
  integer, OPTIONAL:: ilist2(*),ilist3(*)
  RTYPE, OPTIONAL:: flist(*),flist2(*),flist3(*)
  logical up
!
  do i=1,ldim
    ii = i
    vali = ilist(ii)
    do j=i+1,ldim
      valj = ilist(j)
      if (up.EQV..true.) then
        if (valj.lt.vali) then
          ii = j
          vali = ilist(ii)
        end if
      else
        if (valj.gt.vali) then
          ii = j
          vali = ilist(ii)
        end if
      end if
    end do
    if (ii.ne.i) then
      dummi = ilist(i)
      ilist(i) = ilist(ii)
      ilist(ii) = dummi
      if (present(ilist2).EQV..true.) then
        dummi = ilist2(i)
        ilist2(i) = ilist2(ii)
        ilist2(ii) = dummi
      end if
      if (present(ilist3).EQV..true.) then
        dummi = ilist3(i)
        ilist3(i) = ilist3(ii)
        ilist3(ii) = dummi
      end if
      if (present(flist).EQV..true.) then
        dummy = flist(i)
        flist(i) = flist(ii)
        flist(ii) = dummy
      end if
      if (present(flist2).EQV..true.) then
        dummy = flist2(i)
        flist2(i) = flist2(ii)
        flist2(ii) = dummy
      end if
      if (present(flist3).EQV..true.) then
        dummy = flist3(i)
        flist3(i) = flist3(ii)
        flist3(ii) = dummy
      end if
    end if
  end do
!
end
!
!-------------------------------------------------------------------
!
! sort according to a float with up to three optional integer list 
! and two additional float-list arguments to be sorted by the same
! index
!
subroutine fsort(ldim,up,flist,ilist,ilist2,ilist3,flist2,flist3)
!
  implicit none
!
  integer i,ii,j,ldim,dummi
  RTYPE flist(*),dummy,vali,valj
  integer, OPTIONAL:: ilist(*),ilist2(*),ilist3(*)
  RTYPE, OPTIONAL:: flist2(*),flist3(*)
  logical up
!
  do i=1,ldim
    ii = i
    vali = flist(ii)
    do j=i+1,ldim
      valj = flist(j)
      if (up.EQV..true.) then
        if (valj.lt.vali) then
          ii = j
          vali = flist(ii)
        end if
      else
        if (valj.gt.vali) then
          ii = j
          vali = flist(ii)
        end if
      end if
    end do
    if (ii.ne.i) then
      dummy = flist(i)
      flist(i) = flist(ii)
      flist(ii) = dummy
      if (present(ilist).EQV..true.) then
        dummi = ilist(i)
        ilist(i) = ilist(ii)
        ilist(ii) = dummi
      end if
      if (present(ilist2).EQV..true.) then
        dummi = ilist2(i)
        ilist2(i) = ilist2(ii)
        ilist2(ii) = dummi
      end if
      if (present(ilist3).EQV..true.) then
        dummi = ilist3(i)
        ilist3(i) = ilist3(ii)
        ilist3(ii) = dummi
      end if
      if (present(flist2).EQV..true.) then
        dummy = flist2(i)
        flist2(i) = flist2(ii)
        flist2(ii) = dummy
      end if
      if (present(flist3).EQV..true.) then
        dummy = flist3(i)
        flist3(i) = flist3(ii)
        flist3(ii) = dummy
      end if
    end if
  end do
!
end
!
!-------------------------------------------------------------------
!
! binary search for a target interval confinement in a sorted array of integers
! equivalence is allowed with left interval bound and match for series of identical values is rightmost
!
subroutine ibinary_search(arsz,arvec,keyv,ifound)
!
  implicit none
!
  integer arsz,ifound,k,ivsz
  integer keyv,arvec(*)
  logical notdone
!
  notdone = .true.
  if (keyv.lt.arvec(1)) then
    ifound = 0
    return
  else if (keyv.ge.arvec(arsz)) then
    ifound = arsz+1
    return
  end if
!
  k = 1
  ivsz = max(1,nint(dble(arsz)/2.0))
  do while (k.lt.arsz)
    if (arvec(min(arsz,k+ivsz)).gt.keyv) then
!     stay in lower half
    else
      k = min(arsz,k+ivsz)
    end if
    if (ivsz.eq.1) exit
    ivsz = max(1,nint(dble(ivsz)/2.0))
  end do
!
  ifound = k
!
end
!
!-------------------------------------------------------------------
!
! binary search for a target interval confinement in a sorted array of RTYPEs
! equivalence is allowed with left interval bound and match for series of identical values is rightmost
!
subroutine fbinary_search(arsz,arvec,keyv,ifound)
!
  implicit none
!
  integer arsz,ifound,k,ivsz
  RTYPE keyv,arvec(*)
  logical notdone
!
  notdone = .true.
  if (keyv.lt.arvec(1)) then
    ifound = 0
    return
  else if (keyv.ge.arvec(arsz)) then
    ifound = arsz+1
    return
  end if
!
  k = 1
  ivsz = max(1,nint(dble(arsz)/2.0))
  do while (k.lt.arsz)
    if (arvec(min(arsz,k+ivsz)).gt.keyv) then
!     stay in lower half
    else
      k = min(arsz,k+ivsz)
    end if
    if (ivsz.eq.1) exit
    ivsz = max(1,nint(dble(ivsz)/2.0))
  end do
!
  ifound = k
!
end
!
!-------------------------------------------------------------------
!
! binary search for a target interval confinement in a sorted array of sp floats
! equivalence is allowed with left interval bound and match for series of identical values is rightmost
!
subroutine ffbinary_search(arsz,arvec,keyv,ifound)
!
  implicit none
!
  integer arsz,ifound,k,ivsz
  real(KIND=4) keyv,arvec(*)
  logical notdone
!
  notdone = .true.
  if (keyv.lt.arvec(1)) then
    ifound = 0
    return
  else if (keyv.ge.arvec(arsz)) then
    ifound = arsz+1
    return
  end if
!
  k = 1
  ivsz = max(1,nint(dble(arsz)/2.0))
  do while (k.lt.arsz)
    if (arvec(min(arsz,k+ivsz)).gt.keyv) then
!     stay in lower half
    else
      k = min(arsz,k+ivsz)
    end if
    if (ivsz.eq.1) exit
    ivsz = max(1,nint(dble(ivsz)/2.0))
  end do
!
  ifound = k
!
end
!
!-----------------------------------------------------------------------
!
! the PRNG from the references below (ideally, multiple options
! to be offered here)
!
! literature references:
!
! P. L'Ecuyer, Communications of the ACM, 31, 742-774 (1988)
!
! W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P.
! Flannery, Numerical Recipes (Fortran), 2nd Ed., Cambridge
! University Press, 1992, Section 7-1
!
subroutine init_campprng(which)
!
  use iounit
  use interfaces
  use keys
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
#ifdef ENABLE_MPI
  use mpi
#endif
!
  implicit none
!
  integer, INTENT(IN):: which
!
  integer(KIND=8) dts(8)
  character(12) rc(3)
#ifdef ENABLE_MPI
  integer ierr,mrk
#endif
  integer cont,iomess,i,j,k
  character(MAXKWLEN) seedstr
  character(MAXKEYLEN) kline,readfrom
!
  rndstat%prng = which
!
! see L'Ecuyer
  if (which.eq.1) then
!   set parameters
    rndstat%iparams(1) = 141803398
    rndstat%iparams(2) = 32
    rndstat%iparams(3) = 2147483563
    rndstat%iparams(4) = rndstat%iparams(3) - 1
    rndstat%iparams(5) = 40014
    rndstat%iparams(6) = 53668
    rndstat%iparams(7) = 12211
!   the next 5 are not used in initialization
    rndstat%iparams(8) = 2147483399
    rndstat%iparams(9) = 40692
    rndstat%iparams(10)= 52774
    rndstat%iparams(11)= 3791
    rndstat%iparams(12)= 1+rndstat%iparams(4)/rndstat%iparams(2)
    rndstat%rparams(1) = 1.0d0/(1.0*rndstat%iparams(3))
    allocate(rndstat%itable(rndstat%iparams(2)))
    rndstat%itable(:) = 0
!
!   initial values for variables
    rndstat%ivars(1) = rndstat%iparams(1)

    call Date_and_Time(rc(1),rc(2),rc(3),dts)
    dts(1) = mod(int(dts(1)),10)-1
    rndstat%ivars(1) = rndstat%ivars(1) + 32140800*dts(1) + 2678400*(dts(2)-1)
    rndstat%ivars(1) = rndstat%ivars(1) + 86400*(dts(3)-1) + 3600*dts(5)
#ifdef ENABLE_MPI
    call MPI_COMM_RANK(MPI_COMM_WORLD,mrk,ierr)
    rndstat%ivars(1) = rndstat%ivars(1) + 60*dts(6) + dts(7) + 29*mrk
#else
    rndstat%ivars(1) = rndstat%ivars(1) + 60*dts(6) + dts(7) + dts(8)
#endif
!
!   provide further increment by process ID
    rndstat%ivars(1) = abs(rndstat%ivars(1) + handwrapped_getpid())
!
!   get a user-specified seed 
    do i=1,nkey
      cont = 1
      do j=1,size(key(i)%line)
        kline(j:j) = key(i)%line(j)
      end do
      do j=size(key(i)%line)+1,MAXKEYLEN
        kline(j:j) = ' '
      end do
      call extract_str(kline,seedstr,cont)
      call toupper(seedstr)
      if (seedstr(1:11) .eq. 'RANDOMSEED ') then
        readfrom = kline(cont:MAXKEYLEN)
        read (readfrom,*,iostat=iomess) rndstat%ivars(1)
        if (iomess.eq.IOSTAT_END) then
          exit
        else if (iomess.eq.2) then
          write(ilog,*) 'Fatal. I/O error while attempting to read random seed.'
          call fexit()
        end if
        rndstat%ivars(1) = max(1,rndstat%ivars(1))
      end if
    end do
!
    write(ilog,*)
    write(ilog,*) 'Initialized PRNG with Seed of ',rndstat%ivars(1)
    write(ilog,*)
!
!   initialize shuffling table
    rndstat%ivars(2) = rndstat%ivars(1)
    do i=rndstat%iparams(2)+8,1,-1
      k = rndstat%ivars(1) / rndstat%iparams(6)
      rndstat%ivars(1) = rndstat%iparams(5) * (rndstat%ivars(1)-k*rndstat%iparams(6)) - k*rndstat%iparams(7)
      if (rndstat%ivars(1).lt.0)  rndstat%ivars(1) = rndstat%ivars(1) + rndstat%iparams(3)
      if (i.le.rndstat%iparams(2))  rndstat%itable(i) = rndstat%ivars(1)
    end do
    rndstat%ivars(3) = rndstat%itable(1)
  else
    write(ilog,*) 'Fatal. Encountered unsupported PRNG in init_campprng(...). This is most certainly a bug.'
    call fexit()
  end if
!
end
!
!------------------------------------------------------------------------------------------------------------
!
function random()
!
  use keys, ONLY: rndstat
!
  implicit none
!
  RTYPE random
  integer i,k
!
  k = rndstat%ivars(1)/rndstat%iparams(6)
  rndstat%ivars(1) = rndstat%iparams(5)*(rndstat%ivars(1)-k*rndstat%iparams(6)) - k*rndstat%iparams(7)
  if (rndstat%ivars(1).lt.0)  rndstat%ivars(1) = rndstat%ivars(1) + rndstat%iparams(3)
  k = rndstat%ivars(2)/rndstat%iparams(10)
  rndstat%ivars(2) = rndstat%iparams(9)*(rndstat%ivars(2)-k*rndstat%iparams(10)) - k*rndstat%iparams(11)
  if (rndstat%ivars(2).lt.0)  rndstat%ivars(2) = rndstat%ivars(2) + rndstat%iparams(8)
  i = 1 + rndstat%ivars(3)/rndstat%iparams(12)
  rndstat%ivars(3) = rndstat%itable(i) - rndstat%ivars(2)
  rndstat%itable(i) = rndstat%ivars(1)
  if (rndstat%ivars(3).lt.1)  rndstat%ivars(3) = rndstat%ivars(3) + rndstat%iparams(4)
  random = rndstat%rparams(1)*rndstat%ivars(3)
!
end
!
!-------------------------------------------------------------------------------------------------------
!
function random2()
!
  use iounit
  use keys
  use interfaces
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
#ifdef ENABLE_MPI
  use mpi
#endif
!
  implicit none
!
  integer im1,ia1,iq1,ir1
  integer im2,ia2,iq2,ir2
  integer big,ntable
  integer imm1,ndiv,j
  RTYPE factor,random2
  parameter (im1=2147483563)
  parameter (ia1=40014)
  parameter (iq1=53668)
  parameter (ir1=12211)
  parameter (im2=2147483399)
  parameter (ia2=40692)
  parameter (iq2=52774)
  parameter (ir2=3791)
  parameter (big=141803398)
  parameter (ntable=32)
  parameter (imm1=im1-1)
  parameter (ndiv=1+imm1/ntable)
  parameter (factor=1.0d0/im1)
  integer i,k,seed,seed2
  integer iy,itable(ntable)
  integer(KIND=8) dts(8)
  character(12) rc(3)
#ifdef ENABLE_MPI
  integer ierr,mrk
#endif
  integer cont,iomess
  character(MAXKWLEN) seedstr
  character(MAXKEYLEN) kline,readfrom
  logical firstcall
  save firstcall,seed,seed2,iy,itable
  data firstcall /.true./
!
!
! the default seed is set to a large number and then time-incremented
!
  if (firstcall.EQV..true.) then
    firstcall = .false.
    seed = big

    call Date_and_Time(rc(1),rc(2),rc(3),dts)
    dts(1) = mod(int(dts(1)),10)-1
    seed = seed + 32140800*dts(1) + 2678400*(dts(2)-1)
    seed = seed + 86400*(dts(3)-1) + 3600*dts(5)
#ifdef ENABLE_MPI
    call MPI_COMM_RANK(MPI_COMM_WORLD,mrk,ierr)
    seed = seed + 60*dts(6) + dts(7) + 29*mrk
#else
    seed = seed + 60*dts(6) + dts(7) + dts(8)
#endif
!
!   provide further increment by process ID
    seed = abs(seed + handwrapped_getpid())
!
!   get a user-specified seed 
    do i=1,nkey
      cont = 1
      do j=1,size(key(i)%line)
        kline(j:j) = key(i)%line(j)
      end do
      do j=size(key(i)%line)+1,MAXKEYLEN
        kline(j:j) = ' '
      end do
      call extract_str(kline,seedstr,cont)
      call toupper(seedstr)
      if (seedstr(1:11) .eq. 'RANDOMSEED ') then
        readfrom = kline(cont:MAXKEYLEN)
        read (readfrom,*,iostat=iomess) seed
        if (iomess.eq.IOSTAT_END) then
          exit
        else if (iomess.eq.2) then
          write(ilog,*) 'Fatal. I/O error while attempting to read r&
 &andom seed.'
          call fexit()
        end if
        seed = max(1,seed)
      end if
    end do
!
    write(ilog,*)
    write(ilog,*) 'Initialized PRNG with Seed of ',seed
    write(ilog,*)
!
!   initialize shuffling table
!
    seed2 = seed
    do i=ntable+8,1,-1
      k = seed / iq1
      seed = ia1 * (seed-k*iq1) - k*ir1
      if (seed.lt.0)  seed = seed + im1
      if (i.le.ntable)  itable(i) = seed
    end do
    iy = itable(1)
  end if
!
! get a new random number value each call
!
  k = seed/iq1
  seed = ia1*(seed-k*iq1) - k*ir1
  if (seed.lt.0)  seed = seed + im1
  k = seed2/iq2
  seed2 = ia2*(seed2-k*iq2) - k*ir2
  if (seed2.lt.0)  seed2 = seed2 + im2
  i = 1 + iy/ndiv
  iy = itable(i) - seed2
  itable(i) = seed
  if (iy.lt.1)  iy = iy + imm1
  random2 = factor*iy
!
end
!
!----------------------------------------------------------------------
!
! this routine must only be called after initialization of PRNG has happened
!
subroutine get_nrandoms(nr,rvec)
!
  use keys, ONLY: rndstat
!
  integer, INTENT(IN):: nr
!
  RTYPE rvec(nr)
  integer i,j,k
!
  do j=1,nr
    k = rndstat%ivars(1)/rndstat%iparams(6)
    rndstat%ivars(1) = rndstat%iparams(5)*(rndstat%ivars(1)-k*rndstat%iparams(6)) - k*rndstat%iparams(7)
    if (rndstat%ivars(1).lt.0)  rndstat%ivars(1) = rndstat%ivars(1) + rndstat%iparams(3)
    k = rndstat%ivars(2)/rndstat%iparams(10)
    rndstat%ivars(2) = rndstat%iparams(9)*(rndstat%ivars(2)-k*rndstat%iparams(10)) - k*rndstat%iparams(11)
    if (rndstat%ivars(2).lt.0)  rndstat%ivars(2) = rndstat%ivars(2) + rndstat%iparams(8)
    i = 1 + rndstat%ivars(3)/rndstat%iparams(12)
    rndstat%ivars(3) = rndstat%itable(i) - rndstat%ivars(2)
    rndstat%itable(i) = rndstat%ivars(1)
    if (rndstat%ivars(3).lt.1)  rndstat%ivars(3) = rndstat%ivars(3) + rndstat%iparams(4)
    rvec(j) = rndstat%rparams(1)*rndstat%ivars(3)
  end do 
!
end
!
!-----------------------------------------------------------------------
!
! get a normally instead of uniformely distributed random number
! std. dev. is 1.0, and mean is 0.0
! since values are generated in pairwise fashion, overhang is
! stored
!
function normal()
!
  use iounit
!
  implicit none
!
  RTYPE random,f1,f2,fsq,ff
  RTYPE normal,rn1
  logical getnew
!
  save getnew,rn1
  data getnew /.true./
!
! generate new pair
!
  if (getnew) then
    f1 = 2.0*random() - 1.0
    f2 = 2.0*random() - 1.0
    fsq = f1**2 + f2**2
    do while (fsq.ge.1.0)
      f1 = 2.0*random() - 1.0
      f2 = 2.0*random() - 1.0
      fsq = f1**2 + f2**2
    end do
    ff = sqrt(-2.0*log(fsq)/fsq)
    rn1 = f1*ff
    normal = f2*ff
    getnew = .false.
!
! spit out second value from last call and reset flag
!
  else
    normal = rn1
    getnew = .true.
  end if
!
end
!
!-----------------------------------------------------------------------
!
! random variates of the gamma distribution (log-version) according to:
! Marsaglia, G., and W. W. Tsang. "A Simple Method for Generating Gamma Variables."
! ACM Transactions on Mathematical Software. Vol. 26, 2000, pp. 363–372. 
!
function rndgamma(alph)
!
  RTYPE random,normal,xn,deh,vau,larger,smaller,alph,rndgamma
!
  larger = 0.0
  smaller = 1.0
  do while (smaller > larger)
    xn = normal()
    deh = alph - (1./3.)
    vau = (1+xn/sqrt(9.0*deh))**3
    if (vau.le.0.0) cycle
    larger = 0.5*xn*xn + deh - deh*vau + deh*log(vau)
    smaller = log(random())
  end do
!
  rndgamma = deh*vau
!
end
!
!------------------------------------------------------------------------
!
! this subroutine gets a unit vector in random direction
! (== uniformly sampled, random direction)
!
subroutine ranvec(vec)
!
  use iounit
!
  implicit none
!
  RTYPE random,vec(3),f1,f2,fsq
!
! first two components have to be proper
!
  fsq = 2.0
  do while (fsq.ge.1.0)
    f1 = 2.0*random () - 1.0
    f2 = 2.0*random () - 1.0
    fsq = f1**2 + f2**2
  end do
!
! third component is then straightforward
!
  vec(3) = 1.0 - 2.0*fsq
  fsq = 2.0*sqrt(1.0 - fsq)
  vec(1) = fsq*f1
  vec(2) = fsq*f2
!
end
!
!-----------------------------------------------------------------------
!
! a simple helper to resize vectors through a common interface (resize_vec)
! this is the integer variant
!
subroutine resize_ivec(oldsz,iv,newsz,preserve)
!
  integer, INTENT(IN):: oldsz,newsz
  integer, INTENT(INOUT), ALLOCATABLE:: iv(:)
  logical, INTENT(IN):: preserve
!
  integer bu(max(1,oldsz))
!
  if ((oldsz.gt.0).AND.(preserve.EQV..true.)) bu(1:oldsz) = iv(1:oldsz)
!
  if (allocated(iv).EQV..true.) deallocate(iv)
  if (newsz.gt.0) then
    allocate(iv(newsz))
    if (preserve.EQV..true.) iv(1:min(newsz,oldsz)) = bu(1:min(newsz,oldsz))
  end if
!
end
!
!-----------------------------------------------------------------------
!
! a simple helper to resize vectors through a common interface (resize_vec)
! this is the RTYPE variant
!
subroutine resize_fvec(oldsz,iv,newsz,preserve)
!
  integer, INTENT(IN):: oldsz,newsz
  RTYPE, INTENT(INOUT), ALLOCATABLE:: iv(:)
  logical, INTENT(IN):: preserve
!
  RTYPE bu(max(1,oldsz))
!
  if ((oldsz.gt.0).AND.(preserve.EQV..true.)) bu(1:oldsz) = iv(1:oldsz)
!
  if (allocated(iv).EQV..true.) deallocate(iv)
  if (newsz.gt.0) then
    allocate(iv(newsz))
    if (preserve.EQV..true.) iv(1:min(newsz,oldsz)) = bu(1:min(newsz,oldsz))
  end if
!
end  
!
!-----------------------------------------------------------------------
!
! a simple helper to resize vectors through a common interface (resize_vec)
! this is the force sp-fp (real(KIND=4)) variant
!
subroutine resize_ffvec(oldsz,iv,newsz,preserve)
!
  integer, INTENT(IN):: oldsz,newsz
  real(KIND=4), INTENT(INOUT), ALLOCATABLE:: iv(:)
  logical, INTENT(IN):: preserve
!
  real(KIND=4) bu(max(1,oldsz))
!
  if ((oldsz.gt.0).AND.(preserve.EQV..true.)) bu(1:oldsz) = iv(1:oldsz)
!
  if (allocated(iv).EQV..true.) deallocate(iv)
  if (newsz.gt.0) then
    allocate(iv(newsz))
    if (preserve.EQV..true.) iv(1:min(newsz,oldsz)) = bu(1:min(newsz,oldsz))
  end if
!
end  
!
!----------------------------------------------------------------------------------------------------
!
#ifdef ENABLE_LINUX
!
subroutine check_memory_footprint()
! 
  use interfaces
  use iounit, ONLY: ilog
  use, INTRINSIC:: ISO_C_BINDING
!  use IFPORT
!
  implicit none
!
  integer k,ck,myp,s1,s2,s3,s4
!  integer(C_INT) handwrapped_getpid
  character(1024) cmdl
  character(MAXSTRLEN) fnn
  character(10) fx
  logical kl(4)
!
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
  write(ilog,*)
  write(ilog,*) '-------------------- MEMORY INQUIRY BY CAMPARI --------------------------'

  cmdl(:) = ' '
  cmdl="cat /proc/meminfo | awk -vk=0 '{if ($1 == "//'"MemTotal:"'//") {k=k+$2; ff = 0; if (toupper($3) == "//'"KB"'//") &
 &{ff=1.0e-6} else if (toupper($3) == "//'"MB"'//") {ff=1.0e-3} else if (toupper($3) == "//'"GB"'//") {ff=1.0}; &
 &if (ff == 0) {printf("//'"'//" The total memory on this machine is %10.4g %s.\n"//'"'//",k,toupper($3))} else &
 &{printf("//'"'//" The total memory on this machine is %10.4g GB.\n"//'"'//",ff*k)}}}'"
  call strlims(cmdl,s3,s4)
  call EXECUTE_COMMAND_LINE(cmdl(s3:s4),WAIT=.true.,EXITSTAT=k,CMDSTAT=ck)
!  k = system(cmdl(s3:s4)) ! old Intel compiler extension
  cmdl(:) = ' '
  cmdl="cat /proc/meminfo | awk -vk=0 '{if (($1 == "//'"MemFree:"'//") || ($1 == "//'"Inactive:"'//")) {k=k+$2; print k,$3}}' | &
 &tail -n1 | awk '{ff = 0; if (toupper($2) == "//'"KB"'//") &
 &{ff=1.0e-6} else if (toupper($2) == "//'"MB"'//") {ff=1.0e-3} else if (toupper($2) == "//'"GB"'//") {ff=1.0}; &
 &if (ff == 0) {printf("//'"'//" The available memory on this machine is %10.4g %s.\n"//'"'//",$1,toupper($2))} else &
 &{printf("//'"'//" The available memory on this machine is %10.4g GB.\n"//'"'//",ff*$1)}}'"
  call strlims(cmdl,s3,s4)
  call EXECUTE_COMMAND_LINE(cmdl(s3:s4),WAIT=.true.,EXITSTAT=k,CMDSTAT=ck)
!  k = system(cmdl(s3:s4)) ! old Intel compiler extension
  myp = handwrapped_getpid()
  ck = 0
  fx(:) = ' '
  call int2str(myp,fx,ck)
  call strlims(fx,s1,s2)
  fnn='/proc/'//fx(s1:s2)//'/status'
  call strlims(fnn,s3,s4)
  cmdl="cat "//fnn(s3:s4)//" | awk -v nn="//fx(s1:s2)//" '{if ($1 == "//'"VmRSS:"'//") {k=$2; ff = 0; if (toupper($3) == &
 &"//'"KB"'//") {ff=1.0e-6} else if (toupper($3) == "//'"MB"'//") {ff=1.0e-3} else if (toupper($3) == "//'"GB"'//") {ff=1.0}; &
 &if (ff == 0) {printf("//'"'//" CAMPARI (PID: %d) currently uses %10.4g %s.\n"//'"'//",nn,k,toupper($3))} else &
 &{printf("//'"'//" CAMPARI (PID: %d) currently uses %10.4g GB.\n"//'"'//",nn,ff*k)}}}'"
  call strlims(cmdl,s3,s4)
  call EXECUTE_COMMAND_LINE(cmdl(s3:s4),WAIT=.true.,EXITSTAT=k,CMDSTAT=ck)
!  k = system(cmdl(s3:s4)) ! old Intel compiler extension
!  ck = GETLASTERRORQQ()
!  write(*,*) ENOENT,ENOEXEC,ENOMEM,E2BIG
!  write(*,*) k,ck
  write(ilog,*) '-----------------END OF MEMORY INQUIRY BY CAMPARI -----------------------'
  write(ilog,*)
#ifdef ENABLE_THREADS
!$OMP END MASTER
#endif

!
end
!
#endif
! 
