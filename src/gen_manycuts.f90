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
! MAIN AUTHOR:   Andreas Vitalis                                           !
! CONTRIBUTIONS: Nicolas Bloechliger                                       !
! WRAPPER: Davide Garolini                                                 !
!                                                                          !
!--------------------------------------------------------------------------!

!-------------------------------------------------------------------------------
!
subroutine gen_manycuts(n_snaps,start,ntbrks2,&
  pwidth,setis,distv,invvec,ivec2,trbrkslst,mute_in)
  use m_hw_fprintf
  use gutenberg
  use, INTRINSIC :: ISO_C_BINDING
  implicit none
!
  integer, INTENT(IN) :: n_snaps, start
  integer, INTENT(IN) :: invvec(n_snaps+2), ivec2(n_snaps), setis(n_snaps)
!list of breaks (must be passed at least of size 1, but n_breaks can be zero)
  integer, INTENT(IN) :: ntbrks2 !number of breaks == 0
  integer, INTENT(IN) :: trbrkslst(max(1,ntbrks2)) !one value == 0

  real, INTENT(IN) :: distv(n_snaps) !distance list
!
!
  integer tmat(3,3),vv1(3),vv2(3),vv3(3)
  integer pvo(5),pvn(5),i,j,k,pwidth,ii,jj,kk,ll,tl2
  ! tl2
!
  integer, ALLOCATABLE :: cutv(:)
  logical, ALLOCATABLE :: ibrkx(:)
!
  character(12) nod2
  ! character(250) fn
  logical exists
  integer iu, freeunit !this is the return value of the function that looks for new units to use

  integer(c_int) :: retv        ! (C integer for return value)
  integer(c_int) :: newline_sel ! (C integer for argument whether to print \n)
  integer(c_int) :: i1,i3,setis1,cutv1,ivec21,&
  invvec1,cutv12,vv11,vv21,vv31,vv12,vv22,vv32,vv13,vv23,vv33,ii1,jj1,kk1
  real(c_double) :: distv1
  integer i2
!
  character(c_char), TARGET :: fn(0:100) ! C character arrays (for const char*)
  ! character(5000), TARGET :: ff3                    ! Fortran string
  character(250) :: fcfd
  logical mute_in
  mute = mute_in

  call spr('------------------------------------------------------------')
  call spr("Generating annotations...")

  tl2 = 12
  call int2str(start,nod2,tl2)
  ! for the filename
  fcfd = 'REPIX_'//nod2(1:tl2)//'.dat'
  inquire(file=fcfd(1:22),exist=exists)
  if(exists) then
    iu = freeunit()
    open(unit=iu, file=fcfd(1:22), status='old')
    close(unit=iu, status='delete')
  end if

  do i2=1,LEN(TRIM(fcfd))
    fn(i2-1) = fcfd(i2:i2)
  end do
  fn(LEN(TRIM(fcfd))) = C_NULL_CHAR
  retv = 0

  allocate(cutv(n_snaps))
  allocate(ibrkx(n_snaps+1))
  ibrkx(:) = .false.
  do i=1,ntbrks2
    ibrkx(trbrkslst(i)+1) = .true.
  end do
!
  cutv(1) = 2
  if (ibrkx(n_snaps+1).EQV..true.) cutv(1) = 1
  do i=2,n_snaps
    kk = setis(i) + 1
    if ((invvec(kk+1).lt.i).AND.(invvec(kk-1).lt.i)) then
      cutv(i) = cutv(i-1) - 2
      if (ibrkx(kk-1).EQV..true.) cutv(i) = cutv(i) + 1
      if (ibrkx(kk).EQV..true.) cutv(i) = cutv(i) + 1
    else if ((invvec(kk+1).gt.i).AND.(invvec(kk-1).gt.i)) then
      cutv(i) = cutv(i-1) + 2
      if (ibrkx(kk-1).EQV..true.) cutv(i) = cutv(i) - 1
      if (ibrkx(kk).EQV..true.) cutv(i) = cutv(i) - 1
    else
      cutv(i) = cutv(i-1)
      if (ibrkx(kk-1).EQV..true.) then
        if (invvec(setis(i)).gt.i) then
          cutv(i) = cutv(i) - 1
        else
          cutv(i) = cutv(i) + 1
        end if
      end if
      if (ibrkx(kk).EQV..true.) then
        if (invvec(setis(i)+2).gt.i) then
          cutv(i) = cutv(i) - 1
        else
          cutv(i) = cutv(i) + 1
        end if
      end if
    end if
  end do


!----------------------------------------------------------------------
! This is the BAD local cut (3state mfpt)
!----------------------------------------------------------------------
  tmat(:,:) = 0
  tmat(3,3) = n_snaps+1-ntbrks2
  do i=1,pwidth
    kk = setis(i) + 1
    pvo(2) = 3
    if (invvec(kk+1).gt.i) then
      pvo(3) = 3
    else
      pvo(3) = 2
    end if
    if (invvec(kk-1).gt.i) then
      pvo(1) = 3
    else
      pvo(1) = 2
    end if
    pvn(1:3) = pvo(1:3)
    pvn(2) = 2
    do k=1,2
      if (ibrkx(kk+k-2).EQV..true.) cycle
      tmat(pvo(k),pvo(k+1)) = tmat(pvo(k),pvo(k+1)) - 1
      tmat(pvn(k),pvn(k+1)) = tmat(pvn(k),pvn(k+1)) + 1
    end do
  end do
  do i=1,n_snaps
    if ((i.gt.pwidth).AND.((n_snaps-i).ge.pwidth)) then
      if (((abs(setis(i+pwidth)-setis(i-pwidth)).le.1).AND.(abs(setis(i)-setis(i-pwidth)).le.1)).OR.&
 &        ((abs(setis(i+pwidth)-setis(i-pwidth)).le.1).AND.(abs(setis(i)-setis(i+pwidth)).le.1)).OR.&
 &        ((abs(setis(i)-setis(i+pwidth)).le.1).AND.(abs(setis(i)-setis(i-pwidth)).le.1))) then
!       this should be excessively rare
        j = min(setis(i+pwidth),setis(i-pwidth))
        j = min(j,setis(i))
        pvo(1:5) = 3
        if (invvec(j).ge.i) then
          if (abs(invvec(j)-i).lt.pwidth) pvo(1) = 2
        else
          if (abs(invvec(j)-i).le.pwidth) pvo(1) = 1
        end if
        if (invvec(j+4).ge.i) then
          if (abs(invvec(j+4)-i).lt.pwidth) pvo(5) = 2
        else
          if (abs(invvec(j+4)-i).le.pwidth) pvo(5) = 1
        end if
        pvn(1:5) = pvo(1:5)
        do k=j,j+2
          if (k.eq.setis(i)) then
            pvo(k-j+2) = 2
            pvn(k-j+2) = 1
          else if (k.eq.setis(i-pwidth)) then
            pvo(k-j+2) = 1
            pvn(k-j+2) = 3
          else if (k.eq.setis(i+pwidth)) then
            pvo(k-j+2) = 3
            pvn(k-j+2) = 2
          end if
        end do
        do k=1,4
          if (ibrkx(j+k-1).EQV..true.) cycle
          tmat(pvo(k),pvo(k+1)) = tmat(pvo(k),pvo(k+1)) - 1
          tmat(pvn(k),pvn(k+1)) = tmat(pvn(k),pvn(k+1)) + 1
        end do
      else if ((abs(setis(i)-setis(i+pwidth)).le.1).OR.(abs(setis(i)-setis(i-pwidth)).le.1).OR.&
 &(abs(setis(i+pwidth)-setis(i-pwidth)).le.1)) then
        if (abs(setis(i+pwidth)-setis(i-pwidth)).le.1) then
          ll = 0
          ii = min(setis(i+pwidth),setis(i-pwidth)) + 1
          kk = max(setis(i+pwidth),setis(i-pwidth)) + 1
          if (setis(i+pwidth).gt.setis(i-pwidth)) then
            pvo(2) = 1
            pvo(3) = 3
          else
            pvo(2) = 3
            pvo(3) = 1
          end if
          pvo(4) = 3
          if (invvec(kk+1).ge.i) then
            if (abs(invvec(kk+1)-i).lt.pwidth) pvo(4) = 2
          else
            if (abs(invvec(kk+1)-i).le.pwidth) pvo(4) = 1
          end if
          pvo(1) = 3
          if (invvec(ii-1).ge.i) then
            if (abs(invvec(ii-1)-i).lt.pwidth) pvo(1) = 2
          else
            if (abs(invvec(ii-1)-i).le.pwidth) pvo(1) = 1
          end if
          pvn(1:4) = pvo(1:4)
          if (setis(i+pwidth).gt.setis(i-pwidth)) then
            pvn(2) = 3
            pvn(3) = 2
          else
            pvn(2) = 2
            pvn(3) = 3
          end if
          do k=1,3
            if (ibrkx(ii+k-2).EQV..true.) cycle
            tmat(pvo(k),pvo(k+1)) = tmat(pvo(k),pvo(k+1)) - 1
            tmat(pvn(k),pvn(k+1)) = tmat(pvn(k),pvn(k+1)) + 1
          end do
        else if (abs(setis(i)-setis(i-pwidth)).le.1) then
          ll = 1
          ii = min(setis(i),setis(i-pwidth)) + 1
          kk = max(setis(i),setis(i-pwidth)) + 1
          if (setis(i).gt.setis(i-pwidth)) then
            pvo(2) = 1
            pvo(3) = 2
          else
            pvo(2) = 2
            pvo(3) = 1
          end if
          pvo(4) = 3
          if (invvec(kk+1).gt.i) then
            if (abs(invvec(kk+1)-i).lt.pwidth) pvo(4) = 2
          else
            if (abs(invvec(kk+1)-i).le.pwidth) pvo(4) = 1
          end if
          pvo(1) = 3
          if (invvec(ii-1).gt.i) then
            if (abs(invvec(ii-1)-i).lt.pwidth) pvo(1) = 2
          else
            if (abs(invvec(ii-1)-i).le.pwidth) pvo(1) = 1
          end if
          pvn(1:4) = pvo(1:4)
          if (setis(i).gt.setis(i-pwidth)) then
            pvn(2) = 3
            pvn(3) = 1
          else
            pvn(2) = 1
            pvn(3) = 3
          end if
          do k=1,3
            if (ibrkx(ii+k-2).EQV..true.) cycle
            tmat(pvo(k),pvo(k+1)) = tmat(pvo(k),pvo(k+1)) - 1
            tmat(pvn(k),pvn(k+1)) = tmat(pvn(k),pvn(k+1)) + 1
          end do
        else if (abs(setis(i)-setis(i+pwidth)).le.1) then
          ll = -1
          ii = min(setis(i),setis(i+pwidth)) + 1
          kk = max(setis(i),setis(i+pwidth)) + 1
          if (setis(i).gt.setis(i+pwidth)) then
            pvo(2) = 3
            pvo(3) = 2
          else
            pvo(2) = 2
            pvo(3) = 3
          end if
          pvo(4) = 3
          if (invvec(kk+1).gt.i) then
            if (abs(invvec(kk+1)-i).lt.pwidth) pvo(4) = 2
          else
            if (abs(invvec(kk+1)-i).le.pwidth) pvo(4) = 1
          end if
          pvo(1) = 3
          if (invvec(ii-1).gt.i) then
            if (abs(invvec(ii-1)-i).lt.pwidth) pvo(1) = 2
          else
            if (abs(invvec(ii-1)-i).le.pwidth) pvo(1) = 1
          end if
          pvn(1:4) = pvo(1:4)
          if (setis(i).gt.setis(i+pwidth)) then
            pvn(2) = 2
            pvn(3) = 1
          else
            pvn(2) = 1
            pvn(3) = 2
          end if
          do k=1,3
            if (ibrkx(ii+k-2).EQV..true.) cycle
            tmat(pvo(k),pvo(k+1)) = tmat(pvo(k),pvo(k+1)) - 1
            tmat(pvn(k),pvn(k+1)) = tmat(pvn(k),pvn(k+1)) + 1
          end do
        end if
        do j=ll,ll
          kk = setis(i+j*pwidth) + 1
          pvo(2) = j+2
          pvo(3) = 3
          if (invvec(kk+1).gt.i) then
            if (abs(invvec(kk+1)-i).lt.pwidth) pvo(3) = 2
          else
            if (abs(invvec(kk+1)-i).le.pwidth) pvo(3) = 1
          end if
          pvo(1) = 3
          if (invvec(kk-1).gt.i) then
            if (abs(invvec(kk-1)-i).lt.pwidth) pvo(1) = 2
          else
            if (abs(invvec(kk-1)-i).le.pwidth) pvo(1) = 1
          end if
          pvn(1:3) = pvo(1:3)
          pvn(2) = j+1
          if (pvn(2).eq.0) pvn(2) = 3
          do k=1,2
            if (ibrkx(kk+k-2).EQV..true.) cycle
            tmat(pvo(k),pvo(k+1)) = tmat(pvo(k),pvo(k+1)) - 1
            tmat(pvn(k),pvn(k+1)) = tmat(pvn(k),pvn(k+1)) + 1
          end do
        end do
      else
        do j=-1,1
          kk = setis(i+j*pwidth) + 1
          pvo(2) = j+2
          pvo(3) = 3
          if (invvec(kk+1).gt.i) then
            if (abs(invvec(kk+1)-i).lt.pwidth) pvo(3) = 2
          else
            if (abs(invvec(kk+1)-i).le.pwidth) pvo(3) = 1
          end if
          pvo(1) = 3
          if (invvec(kk-1).gt.i) then
            if (abs(invvec(kk-1)-i).lt.pwidth) pvo(1) = 2
          else
            if (abs(invvec(kk-1)-i).le.pwidth) pvo(1) = 1
          end if
          pvn(1:3) = pvo(1:3)
          pvn(2) = j+1
          if (pvn(2).eq.0) pvn(2) = 3
          do k=1,2
            if (ibrkx(kk+k-2).EQV..true.) cycle
            tmat(pvo(k),pvo(k+1)) = tmat(pvo(k),pvo(k+1)) - 1
            tmat(pvn(k),pvn(k+1)) = tmat(pvn(k),pvn(k+1)) + 1
          end do
        end do
      end if
    else if (i.le.pwidth) then
      if (abs(setis(i)-setis(i+pwidth)).le.1) then
        ii = min(setis(i),setis(i+pwidth)) + 1
        kk = max(setis(i),setis(i+pwidth)) + 1
        if (setis(i).gt.setis(i+pwidth)) then
          pvo(2) = 3
          pvo(3) = 2
        else
          pvo(2) = 2
          pvo(3) = 3
        end if
        pvo(4) = 3
        if (invvec(kk+1).gt.i) then
          if (abs(invvec(kk+1)-i).lt.pwidth) pvo(4) = 2
        else
          if (abs(invvec(kk+1)-i).le.pwidth) pvo(4) = 1
        end if
        pvo(1) = 3
        if (invvec(ii-1).gt.i) then
          if (abs(invvec(ii-1)-i).lt.pwidth) pvo(1) = 2
        else
          if (abs(invvec(ii-1)-i).le.pwidth) pvo(1) = 1
        end if
        pvn(1:4) = pvo(1:4)
        if (setis(i).gt.setis(i+pwidth)) then
          pvn(2) = 2
          pvn(3) = 1
        else
          pvn(2) = 1
          pvn(3) = 2
        end if
        do k=1,3
          if (ibrkx(ii+k-2).EQV..true.) cycle
          tmat(pvo(k),pvo(k+1)) = tmat(pvo(k),pvo(k+1)) - 1
          tmat(pvn(k),pvn(k+1)) = tmat(pvn(k),pvn(k+1)) + 1
        end do
      else
        do j=0,1
          kk = setis(i+j*pwidth) + 1
          pvo(2) = j+2
          pvo(3) = 3
          if (invvec(kk+1).gt.i) then
            if (abs(invvec(kk+1)-i).lt.pwidth) pvo(3) = 2
          else
            if (abs(invvec(kk+1)-i).le.pwidth) pvo(3) = 1
          end if
          pvo(1) = 3
          if (invvec(kk-1).gt.i) then
            if (abs(invvec(kk-1)-i).lt.pwidth) pvo(1) = 2
          else
            if (abs(invvec(kk-1)-i).le.pwidth) pvo(1) = 1
          end if
          pvn(1:3) = pvo(1:3)
          pvn(2) = j+1
          if (pvn(2).eq.0) pvn(2) = 3
          do k=1,2
            if (ibrkx(kk+k-2).EQV..true.) cycle
            tmat(pvo(k),pvo(k+1)) = tmat(pvo(k),pvo(k+1)) - 1
            tmat(pvn(k),pvn(k+1)) = tmat(pvn(k),pvn(k+1)) + 1
          end do
        end do
      end if
    else
      if (abs(setis(i)-setis(i-pwidth)).le.1) then
        ii = min(setis(i),setis(i-pwidth)) + 1
        kk = max(setis(i),setis(i-pwidth)) + 1
        if (setis(i).gt.setis(i-pwidth)) then
          pvo(2) = 1
          pvo(3) = 2
        else
          pvo(2) = 2
          pvo(3) = 1
        end if
        pvo(4) = 3
        if (invvec(kk+1).gt.i) then
          if (abs(invvec(kk+1)-i).lt.pwidth) pvo(4) = 2
        else
          if (abs(invvec(kk+1)-i).le.pwidth) pvo(4) = 1
        end if
        pvo(1) = 3
        if (invvec(ii-1).gt.i) then
          if (abs(invvec(ii-1)-i).lt.pwidth) pvo(1) = 2
        else
          if (abs(invvec(ii-1)-i).le.pwidth) pvo(1) = 1
        end if
        pvn(1:4) = pvo(1:4)
        if (setis(i).gt.setis(i-pwidth)) then
          pvn(2) = 3
          pvn(3) = 1
        else
          pvn(2) = 1
          pvn(3) = 3
        end if
        do k=1,3
          if (ibrkx(ii+k-2).EQV..true.) cycle
          tmat(pvo(k),pvo(k+1)) = tmat(pvo(k),pvo(k+1)) - 1
          tmat(pvn(k),pvn(k+1)) = tmat(pvn(k),pvn(k+1)) + 1
        end do
      else
        do j=-1,0
          kk = setis(i+j*pwidth) + 1
          pvo(2) = j+2
          pvo(3) = 3
          if (invvec(kk+1).gt.i) then
            if (abs(invvec(kk+1)-i).lt.pwidth) pvo(3) = 2
          else
            if (abs(invvec(kk+1)-i).le.pwidth) pvo(3) = 1
          end if
          pvo(1) = 3
          if (invvec(kk-1).gt.i) then
            if (abs(invvec(kk-1)-i).lt.pwidth) pvo(1) = 2
          else
            if (abs(invvec(kk-1)-i).le.pwidth) pvo(1) = 1
          end if
          pvn(1:3) = pvo(1:3)
          pvn(2) = j+1
          if (pvn(2).eq.0) pvn(2) = 3
          do k=1,2
            if (ibrkx(kk+k-2).EQV..true.) cycle
            tmat(pvo(k),pvo(k+1)) = tmat(pvo(k),pvo(k+1)) - 1
            tmat(pvn(k),pvn(k+1)) = tmat(pvn(k),pvn(k+1)) + 1
          end do
        end do
      end if
    end if
    ii = min(i,pwidth)
    jj = min(n_snaps-i,pwidth)
    kk = max(0,n_snaps - ii - jj)
    vv1(:) = tmat(1,:)
    vv2(:) = tmat(2,:)
    vv3(:) = tmat(3,:)
    ! integer(c_int) :: i1,i3,setis1,cutv1,ivec21,&
    ! invvec1,cutv1,vv11,vv21,vv31,ii1,jj1,kk1
    ! real(c_double) :: distv1
    ! write(ff3,666) i,n_snaps-i,setis(i),cutv(i),distv(i),ivec2(i),&
    ! invvec(ivec2(i)+1),cutv(invvec(ivec2(i)+1)),vv1,vv2,vv3,ii,jj,kk
    ! transcribe Fortran string into C character array (ff3 to f2)
    ! do i2=1,LEN(TRIM(ff3))
    !   f2(i2-1) = ff3(i2:i2)
    ! end do
    ! !
    ! ! terminate C const char appropriately
    ! f2(LEN(TRIM(ff3))) = C_NULL_CHAR

    i1 = i
    i3 = n_snaps-i
    setis1 = setis(i)
    cutv1 = cutv(i)
    distv1 = distv(i)
    ivec21 = ivec2(i)
    invvec1 = invvec(ivec2(i)+1)
    cutv12 = cutv(invvec(ivec2(i)+1))
    vv11 = vv1(1)
    vv21 = vv2(1)
    vv31 = vv3(1)
    vv12 = vv1(2)
    vv22 = vv2(2)
    vv32 = vv3(2)
    vv13 = vv1(3)
    vv23 = vv2(3)
    vv33 = vv3(3)
    ii1 = ii
    jj1 = jj
    kk1 = kk


    ! invoke fprintf wrapper: the return value is that from the fprintf call
    newline_sel = 1
    retv = handwrapped_fprint(fn,i1,i3,setis1,cutv1,distv1,ivec21,&
    invvec1,cutv12,vv11,vv21,vv31,vv12,vv22,vv32,vv13,vv23,vv33,ii1,jj1,kk1,&
    newline_sel)

  end do
 ! 666 format(4(i10,1x),g12.5,1x,1000(i10,1x))
!
  if(retv.ne.0) then
    call spr("File written successfully...")
  else
    call spr("The file writing was not completely successful...")
  end if
  deallocate(cutv)
  deallocate(ibrkx)
!
  ! close(unit=iu)
  call spr("...done")
  call spr('------------------------------------------------------------')!
end
