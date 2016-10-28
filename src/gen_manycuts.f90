
!-------------------------------------------------------------------------------
!
subroutine gen_manycuts(n_snaps,start,ntbrks2,&
  pwidth,setis,distv,invvec,ivec2,trbrkslst)
!
  implicit none
!
  integer, INTENT(IN) :: n_snaps, start
  integer, INTENT(IN) :: invvec(n_snaps+2), ivec2(n_snaps), setis(n_snaps)
!list of breaks (must be passed at least of size 1, but n_breaks can be zero)
  integer, INTENT(IN) :: ntbrks2 !number of breaks == 0
  integer, INTENT(IN) :: trbrkslst(max(1,ntbrks2)) !one value == 0

  real(KIND=4), INTENT(IN) :: distv(n_snaps) !distance list
!
!
  integer tmat(3,3),vv1(3),vv2(3),vv3(3)
  integer pvo(5),pvn(5),iu,i,j,k,pwidth,ii,jj,kk,ll,tl2
  integer freeunit !this is the return value of the function that looks for new units to use
!
  integer, ALLOCATABLE:: cutv(:)
  logical, ALLOCATABLE:: ibrkx(:)
!
  character(12) nod2
  character(250) fn
  logical exists

  write(*,*) "Generating annotations..."

  tl2 = 12
  call int2str(start,nod2,tl2)
  fn = 'REPIX_'//nod2(1:tl2)//'.dat'
  inquire(file=fn(1:22),exist=exists)
  if(exists) then
    iu = freeunit()
    open(unit=iu,file=fn(1:22),status='old',position='append')
    close(unit=iu,status='delete')
  end if
  iu=freeunit()
  open(unit=iu,file=fn(1:22),status='new')
!
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
    write(iu,666) i,n_snaps-i,setis(i),cutv(i),distv(i),ivec2(i),invvec(ivec2(i)+1),cutv(invvec(ivec2(i)+1)),vv1,vv2,vv3,ii,jj,kk
  end do
 666 format(4(i10,1x),g12.5,1x,1000(i10,1x))
!
  deallocate(cutv)
  deallocate(ibrkx)
!
  close(unit=iu)
  write(*,*) "DONE"
!
end

!---------------------------------------------------------------------------------------
!
subroutine int2str(ii,string,strsz)
  implicit none
!
  integer mx(0:16),leng,i,ii,strsz,minsz,cpy,sumi
  character it(0:9)
  character(*) string
  logical just_r
  data it /'0','1','2','3','4','5','6','7','8','9'/
!
  if (strsz.le.0) then
    just_r = .true.
    strsz = 1
  else
    just_r = .false.
  end if
  minsz = strsz
  leng = len(string)
!
! invert sign if necessary (technically this routines only
! processes unsigned integers)
  if (ii .lt. 0) then
    i = -ii
  end if
!
! this is standard: the mx(...) are integers of course
!
  sumi = 0
  do i=16,0,-1
    mx(i) = (ii-sumi)/(10.0**i)
    if (i.gt.0) then
      sumi = sumi + (10.0**i) * mx(i)
    end if
  end do
!
! now we have the digits of our string in integer form (0 through 9)
!
! check for integer-size
!
  if ((mx(10).gt.0).OR.(mx(11).gt.0).OR.(mx(12).gt.0).OR.&
 &    (mx(13).gt.0).OR.&
 &    (mx(14).gt.0).OR.(mx(15).gt.0).OR.(mx(16).gt.0)) then
    write(*,*) 'Warning. Possibly bad result from int2str(...) (inte&
 &ger overflow).'
  end if
!
! find the final string length
!
  strsz = 1
  do i=16,1,-1
    if (mx(i).ne.0) then
      strsz = i+1
      exit
    end if
  end do
!
! correct if too large or small
  strsz = min(strsz,leng)
  strsz = max(strsz,minsz)
!
  if (strsz.gt.17) then
    write(*,*) 'Warning. Possibly bad result from int2str(...) (digi&
 &t overflow).'
  end if
!
! convert individual digits to a string of numeric characters
!
  do i=1,strsz
    string(i:i) = it(mx(strsz-i))
  end do
!
! left/right-justification
!
  if (just_r.EQV..true.) then
    do i=strsz,1,-1
      cpy = leng-strsz+i
      string(cpy:cpy) = string(i:i)
    end do
    do i=1,leng-strsz
      string(i:i) = ' '
    end do
  else
    do i = strsz+1,leng
      string(i:i) = ' '
    end do
  end if
!
end
!
