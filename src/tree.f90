! ELIMINATE ALST, REINSERT DISTV PASS
!
!------------------------------------------------------------------------------------------------------
!
! using PRIM's algorithm, this routine will trace the MST adjacency list to derive the progress index
! per se (along with the minimal distance to the current set (distv), the inverse map (invvec)
! and the parent/source vector (iv2))
!
subroutine gen_progind_from_adjlst(n_snaps, starter, mnb, alnbs, alst, &
  aldis, progind, distv, invvec, iv2)
!
  implicit none
!
  integer, INTENT(IN):: n_snaps,mnb,starter,alnbs(n_snaps),alst(n_snaps,mnb)
  real(KIND=4), INTENT(IN):: aldis(n_snaps,mnb)
!
  integer, INTENT(OUT):: invvec(n_snaps+2),progind(n_snaps),iv2(n_snaps)
  real(KIND=4), INTENT(OUT):: distv(n_snaps)
!
  integer heapsize,lprogind,j
  logical, ALLOCATABLE:: added(:), inprogind(:)
  integer, ALLOCATABLE:: heap(:),hsource(:) !dynamic dimensions
  real(KIND=4), ALLOCATABLE:: key(:)
  write(ilog,*) "Generating progress index..."
  allocate(inprogind(n_snaps)) !the maximum length of these are the number of snapshot
  allocate(added(n_snaps))
  allocate(key(n_snaps))
  allocate(heap(n_snaps))
  allocate(hsource(n_snaps))
  added(:) = .false.
  inprogind(:) = .false.
  key(:) = 0.0
  hsource(:) = 0
! add first snapshot to progress index:
  if (starter.le.n_snaps) then
    progind(1) = starter
  else
    progind(1) = 1
  end if
  lprogind = 1 !this is the position in the list of progind (it is the order)
  added(progind(lprogind)) = .true.
  inprogind(progind(lprogind)) = .true. !what is the difference with the above?
  invvec(progind(lprogind)+1) = 1 !inverse map
  iv2(lprogind) = progind(lprogind) !parent/source vector
!
! build heap:
  heapsize = alnbs(progind(lprogind)) !number of connection of a certain vertex
  heap(1:heapsize) = alst(progind(lprogind),1:heapsize) !the specific connections indexes
  key(1:heapsize) = aldis(progind(lprogind),1:heapsize) !key is the value of the distances
  hsource(1:heapsize) = progind(lprogind) !value of the prog index that is generating the connections
  call hbuild(heapsize,heap(1:heapsize),key(1:heapsize),hsource(1:heapsize))
  do j=1,alnbs(progind(lprogind))
    added(alst(progind(lprogind),j)) = .true.
  end do
  do while (lprogind.lt.n_snaps)
    do j=1,alnbs(progind(lprogind))
!! add neighbors of last snapshot in progress index to heap
      if (added(alst(progind(lprogind),j)).EQV..false.) then
        call hinsert(heapsize,heap(1:(heapsize+1)),key(1:(heapsize+1)),hsource(1:(heapsize+1)),&
 &                   alst(progind(lprogind),j),aldis(progind(lprogind),j),progind(lprogind))
        added(alst(progind(lprogind),j)) = .true.
      end if
    end do
! append next snapshot to progind():
    lprogind = lprogind+1
    progind(lprogind) = heap(1)
    iv2(lprogind) = hsource(1)
    distv(lprogind) = key(1)
    invvec(heap(1)+1) = lprogind
    inprogind(progind(lprogind)) = .true.
! remove added snapshot from heap:
    call hremovemin(heapsize,heap(1:heapsize),key(1:heapsize),hsource(1:heapsize))
  end do
!
  do j=1,n_snaps
    if (distv(j).lt.0.0) then
      if (distv(j).ge.-10.0*TINY(distv(j))) then
        distv(j) = 0.0
      else
        distv(j) = -1.0/distv(j)
      end if
    end if
  end do
!
  deallocate(key)
  deallocate(added)
  deallocate(inprogind)
  deallocate(heap)
  deallocate(hsource)
  write(ilog,*) "DONE"
!
end
!---------------------------------------------------------------------------------------
!
subroutine heapify(n,heap,key,source,a)
!this is the input:
!heapsize, !number of connections
!heap(1:heapsize), !specific connections ids
!key(1:heapsize), !distances of those
!hsource(1:heapsize) !prog vertex id
!
!each variable was already defined in the upper function. Here they ara again defined except for
!do a=hparent(n),1,-1
!  call heapify(n,heap,key,source,a)
!end do
  implicit none
!
  integer n,a,i,x,hleft,hright
  integer heap(n),source(n)
  real(KIND=4) key(n)
  logical isheap
!
  isheap = .false.
  i = a
  do while(isheap.EQV..false.)
    x = i
    if (hleft(i).le.n) then
      if (key(hleft(i)).lt.key(x)) then
        x = hleft(i)
      end if
    end if
    if (hright(i).le.n) then !no double logic? &&??
      if (key(hright(i)).lt.key(x)) then
        x = hright(i)
      end if
    end if
    if (x.eq.i) then
      isheap = .true.
    else
      call hswap(n,heap,key,source,i,x)
      i = x
    end if
  end do
end
!
!-----------------------------------------------------------------------------------
!
subroutine hswap(n,heap,key,source,i,x)
!
  implicit none
!
  integer n,i,x
  integer heap(n),source(n)
  real(KIND=4) key(n)
  integer temp,tempsource
  real(KIND=4) tempkey
!
  temp = heap(i)
  tempkey = key(i)
  tempsource = source(i)
  heap(i) = heap(x)
  source(i) = source(x)
  key(i) = key(x)
  heap(x) = temp
  key(x) = tempkey
  source(x) = tempsource
end
!
!-----------------------------------------------------------------------------------
!
subroutine hbuild(n,heap,key,source)
!this is the input:
!heapsize, !number of connections
!heap(1:heapsize), !specific connections ids
!key(1:heapsize), !distances of those
!hsource(1:heapsize) !prog vertex id
  implicit none
!
  integer n,hparent,a
  integer heap(n),source(n)
  real(KIND=4) key(n)
!
  do a=hparent(n),1,-1
    call heapify(n,heap,key,source,a)
  end do
end
!
!-----------------------------------------------------------------------------------
!
function hparent(i)
!
  implicit none
!
  integer hparent,i
!
  if (i.eq.1) then
!    write(ilog,*) 'NICO: Bug! Parent of root of heap requested.'
!    call fexit()
  end if
  hparent = i/2
!
end
!
!-----------------------------------------------------------------------------------
!
function hleft(i)
!
  implicit none
!
  integer hleft,i
!
  hleft = 2*i
end
!
!-----------------------------------------------------------------------------------
!
function hright(i)
!
  implicit none
!
  integer hright,i
!
  hright = 2*i+1
end
!
!-----------------------------------------------------------------------------------
!
subroutine hinsert(n,heap,key,source,new,newkey,newsource)
!
  implicit none
!
  integer n,new,hparent,i,newsource
  integer heap(n+1),source(n+1)
  real(KIND=4) key(n+1)
  real(KIND=4) newkey
!
  n = n+1
  heap(n) = new
  key(n) = newkey
  source(n) = newsource
  i = n
  do while (i.gt.1)
    if (key(i).lt.key(hparent(i))) then
      call hswap(n,heap,key,source,i,hparent(i))
      i = hparent(i)
    else
      exit
    end if
  end do
end
!
!-----------------------------------------------------------------------------------
!
subroutine hremovemin(n,heap,key,source)
!
  implicit none
!
  integer n
  integer heap(n),source(n)
  real(KIND=4) key(n)
!
  heap(1) = heap(n)
  key(1) = key(n)
  source(1) = source(n)
  n = n-1
  call heapify(n,heap,key,source,1)
!
end
!
!---------------------------------------------------------------------------------------
!
subroutine int2str(ii,string,strsz)
!
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
    write(ilog,*) 'Warning. Possibly bad result from int2str(...) (inte&
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
    write(ilog,*) 'Warning. Possibly bad result from int2str(...) (digi&
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
!-----------------------------------------------------------------------
!
! a function which provides an empty file handle
! it just loops up until it finds an open one (tests explicitly)
!
function freeunit()
!
  implicit none
!
  integer freeunit
  logical inuse
!
! try each logical unit until an unopened one is found
!
  freeunit = 0
  inuse = .true.
  do while (inuse)
    freeunit = freeunit + 1
    if ((freeunit.ne.5).AND.(freeunit.ne.6)) then
      inquire (unit=freeunit,opened=inuse)
    end if
  end do
!
end
!---------------------------------------------------------------------------------------------------
!
!
subroutine contract_mst(n_snaps,mnb,alnbs,alst,aldis,nrnds,istats)
!
  implicit none
!
  integer, INTENT(IN):: n_snaps,mnb,nrnds,alnbs(n_snaps),alst(n_snaps,mnb)
  real(KIND=4), INTENT(INOUT):: aldis(n_snaps,mnb)
  integer, INTENT(OUT):: istats(2)
!
  integer i,j,kk,thej,ll
!
  logical, ALLOCATABLE:: terminal(:)
!
  allocate(terminal(n_snaps))
  terminal(:) = .false.
!
  istats(2) = 0
  do i=1,n_snaps
    if (alnbs(i).eq.1) then
      terminal(i) = .true.
      istats(2) = istats(2) + 1
    end if
  end do
  write(ilog,*) nrnds,n_snaps,istats(2)
!
  istats(1) = 0
  do kk=1,nrnds
    do i=1,n_snaps
      if (terminal(i).EQV..true.) then
        do j=1,alnbs(i)
          if (aldis(i,j).ge.0.0) then
            thej = alst(i,j)
            if (aldis(i,j).eq.0.0) then
              aldis(i,j) = -10.0*TINY(aldis(i,j))
            else
              aldis(i,j) = -1.0/aldis(i,j)
            end if
            istats(1) = istats(1) + 1
            exit
          end if
        end do
        do ll=1,alnbs(thej)
          if (alst(thej,ll).eq.i) then
            aldis(thej,ll) = aldis(i,j)
            exit
          end if
        end do
      end if
    end do
!   reset terminal status
    terminal(:) = .false.
    istats(2) = 0
    do i=1,n_snaps
      thej = 0
      do j=1,alnbs(i)
        if (aldis(i,j).ge.0.0) thej = thej + 1
      end do
      if (thej.eq.1) then
        terminal(i) = .true.
        istats(2) = istats(2) + 1
      end if
    end do
    if (istats(2).eq.0) exit
  end do
!
  deallocate(terminal)
!
end
!
!-------------------------------------------------------------------------------
!
subroutine gen_manycuts(n_snaps,start,ntbrks2,pwidth,setis,distv,invvec,ivec2,trbrkslst)
!
  implicit none
!
  integer, INTENT(IN):: n_snaps,start,ntbrks2
  integer, INTENT(IN):: invvec(n_snaps+2),ivec2(n_snaps),setis(n_snaps),trbrkslst(max(1,ntbrks2))
  real(KIND=4), INTENT(IN):: distv(n_snaps)
!
!
  integer tmat(3,3),vv1(3),vv2(3),vv3(3)
  integer pvo(5),pvn(5),iu,i,j,k,pwidth,ii,jj,kk,ll,tl2,freeunit
!
  integer, ALLOCATABLE:: cutv(:)
  logical, ALLOCATABLE:: ibrkx(:)
!
  character(12) nod2
  character(250) fn
  logical exists

  write(ilog,*) "Generating annotations..."

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
!
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
  write(ilog,*) "DONE"
!
end
!
!----------------------------------------------------------------------------------------------------------------------------------!
