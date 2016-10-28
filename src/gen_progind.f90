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
  write(*,*) "Generating progress index..."
  allocate(inprogind(n_snaps)) !the maximum length of these are the number of snapshot
  allocate(added(n_snaps))
  allocate(key(n_snaps))
  allocate(heap(n_snaps))
  allocate(hsource(n_snaps))
  added(:) = .false.
  inprogind(:) = .false.
  distv(:) = 0.0
  key(:) = 0.0
  hsource(:) = 0
  !is this the boundary condition?
  invvec(1) = 3*n_snaps !-> boundary conditions?
  invvec(n_snaps+2) = 6*n_snaps
! add first snapshot to progress index:
  if (starter.le.n_snaps) then
    progind(1) = starter
  else
    progind(1) = 1
    write(*,*) 'Warning. The snapshot index requested is not available (there &
    &are ',n_snaps,' snapshots in memory). &
    &Using first one instead.'
  end if
  lprogind = 1 !this is the position in the list of progind (it is the order)
  distv(lprogind) = 0.0
  added(progind(lprogind)) = .true.
  inprogind(progind(lprogind)) = .true. !what is the difference with the above?
  invvec(progind(lprogind)+1) = 1 !inverse map
  iv2(lprogind) = progind(lprogind) !parent/source vector
!
! build heap:
  ! inizialization with the first point of the progress index (1 or starter)
  heapsize = alnbs(progind(lprogind)) !number of connection of the first vertex
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
  write(*,*) "DONE"
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
!    write(*,*) 'NICO: Bug! Parent of root of heap requested.'
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
!
!----------------------------------------------------------------------------------------------------------------------------------!
