!------------------------------------------------------------------------------------------------------
!
! using PRIM's algorithm, this routine will trace the MST adjacency list to derive the progress index
! per se (along with the minimal distance to the current set (distv), the inverse map (invvec)
! and the parent/source vector (iv2))
!
subroutine gen_progind_from_adjlst_r(n_snaps_in, starter,  progind, distv, invvec, iv2)
  ! mnb, alnbs, alst, aldis,
  use m_variables_gen
  use m_mst_dumping
  implicit none
!
  integer, INTENT(IN):: n_snaps_in, starter
!
  integer, INTENT(OUT):: invvec(n_snaps+2),progind(n_snaps),iv2(n_snaps)
  real, INTENT(OUT):: distv(n_snaps)
!
  integer heapsize,lprogind,j,i
  logical, ALLOCATABLE:: added(:), inprogind(:)
  integer, ALLOCATABLE:: heap(:),hsource(:) !dynamic dimensions
  real, ALLOCATABLE:: key(:)
  real t4,t3 !timing variables
  allocate(inprogind(n_snaps)) !the maximum length of these are the number of snapshot
  allocate(added(n_snaps))
  allocate(key(n_snaps))
  allocate(heap(n_snaps))
  allocate(hsource(n_snaps))
  n_snaps = n_snaps_in
  ilog = 6
  write(ilog,*) '---------------------------------------------------------------------'
  write(ilog,*) "Generating progress index from netcdf file..."
  write(ilog,*)

  write(ilog,*) 'Reading mst from working directory file "MST_DUMPLING.nc"...'
  call CPU_time(t3)
  call read_nbl_nc()
  call CPU_time(t4)
  write(ilog,*) '...reading completed after ',t4-t3, ' [s]'
  write(ilog,*)

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
    write(ilog,*) 'Warning. The snapshot index requested is not available (there &
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
  heapsize = approxmst(progind(lprogind))%deg !number of connection of the first vertex
  heap(1:heapsize) = approxmst(progind(lprogind))%adj(1:heapsize) !the specific connections indexes
  key(1:heapsize) = approxmst(progind(lprogind))%dist(1:heapsize) !key is the value of the distances
  hsource(1:heapsize) = progind(lprogind) !value of the prog index that is generating the connections
  call hbuild(heapsize,heap(1:heapsize),key(1:heapsize),hsource(1:heapsize))
  do j=1,approxmst(progind(lprogind))%deg
    added(approxmst(progind(lprogind))%adj(j)) = .true.
  end do
  do while (lprogind.lt.n_snaps)
    do j=1,approxmst(progind(lprogind))%deg
!! add neighbors of last snapshot in progress index to heap
      if (added(approxmst(progind(lprogind))%adj(j)).EQV..false.) then
        call hinsert(heapsize,heap(1:(heapsize+1)),key(1:(heapsize+1)),hsource(1:(heapsize+1)),&
 &                   approxmst(progind(lprogind))%adj(j),approxmst(progind(lprogind))%dist(j),progind(lprogind))
        added(approxmst(progind(lprogind))%adj(j)) = .true.
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
  if (allocated(approxmst).EQV..true.) then
    do i=1,size(approxmst)
      if (allocated(approxmst(i)%adj).EQV..true.) deallocate(approxmst(i)%adj)
      if (allocated(approxmst(i)%dist).EQV..true.) deallocate(approxmst(i)%dist)
    end do
    deallocate(approxmst)
  end if
  write(ilog,*)
  write(ilog,*) "...progrex index generated"
  write(ilog,*) '---------------------------------------------------------------------'
!
end

! please note that here I tried to copy the other things in gen_progind but it
! was duplicating the name. There could be a missing link to that functions though
