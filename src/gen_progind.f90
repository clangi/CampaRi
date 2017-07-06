! ELIMINATE ALST, REINSERT DISTV PASS
!
!------------------------------------------------------------------------------------------------------
!
! using PRIM's algorithm, this routine will trace the MST adjacency list to derive the progress index
! per se (along with the minimal distance to the current set (distv), the inverse map (invvec)
! and the parent/source vector (iv2))
!
subroutine gen_progind_from_adjlst(n_snaps_in, starting_snap, max_degree, &
  alnbs, alst, aldis, & !minimum spanning tree in input (not necessary in netcdf)
  progind, distv, &
  dfffo, invvec, iv2, mute_in, use_tree_from_r)

  ! added that must be modified in R:
  ! - dfffo
  ! - mute_in
  ! - use_tree_from_r
  ! - n_snaps_in
  ! TODO: check if the file is available and use it or R (if both inserted)
  ! no probably better to keep using the r data -> better checks in R

  use m_variables_gen
#ifdef LINK_NETCDF
  use m_mst_dumping
#endif
  use gutenberg
  implicit none
!
  integer, INTENT(IN) :: n_snaps_in, dfffo, max_degree, starting_snap
  integer, INTENT(IN) :: alnbs(dfffo), alst(dfffo, max_degree)
  real, INTENT(IN) :: aldis(dfffo,max_degree)
  logical, INTENT(IN) :: mute_in
!
  integer, INTENT(OUT) :: invvec(n_snaps_in+2),progind(n_snaps_in),iv2(n_snaps_in)
  real, INTENT(OUT) :: distv(n_snaps_in)
!
  integer heapsize, lprogind, j, i
  real t4, t3 !timing variables
  logical, ALLOCATABLE :: added(:), inprogind(:)
  integer, ALLOCATABLE :: heap(:),hsource(:) !dynamic dimensions
  real, ALLOCATABLE :: key(:)

  mute = mute_in
  n_snaps = n_snaps_in

  ! CHEKS IF THE USER HAS WHAT HE WANTS - NETCDF
  ! ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---
#ifdef LINK_NETCDF
  if(use_tree_from_r) then
    call sl()
    call spr('------------------------------------------------------------------')
    call spr('ATTENTION: Even if CampaRi was installed using netcdf support, you selected &
    to use the R data management system. Both options will be followed at the same time.')
    call spr('------------------------------------------------------------------')
  end if
#else
  if(.not.use_tree_from_r) then
    call sl()
    call spr('------------------------------------------------------------------')
    call spr('ATTENTION: Even if CampaRi was installed without the netcdf support, &
    the user tried to use the netcdf dumping functionality. This run will follow the &
    usual flow without netcdf. If you want to use it install the full version of the package.')
    call spr('------------------------------------------------------------------')
  end if
#endif
  ! ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---

  call spr('---------------------------------------------------------------------')
#ifdef LINK_NETCDF
  call spr("Generating progress index from netcdf file...")
  call sl()
  call spr('Reading mst from working directory file "MST_DUMPLING.nc"...')
  call CPU_time(t3)
  call read_nbl_nc()
  call CPU_time(t4)
  call srpr('...reading completed after (s) ',t4-t3)
  call sl()
#else
  call spr("Generating progress index from inserted minimum spanning tree (no netcdf)...")
#endif

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
  if (starting_snap.le.n_snaps) then
    progind(1) = starting_snap
  else
    progind(1) = 1
    call spr('Warning. The snapshot index requested as starting point is not available. &
    Using first one instead.')
  end if
  lprogind = 1 !this is the position in the list of progind (it is the order)
  distv(lprogind) = 0.0
  added(progind(lprogind)) = .true.
  inprogind(progind(lprogind)) = .true. !what is the difference with the above?
  invvec(progind(lprogind)+1) = 1 !inverse map
  iv2(lprogind) = progind(lprogind) !parent/source vector
!
! build heap:
  ! inizialization with the first point of the progress index (1 or starting_snap)
#ifdef LINK_NETCDF
  heapsize = approxmst(progind(lprogind))%deg !number of connection of the first vertex
  heap(1:heapsize) = approxmst(progind(lprogind))%adj(1:heapsize) !the specific connections indexes
  key(1:heapsize) = approxmst(progind(lprogind))%dist(1:heapsize) !key is the value of the distances
#else
  heapsize = alnbs(progind(lprogind)) !number of connection of the first vertex
  heap(1:heapsize) = alst(progind(lprogind),1:heapsize) !the specific connections indexes
  key(1:heapsize) = aldis(progind(lprogind),1:heapsize) !key is the value of the distances
#endif


  hsource(1:heapsize) = progind(lprogind) !value of the prog index that is generating the connections
  call hbuild(heapsize,heap(1:heapsize),key(1:heapsize),hsource(1:heapsize))


#ifdef LINK_NETCDF
  do j=1,approxmst(progind(lprogind))%deg
    added(approxmst(progind(lprogind))%adj(j)) = .true.
  end do
#else
  do j=1,alnbs(progind(lprogind))
    added(alst(progind(lprogind),j)) = .true.
  end do
#endif


  do while (lprogind.lt.n_snaps)

#ifdef LINK_NETCDF
    do j=1,approxmst(progind(lprogind))%deg
    !! add neighbors of last snapshot in progress index to heap
      if (added(approxmst(progind(lprogind))%adj(j)).EQV..false.) then
        call hinsert(heapsize,heap(1:(heapsize+1)),key(1:(heapsize+1)),hsource(1:(heapsize+1)),&
    &                   approxmst(progind(lprogind))%adj(j),approxmst(progind(lprogind))%dist(j),progind(lprogind))
        added(approxmst(progind(lprogind))%adj(j)) = .true.
      end if
    end do
#else
    do j=1,alnbs(progind(lprogind))
!! add neighbors of last snapshot in progress index to heap
      if (added(alst(progind(lprogind),j)).EQV..false.) then
        call hinsert(heapsize,heap(1:(heapsize+1)),key(1:(heapsize+1)),hsource(1:(heapsize+1)),&
 &                   alst(progind(lprogind),j),aldis(progind(lprogind),j),progind(lprogind))
        added(alst(progind(lprogind),j)) = .true.
      end if
    end do
#endif
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
#ifdef LINK_NETCDF
  if (allocated(approxmst).EQV..true.) then
    do i=1,size(approxmst)
      if (allocated(approxmst(i)%adj).EQV..true.) deallocate(approxmst(i)%adj)
      if (allocated(approxmst(i)%dist).EQV..true.) deallocate(approxmst(i)%dist)
    end do
    deallocate(approxmst)
  end if
#endif
  call sl()
  call spr("...progrex index generated")
  call spr('---------------------------------------------------------------------')
!
!
end
!
!-------------------------------------------------------------------------------
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
  real key(n)
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
  real key(n)
  integer temp,tempsource
  real tempkey
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
  real key(n)
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
  real key(n+1)
  real newkey
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
  real key(n)
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
