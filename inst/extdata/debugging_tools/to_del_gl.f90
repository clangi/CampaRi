program to_delexe

  implicit none

  integer :: dist = 5
  integer, parameter :: n_snaps = 100
  integer, parameter :: n_xyz = 10
  real INPUT_trj(n_snaps, n_xyz)
  real :: mahala_mat(n_xyz, n_xyz) = 0
  ! real outp(n_snaps,n_snaps)
  real radius
  real hardcut
  logical :: do_mst = .true.
  logical :: birch_it = .false.
  logical :: normalize_it = .false.
  logical :: do_netcdf = .false.
  logical :: mute_it = .false.
  logical :: return_tree = .true.

  integer :: a_deg(n_snaps) = 0
  integer :: a_ix(n_snaps,n_snaps) = 0
  real :: a_dis(n_snaps,n_snaps) = 0
  integer :: max_deg = 0

  integer :: o_invec(n_snaps+2) = 0
  integer :: o_iv2(n_snaps) = 0
  integer :: trbrkslst2(1) = 0
  integer :: o_progind(n_snaps) = 0
  real :: o_distv(n_snaps) = 0
  integer :: starting_sn = 1
  integer i

  !sst
  real :: rootmax_radius
  integer :: tree_height = 5
  integer :: n_search_attempts

  ! forall(i = 1:11) dist(i) = 0


  call RANDOM_NUMBER(INPUT_trj)
  rootmax_radius = (sum(INPUT_trj)/(n_snaps*n_xyz))*(4.0/4.0)
  n_search_attempts = floor(size(INPUT_trj(:,1))/100.0)
  radius = rootmax_radius/tree_height
  if(.not.birch_it) radius = huge(radius)
  hardcut = HUGE(hardcut)

  do_netcdf = .false.
  if(do_netcdf) then
    print *, 'INIT NETCDF DUMPING TEST'
    ! call generate_neighbour_list_w( &
    ! INPUT_trj, n_xyz, n_snaps, radius, hardcut, & !input
    ! dist, dis_wei_in, birch_it, & !algorithm details
    ! rootmax_radius, tree_height, n_search_attempts, & !sst details
    ! normalize_it, logging, make_it_verbose) !modes
    ! print *, ""
    ! print *, ""
  else
    print *,"INITIALIZATION BEFORE====maxmin:::",maxval(INPUT_trj),minval(INPUT_trj)
    call generate_neighbour_list( &
    INPUT_trj, n_xyz, n_snaps, n_snaps, radius, hardcut, & !input
    a_deg, a_ix, a_dis, max_deg, & !output
    dist, mahala_mat, birch_it, do_mst, & !algorithm details
    rootmax_radius, tree_height, n_search_attempts, & !sst details
    normalize_it, return_tree, mute_it) !modes
    print *, ""
    print *, "a_deg(1:50)",a_deg(1:50)
    print *, ""
    print *, "a_ix(1,1:100)",a_ix(1,1:100)
    print *, ""
    print *, "a_dis(1,1:100)",a_dis(1,1:100)
    print *, ""
    print *, "max_deg", max_deg
    print *, ""
    print *, ""
    print *, '----------------------------------------------------------------'
    print *, '----------------------------------------------------------------'
    call gen_progind_from_adjlst(n_snaps, starting_sn, &
    max_deg, a_deg, a_ix, a_dis, &
    o_progind, o_distv, n_snaps, o_invec, o_iv2,&
    mute_it, return_tree)
    print *, ""
    print *, o_progind(1:10)
    print *, ""
    print *, o_distv(1:10)
    print *, ""
    print *, ""
    call gen_manycuts(n_snaps, starting_sn,0,50,&
    o_progind,o_distv,o_invec,o_iv2,trbrkslst2, mute_it)

    print *,o_invec(1:10)
  endif
end program to_delexe
