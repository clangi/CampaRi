program trial

  implicit none

  integer :: distances(11)
  real :: dis_wei_in(11)
  integer, parameter :: n_snapshots = 6000
  integer, parameter :: xyz_3coo = 42
  real INPUT_trj(n_snapshots,xyz_3coo)
  integer :: d_meth_i = 13
  integer :: starter1

  ! real outp(n_snapshots,n_snapshots)
  real rad
  real inter_rad
  logical mst_iin
  logical bir_in
  logical logging
  logical normalize_it
  logical netcdfT
  logical make_it_verbose
  integer max_deg

  integer a_deg(n_snapshots)
  integer a_ix(n_snapshots,n_snapshots)
  real a_dis(n_snapshots,n_snapshots)

  integer o_invec(n_snapshots+2)
  integer o_iv2(n_snapshots)
  integer trbrkslst2(1)
  integer o_progind(n_snapshots)
  real o_distv(n_snapshots)
  integer i

  !sst
  real :: rootmax_rad_in1
  integer :: tree_height_in1
  integer :: n_search_attempts_in1

  forall(i = 1:11) distances(i) = 0
  distances(1) = 5
  ! distances(2) = 11
  ! distances(1) = 11
  forall(i = 1:11) dis_wei_in(i) = 0
  ! dis_wei_in(1) = 1
  ! dis_wei_in(2) = 0.001
  ! rad = HUGE(rad)

  mst_iin = .true.
  bir_in = .false.
  logging = .false.
  normalize_it = .true.
  make_it_verbose = .true.

  call RANDOM_NUMBER(INPUT_trj)
  INPUT_trj(500:1000,:) = INPUT_trj(500:1000,:) + 5
  ! INPUT_trj(1500:5000,:) = INPUT_trj(1500:5000,:) + 11
  rootmax_rad_in1 = (sum(INPUT_trj)/(n_snapshots*xyz_3coo))*(4.0/4.0)
  ! rootmax_rad_in1 = 7.5
  tree_height_in1 = 5
  n_search_attempts_in1 = floor(size(INPUT_trj(:,1))/100.0)
  rad = rootmax_rad_in1/tree_height_in1
  if(.not.bir_in) rad = huge(rad)
  inter_rad = HUGE(inter_rad)
  a_deg = 0
  a_ix(:,:) = 0
  a_dis(:,:) = 0.0
  max_deg = 0

! print *, n_search_attempts_in1, "ADSAJFKSFLSDGJDG"
  trbrkslst2 = 0
  o_invec = 0
  o_iv2 = 0

  o_progind = 0
  o_distv = 0

  netcdfT = .true.
  if(netcdfT) then
    print *, 'INIT NETCDF DUMPING TEST'
    call generate_neighbour_list_w( &
    INPUT_trj, xyz_3coo, n_snapshots, rad, inter_rad, & !input
    distances, dis_wei_in, bir_in, & !algorithm details
    rootmax_rad_in1, tree_height_in1, n_search_attempts_in1, & !sst details
    normalize_it, logging, make_it_verbose) !modes
    print *, ""
    ! print *, "a_deg(1:50)",a_deg(1:50)
    ! print *, ""
    ! print *, "a_ix(1,1:100)",a_ix(1,1:100)
    ! print *, ""
    ! print *, "a_dis(1,1:100)",a_dis(1,1:100)
    ! print *, ""
    ! print *, "max_deg", max_deg
    print *, ""
    print *, ""

    starter1 = 1
    call gen_progind_from_adjlst_r(n_snapshots, starter1, &
    o_progind, o_distv, o_invec, o_iv2)
    print *, ""
    print *, o_progind(1:10)
    print *, ""
    print *, o_distv(1:10)
    print *, ""
    print *, ""
    call gen_manycuts(n_snapshots, starter1,0,50,&
    o_progind,o_distv,o_invec,o_iv2,trbrkslst2)

    print *,o_invec(1:10)
  else
    ! print *, INPUT_trj
    print *,"INITIALIZATION BEFORE====maxmin:::",maxval(INPUT_trj),minval(INPUT_trj)
    call generate_neighbour_list( &
    INPUT_trj, xyz_3coo, n_snapshots, rad, inter_rad, & !input
    a_deg, a_ix, a_dis, max_deg, & !output
    distances, dis_wei_in, bir_in, mst_iin, & !algorithm details
    rootmax_rad_in1, tree_height_in1, n_search_attempts_in1, & !sst details
    normalize_it, logging, make_it_verbose) !modes
    print *, ""
    ! print *, "a_deg(1:50)",a_deg(1:50)
    ! print *, ""
    ! print *, "a_ix(1,1:100)",a_ix(1,1:100)
    ! print *, ""
    ! print *, "a_dis(1,1:100)",a_dis(1,1:100)
    ! print *, ""
    ! print *, "max_deg", max_deg
    ! print *, ""
    ! print *, ""
    ! print *,"INITIALIZATION BEFORE====maxmin:::",maxval(INPUT_trj),minval(INPUT_trj)
    ! call generate_neighbour_list( &
    ! INPUT_trj, xyz_3coo, n_snapshots, rad, inter_rad, & !input
    ! a_deg, a_ix, a_dis, max_deg, & !output
    ! distances, dis_wei_in, bir_in, mst_iin, & !algorithm details
    ! rootmax_rad_in1, tree_height_in1, n_search_attempts_in1, & !sst details
    ! normalize_it, logging, make_it_verbose) !modes
    ! print *, ""
    ! print *, "a_deg(1:50)",a_deg(1:50)
    ! print *, ""
    ! print *, "a_ix(1,1:100)",a_ix(1,1:100)
    ! print *, ""
    ! print *, "a_dis(1,1:100)",a_dis(1,1:100)
    ! print *, ""
    ! print *, "max_deg", max_deg
    ! print *, ""
    ! print *, ""
    !     !

    starter1 = 1
    call gen_progind_from_adjlst(n_snapshots, starter1, &
    max_deg, a_deg, a_ix, a_dis, &
    o_progind, o_distv, o_invec, o_iv2)
    print *, ""
    print *, o_progind(1:10)
    print *, ""
    print *, o_distv(1:10)
    print *, ""
    print *, ""
    call gen_manycuts(n_snapshots, starter1,0,50,&
    o_progind,o_distv,o_invec,o_iv2,trbrkslst2)

    print *,o_invec(1:10)
  endif
end program trial
