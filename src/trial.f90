program trial

  implicit none

  integer :: distances(11)
  real :: dis_wei_in(11)
  integer, parameter :: n_snapshots = 1000
  integer, parameter :: xyz_3coo = 42
  real INPUT_trj(n_snapshots,xyz_3coo)
  integer :: d_meth_i = 13

  ! real outp(n_snapshots,n_snapshots)
  real rad
  real inter_rad
  logical mst_iin
  logical bir_in
  logical logging
  logical normalize_it
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

  inter_rad = HUGE(inter_rad)
  mst_iin = .true.
  bir_in = .false.
  logging = .false.
  normalize_it = .true.
  make_it_verbose = .true.

  call RANDOM_NUMBER(INPUT_trj)
  rootmax_rad_in1 = (sum(INPUT_trj)/(n_snapshots*xyz_3coo))*(10.0/4.0)
  tree_height_in1 = 5
  n_search_attempts_in1 = floor(size(INPUT_trj)/10.0)
  rad = rootmax_rad_in1/tree_height_in1
  if(.not.bir_in) rad = huge(rad)
  a_deg = 0
  a_ix = 0
  a_dis = 0.0
  max_deg = 0

  trbrkslst2 = 0
  o_invec = 0
  o_iv2 = 0

  o_progind = 0
  o_distv = 0
  ! print *, INPUT_trj
  print *,"INITIALIZATION BEFORE====max:::::::::::::::::::::::::::::::",maxval(INPUT_trj)
  call generate_neighbour_list( &
  INPUT_trj, xyz_3coo, n_snapshots, rad, inter_rad, & !input
  a_deg, a_ix, a_dis, max_deg, & !output
  distances, dis_wei_in, bir_in, mst_iin, & !algorithm details
  rootmax_rad_in1, tree_height_in1, n_search_attempts_in1, & !sst details
  d_meth_i, normalize_it, logging, make_it_verbose) !modes
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

  ! call gen_progind_from_adjlst(n_snapshots, 10, &
  ! max_deg, a_deg, a_ix, &
  ! a_dis, o_progind, o_distv, o_invec, o_iv2)
  ! print *, ""
  ! print *, o_progind(1:10)
  ! print *, ""
  ! print *, o_distv(1:10)
  ! print *, ""
  ! print *, ""
  ! call gen_manycuts(n_snapshots,10,0,50,&
  ! o_progind,o_distv,o_invec,o_iv2,trbrkslst2)
  !
  ! print *,o_invec(1:10)

end program trial
