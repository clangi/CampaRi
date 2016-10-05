program trial

  implicit none


  integer :: distances(11)
  real :: dis_wei_in(11)
  integer, parameter :: dim1 = 1000
  integer, parameter :: dim2 = 42
  real inp(dim1,dim2)
  integer :: d_meth_i = 13

  real outp(dim1,dim1)
  real r1,r2
  logical mst_in
  integer max_deg
  integer i

  integer a_deg(dim1)
  integer a_ix(dim1,dim1)
  real a_dis(dim1,dim1)

  integer o_invec(dim1+2)
  integer o_iv2(dim1)
  integer trbrkslst2(1)
  integer o_progind(dim1)
  real o_distv(dim1)


  forall(i = 1:11) distances(i) = 0
  distances(1) = 5
  distances(2) = 11
  ! distances(1) = 11
  forall(i = 1:11) dis_wei_in(i) = 0
  ! dis_wei_in(1) = 1
  ! dis_wei_in(2) = 0.001
  r1 = 1000
  r2 = 1050
  mst_in = .true.

  call RANDOM_NUMBER(inp)
  a_deg = 0
  a_ix = 0
  a_dis = 0.0
  max_deg = 0

  trbrkslst2 = 0
  o_invec = 0
  o_iv2 = 0

  o_progind = 0
  o_distv = 0
  ! print *, inp

  call generate_neighbour_list( &
  inp, dim2, dim1, r1, r2, & !input
  a_deg, a_ix, a_dis, max_deg, & !output
  distances,dis_wei_in, mst_in, d_meth_i,.true.,.true.) !modes
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

  ! call gen_progind_from_adjlst(dim1, 10, &
  ! max_deg, a_deg, a_ix, &
  ! a_dis, o_progind, o_distv, o_invec, o_iv2)
  ! print *, ""
  ! print *, o_progind(1:10)
  ! print *, ""
  ! print *, o_distv(1:10)
  ! print *, ""
  ! print *, ""
  ! call gen_manycuts(dim1,10,0,50,&
  ! o_progind,o_distv,o_invec,o_iv2,trbrkslst2)
  !
end program trial
