program trial

 implicit none

 real inp(100,42)
 integer :: d5 = 5
 integer :: dim1 = 100
 integer :: dim2 = 42
 integer :: d_meth_i = 11
 real outp(100,100)
 real r1,r2
 logical mst_in
 integer max_deg
 integer i

 integer a_deg(100)
 integer a_ix(100,100)
 real a_dis(100,100)

integer o_invec(102)
integer o_iv2(100)
integer trbrkslst2(1)
integer o_progind(100)
real o_distv(100)

 r1 = 1.0
 r2 = 1.5
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
  d5, mst_in, d_meth_i, .true.) !modes
 print *, ""
 print *, a_deg
 print *, ""
 print *, a_ix(1,:)
 print *, ""
 print *, a_dis(1,:)
 print *, ""
 print *, max_deg
 print *, ""
 print *, ""

 call gen_progind_from_adjlst(dim1, 10, &
 max_deg, a_deg, a_ix, &
 a_dis, o_progind, o_distv, o_invec, o_iv2)
 print *, ""
 print *, o_progind(1:10)
 print *, ""
 print *, o_distv(1:10)
 print *, ""
 print *, ""
 call gen_manycuts(dim1,10,0,50,&
 o_progind,o_distv,o_invec,o_iv2,trbrkslst2)



end program trial
