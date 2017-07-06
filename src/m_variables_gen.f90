module m_variables_gen
  integer dis_method
  logical normalize_dis !normalizing distances
  integer :: ilog = 6  !code reference to logging file
  logical birch
  logical do_mst
  integer n_xyz ! allocation sizes of length of a snapshot array (xyz*3 in 3D)
  integer n_snaps ! global number of snapshots. It is necessary as alsz roof
  integer rand_seed !counts of the reliability of the random generation
  logical firstcall !keeps the seeding under control
  integer seed, seed2, iy, itable(32)!usually itable has ntable dim but it is 32
end module m_variables_gen
