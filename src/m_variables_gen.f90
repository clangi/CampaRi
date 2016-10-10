module m_variables_gen
  integer dis_method(11)
  real dis_weight(11) !distance weight when averaged
  integer tmp_dis_method
  integer n_dis_method ! number of distances to be averaged
  logical normalize_dis !normalizing distances
  logical birch
  integer n_xyz ! allocation sizes of length of a snapshot array (xyz*3 in 3D)
  integer n_snaps ! global number of snapshots. It is necessary as alsz roof
  logical ver !global verbose var
  logical superver ! dev flag for superverbose output
  integer :: ilog = 378 !code reference to logging file
end module m_variables_gen
