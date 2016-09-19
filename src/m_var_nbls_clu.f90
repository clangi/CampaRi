module m_var_nbls_clu
  integer dis_method
  integer n_xyz ! allocation sizes of length of a snapshot array (xyz*3 in 3D)
  integer n_snaps ! global number of snapshots. It is necessary as alsz roof
  logical ver !global verbose var
  logical superver ! dev flag for superverbose output
end module m_var_nbls_clu
