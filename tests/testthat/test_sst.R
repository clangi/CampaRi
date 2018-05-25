context('sst')

test_that('Building a short spanning tree', {
  
  silent <- T
  plt_stff <- !silent
  if(!silent) require(testthat)
  
  expect_error(adjl <- mst_from_trj(trj = matrix(rnorm(1000), nrow = 100, ncol = 10), dump_to_netcdf = FALSE, birch_clu = T, mute_fortran = silent), NA)
  expect_error(ret <- gen_progindex(adjl, snap_start = 21), NA)
  expect_error(ret2 <- gen_annotation(ret, snap_start = 21), NA)
  # plottt <- sapphire_plot(sap_file = 'REPIX_000000000021.dat', title = "CAMPARI WRAPPER - MST", return_plot = TRUE)
  
  expect_error(mst_from_trj(trj = matrix(rnorm(1000), nrow = 100, ncol = 10), dump_to_netcdf = TRUE, 
               normalize_d = FALSE, clu_radius = 10, clu_hardcut = 100, min_span_tree = TRUE, birch_clu = T), NA)
  expect_true(file.exists("MST_DUMPLING.nc"))
  expect_error(ret <- gen_progindex(nsnaps = 100, read_from_netcdf = T, snap_start = 21), NA)
  expect_error(ret2 <- gen_annotation(ret, snap_start = 21, local_cut_width = 20), NA)
  # plottt <- sapphire_plot(sap_file = 'REPIX_000000000021.dat', title = "CAMPARI WRAPPER - MST", return_plot = TRUE)

  if(file.exists('MST_DUMPLING.nc')) file.remove('MST_DUMPLING.nc')
  if(file.exists('REPIX_000000000021.dat')) file.remove('REPIX_000000000021.dat')
})
