context('sst')


test_that('Building a short spanning tree', {
  adjl <- mst_from_trj(trj = matrix(rnorm(1000), nrow = 100, ncol = 10), dump_to_netcdf = FALSE, birch_clu = T)
  expect_true(!is.null(adjl))
  ret <- gen_progindex(adjl, snap_start = 21)
  expect_true(!is.null(ret))
  ret2 <- gen_annotation(ret, snap_start = 21)
  expect_true(!is.null(ret2))
  plottt <- sapphire_plot(sap_file = 'REPIX_000000000021.dat', title = "CAMPARI WRAPPER - MST", return_plot = TRUE)
  expect_true(!is.null(plottt))
  mst_from_trj(trj = matrix(rnorm(1000), nrow = 100, ncol = 10), dump_to_netcdf = TRUE, 
               normalize_d = FALSE, clu_radius = 10, clu_hardcut = 100, min_span_tree = TRUE, birch_clu = T)
  expect_true(file.exists("MST_DUMPLING.nc"))
  ret <- gen_progindex(nsnaps = 100, read_from_netcdf = T, snap_start = 21)
  expect_true(!is.null(ret))
  ret2 <- gen_annotation(ret, snap_start = 21, local_cut_width = 20)
  expect_true(!is.null(ret2))
  plottt <- sapphire_plot(sap_file = 'REPIX_000000000021.dat', title = "CAMPARI WRAPPER - MST", return_plot = TRUE)
  expect_true(!is.null(plottt))
})
