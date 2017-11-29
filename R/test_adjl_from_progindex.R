context('adjl_from_progindex')

test_that('Test adjl_from_progindex', {
  adjl <- mst_from_trj(trj = matrix(rnorm(1000), nrow = 100, ncol = 10), dump_to_netcdf = FALSE)
  expect_true(!is.null(adjl))
  adjl2 <- contract_mst(adjl, n_fold = 2)
  expect_true(!is.null(adjl2))
  ret <- gen_progindex(adjl, snap_start = 21)
  expect_true(!is.null(ret))
  ret2 <- gen_annotation(ret, snap_start = 21)
  expect_error(a <- adjl_from_progindex(prog_index_file = 'REPIX_000000000021.dat'), NA)
  
  if(file.exists('MST_DUMPLING.nc')) file.remove('MST_DUMPLING.nc')
  if(file.exists('REPIX_000000000021.dat')) file.remove('REPIX_000000000021.dat')
  })



