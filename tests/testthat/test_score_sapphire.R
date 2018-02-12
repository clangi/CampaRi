context('score_sapphire')

test_that('scoring sapphire plots', {
  # require(testthat)
  silent <- T
  plt_stff <- T
  # thedata <- matrix(c(rnorm(10000, sd = 1), rnorm(10000, mean = 10, sd = 1), rnorm(10000, mean = 20, sd = 1)), nrow = 3000, ncol = 10)
  thedata <- matrix(c(rnorm(10000, sd = 1), rnorm(10000, mean = 5, sd = 1), rnorm(10000, mean = 10, sd = 1)), nrow = 3000, ncol = 10)
  # thedata <- matrix(c(rnorm(1000, sd = 1), rnorm(1000, mean = 5, sd = 1), rnorm(1000, mean = 10, sd = 1)), nrow = 300, ncol = 10)
  if(plt_stff) plot(c(thedata))
  
  expect_error(adjl <- mst_from_trj(trj = thedata, dump_to_netcdf = FALSE, mute_fortran = silent), NA)
  expect_error(ret <- gen_progindex(adjl, snap_start = 21, mute = silent), NA)
  expect_error(ret2 <- gen_annotation(ret, snap_start = 21, mute = silent), NA)

  the_sap <- 'REPIX_000000000021.dat'
  if(plt_stff) sapphire_plot(sap_file = the_sap, ann_trace = ann, timeline = T); cat(round(sqrt(nrow(thedata)*10)))
  nbin <- round(sqrt(nrow(thedata)*10))
  nbin <- 20
  if(plt_stff) basins_recognition(the_sap, nx = nbin, match = F, plot = T, out.file = F, new.dev = F)
  
  ann <- c(rep(1, 1000), rep(2, 1000), rep(3, 1000))
  # debugonce(score_sapphire)
  # score_sapphire(the_sap = the_sap, ann = ann)
  expect_error(qres <- score_sapphire(the_sap = the_sap, ann = ann, silent = sileit), NA)
  # print(qres)
  
  thedata <- matrix(c(rnorm(10000, sd = 1), rnorm(10000, mean = 5, sd = 1), rnorm(10000, mean = 10, sd = 1)), nrow = 3000, ncol = 10)
  expect_error(adjl <- mst_from_trj(trj = thedata, dump_to_netcdf = FALSE), NA)
  expect_error(ret <- gen_progindex(adjl, snap_start = 21), NA)
  expect_error(ret2 <- gen_annotation(ret, snap_start = 21), NA)
  the_sap <- 'REPIX_000000000021.dat'
  # sapphire_plot(sap_file = the_sap, ann_trace = ann)
  ann <- c(rep(1, 1000), rep(2, 1000), rep(3, 1000))
  expect_error(qres <- score_sapphire(the_sap = the_sap, ann = ann, silent = sileit), NA)
  # qres
  
  if(file.exists('MST_DUMPLING.nc')) file.remove('MST_DUMPLING.nc')
  if(file.exists('REPIX_000000000021.dat')) file.remove('REPIX_000000000021.dat')
})
