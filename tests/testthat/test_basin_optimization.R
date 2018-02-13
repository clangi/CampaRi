context('basin_optimization')

test_that('optimize sapphire plot clusters using SBR', {
  # require(testthat)
  silent <- F
  plt_stff <- T
  # thedata <- matrix(c(rnorm(10000, sd = 1), rnorm(10000, mean = 10, sd = 1), rnorm(10000, mean = 20, sd = 1)), nrow = 3000, ncol = 10)
  thedata <- matrix(c(rnorm(10000, sd = 1), rnorm(10000, mean = 5, sd = 1), rnorm(10000, mean = 10, sd = 1)), nrow = 3000, ncol = 10)
  # thedata <- matrix(c(rnorm(1000, sd = 1), rnorm(1000, mean = 5, sd = 1), rnorm(1000, mean = 10, sd = 1)), nrow = 300, ncol = 10)
  if(plt_stff) plot(c(thedata))
  
  expect_error(adjl <- mst_from_trj(trj = thedata, dump_to_netcdf = FALSE, mute_fortran = silent), NA)
  expect_error(ret <- gen_progindex(adjl, snap_start = 21, mute = silent), NA)
  expect_error(ret2 <- gen_annotation(ret, snap_start = 21, mute = silent), NA)
  
  ann <- c(rep(1, 1000), rep(2, 1000), rep(3, 1000))
  the_sap <- 'REPIX_000000000021.dat'
  if(plt_stff) sapphire_plot(sap_file = the_sap, ann_trace = ann, timeline = T); cat(round(sqrt(nrow(thedata)*10)))
  nbin <- round(sqrt(nrow(thedata)*10))
  if(plt_stff) basins_recognition(data = the_sap, nx = nbin, match = F, plot = T, out.file = F, new.dev = F)
  # if(plt_stff) basins_recognition(the_sap, nx = nbin, match = T, plot = T, out.file = F, new.dev = F)
  expect_error(optimal_bas <- CampaRi::basin_optimization(the_sap = the_sap,  how_fine_search = 100, number_of_clusters = 3, force_matching = T, silent = silent), NA)
  expect_error(optimal_bas <- CampaRi::basin_optimization(the_sap = the_sap,  how_fine_search = 100, force_matching = T, silent = silent))
  # optimal_bas$Optimal_nbins
  # optimal_bas$bas
  
  if(file.exists('MST_DUMPLING.nc')) file.remove('MST_DUMPLING.nc')
  if(file.exists('REPIX_000000000021.dat')) file.remove('REPIX_000000000021.dat')
})

# require(microbenchmark); require(ggplot2); v <- 1:10e7
# a <- microbenchmark("bs" = .BiSearch(v, 432), "classic" = (432 %in% v), times = 10)
# autoplot(a)