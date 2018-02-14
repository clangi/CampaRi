context('basin_optimization')

test_that('optimize sapphire plot clusters using SBR', {
  
  # CREATION OF THE DATASET
  silent <- F
  plt_stff <- !silent
  if(!silent) require(testthat)
  # stdd <- 0.1; n_dim <- 20; n_snap <- 6000; n_tot <- n_dim*n_snap/4; if(!silent) print(n_tot)
  # thedata <- matrix(c(rnorm(n_tot, sd = stdd), rnorm(n_tot, mean = 5, sd = stdd), rnorm(n_tot, mean = 10, sd = stdd),
  #                     rnorm(n_tot, sd = stdd)), nrow = n_snap, ncol = n_dim)
  # stdd <- 1; n_dim <- 2; n_snap <- 3000; n_tot <- n_dim*n_snap/6; if(!silent) print(n_tot)
  # thedata <- matrix(c(rnorm(n_tot, sd = stdd), rnorm(n_tot, mean = 5, sd = stdd), rnorm(n_tot, mean = 10, sd = stdd),
  #                     rnorm(n_tot, sd = stdd), rnorm(n_tot, mean = 15, sd=stdd), rnorm(n_tot, mean = 5, sd=stdd)), nrow = n_snap, ncol = n_dim)
  stdd <- 0.5; n_dim <- 3; n_snap <- 3000; n_tot <- n_dim*n_snap/3; if(!silent) print(n_tot)
  thedata <- matrix(c(rnorm(n_tot, sd = stdd), rnorm(n_tot/2, mean = 5, sd = stdd), rnorm(n_tot, mean = 10, sd = stdd),
                      rnorm(n_tot/2, mean = 5, sd = stdd)), nrow = n_snap, ncol = n_dim)
  if(plt_stff) plot(c(thedata))
  
  adjl <- mst_from_trj(trj = thedata, normalize_d = T, dump_to_netcdf = FALSE, mute_fortran = silent)
  ret <- gen_progindex(adjl, snap_start = 21, mute = silent)
  ret2 <- gen_annotation(ret, snap_start = 21, mute = silent)
  
  # ann <- c(rep(1, n_snap/4), rep(2, n_snap/4), rep(3, n_snap/4),  rep(1, n_snap/4))
  # ann <- c(rep(1, n_snap/6), rep(2, n_snap/6), rep(3, n_snap/6),  rep(1, n_snap/6),  rep(4, n_snap/6),  rep(2, n_snap/6))
  # ann <- c(rep(1, n_snap/3), rep(2, n_snap/3), rep(3, n_snap/3))
  ann <- c(rep(1, n_snap/3), rep(2, n_snap/6), rep(3, n_snap/3), rep(2, n_snap/6))
  the_sap <- 'REPIX_000000000021.dat'
  # the_sap <- 'PROGIDX_000000000021.dat'
  if(plt_stff) sapphire_plot(sap_file = the_sap, ann_trace = ann,only_timeline=F, timeline = T)
  
  # BASINS OPTIMIZATION
  nbin <- round(sqrt(nrow(thedata)*10)); cat(round(sqrt(nrow(thedata)*10)))
  if(plt_stff) basins_recognition(data = the_sap, nx = nbin, match = F, plot = T, out.file = F, new.dev = F)
  # if(plt_stff) basins_recognition(the_sap, nx = nbin, match = T, plot = T, out.file = F, new.dev = F)
  expect_error(optimal_bas <- CampaRi::basin_optimization(the_sap = the_sap,  how_fine_search = 100, number_of_clusters = 3, force_matching = T, silent = silent), NA)
  expect_error(optimal_bas <- CampaRi::basin_optimization(the_sap = the_sap,  how_fine_search = 100, force_matching = T, silent = silent))
  # optimal_bas$Optimal_nbins
  # optimal_bas$bas
  
  # now testing for internal optimization
  expect_error(optimal_bas <- CampaRi::basin_optimization(the_sap = the_sap,  how_fine_search = 100, number_of_clusters = 3, force_matching = T, silent = silent), NA)
  
  if(file.exists('MST_DUMPLING.nc')) file.remove('MST_DUMPLING.nc')
  if(file.exists('REPIX_000000000021.dat')) file.remove('REPIX_000000000021.dat')
})

# require(microbenchmark); require(ggplot2); v <- 1:10e7
# a <- microbenchmark("bs" = .BiSearch(v, 432), "classic" = (432 %in% v), times = 10)
# autoplot(a)