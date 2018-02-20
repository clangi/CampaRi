context('basin_optimization')

test_that('optimize sapphire plot clusters using SBR', {
  
  # CREATION OF THE DATASET
  silent <- T
  plt_stff <- !silent
  if(!silent) {require(testthat); require(CampaRi)} 
  # data generation is now already made and put in the package
  # myNormalize <- function(x) return((x-min(x))/(max(x)-min(x)))
  # stdd <- 0.1; n_dim <- 20; n_snap <- 6000; n_tot <- n_dim*n_snap/4; if(!silent) print(n_tot)
  # thedata <- matrix(c(rnorm(n_tot, sd = stdd), rnorm(n_tot, mean = 5, sd = stdd), rnorm(n_tot, mean = 10, sd = stdd),
  #                     rnorm(n_tot, sd = stdd)), nrow = n_snap, ncol = n_dim)
  # stdd <- 1; n_dim <- 2; n_snap <- 3000; n_tot <- n_dim*n_snap/6; if(!silent) print(n_tot)
  # thedata <- matrix(c(rnorm(n_tot, sd = stdd), rnorm(n_tot, mean = 5, sd = stdd), rnorm(n_tot, mean = 10, sd = stdd),
  #                     rnorm(n_tot, sd = stdd), rnorm(n_tot, mean = 15, sd=stdd), rnorm(n_tot, mean = 5, sd=stdd)), nrow = n_snap, ncol = n_dim)
  stdd <- 3; n_dim <- 10; n_snap <- 3000; n_tot <- n_dim*n_snap/3; if(!silent) print(n_tot)
  # thedata <- rbind(matrix(rnorm(n_tot/2, mean = 0, sd = stdd), nrow = n_snap/6, ncol = n_dim),
  #                  matrix(rnorm(n_tot, mean = 5, sd = stdd), nrow = n_snap/3, ncol = n_dim),
  #                  matrix(rnorm(n_tot/2, mean = 0, sd = stdd), nrow = n_snap/6, ncol = n_dim),
  #                  matrix(rnorm(n_tot/2, mean = 10, sd = stdd), nrow = n_snap/3, ncol = n_dim))
  # thedata <- myNormalize(thedata)
  # if(plt_stff) plot(apply(thedata, 1, mean))
  # if(plt_stff) plot(apply(thedata, 2, mean)) # this must be similar! It is similar for each dimension
  # if(plt_stff) plot(c(thedata))
  
  # adjl <- mst_from_trj(trj = thedata, normalize_d = T, dump_to_netcdf = FALSE, mute_fortran = silent)
  # ret <- gen_progindex(adjl, snap_start = 21, mute = silent)
  # ret2 <- gen_annotation(ret, snap_start = 21, mute = silent)
  
  # ann <- c(rep(1, n_snap/4), rep(2, n_snap/4), rep(3, n_snap/4),  rep(1, n_snap/4))
  # ann <- c(rep(1, n_snap/6), rep(2, n_snap/6), rep(3, n_snap/6),  rep(1, n_snap/6),  rep(4, n_snap/6),  rep(2, n_snap/6))
  # ann <- c(rep(1, n_snap/3), rep(2, n_snap/3), rep(3, n_snap/3))
  ann <- c(rep(1, n_snap/6), rep(2, n_snap/3), rep(1, n_snap/6), rep(3, n_snap/3))
  file.pi <- system.file("extdata", "REPIX_000000000021.dat", package = "CampaRi")
  if(plt_stff) sapphire_plot(sap_file = file.pi, ann_trace = ann, only_timeline=F, timeline = T)
  
  # BASINS OPTIMIZATION
  nbin <- round(sqrt(n_snap*10)); if(!silent) print(nbin)
  if(plt_stff) basins_recognition(data = file.pi, nx = nbin, match = F, plot = T, out.file = F, new.dev = F)
  # if(plt_stff) basins_recognition(file.pi, nx = nbin, match = T, plot = T, out.file = F, new.dev = F)
  expect_error(optimal_bas <- CampaRi::basin_optimization(the_sap = file.pi,  how_fine_search = 100, number_of_clusters = 3, force_matching = F, silent = F, 
                                                          plot_basin_identification = plt_stff), NA)
  expect_error(optimal_bas <- CampaRi::basin_optimization(the_sap = file.pi,  how_fine_search = 100, force_matching = T, silent = silent))
  # optimal_bas$Optimal_nbins
  # optimal_bas$bas
  
  # now testing for internal optimization
  expect_error(optimal_bas <- CampaRi::basin_optimization(the_sap = file.pi,  how_fine_search = 10, basin_optimization_method = "MI_barrier_weighting", force_matching = T, 
                                                          plot_basin_identification = plt_stff, silent = silent), NA)
  expect_error(optimal_bas <- CampaRi::basin_optimization(the_sap = file.pi,  how_fine_search = 10, basin_optimization_method = "MI_barrier_weighting",
                                                          number_of_clusters = 3, force_matching = T, 
                                                          plot_basin_identification = plt_stff, silent = silent), NA)
  
  # if(file.exists('MST_DUMPLING.nc')) file.remove('MST_DUMPLING.nc')
  # if(file.exists('REPIX_000000000021.dat')) file.remove('REPIX_000000000021.dat')
})

# require(microbenchmark); require(ggplot2); v <- 1:10e7
# a <- microbenchmark("bs" = .BiSearch(v, 432), "classic" = (432 %in% v), times = 10)
# autoplot(a)