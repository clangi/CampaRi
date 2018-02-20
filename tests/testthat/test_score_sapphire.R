context('score_sapphire')

test_that('scoring sapphire plots', {
  silent <- T
  plt_stff <- !silent
  if(!silent) require(testthat)
  # thedata <- matrix(c(rnorm(10000, sd = 1), rnorm(10000, mean = 5, sd = 1), rnorm(10000, mean = 10, sd = 1)), nrow = 3000, ncol = 10)
  # if(plt_stff) plot(c(thedata))
  
  # expect_error(adjl <- mst_from_trj(trj = thedata, dump_to_netcdf = FALSE, mute_fortran = silent), NA)
  # expect_error(ret <- gen_progindex(adjl, snap_start = 21, mute = silent), NA)
  # expect_error(ret2 <- gen_annotation(ret, snap_start = 21, mute = silent), NA)

  stdd <- 3; n_dim <- 10; n_snap <- 3000; n_tot <- n_dim*n_snap/3; if(!silent) print(n_tot)
  ann <- c(rep(1, n_snap/6), rep(2, n_snap/3), rep(1, n_snap/6), rep(3, n_snap/3))
  file.pi <- system.file("extdata", "REPIX_000000000021.dat", package = "CampaRi")
  if(plt_stff) sapphire_plot(sap_file = file.pi, ann_trace = ann, timeline = T)
  # cat(round(sqrt(nrow(thedata)*10)))
  nbin <- round(sqrt(3000*10))
  # nbin <- 20
  if(plt_stff) basins_recognition(file.pi, nx = nbin, match = F, plot = T, out.file = F, new.dev = F)

  # blind run with multiple barriers unbound and with merging (in this way it is supervised!)
  expect_error(qres <- score_sapphire(the_sap = file.pi, ann = ann, plot_basin_identification = plt_stff, plot_pred_true_resume = plt_stff,
                                      merge_clusters = T, silent = silent), NA)
  
  # blind run without optimization and 173 bins and force_matching (no barrier will be found)
  expect_error(qres <- score_sapphire(the_sap = file.pi, ann = ann, force_matching = T, plot_basin_identification = plt_stff, 
                                      plot_pred_true_resume = plt_stff, merge_clusters = T, silent = silent), NA)
  
  # fixed number of clusters
  expect_error(qres <- score_sapphire(the_sap = file.pi, ann = ann, silent = silent, basin_optimization = T, plot_pred_true_resume = plt_stff,
                                      number_of_clusters = 3, plot_basin_identification = plt_stff), NA)
  
  # run with force_matching
  expect_error(qres <- score_sapphire(the_sap = file.pi, ann = ann, silent = silent, basin_optimization = T, plot_pred_true_resume = plt_stff,
                                      force_matching = T, number_of_clusters = 3, plot_basin_identification = plt_stff), NA)
  
  # MI - ratio based optimization without knowing the number of clusters - !long
  # expect_error(qres <- score_sapphire(the_sap = file.pi, ann = ann, silent = silent, basin_optimization = T, plot_pred_true_resume = plt_stff,
  #                                     force_matching = T,  basin_optimization_method = "MI_barrier_weighting", plot_basin_identification = plt_stff), NA)
  
  # not silent run
  expect_error(qres <- score_sapphire(the_sap = file.pi, ann = ann, silent = F, basin_optimization = T, plot_pred_true_resume = plt_stff,
                                      force_matching = T, number_of_clusters = 3, plot_basin_identification = plt_stff), NA)
  
  # if(file.exists('MST_DUMPLING.nc')) file.remove('MST_DUMPLING.nc')
  # if(file.exists('REPIX_000000000021.dat')) file.remove('REPIX_000000000021.dat')
})
