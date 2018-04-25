context('basin_optimization')

test_that('optimize sapphire plot clusters using SBR', {
  
  # CREATION OF THE DATASET
  silent <- T
  plt_stff <- !silent
  if(!silent) {require(testthat); require(CampaRi)} 
  stdd <- 3; n_dim <- 10; n_snap <- 3000; n_tot <- n_dim*n_snap/3; if(!silent) print(n_tot)
  ann <- c(rep(1, n_snap/6), rep(2, n_snap/3), rep(1, n_snap/6), rep(3, n_snap/3))
  file.pi <- system.file("extdata", "REPIX_000000000021.dat", package = "CampaRi")
  if(plt_stff) sapphire_plot(sap_file = file.pi, ann_trace = ann, only_timeline=F, timeline = T)
  
  # BASINS OPTIMIZATION
  nbin <- round(sqrt(n_snap*10)); if(!silent) print(nbin)
  if(plt_stff) basins_recognition(data = file.pi, nx = nbin, match = F, plot = T, out.file = F, new.dev = F)
  if(plt_stff) basins_recognition(data = file.pi, nx = nbin, match = T, plot = T, out.file = F, new.dev = F)
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
  
  # ------------------------------------------------------- neuro tests
  # the following tests have been executed only 
  # extensive testing on real cases. It needs not to be
  # executed at all
  do_it <- FALSE
  
  expect_error(optimal_bas <- CampaRi::basin_optimization(the_sap = file.pi,  how_fine_search = 10, basin_optimization_method = "MI_barrier_weighting",
                                                          force_matching = T, number_of_clusters = 3, denat_opt = 'process_subtraction',
                                                          plot_basin_identification = plt_stff, silent = silent), NA)
  
  if(do_it){
    ncl <- 4
    
    # fpi_empty
    expect_error(optimal_bas <- CampaRi::basin_optimization(the_sap = fpi_empty,  how_fine_search = 10, basin_optimization_method = "MI_barrier_weighting",
                                                            force_matching = T, number_of_clusters = 3, denat_opt = 'process_subtraction',
                                                            plot_basin_identification = plt_stff, silent = silent), NA)
    
    
    # fpi_best - explicative  
    optimal_bas <- CampaRi::basin_optimization(the_sap = fpi_best,                                   # PI data
                                               how_fine_search = 20,                                 # number of bins searched between nbins_x_min and nbins_x_max
                                               nbins_x_min = 7,                                      # minimum number of bins (range)
                                               nbins_x_max = 200,                                    # maximum number of bins (range)
                                               basin_optimization_method = "MI_barrier_weighting",   # ranking method
                                               cl.stat.MI_comb = 'kin_MI',                           # Final score combination 
                                               cl.stat.nUni = c(5,10,15,20,25,30,40,50,60),          # MI curve splits 
                                               force_matching = T,                                   # match kin and tempann for the barriers
                                               number_of_clusters = ncl,                             # number of cluster approach. NULL in this case it must be
                                               denat_opt = 'process_subtraction',                    # kin ann is corrected for parabolic artefacts
                                               cl.stat.denat.MI = 7,                                 # if a number also the MI curve is corrected
                                               plot_basin_identification = plt_stff,                 # final plot?
                                               # dbg_basin_optimization = T,                          # debug?
                                               silent = silent)                                      # silent?
    
    # fpi_best  
    optimal_bas <- CampaRi::basin_optimization(the_sap = fpi_best,  how_fine_search = 20, nbins_x_min = 50, nbins_x_max = 200,
                                               basin_optimization_method = "MI_barrier_weighting",  cl.stat.MI_comb = 'kin_MI',
                                               cl.stat.nUni = c(5,10,15,20,25),
                                               force_matching = T, number_of_clusters = ncl, denat_opt = 'process_subtraction', 
                                               cl.stat.denat.MI = -1, # it is dampening the score of the borders (parabolic artefacts)
                                               plot_basin_identification = plt_stff, silent = silent)
    # fpi_ave  
    optimal_bas <- CampaRi::basin_optimization(the_sap = fpi_ave,  how_fine_search = 10, basin_optimization_method = "MI_barrier_weighting", 
                                               force_matching = T, number_of_clusters = 4, nbins_x_min = 7, nbins_x_max = 200,
                                               denat_opt = 'process_subtraction', cl.stat.MI_comb = 'kin_MI',
                                               cl.stat.nUni = c(5,10,15,20,25,30,40,50,60), cl.stat.denat.MI = NULL,
                                               plot_basin_identification = plt_stff, silent = silent, dbg_basin_optimization =F)
    # fpi_worst  
    optimal_bas <- CampaRi::basin_optimization(the_sap = fpi_worst,  how_fine_search = 20, basin_optimization_method = "MI_barrier_weighting", 
                                               force_matching = T, number_of_clusters = 4, nbins_x_min = 7, nbins_x_max = 200,
                                               denat_opt = 'process_subtraction', cl.stat.MI_comb = 'kin_MI',
                                               cl.stat.nUni = seq(20,200,10), cl.stat.denat.MI = NULL,
                                               plot_basin_identification = plt_stff, silent = silent, dbg_basin_optimization = F)
  }
})
# trials for bisearch vs classic binary search (%in%)
# require(microbenchmark); require(ggplot2); v <- 1:10e7
# a <- microbenchmark("bs" = .BiSearch(v, 432), "classic" = (432 %in% v), times = 10)
# autoplot(a)