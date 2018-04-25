context('score_sapphire')

test_that('scoring sapphire plots', {
  silent <- F
  plt_stff <- !silent
  if(!silent) {require(testthat); require(CampaRi)} 
  
  stdd <- 3; n_dim <- 10; n_snap <- 3000; n_tot <- n_dim*n_snap/3; if(!silent) print(n_tot)
  ann <- c(rep(1, n_snap/6), rep(2, n_snap/3), rep(1, n_snap/6), rep(3, n_snap/3))
  file.pi <- system.file("extdata", "REPIX_000000000021.dat", package = "CampaRi")
  if(plt_stff) sapphire_plot(sap_file = file.pi, ann_trace = ann, only_timeline=F, timeline = T)

  # using directly the basin recognition with a certain nbins <- 50
  expect_error(qres <- score_sapphire(the_sap = file.pi, ann = ann, 
                                      force_matching = T, nbins_x = 50,
                                      plot_basin_identification = plt_stff, 
                                      silent = silent), NA)
  expect_error(qres <- score_sapphire(the_sap = file.pi, ann = ann, basin_optimization = T,
                                      force_matching = T, nbins_x = 50, how_fine_search = 10,
                                      plot_basin_identification = plt_stff, 
                                      silent = silent), NA)
  expect_error(qres <- score_sapphire(the_sap = file.pi, ann = ann, basin_optimization = T,
                                      force_matching = F,
                                      plot_basin_identification = plt_stff, 
                                      cl.stat.MI_comb='MI', # can be an option
                                      nbins_x_max = 30, 
                                      # how_fine_search = 50,
                                      how_fine_search = 5, merge_clusters = T,
                                      denat_opt = 'process_subtraction',  
                                      # number_of_clusters = 5,
                                      silent = silent), NA)
  
  
  optimal_bas <- CampaRi::basin_optimization(the_sap = file.pi,                                    # PI data
                                             how_fine_search = 5,                                  # number of bins searched between nbins_x_min and nbins_x_max
                                             nbins_x_min = 20,                                     # minimum number of bins (range)
                                             nbins_x_max = 30,                                     # maximum number of bins (range)
                                             basin_optimization_method = "MI_barrier_weighting",   # ranking method
                                             cl.stat.MI_comb = 'kin_MI',                           # Final score combination 
                                             cl.stat.nUni = c(5,10,15,20,25,30,40,50,60),          # MI curve splits 
                                             force_matching = T,                                   # match kin and tempann for the barriers
                                             number_of_clusters = 3,                               # number of cluster approach. NULL in this case it must be
                                             denat_opt = 'poly_interpolation',                     # kin ann is corrected for parabolic artefacts
                                             cl.stat.denat.MI = 7,                                 # if a number also the MI curve is corrected
                                             plot_basin_identification = plt_stff,                 # final plot?
                                             dbg_basin_optimization = F,                           # debug?
                                             silent = silent)                                      # silent?  
  
  
  expect_error(qres <- score_sapphire(the_sap = file.pi, ann = ann, 
                                      force_matching = T, nbins_x = 50, basin_optimization = optimal_bas$bas,
                                      plot_basin_identification = plt_stff, 
                                      silent = silent))
  expect_error(qres <- score_sapphire(ann = ann,
                                      force_matching = T, nbins_x = 50, basin_optimization = optimal_bas$bas,
                                      plot_basin_identification = plt_stff,
                                      silent = silent)) # error! to fix tomorrow

  
  # new tests using the nSBR function
  expect_error(a1 <- CampaRi::nSBR(data = file.pi, n.cluster = 3, 
                                  comb_met = c('MIC'),
                                  unif.splits = seq(5, 100, 8),  
                                  pk_span = 500, ny = 50, plot = T, 
                                  silent = silent, dbg_nSBR = F, return_plot = T), NA)
  expect_error(ahahscore <- CampaRi::score_sapphire(the_sap = file.pi, ann = ann, manual_barriers = a1$barriers[1:2], silent = silent), NA)
})
