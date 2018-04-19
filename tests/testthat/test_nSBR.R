context('nSBR')

test_that('new trials for SBR', {
  
  # CREATION OF THE DATASET
  silent <- F
  plt_stff <- !silent
  if(!silent) {require(testthat); require(CampaRi)} 
  stdd <- 3; n_dim <- 10; n_snap <- 3000; n_tot <- n_dim*n_snap/3; if(!silent) print(n_tot)
  ann <- c(rep(1, n_snap/6), rep(2, n_snap/3), rep(1, n_snap/6), rep(3, n_snap/3))
  file.pi <- system.file("extdata", "REPIX_000000000021.dat", package = "CampaRi")
  if(plt_stff) sapphire_plot(sap_file = file.pi, ann_trace = ann, only_timeline=F, timeline = T)
  
  # BASINS OPTIMIZATION
  nbin <- round(sqrt(n_snap*10)); if(!silent) print(nbin)
  expect_error(optimal_bas <- CampaRi::nSBR(data = file.pi, nx = 100, ny.aut = T, plot = T, silent = silent, dbg_nSBR = T), NA)
  
  
  # ------------------------------------------------------- neuro tests
  # the following tests have been executed only 
  # extensive testing on real cases. It needs not to be
  # executed at all
  do_it <- FALSE
  if(do_it){
    fpi_empty <- data.table::fread('/disk2a/dgarolini/neuro/simulated_data/analysis/ISI_networks_end17/PROGIDX_000000000100.dat', data.table = F)
    fpi_best  <- data.table::fread('/disk2a/dgarolini/neuro/simulated_data/analysis/ISI_networks_end17/PROGIDX_000000000102.dat', data.table = F)
    fpi_ave   <- data.table::fread('/disk2a/dgarolini/neuro/simulated_data/analysis/ISI_networks_end17/PROGIDX_000000001102.dat', data.table = F)
    fpi_worst <- data.table::fread('/disk2a/dgarolini/neuro/simulated_data/analysis/ISI_networks_end17/PROGIDX_000000001281.dat', data.table = F)
  }else{
    fpi_empty <- fpi_best <- fpi_ave <- fpi_worst <- file.pi
  }
  
  expect_error(optimal_bas <- CampaRi::nSBR(the_sap = file.pi,  how_fine_search = 10, nSBR_method = "MI_barrier_weighting",
                                                          force_matching = T, number_of_clusters = 3, denat_opt = 'process_subtraction',
                                                          plot_basin_identification = plt_stff, silent = silent), NA)
  
  
  if(do_it){
    ncl <- NULL
    ncl <- 4
    
    # fpi_empty
    expect_error(optimal_bas <- CampaRi::nSBR(the_sap = fpi_empty,  how_fine_search = 10, nSBR_method = "MI_barrier_weighting",
                                                            force_matching = T, number_of_clusters = 3, denat_opt = 'process_subtraction',
                                                            plot_basin_identification = plt_stff, silent = silent), NA)
    
  }
})
