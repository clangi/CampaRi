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
  # print(qres)
  
  thedata <- matrix(c(rnorm(10000, sd = 1), rnorm(10000, mean = 5, sd = 1), rnorm(10000, mean = 10, sd = 1)), nrow = 3000, ncol = 10)
  expect_error(adjl <- mst_from_trj(trj = thedata, dump_to_netcdf = FALSE), NA)
  expect_error(ret <- gen_progindex(adjl, snap_start = 21), NA)
  expect_error(ret2 <- gen_annotation(ret, snap_start = 21), NA)
  the_sap <- 'REPIX_000000000021.dat'
  # sapphire_plot(sap_file = the_sap, ann_trace = annotate_it)
  annotate_it <- c(rep(1, 1000), rep(2, 1000), rep(3, 1000))
  expect_error(qres <- score_sapphire(the_sap = the_sap, ann = annotate_it, silent = silent), NA)
  # qres
  
  if(file.exists('MST_DUMPLING.nc')) file.remove('MST_DUMPLING.nc')
  if(file.exists('REPIX_000000000021.dat')) file.remove('REPIX_000000000021.dat')
})
