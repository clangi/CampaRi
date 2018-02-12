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
  
  the_sap <- 'REPIX_000000000021.dat'
  if(plt_stff) sapphire_plot(sap_file = the_sap, ann_trace = ann, timeline = T); cat(round(sqrt(nrow(thedata)*10)))
  nbin <- round(sqrt(nrow(thedata)*10)); nbin
  nbin <- 20
  if(plt_stff) basins_recognition(the_sap, nx = nbin, match = F, plot = T, out.file = F, new.dev = F)
  # if(plt_stff) basins_recognition(the_sap, nx = nbin, match = T, plot = T, out.file = F, new.dev = F)
  CampaRi::basin_optimization()
  
  
  abba <- unique(round(seq(2, nbinsxy, length.out = 100)))
  .BiSBR(st = st, ncl_found = 1, ncl_teo = 3, start.idx = 1, end.idx = length(abba), lin_scale =  abba, force_matching = F, silent = F)
  .BiSBR(st = st, ncl_found = 1, ncl_teo = 3, start.idx = 1, end.idx = length(abba), lin_scale =  abba, force_matching = T, silent = F)
  # testing for minimum number of bins without crashing
  CampaRi::basins_recognition(st[1:1000,], nx = 4, new.dev = F, out.file = F, plot = T)
  CampaRi::basins_recognition(st, nx = 12, new.dev = F, out.file = F, plot = T)
  CampaRi::basins_recognition(st, nx = 88, new.dev = F, out.file = F, plot = T, match = T)
  CampaRi::basins_recognition(st, nx = 97, new.dev = F, out.file = F, plot = T, match = T)
  CampaRi::basins_recognition(st, nx = 343, new.dev = F, out.file = F, plot = T, match = T)
  alu <- CampaRi::basins_recognition(st, nx = 905, new.dev = F, out.file = F, plot = T, match = T)
  
  if(file.exists('MST_DUMPLING.nc')) file.remove('MST_DUMPLING.nc')
  if(file.exists('REPIX_000000000021.dat')) file.remove('REPIX_000000000021.dat')
})


# require(microbenchmark); require(ggplot2); v <- 1:10e7
# a <- microbenchmark("bs" = .BiSearch(v, 432), "classic" = (432 %in% v), times = 10)
# autoplot(a)