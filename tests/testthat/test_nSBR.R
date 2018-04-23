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
  expect_error(optimal_bas <- CampaRi::nSBR(data = file.pi, ny = 30, n.cluster = 4, plot = T, silent = silent, dbg_nSBR = F), NA)
  
  
  # ------------------------------------------------------- neuro tests
  # the following tests have been executed only 
  # extensive testing on real cases. It needs not to be
  # executed at all
  do_it <- FALSE
  if(do_it){
    fpi_empty <- data.table::fread('/disk2a/dgarolini/neuro/simulated_data/analysis/ISI_networks_end17/PROGIDX_000000000100.dat', data.table = F)
    fpi_best  <- data.table::fread('/disk2a/dgarolini/neuro/simulated_data/analysis/ISI_networks_end17/PROGIDX_000000000102.dat', data.table = F)
    fpi_ave1   <- data.table::fread('/disk2a/dgarolini/neuro/simulated_data/analysis/ISI_networks_end17/PROGIDX_000000000221.dat', data.table = F)
    fpi_ave2   <- data.table::fread('/disk2a/dgarolini/neuro/simulated_data/analysis/ISI_networks_end17/PROGIDX_000000000261.dat', data.table = F)
    fpi_ave3   <- data.table::fread('/disk2a/dgarolini/neuro/simulated_data/analysis/ISI_networks_end17/PROGIDX_000000001102.dat', data.table = F)
    fpi_ave4   <- data.table::fread('/disk2a/dgarolini/neuro/simulated_data/analysis/ISI_networks_end17/PROGIDX_000000001231.dat', data.table = F)
    fpi_worst <- data.table::fread('/disk2a/dgarolini/neuro/simulated_data/analysis/ISI_networks_end17/PROGIDX_000000001281.dat', data.table = F)
  }else{
    fpi_empty <- fpi_best <- fpi_ave <- fpi_worst <- file.pi
  }
  sapphire_plot(sap_table = fpi_worst, timeline = T, sub_sampling_factor = 10)
  wrapperone <- function(x,...){ invisible(CampaRi::nSBR(data = x, n.cluster = 4, 
                                                     comb_met = c('MIC', 'MIC', 'kin'),
                                                     unif.splits = seq(5, 100, 8),  
                                                     pk_span = 5000, ny = 50, plot = T, 
                                                     silent = F, dbg_nSBR = F, return_plot = F,...))} 
  expect_error(ob1 <- wrapperone(fpi_best, data.out.it = F), NA); ob1
  expect_error(a1 <- CampaRi::nSBR(data = file.pi, n.cluster = 3, 
                                   comb_met = c('MIC'),
                                   unif.splits = seq(5, 100, 8),  
                                   pk_span = 500, ny = 50, plot = T, 
                                   silent = F, dbg_nSBR = F, return_plot = F), NA)
  expect_error(ahahscore <- CampaRi::score_sapphire(the_sap = file.pi, ann = ann, manual_barriers = a1$barriers[1:2]), NA)
  
  
  if(do_it){
    # expect_error(ahahscore <- CampaRi::score_sapphire(the_sap = fpi_best, ann = rep(1:4, each = 20000), manual_barriers = ob1$barriers[1:3]))
    expect_error(ob21 <- wrapperone(fpi_ave1), NA); ob21
    expect_error(ob22 <- wrapperone(fpi_ave2), NA); ob22
    expect_error(ob23 <- wrapperone(fpi_ave3), NA); ob23
    expect_error(ob24 <- wrapperone(fpi_ave4), NA); ob24
    expect_error(ob3 <- wrapperone(fpi_worst), NA); ob3
    
    library(cowplot)
    cowplot::plot_grid(ob2, ob1, nrow = 2, labels = c('A', 'B'))
    cowplot::plot_grid(ob2, ob1, ob3, nrow = 3, labels = c('A', 'B', 'C'))
    cwpl <- cowplot::plot_grid(ob1, ob21, ob22, ob23, ob24, ob3, nrow = 2, labels = c('A', 'B', 'C', 'D', 'E', 'F'))
    
    ggsave('p1.png', width = 6, height = 8)
    ggsave('p2.png', width = 6, height = 10)
    ggsave(plot = cwpl, filename = 'p3.png', width = 13, height = 9)
    ncl <- NULL
    ncl <- 4
    
    # fpi_empty
    expect_error(optimal_bas <- CampaRi::nSBR(the_sap = fpi_empty,  how_fine_search = 10, nSBR_method = "MI_barrier_weighting",
                                                            force_matching = T, number_of_clusters = 3, denat_opt = 'process_subtraction',
                                                            plot_basin_identification = plt_stff, silent = silent), NA)
    
    # just curiosity
    df.test <- data.frame(a = rnorm(n = 1000, mean = 0, sd = 1), b = rnorm(n = 1000, mean = 0, sd = 1), c = rnorm(n = 1000, mean = 6, sd = 1))
    ggplot(data = df.test) + theme_classic() +
      geom_freqpoly(mapping = aes(a), binwidth = 0.4, col = 'darkred') +
      geom_freqpoly(mapping = aes(b), binwidth = 0.4, col = 'darkblue') +
      geom_freqpoly(mapping = aes(c), binwidth = 0.4) 
    
    .SHEN.hist(df.test$a, 100) + .SHEN.hist(df.test$b, 100)
    .SHEN.hist(df.test$a, 100) + .SHEN.hist(df.test$c, 100)
    .SHEN.hist(c(df.test$a+2, df.test$b), 200)
    .SHEN.hist(c(df.test$a, df.test$c), 200)

    # tring to make a sense if there is 
    slide <- seq(-100, 100, 5)
    binning <- 300
    multivar.SHEN <- sapply(X = slide, FUN = function(x) .SHEN.hist(c(df.test$a, df.test$b + x), binning))
    sum.SHEN <- sapply(X = slide, FUN = function(x) .SHEN.hist(df.test$a, binning) + .SHEN.hist(df.test$b + x, binning))
    df.t <- data.frame('diff' = .normalize(sum.SHEN - multivar.SHEN), 'sSHEN' = .normalize(sum.SHEN), 'mSHEN' = .normalize(multivar.SHEN))
    ggplot(data = df.t, mapping = aes(x = slide)) + theme_classic() +
      geom_line(mapping = aes(y = diff), col = 'darkred') + 
      geom_line(mapping = aes(y = sSHEN), linetype = 'dotted') + 
      geom_line(mapping = aes(y = mSHEN)) + 
      geom_point(x = which.max(df.t$diff), y = max(df.t$diff)) + 
      geom_point(x = 0, y = 0, col = 'red')
    
    # entropy as a function of the binning
    binning <- seq(1, 1000, 10)
    a.shen.line <- sapply(X = binning, FUN = function(x) .SHEN.hist(df.test$a, x))
    b.shen.line <- sapply(X = binning, FUN = function(x) .SHEN.hist(df.test$b, x))
    c.shen.line <- sapply(X = binning, FUN = function(x) .SHEN.hist(df.test$c, x))
    df.shen.line <- cbind(a.shen.line, b.shen.line, c.shen.line)  
    ggplot(data = df.shen.line, mapping = aes(x = binning)) + theme_classic() +
      geom_line(mapping = aes(y = a.shen.line), col = 'darkred') +
      geom_line(mapping = aes(y = b.shen.line), col = 'darkblue') +
      geom_line(mapping = aes(y = c.shen.line)) 
  }
  
})
.SHEN.hist <- function(x, brks){
  hist.int <- hist(x, breaks = brks, plot=FALSE)
  dens <- hist.int$density
  return(.myShEn(dens))
}
.normalize <- function(x, xmax = NULL, xmin = NULL) {
  # if(.isSingleElement(x)) return(1)
  if(is.null(xmax)) xmax <- max(x)
  if(is.null(xmin)) xmin <- min(x)
  if(xmin == xmax) return(x/xmax)
  else return((x*1.0 - xmin)/(xmax - xmin))
}
