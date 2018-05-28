context('nSBR')

test_that('new trials for SBR', {
  
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
  expect_error(optimal_bas <- CampaRi::nSBR(data = file.pi, ny = 30, n.cluster = 4, plot = T, silent = silent, dbg_nSBR = F), NA)
  
  
  # ------------------------------------------------------- classic test
  expect_error(a1 <- CampaRi::nSBR(data = file.pi, n.cluster = 3, 
                                   comb_met = c('MIC_kin', 'MIC', 'MAS', 'MEV', 'MCN', 'MICR2', 'MI', 'kin', 'diff', 'convDiff', 'conv'),
                                   # comb_met = c('MIC_kin', 'MIC', 'MI', 'kin'),
                                   # comb_met = c('MIC', 'MICR2', 'diff'),
                                   unif.splits = seq(5, 100, 8),  
                                   pk_span = 500, ny = 50, plot = T, 
                                   silent = silent, dbg_nSBR = F, return_plot = F), NA)
  
  # ------------------------------------------------------- random test
  expect_error(a1 <- CampaRi::nSBR(data = file.pi, n.cluster = 3, 
                                   comb_met = c('MIC_kin'),
                                   unif.splits = seq(5, 100, 8),  
                                   pk_span = 500, ny = 50, plot = T, random_picks = 100, 
                                   silent = silent, dbg_nSBR = F, return_plot = F))
  expect_error(a1 <- CampaRi::nSBR(data = file.pi, n.cluster = 3, 
                                   comb_met = c('MIC_kin'),
                                   unif.splits = seq(5, 100, 8),  
                                   pk_span = 500, ny = 50, plot = T, random_picks = 100, ann = ann,
                                   silent = silent, dbg_nSBR = F, return_plot = F), NA)
  # ------------------------------------------------------- shuffles test
  expect_error(a1 <- CampaRi::nSBR(data = file.pi, n.cluster = 3, shuffles = T, 
                                   comb_met = c('MIC'),
                                   unif.splits = seq(5, 100, 8),  
                                   pk_span = 500, ny = 20, plot = T, random_picks = 100, ann = ann,
                                   silent = silent, dbg_nSBR = F, return_plot = F), NA)
  
  
  ######################### evaluating the fluctuation and randomicity of the score ##########################
  do_it <- FALSE
  if(do_it){
    # just curiosity - entropy considerations
    wrapperone <- function(x,...){ invisible(CampaRi::nSBR(data = x, n.cluster = 4, 
                                                       comb_met = c('MIC', 'MIC', 'kin'),
                                                       unif.splits = seq(5, 100, 8),  
                                                       pk_span = 5000, ny = 50, plot = T, 
                                                       silent = silent, dbg_nSBR = F, return_plot = F,...))} 
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





