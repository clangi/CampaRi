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
  # ------------------------------------------------------- force_correct_ncl test
  expect_warning(a1 <- CampaRi::nSBR(data = file.pi, n.cluster = 15, shuffles = T, 
                                   comb_met = c('MIC'), force_correct_ncl = T,
                                   unif.splits = seq(5, 100, 8),  
                                   pk_span = 500, ny = 20, plot = T, random_picks = 100, ann = ann,
                                   silent = F, dbg_nSBR = F, return_plot = F))
  # expect_warning(a1 <- CampaRi::nSBR(data = file.pi, n.cluster = 15, shuffles = T, 
  #                                    comb_met = c('MIC'), force_correct_ncl = F,
  #                                    unif.splits = seq(5, 100, 8),  
  #                                    pk_span = 500, ny = 20, plot = T, random_picks = 100, ann = ann,
  #                                    silent = F, dbg_nSBR = F, return_plot = F))
  
  do_it <- FALSE
  if(do_it){
    require(ggplot2)
  ######################### evaluating an automatic way to select the number of cluster -> grid search ##########################
    
    
    # initialization of the vars 
    ncl.vec <- 2:6
    sco.mat <- data.frame(array(NA, dim = c(length(ncl.vec), length(ncl.vec))))
    sco.rnd.mat <- data.frame(array(NA, dim = c(length(ncl.vec), length(ncl.vec))))
    h.split.mat <- data.frame(array(NA, dim = c(length(ncl.vec), length(ncl.vec))))
    span.mat <- data.frame(array(NA, dim = c(length(ncl.vec), length(ncl.vec))))
    ny.mat <- data.frame(array(NA, dim = c(length(ncl.vec), length(ncl.vec))))
    rownames(sco.mat) <- rownames(span.mat) <- rownames(sco.rnd.mat) <- rownames(h.split.mat) <- rownames(ny.mat) <- paste0('ncl.', ncl.vec, '.i')
    colnames(sco.mat) <- colnames(span.mat) <- colnames(sco.rnd.mat) <- colnames(ny.mat) <- colnames(h.split.mat) <- paste0('ncl.', ncl.vec, '.j')
    
    # main loop over clusters I ----------------------------------------------------------------------------------
    for(ncl.i in ncl.vec){
      cat('NCI: Doing', ncl.i, '(', which(ncl.vec == ncl.i), '/', length(ncl.vec), ')\n')
      
      # first call for barriers
      temp1 <- repeat.until.cl(ncl.x = ncl.i,  max_rounds = 15, plot = F,
                               h.splits.t = 10, span.t = 500, ny.t = 30, silent = T)
      
      if(temp1$found) cat('Found the following barriers:', temp1$tmp$barriers, '\n')
      else stop('No correct number!!')
      
      # extraction of the annotation from barriers (numbering not relevant)
      ref.vec <- vec.from.barriers(bar.vec = temp1$tmp$barriers[1:(ncl.i-1)])
      
      # internal loop over clusters J ---------------------------------------------------------------------------
      for(ncl.j in ncl.vec){
        cat('NCJ: Doing', ncl.j, '(', which(ncl.vec == ncl.j), '/', length(ncl.vec), ')\n')
        
        # second call for barriers
        temp2 <- repeat.until.cl(ncl.x = ncl.j,  max_rounds = 15, plot = F,
                                 h.splits.t = 10, span.t = 500, ny.t = 30, silent = T)
        
        if(temp2$found) cat('Found the following barriers:', temp2$tmp$barriers, '\n')
        else stop('No correct number!!')
        
        # scoring
        ref_sc <- CampaRi::score_sapphire(the_sap = file.pi, ann = ref.vec, 
                                          manual_barriers = temp2$tmp$barriers[1:(ncl.j-1)], silent = T)
        
        # final assignment
        h.split.mat[paste0('ncl.', ncl.i, '.i'), paste0('ncl.', ncl.j, '.j')] <- temp2$h.splits.t
        span.mat[paste0('ncl.', ncl.i, '.i'), paste0('ncl.', ncl.j, '.j')] <- temp2$span.t
        ny.mat[paste0('ncl.', ncl.i, '.i'), paste0('ncl.', ncl.j, '.j')] <- temp2$ny.t
        sco.mat[paste0('ncl.', ncl.i, '.i'), paste0('ncl.', ncl.j, '.j')] <- ref_sc$score.out
        # sco.rnd.mat[paste0('ncl.', ncl.i), paste0('ncl.', ncl.j)] <- res2
      }
    }
    # apply(sco.mat, 2, mean) 
    # apply(sco.mat, 1, mean) 
    
    # Plotting results for visual inspection
    ggplot_tiles(mat = sco.mat, labx = 'Clusters', laby = 'Clusters', legtit = 'ARI') + geom_text(label = round(unlist(c(sco.mat)),2), col = 'white')
    ggplot_tiles(mat = h.split.mat, labx = 'Clusters', laby = 'Clusters', legtit = 'ARI') + geom_text(label = round(unlist(c(h.split.mat)),2), col = 'white')
    ggplot_tiles(mat = span.mat, labx = 'Clusters', laby = 'Clusters', legtit = 'ARI') + geom_text(label = round(unlist(c(span.mat)),2), col = 'white')
    ggplot_tiles(mat = ny.mat, labx = 'Clusters', laby = 'Clusters', legtit = 'ARI') + geom_text(label = round(unlist(c(ny.mat)),2), col = 'white')
    
    # best scenario / selection of the number of clusters
    CampaRi::nSBR(data = file.pi, n.cluster = 5,
                  comb_met = c('MIC'),
                  unif.splits = unique(round(seq(5, 100, length.out = 40))),  
                  pk_span = 350, ny = 60, plot = T,
                  silent = F, dbg_nSBR = F, return_plot = F)
    
    # --------------------------- functions ---------------------------------
    
    # tile plot from matrix
    ggplot_tiles <- function(mat, labx = '', laby = '', legtit = '', normalize = FALSE, return_plot = FALSE){
      
      stopifnot(is.character(labx), is.character(laby), is.character(legtit))
      stopifnot(is.logical(normalize), is.logical(return_plot))
      
      melted <- reshape2::melt(as.matrix(mat))
      if(normalize) melted[,3] <- .normalize(melted[,3])
      
      itplots <- ggplot(melted, aes(x = Var2, y = Var1)) + theme_minimal() +
                      geom_raster(aes(fill=value)) + 
        scale_fill_gradientn(legtit, colours = c('black', 'darkred')) + 
        xlab(labx) + ylab(laby) 
      
      if(return_plot) return(itplots)
      else print(itplots)
    }
    
    # create annotation vector from barrier vector
    vec.from.barriers <- function(bar.vec, end.point = 3000){
      
      stopifnot(all(end.point > bar.vec))
      
      bar.vec <- sort(c(bar.vec, end.point))
      bar.vec <- c(bar.vec[1], diff(bar.vec))
      
      return(c(unlist(sapply(1:length(bar.vec), function(x) rep(x, bar.vec[x])))))
    }
    
    # repear until the right number of clusters is found
    repeat.until.cl <- function(ncl.x, max_rounds = 10, h.splits.t = 20,
                                span.t = 500, ny.t = 30, silent = F, plot = F){
      round <- 1
      repeat {
        sspplliitts <- unique(round(seq(5, 100, length.out = h.splits.t)))
        tmp <- CampaRi::nSBR(data = file.pi, n.cluster = ncl.x,
                                                  comb_met = c('MIC'),
                                                  unif.splits = sspplliitts,  
                                                  pk_span = span.t, ny = ny.t, plot = plot,
                                                  silent = silent, dbg_nSBR = F, return_plot = F)
        if(length(tmp$barriers) > (ncl.x - 2)){
          if(!silent) cat('Found', length(tmp$barriers), 'barriers over', ncl.x - 1, '\n')
          fnd <- TRUE
          break
        }else{
          if(!silent) cat('Found', length(tmp$barriers), 'barriers instead of', ncl.x - 1, '\n')
          h.splits.t <- h.splits.t + 10
          span.t <- max(100, span.t - 50)
          ny.t <- ny.t + 10
          round <- round + 1
          if(round > max_rounds){
            fnd <- FALSE
            break
          }
        }
      }
      return(list('found' = fnd, 'tmp' = tmp, 'h.splits.t' = h.splits.t, 'span.t' = span.t, 'ny.t' = ny.t, 'round' = round))
    }
  ######################### evaluating the fluctuation and randomicity of the score ##########################
    # just curiosity - entropy considerations
    df.test <- data.frame(a = rnorm(n = 1000, mean = 0, sd = 1), b = rnorm(n = 1000, mean = 0, sd = 1), c = rnorm(n = 1000, mean = 6, sd = 1))
    ggplot(data = df.test) + theme_classic() +
      geom_freqpoly(mapping = aes(a), binwidth = 0.4, col = 'darkred') +
      geom_freqpoly(mapping = aes(b), binwidth = 0.4, col = 'darkblue') +
      geom_freqpoly(mapping = aes(c), binwidth = 0.4) 
    
    .SHEN.hist(df.test$a, 100); .SHEN.hist(df.test$b, 100)
    .SHEN.hist(df.test$a, 100) + .SHEN.hist(df.test$b, 100)
    .SHEN.hist(df.test$a, 100) + .SHEN.hist(df.test$c, 100)
    .SHEN.hist(c(df.test$a, df.test$b), 200)
    .SHEN.hist(c(df.test$a, df.test$c), 200)

    # tring to make a sense if there is 
    slide <- seq(-100, 100, 5)
    binning <- 300
    multivar.SHEN <- sapply(X = slide, FUN = function(x) .SHEN.hist(c(df.test$a, df.test$b + x), binning))
    sum.SHEN <- sapply(X = slide, FUN = function(x) .SHEN.hist(df.test$a, binning) + .SHEN.hist(df.test$b + x, binning))
    df.t <- data.frame('diff' = .normalize(sum.SHEN - multivar.SHEN), 'sSHEN' = .normalize(sum.SHEN), 'mSHEN' = .normalize(multivar.SHEN))
    ggplot(data = df.t, mapping = aes(x = slide)) + theme_classic() +
      geom_line(mapping = aes(y = diff, col = 'diff')) + 
      geom_line(mapping = aes(y = sSHEN, col = 'Entropy sum'), linetype = 'dotted', size = 1) + 
      geom_line(mapping = aes(y = mSHEN, col = 'Dist sum')) + 
      geom_point(x = which.max(df.t$diff), y = max(df.t$diff)) + 
      geom_point(x = 0, y = 0, col = 'red')
    
    # entropy as a function of the binning
    binning <- seq(1, 1000, 10)
    a.shen.line <- sapply(X = binning, FUN = function(x) .SHEN.hist(df.test$a, x))
    b.shen.line <- sapply(X = binning, FUN = function(x) .SHEN.hist(df.test$b, x))
    c.shen.line <- sapply(X = binning, FUN = function(x) .SHEN.hist(df.test$c, x))
    df.shen.line <- data.frame(a.shen.line, b.shen.line, c.shen.line)  
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
.myShEn <- function(x){
  x2 <- replace(x = x, list = which(x == 0), values = 1)
  return(-sum(x2*log(x2)))
}
.normalize <- function(x, xmax = NULL, xmin = NULL) {
  # if(.isSingleElement(x)) return(1)
  if(is.null(xmax)) xmax <- max(x)
  if(is.null(xmin)) xmin <- min(x)
  if(xmin == xmax) return(x/xmax)
  else return((x*1.0 - xmin)/(xmax - xmin))
}





