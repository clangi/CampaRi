context('generate_network')

test_that('pre-processing with network inference', {
  test_plotting <- F
  my_libs <- F
  silent <- F
  if(my_libs) library(TSA)
  trj <- as.data.frame(matrix(rnorm(1000), nrow = 100, ncol = 10))
  expect_true(!is.null(trj))
  
  # checking some (expected) errors
  # -----------------------------------
  expect_error(a <- generate_network(trj, window = 20, method = "sadasd"))
  expect_error(a <- generate_network(trj, window = 20, post_processing_method = "asdsads"))
  # testing for distances and postprocessing
  # -----------------------------------
  expect_error(net_joint_tr <- generate_network(trj, method = 'minkowski', post_processing_method = 'svd', window = 12), NA)
  expect_error(net1 <- generate_network(trj, method = 'minkowski', window = 20), NA)
  expect_error(path_net <- generate_network(trj, method = 'minkowski', post_processing_method = 'path_minkowski', window = 15), NA)
  expect_error(SVD_net <- generate_network(trj, method = 'minkowski', post_processing_method = 'svd', window = 20), NA)
  # expect_error(MI_net <- generate_network(trj, method = 'MI', post_processing_method = "SymmetricUncertainty", window = 7), NA)
  
  # testing for tsne
  # -----------------------------------
  expect_error(tsne_net <- generate_network(trj, method = 'minkowski', post_processing_method = 'tsne', window = 20, tsne_pearagd = 3))
  # expect_error(tsne_net <- generate_network(trj, method = 'minkowski', post_processing_method = 'tsne', window = 20, tsne_perplexity = 3), NA)
  expect_error(tsne_net <- generate_network(trj, method = 'none', post_processing_method = 'tsne', window = 20), NA)
  
  # testing multiplication
  # -----------------------------------
  expect_error(a <- generate_network(trj, window = 20, method = "sadasd"))
  
  
  # testing for spectrum transformation
  # -----------------------------------
  trj <- as.data.frame(matrix(NA, nrow = 120, ncol = 2))
  trj[,1] <- sin(1:120)*(1/rep(1:30, 4) + tan(120:1)*(1/100))
  trj[,2] <- cos(1:120)*(1/rep(1:30, 4) + tanh(120:1)*(1/sin(120:1)))
  if(test_plotting) plot(trj[,1], type = 'l')
  if(test_plotting) plot(trj[,2], type = 'l')
  if(my_libs) a <- periodogram(y = trj[ ,1], plot = test_plotting)
  if(my_libs) b <- periodogram(y = trj[ ,2], plot = test_plotting)
  if(test_plotting) plot(y = b$spec, x = b$freq,type = "h") # spec is what it was needed

  expect_error(asd <- generate_network(trj, method = 'fft', window = 10), NA)
  if(my_libs) c <- periodogram(y = c(trj[1:5, 1], trj[1:5, 1]), plot = test_plotting)
  if(my_libs) d <- periodogram(y = c(trj[1:5, 2], trj[1:5, 2]), plot = test_plotting)
  # final check of equality
  if(my_libs) expect_true(all(c$spec == asd[1, 1:5]))
  if(my_libs) expect_true(all(d$spec == asd[1, 6:10]))
  
  if(test_plotting) plot(y = c$spec, x = c$freq, type = 'h')
  if(test_plotting) plot(y = asd[1, 1:5], x = c$freq, type = 'h')
  if(test_plotting) plot(y = d$spec, x = d$freq, type = 'h')
  if(test_plotting) plot(y = asd[1, 6:10], x = c$freq, type = 'h')
  
  # more staff -> range & freq
  # amplitude (range)
  expect_warning(asd <- generate_network(trj, post_processing_method = "amplitude", window = 10))
  expect_error(asd <- generate_network(trj, post_processing_method = "amplitude", window = 10), NA)
  res <- c()
  for(i in 1:ncol(trj))
    res <- c(res, diff(range(c(trj[1:5, i], trj[1:5, i]))))
  expect_true(all(res == asd$trj_out[1, ]))
  # maxfreq
  expect_warning(asd <- generate_network(trj, post_processing_method = "maxfreq", window = 10))
  expect_error(asd <- generate_network(trj, post_processing_method = "maxfreq", window = 10), NA)
  res <- c()
  if(my_libs){
    for(i in 1:ncol(trj)){
      P <- periodogram(c(trj[1:5, i], trj[1:5, i]), plot = F)
      res <- c(res, P$freq[which(P$spec == max(P$spec))])
    }
    expect_true(all(res == asd[1, ]))
  }
  # both
  expect_warning(asd <- generate_network(trj, post_processing_method = "amplitude_maxfreq", window = 10))
  expect_error(asd <- generate_network(trj, post_processing_method = "amplitude_maxfreq", window = 10), NA)
  
  
  # MINE FAMILY
  library(minerva)
  lung <- 2000
  trj <- as.data.frame(matrix(NA, nrow = lung, ncol = 4))
  trj[,1] <- sin(1:lung) + tan(lung:1)*(1/100) + rnorm(lung, mean = 0, sd = 0.1)
  trj[,3] <- sin(1:lung) + tan(lung:1)*(1/100) + rnorm(lung, mean = 0, sd = 0.5)
  trj[,2] <- cos(1:lung) + tanh(lung:1)*(1/sin(lung:1))
  trj[,4] <- cos(1:lung) + tanh(lung:1)*(1/sin(lung:1)) - rnorm(lung, mean = 1, sd = 0.5)
  trj <- rbind(trj, as.data.frame(matrix(rnorm(4000), nrow = 1000, ncol = 4)))
  res <- minerva::mine(x = trj[1:50,]) # real exit values
  a <- res$MIC[lower.tri(x = res$MIC, diag = FALSE)]
  expect_error(b <- generate_network(trj = trj, window = 10, method = 'MIC', silent = silent), NA) # multithreading MIC is super bad
  expect_equal(a, b$trj_out[1,])
  
  # microtest
  if(F){
    library(tictoc)
    mobj <- microbenchmark::microbenchmark(
      mink = generate_network(trj = trj, window = 100, method = 'minkowski', silent = T),
      mink_threads = generate_network(trj = trj, window = 100, method = 'minkowski', silent = T, do_multithreads = T),
      times = 10
    ) # much faster without multithreads! the windows are too small
    print(mobj)
    autoplot(mobj)

    # my problem: MIC is slow TOOO MUCH
    mobj <- microbenchmark::microbenchmark(
      mic = minerva::mine(x = trj[1:500,], n.cores = 1),
      mic_threads = minerva::mine(x = trj[1:500,], n.cores = 7),
      times = 10
    ) # much worse with threads
    mobj <- microbenchmark::microbenchmark(
      mic = minerva::mine(x = trj[1:500,], n.cores = 1),
      mic_less_alpha = minerva::mine(x = trj[1:500,], alpha = 0.4),
      times = 10
    )
    print(mobj)
    autoplot(mobj)
    
    # Analysing alpha factor in face of different number of rows in inference data
    # ------------------------------------------------------------------------------------------------------------------
    
    # alpha seems a promising parameters to be reduced. Analysing how the mic changes with different alpha
    alphass <- seq(0.1, 0.6, length.out = 30)
    SIZES <- round(seq(10, 200, length.out = 20))
    count <- 0; the_minerva_out <- list()
    
    # main loop
    for(size in SIZES){
      count <- count + 1
      the_minerva_out[[count]] <- lapply(alphass, FUN = function(x) {
        tic()
        out <- minerva::mine(x = trj[1:size,], alpha = x)
        out.t <- toc(quiet = T)
        return(list('m.out' = out, 't.out' = out.t))
        })
    }
    
    mic.alph.times <- unlist(lapply(the_minerva_out, FUN = function(y) sapply(y, function(x) x$t.out$toc - x$t.out$tic)))
    ref_mic <- the_minerva_out[[length(SIZES)]][[length(alphass)]]$m.out$MIC
    mic.vals.diff <- unlist(lapply(the_minerva_out, FUN = function(y) sapply(y, function(x)
      return(sum(abs(x$m.out$MIC-ref_mic))/ (nrow(ref_mic)^2)))))
    
    library(RColorBrewer); RColorBrewer::display.brewer.all()
    pb <- round(df.pl[300, ], 3); pmy <- round(df.pl[290, ], 3)
    df.pl <- cbind('mt' = unlist(mic.alph.times), 'mv' = unlist(mic.vals.diff), 'alpha' = rep(alphass, 20), 'sizes' = rep(SIZES, each = 30))
    ggplot(data = df.pl, mapping = aes(x = mt, y = mv, group = as.factor(sizes))) + 
      theme_classic() + geom_point(mapping = aes(fill = alpha), size = 4, shape=21) + geom_line(aes(col = sizes), size = 0.8) + 
      geom_point(x = pmy[1], y = pmy[2], col = 'red') + 
      geom_point(x = pb[1], y = pb[2], col = 'darkred') + 
      annotate('text', x = max(mic.alph.times)/2, y = max(mic.vals.diff)*2/3, 
               label = paste0('A: ', pmy[3], ' | tm-gain: ', pb[1] - pmy[1],
                              '; vm: ', pmy[2]), col = 'red') + 
      ggtitle(label = 'Difference in result with lower values of alpha vs time used', 
            subtitle = 'INDIPENDENT on number of samples. Around 0.4 for best alpha') + 
      xlab('Time (s)') + ylab('Average diff for each value in sim matrix') + 
      scale_fill_gradientn(colors = brewer.pal(9, 'Greens')) + 
      scale_color_gradientn(colors = brewer.pal(9, 'Blues')) 
    
    ggsave('diff_alpha_for_mic_times.png')
    # ------------------------------------------------------------------------------------------------------------------
    
    # continuing la reserch with different parameters
    mobj <- microbenchmark::microbenchmark(
      mica04_approx = minerva::mine(x = trj[1:500,], alpha = 0.4, est = 'mic_approx'),
      mica04 = minerva::mine(x = trj[1:500,], alpha = 0.4),
      times = 100
    ) # not so relevant the approximated proc.
    print(mobj)
    autoplot(mobj)
    mobj <- microbenchmark::microbenchmark(
      mica04_lessC = minerva::mine(x = trj[1:500,], alpha = 0.4, C = 10),
      mica04 = minerva::mine(x = trj[1:500,], alpha = 0.4),
      times = 100
    ) # not so relevant the approximated proc.
    print(mobj)
    autoplot(mobj)
    
    # Analysing C factor in face of different number of rows in inference data
    # ------------------------------------------------------------------------------------------------------------------
    # C seems a promising parameters to be reduced. Analysing how the mic changes with different alpha
    Css <- round(seq(2, 15, length.out = 14))
    # Css <- c(2,2)
    SIZES <- round(seq(10, 300, length.out = 40))
    # SIZES <- c(100, 100)
    count <- 0; the_minerva_out <- list()
    
    # main loop
    for(size in SIZES){
      count <- count + 1
      the_minerva_out[[count]] <- lapply(Css, FUN = function(x) {
        tic()
        out <- minerva::mine(x = trj[1:size,], alpha = 0.6, C = x)
        out.t <- toc(quiet = T)
        return(list('m.out' = out, 't.out' = out.t))
      })
    }
    
    mic.Cs.times <- unlist(lapply(the_minerva_out, FUN = function(y) sapply(y, function(x) x$t.out$toc - x$t.out$tic)))
    ref_mic <- the_minerva_out[[length(SIZES)]][[length(Css)]]$m.out$MIC
    mic.vals.diff <- unlist(lapply(the_minerva_out, FUN = function(y) sapply(y, function(x)
      return(sum(abs(x$m.out$MIC-ref_mic))/ (nrow(ref_mic)^2)))))
    
    library(RColorBrewer); RColorBrewer::display.brewer.all()
    df.pl <- cbind('mt' = mic.Cs.times, 
                   'mv' = mic.vals.diff, 
                   'Css' = rep(Css, length(SIZES)), 
                   'sizes' = rep(SIZES, each = length(Css)))
    
    gg <- ggplot(data = df.pl, mapping = aes(x = mt, y = mv, group = as.factor(sizes))) + 
          theme_classic() + geom_point(mapping = aes(fill = Css), size = 4, shape=21) + geom_line(aes(col = sizes), size = 0.8) + 
          ggtitle(label = 'Difference in result with lower values of Css vs time used', 
                  subtitle = 'INDIPENDENT on number of samples. Around 3 for best Css. alpha at 0.6') + 
          xlab('Time (s)') + ylab('Average diff for each value in sim matrix') + 
          scale_fill_gradientn(colors = brewer.pal(9, 'Greens')) + 
          scale_color_gradientn(colors = brewer.pal(9, 'Blues')) 
    gg
    
    ggsave('diff_Css_for_mic_times.png')
    # ------------------------------------------------------------------------------------------------------------------
    
    # final difference between two methods
    mobj <- microbenchmark::microbenchmark(
      mica04_lessC = minerva::mine(x = trj[1:100,], alpha = 0.4, C = 2),
      mic_std = minerva::mine(x = trj[1:100,]),
      times = 100
    ) # not so relevant the approximated proc.
    print(mobj)
    autoplot(mobj)
    ggsave('finaloptalpha04C2.png')
  }
  
  # SPEED TESTS =====================================================================================================
  if(F){
    
    # THE DATASET
    silent <- F
    plt_stff <- !silent
    if(!silent) {require(microbenchmark); require(testthat); require(CampaRi); library(RcppRoll); library(ggfortify)} 
    ttrj <- system.file("extdata", "NBU.fyc", package = "CampaRi")
    ttrj <- data.table::fread(ttrj, data.table = F)
    
    # simple it is working tests
    a <- struct1(test_data = ttrj, w = 100, m = 'minkowski'); dim(a)
    e <- struct2(test_data = ttrj, w = 100, m = 'minkowski'); dim(e)
    a[1,] == e[1,]
    
    # microbenchmarking it
    mobj <- microbenchmark::microbenchmark(
      forloop10 = struct1(test_data = ttrj, w = 10, m = 'minkowski'),
      applyloop10 = struct2(test_data = ttrj, w = 10, m = 'minkowski'), 
      forloop100 = struct1(test_data = ttrj, w = 100, m = 'minkowski'),
      applyloop100 = struct2(test_data = ttrj, w = 100, m = 'minkowski'), 
      times = 10
    )
    autoplot(mobj)
    
    # comparing a bit of stuff to understand what it is really doing
    a[1,] # this must be equal to the following:
    b <- WGCNA::adjacency(ttrj[1:50,], type = "distance", distOptions = 'method = "minkowski"', distFnc = 'dist')
    b <- b[lower.tri(x = b, diag = FALSE)]
    all(b == a[1,]) # must be true
    
    # lets test wgcna vs distance
    b <- WGCNA::adjacency(ttrj[1:50,], type = "distance", distOptions = 'method = "minkowski"', distFnc = 'dist')
    dd <- dist(t(ttrj[1:50,]), method = 'minkowski') # d = do.call(distFnc, c(list(x = t(datExpr)), distOptions))
    c <- 1 - as.matrix((dd/max(dd, na.rm = TRUE))^2) # this is exaclty what it does! it seems there is not a lot to be faster here
    c <- c[lower.tri(x = c, diag = FALSE)]
    all(c == b)
    
    
    # FUNCTIONS ----------------------------------------------------------------------------------------------------
    # Now we will test some different speeds, i.e. is it possible to avoid WGCNA?
    struct1 <- function(test_data, w=round(nrow(test_data)/100), m='minkowski'){
      trj_out <- matrix(NA, nrow = nrow(test_data), ncol = ((ncol(test_data)-1)*ncol(test_data)/2)) 
      if(((w-1)/2)%%1 == 0) {
        window_r <- window_l <- (w-1)/2
      }else{
        window_r <- floor((w-1)/2)
        window_l <- ceiling((w-1)/2)
      } 
      for(i in 1:dim(test_data)[1]){
        # taking the piece of trj to transform
        if((i - window_l) <= 0) {                  # left border collision
          tmp_trj <- test_data[1:(i + window_r),]
        }else if((window_r + i) > nrow(test_data)){      # right border collision
          tmp_trj <- test_data[(i - window_l):dim(test_data)[1],]
        }else{                                     # usual center piece
          tmp_trj <- test_data[(i - window_l):(i + window_r),]
        } 
        if(i == 1) print(tmp_trj) 
        if(i == 1) print(window_r) 
        if(i == 1) print(window_l) 
        # built_net <- WGCNA::adjacency(tmp_trj, type = wgcna_type, corFnc = 'cor', power = wgcna_power, corOptions = wgcna_corOp)
        distOptions_manual <- paste0("method = '", m, "'")
        built_net <- WGCNA::adjacency(tmp_trj, type = "distance", distOptions = distOptions_manual, distFnc = 'dist') 
        # type = "distance", adjacency = (1-(dist/max(dist))^2)^power.
            
        built_net_fin <- built_net
        trj_out[i,] <- built_net_fin[lower.tri(x = built_net_fin, diag = FALSE)]
      }
      return(trj_out)
    }
    
    
    
    # apply application of the same struct
    struct2 <- function(test_data, w=int(nrow(test_data)/100), m='minkowski'){
      if(((w-1)/2)%%1 == 0) {
        window_r <- window_l <- (w-1)/2
      }else{
        window_r <- floor((w-1)/2)
        window_l <- ceiling((w-1)/2)
      } 
      trj_out <- sapply(X = 1:dim(test_data)[1], FUN = function(i){
        # taking the piece of trj to transform
        if((i - window_l) <= 0) {                  # left border collision
          tmp_trj <- test_data[1:(i + window_r),]
        }else if((window_r + i) > nrow(test_data)){      # right border collision
          tmp_trj <- test_data[(i - window_l):dim(test_data)[1],]
        }else{                                     # usual center piece
          tmp_trj <- test_data[(i - window_l):(i + window_r),]
        } 
        if(i == 1) print(tmp_trj) 
        if(i == 1) print(window_r) 
        if(i == 1) print(window_l) 
        # built_net <- WGCNA::adjacency(tmp_trj, type = wgcna_type, corFnc = 'cor', power = wgcna_power, corOptions = wgcna_corOp)
        distOptions_manual <- paste0("method = '", m, "'")
        built_net <- WGCNA::adjacency(tmp_trj, type = "distance", distOptions = distOptions_manual, distFnc = 'dist') 
        # type = "distance", adjacency = (1-(dist/max(dist))^2)^power.
        
        built_net_fin <- built_net
        return(built_net_fin[lower.tri(x = built_net_fin, diag = FALSE)])
      })
      return(t(trj_out))
    }
    
    
    # Rollit solution trials -> NO GO for Rollit
    ## Not run: 
    x <- matrix(1:160, ncol=4); dim(x)
    
    const_vars <- list(m = "mean(x)")
    var_fun <- "( (x-m) * (x-m) )/(N-1)"
    rolling_var <- rollit( var_fun, const_vars=const_vars )
    
    rolling_var(x, 2)
    rolling_d <- rollit(fun = 'dist(x)') 
    # I need to use dist to do it but in this case is not elastic. No 'complex' r code!
    rolling_d(x, 2)
    
    x <- c(1, 5, 10, 15)
    cbind( rolling_var(x, 2), roll_var(x, 2) )
    
    ## a benchmark
    if( require("microbenchmark") && require("zoo") ) {
      x <- rnorm(1E4)
      microbenchmark(
        rolling_var(x, 100),
        roll_var(x, 100),
        rollapply(x, 100, var),
        times=10
      )
    }
    ## End(Not run)
  }
  # this test needs further packages. They should not be a forced installation
  # expect_error(tsne_net <- generate_network(trj, method = 'MI', post_processing_method = 'tsne', window = 20), NA)
  })
