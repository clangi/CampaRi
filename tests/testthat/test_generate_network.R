context('generate_network')

test_that('pre-processing with network inference', {
  test_plotting <- F
  my_libs <- F
  if(my_libs) library(TSA)
  trj <- data.frame(rnorm(1000), nrow = 100, ncol = 10)
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
  
  # SPEED TESTS =====================================================================================================
  
  # THE DATASET
  silent <- F
  plt_stff <- !silent
  if(!silent) {require(microbenchmark); require(testthat); require(CampaRi); library(RcppRoll)} 
  ttrj <- system.file("extdata", "NBU.fyc", package = "CampaRi")
  ttrj <- data.table::fread(ttrj, data.table = F)
  
  # simple it is working tests
  a <- struct1(test_data = ttrj, w = 50, m = 'minkowski')
  
  
  
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
  
  
  
  
  struct2 <- function(test_data, w=int(nrow(test_data)/100), m='minkowski'){
    
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
  
  
  
  
  
  
  # this test needs further packages. They should not be a forced installation
  # expect_error(tsne_net <- generate_network(trj, method = 'MI', post_processing_method = 'tsne', window = 20), NA)
  })
