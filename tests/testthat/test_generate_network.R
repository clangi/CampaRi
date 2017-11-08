context('generate_network')

test_that('pre-processing with network inference', {
  test_plotting <- F
  my_libs <- F
  if(my_libs) library(TSA)
  trj <- matrix(rnorm(1000), nrow = 100, ncol = 10)
  expect_true(!is.null(trj))
  
  # checking some (expected) errors
  # -----------------------------------
  expect_error(a <- generate_network(trj, window = 20, method = "sadasd"))
  expect_error(a <- generate_network(trj, window = 20, post_processing_method = "asdsads"))
  # testing for distances and postprocessing
  # -----------------------------------
  expect_error(net_joint_tr <- generate_network(trj, method = 'minkowski', post_processing_method = 'svd', window = 18), NA)
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
  trj <- matrix(NA, nrow = 120, ncol = 2)
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
  res <- c()
  for(i in 1:ncol(trj))
    res <- c(res, diff(range(c(trj[1:5, i], trj[1:5, i]))))
  expect_true(all(res == asd[1, ]))
  # maxfreq
  expect_warning(asd <- generate_network(trj, post_processing_method = "maxfreq", window = 10))
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
  
  # this test needs further packages. They should not be a forced installation
  # expect_error(tsne_net <- generate_network(trj, method = 'MI', post_processing_method = 'tsne', window = 20), NA)
  })
