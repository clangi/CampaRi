context('generate_network')

test_that('pre-processing with network inference', {
  test_plotting <- T
  trj <- matrix(rnorm(1000), nrow = 100, ncol = 10)
  expect_true(!is.null(trj))
  
  # checking some (expected) errors
  # -----------------------------------
  expect_that(a <- generate_network(trj, window = 20, method = "sadasd"), throws_error())
  expect_that(a <- generate_network(trj, window = 20, post_processing_method = "asdsads"), throws_error())
  # testing for distances and postprocessing
  # -----------------------------------
  expect_that(net_joint_tr <- generate_network(trj, method = 'minkowski', post_processing_method = 'svd', window = 18), not(throws_error()))
  expect_that(net1 <- generate_network(trj, method = 'minkowski', window = 20), not(throws_error()))
  expect_that(path_net <- generate_network(trj, method = 'minkowski', post_processing_method = 'path_minkowski', window = 15), not(throws_error()))
  expect_that(SVD_net <- generate_network(trj, method = 'minkowski', post_processing_method = 'svd', window = 20), not(throws_error()))
  # expect_that(MI_net <- generate_network(trj, method = 'MI', post_processing_method = "SymmetricUncertainty", window = 7), not(throws_error()))
 
  # testing for tsne
  # -----------------------------------
  expect_that(tsne_net <- generate_network(trj, method = 'minkowski', post_processing_method = 'tsne', window = 20), not(throws_error()))
  expect_that(tsne_net <- generate_network(trj, method = 'none', post_processing_method = 'tsne', window = 20), not(throws_error()))
  
  
  # testing for spectrum transformation
  # -----------------------------------
  trj <- matrix(NA, nrow = 120, ncol = 2)
  trj[,1] <- sin(1:120)*(1/rep(1:30, 4) + tan(120:1)*(1/100))
  trj[,2] <- cos(1:120)*(1/rep(1:30, 4) + tanh(120:1)*(1/sin(120:1)))
  if(test_plotting) plot(trj[,1], type = 'l')
  if(test_plotting) plot(trj[,2], type = 'l')
  a <- periodogram(y = trj[ ,1], plot = test_plotting)
  b <- periodogram(y = trj[ ,2], plot = test_plotting)
  if(test_plotting) plot(y = b$spec, x = b$freq,type = "h") # spec is what it was needed

  expect_that(asd <- generate_network(trj, method = 'fft', window = 10), not(throws_error()))
  c <- periodogram(y = c(trj[1:5, 1], trj[1:5, 1]), plot = test_plotting)
  d <- periodogram(y = c(trj[1:5, 2], trj[1:5, 2]), plot = test_plotting)
  # final check of equality
  expect_true(all(c$spec == asd[1, 1:5]))
  expect_true(all(d$spec == asd[1, 6:10]))
  
  if(test_plotting) plot(y = c$spec, x = c$freq, type = 'h')
  if(test_plotting) plot(y = asd[1, 1:5], x = c$freq, type = 'h')
  if(test_plotting) plot(y = d$spec, x = d$freq, type = 'h')
  if(test_plotting) plot(y = asd[1, 6:10], x = c$freq, type = 'h')
  
  # more staff -> range & freq
  expect_that(asd <- generate_network(trj, method = 'fft', window = 10), not(throws_error()))
  
  
  
  
  # this test needs further packages. They should not be a forced installation
  # expect_that(tsne_net <- generate_network(trj, method = 'MI', post_processing_method = 'tsne', window = 20), not(throws_error()))
  })