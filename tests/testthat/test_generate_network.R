context('generate_network')

test_that('pre-processing with network inference', {
  trj <- matrix(rnorm(1000), nrow = 100, ncol = 10)
  expect_true(!is.null(trj))
  expect_that(net_joint_tr <- generate_network(trj, method = 'minkowski', post_processing_method = 'svd', window = 18), not(throws_error()))
  expect_that(net1 <- generate_network(trj, method = 'minkowski', window = 20), not(throws_error()))
  expect_that(path_net <- generate_network(trj, method = 'minkowski', post_processing_method = 'path_minkowski', window = 15), not(throws_error()))
  expect_that(SVD_net <- generate_network(trj, method = 'minkowski', post_processing_method = 'svd', window = 20), not(throws_error()))
  # expect_that(MI_net <- generate_network(trj, method = 'MI', post_processing_method = "SymmetricUncertainty", window = 7), not(throws_error()))
 })