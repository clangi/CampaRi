context('multiplicate_trj')

test_that('multiplicate trj', {
  trj <- matrix(rnorm(1000), nrow = 100, ncol = 10)
  expect_true(!is.null(trj))
  expect_error(multiplic1 <- multiplicate_trj(trj, window = 15), NA)
  expect_error(multiplic1 <- multiplicate_trj(trj, window = 15, method = 'sincos'), NA)
  expect_error(multiplic1 <- multiplicate_trj(trj, window = 3, spacing_window = 10), NA)
})






