context('multiplicate_trj')

test_that('multiplicate trj', {
  trj <- matrix(rnorm(1000), nrow = 100, ncol = 10)
  expect_true(!is.null(trj))
  expect_that(multiplic1 <- multiplicate_trj(trj, window = 15), not(throws_error()))
  expect_that(multiplic1 <- multiplicate_trj(trj, window = 15, method = 'sincos'), not(throws_error()))
})






