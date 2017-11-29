context('adjl_from_adjmat')

test_that('Test adjl_from_adjmat', {
  expect_error(a <- adjl_from_adjmat(matrix(rnorm(1000), nrow = 10, ncol = 100)), NA)
})

