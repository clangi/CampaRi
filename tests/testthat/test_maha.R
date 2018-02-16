context('mahalanobis distance')

test_that('Using the pipeline for finding the mahalanobis distance', {
  silent <- T
  plt_stff <- !silent
  if(!silent) require(testthat)
  
  # expect_error(adjl <- mst_from_trj(trj = matrix(rnorm(1000), nrow = 100, ncol = 10), mute_fortran = silent), NA)
  # expect_error(adjl <- mst_from_trj(trj = matrix(rnorm(1000), ncol = 10, nrow = 100),
  #                      distance_method = 5, clu_radius = 100, clu_hardcut = 100,
  #                      birch_clu = FALSE, mode = "fortran", mute_fortran = silent), NA)
  expect_error(adjl <- mst_from_trj(trj = matrix(rnorm(1000),ncol=10,nrow=100),
                       distance_method = 5, clu_radius = 0.1,
                       birch_clu = TRUE, mode = "fortran", rootmax_rad = 1.3,
                       tree_height = 5, n_search_attempts = 50, mute_fortran = silent), NA)
  
  # test mahalanobis distance
  n_elem <- 1000
  n_feat <- 3
  shift <- 0.9
  which_feat_to_shift <- 2:3
  
  # random initiatialization
  a <- matrix(rnorm(n_elem*n_feat), nrow = n_feat, ncol = n_elem)
  b <- matrix(rnorm(n_elem*n_feat), nrow = n_feat, ncol = n_elem)
  c <- matrix(rnorm(n_elem*n_feat), nrow = n_feat, ncol = n_elem)
  # shifting
  a[which_feat_to_shift,] <- a[which_feat_to_shift,] + shift
  b[which_feat_to_shift,] <- b[which_feat_to_shift,] + shift
  c[which_feat_to_shift,] <- c[which_feat_to_shift,] - shift
  # normalizing
  a <- a/max(a)
  b <- b/max(b)
  c <- c/max(c)
  tra <- cbind(a,b,c)
  
  # find the distance
  expect_error(Xad <- find_mahalanobis(a, b, c, mute_fortran = F), NA)
  expect_error(adjl <- mst_from_trj(trj = t(tra), distance_method = 12, distance_matrix = Xad, mute_fortran = F), NA)
  
  
  if(file.exists('MST_DUMPLING.nc')) file.remove('MST_DUMPLING.nc')
  if(file.exists('Rplots.pdf')) file.remove('Rplots.pdf')
})
