context('score_sapphire')

test_that('scoring sapphire plots', {
  
  # standard initialization
  silent <- F
  plt_stff <- !silent
  if(!silent) {require(testthat); require(CampaRi)} 
  library(clusterCrit); library(ClusterR)
  
  # Create some spheric data around three distinct centers
  x <- rbind(matrix(rnorm(1000, mean = 0, sd = 0.5), ncol = 2),
             matrix(rnorm(1000, mean = 2, sd = 0.5), ncol = 2),
             matrix(rnorm(1000, mean = 4, sd = 0.5), ncol = 2))
  plot(x[,1], x[,2])
  vals <- vector()
  for (k in 2:6) {
    # Perform the kmeans algorithm
    if(F) {
      cl <- kmeans(x, k)
      part <- cl$cluster
    }else{
      cl <- ClusterR::KMeans_rcpp(x, k, initializer = "kmeans++", max_iters = 200, verbose = F)
      part <- cl$clusters
    }
    vals <- c(vals, as.numeric(intCriteria(x, as.vector(part, mode = 'integer'), "PBM")))
  }
  plot(vals) # standard PBM is finding 3 clusters
  debugonce(gmrq.kfold)
  expect_error(out <- CampaRi::gmrq.kfold(traj = x, gmrq = c(2:3), lag = 1, kfolds = 5, clust.method = "kmeans", clust.param = data.frame('k'=2:4), plot = T), NA)
})
