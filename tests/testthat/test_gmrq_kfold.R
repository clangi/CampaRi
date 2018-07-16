context('score_sapphire')

test_that('scoring sapphire plots', {
  
  # standard initialization
  silent <- T
  plt_stff <- !silent
  if(!silent) {require(testthat); require(CampaRi); require(ggplot2)} 
  file.pi <- system.file("extdata", "REPIX_000000000021.dat", package = "CampaRi")
  # nSBR(file.pi, ny = 20, unif.splits = seq(10,100), return_ordered_predicted = T, dbg_nSBR =F)
  
  
  # Create some spheric data around three distinct centers
  x <- rbind(matrix(rnorm(1000, mean = 0, sd = 0.5), ncol = 2),
             matrix(rnorm(1000, mean = 2, sd = 0.5), ncol = 2),
             matrix(rnorm(1000, mean = 4, sd = 0.5), ncol = 2))
  if(plt_stff) ggplot(data.frame(x = x[,1], y = x[,2]), aes(x=x, y=y)) + theme_classic() + geom_point()
  
  # debugonce(gmrq.kfold)
  params <- data.frame('nx'=c(7,20,50, 70, 100))
  preproc <- list(basename = 't',  FMCSC_CPROGINDMODE=1, #mst
                                   FMCSC_CCOLLECT=1, 
                                   FMCSC_CMODE=4,
                                   FMCSC_CDISTANCE=7, #rmsd without alignment 7 - dihedral distances need a complete analysis (pdb_format dcd pdb etc...) 
                                   # FMCSC_CPROGINDSTART=1, #starting snapshot  # this is a standard now!!!
                                   FMCSC_CMAXRAD=10880, #clustering
                                   FMCSC_CRADIUS=10880,
                                   FMCSC_CCUTOFF=10880,
                                   FMCSC_CPROGINDWIDTH=1000)
  # ------------------------------------------------------ testing kmeans
  expect_error(out <- CampaRi::gmrq.kfold(traj = x, gmrq = c(2:3), lag = 1, kfolds = 5, clust.method = "kmeans", clust.param = data.frame('k'=2:4), 
                                          chunks = F, plot = T), NA)
  # ------------------------------------------------------ testing sbr error - no preproc
  expect_error(out <- CampaRi::gmrq.kfold(traj = x, gmrq = c(2:3), lag = 1, kfolds = 5, clust.method = "sbr", clust.param = params, 
                                          chunks = F, plot = T))
  # ------------------------------------------------------ testing kmeans error
  # expect_error(out <- CampaRi::gmrq.kfold(traj = x, gmrq = c(2:3), lag = 1, kfolds = 5, clust.method = "kmeans", clust.param = data.frame('k'=2:4), 
  #                                         chunks = T, plot = T))
  
  
  # ------------------------------------------------------ testing SBR
  # expect_error(out <- CampaRi::gmrq.kfold(traj = x, gmrq = c(2:3), lag = 1, kfolds = 5, clust.method = "sbr", clust.param = params, dbg_gmrq.kfold = F,
  #                                         chunks = F, preproc = preproc, plot = T), NA)
  # expect_error(out <- CampaRi::gmrq.kfold(traj = x, gmrq = 2, lag = 1, kfolds = 5, clust.method = "sbr", clust.param = params, dbg_gmrq.kfold = F, 
  #                                         chunks = F, preproc = preproc, plot = T), NA)
  
  # ------------------------------------------------------ testing nSBR
  params <- data.frame(comb_met = c('MIC'),
                       unif.splits = I(rep(list(seq(8,150, 4)),4)),  
                       pk_span = 100, ny = 50,
                       n.cluster = 2:5, stringsAsFactors = F)
  # do.call(nSBR, c(data = 'PROGIDX_000000000001_fold5.dat', plot = T, lapply(params, function(x) x[[4]])))
  # expect_error(out <- CampaRi::gmrq.kfold(traj = x, gmrq = 2, lag = 1, kfolds = 5, clust.method = "nsbr", clust.param = params, dbg_gmrq.kfold = F, 
  #                                         chunks = T, preproc = preproc, plot = T, dbg_gmrq.kfold = T), NA)
  
  if(file.exists('FRAMES_NBL.nc')) file.remove('FRAMES_NBL.nc')
  if(file.exists('PROGIDX_000000000021.dat')) file.remove('PROGIDX_000000000021.dat')
  file_to_d <- c(list.files(pattern = 'STRUCT_CLUSTERING*'), list.files(pattern = 'ascii_based_analysis.*'), list.files(pattern = 'PROG'))
  for(i in file_to_d) {if(file.exists(i)) file.remove(i)}
  
  
  if(F){
    library(clusterCrit); library(ClusterR)
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
    profvis::profvis(out <- CampaRi::gmrq.kfold(traj = x, gmrq = c(2:3), lag = 1, kfolds = 5, clust.method = "kmeans", clust.param = data.frame('k'=2:4), plot = F))
  }
})
