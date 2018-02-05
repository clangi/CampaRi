# context('mst-alternatives')
# 
# test_that('Comparing times for building a minimum spanning tree', {
#   dataset <- matrix(rnorm(10000), nrow = 1000, ncol = 10)
#   expect_error(adjl <- mst_from_trj(trj = dataset, distance_method = 1, dump_to_netcdf = TRUE), NA)
# 
#   
#   expect_error(ret <- gen_progindex(nsnaps = 100, read_from_netcdf = T, snap_start = 21), NA)
#   expect_true(!is.null(ret))
#   expect_error(ret2 <- gen_annotation(ret, snap_start = 21, local_cut_width = 10), NA)
#   expect_true(!is.null(ret2))
#   
#   # testing alternatives for mst building
#   library(microbenchmark)
#   library(igraph)
#   library(vegan)
#   mbm <- microbenchmark(
#     'CampaRi1' = mst_from_trj(dataset, distance_method = 1, dump_to_netcdf = TRUE, mute_fortran = T),
#     'CampaRi5' = mst_from_trj(dataset, distance_method = 5, dump_to_netcdf = TRUE, mute_fortran = T),
#     'igraph_euc' = igraph::mst(graph = igraph::graph_from_adjacency_matrix(dist(dataset))),
#     'vegan_euc' = vegan::spantree(d = dist(dataset))
#   )
#   library(ggplot2)
#   autoplot(mbm)  
#   adjl <- mst_from_trj(trj = dataset, distance_method = 1, dump_to_netcdf = F)
#   
#   adjl2 <- vegan::spantree(d = dist(dataset))
#   adjl2 <- igraph::mst(graph = igraph::graph_from_edgelist(dist(dataset)))
#   # as.data.frame(table(adjl2$kid)) #!!!! to count stuff
#   adjl3 <- list(as.data.frame(table(adjl2$kid))[,2], adjl2$kid, adjl2$dist)
#   ret <- gen_progindex(adjl = adjl3, snap_start = 1)
#   
#   
#   if(file.exists('MST_DUMPLING.nc')) file.remove('MST_DUMPLING.nc')
#   if(file.exists('REPIX_000000000021.dat')) file.remove('REPIX_000000000021.dat')
# })
# 
