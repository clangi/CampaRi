#----------------------------------------------
#         WELCOME TO CAMPARI WR TESTING
#
# Davide G.
#----------------------------------------------
# Testing simple wgcna creation pre-processing




# Set and making working directory
#----------------------------------
w_dir <- getwd()
system("mkdir output_to_delete")
setwd("output_to_delete/")
# library(CampaRi)

# butane chain simulated data
# ---------------------------
library(bio3d)
trj <- read.dcd("~/projects/CampaR/CampaRi/inst/extdata/NBU_1250fs.dcd")
options(CampaRi.data_management = "netcdf") # netcdf handling


# variables initialization
mst <- F
birch_height <- 8
cluster_maxrad <- 120.0
cluster_rad <- 60.0

distance_measure <- 1

if(mst) {
  cluster_rad <- cluster_maxrad <-.Machine$integer.max
  mst1_sst2 <- 1
  birch_clustering <- F
}else{
  mst1_sst2 <- 2
  birch_clustering <- T
}
# ------------------------------------------------------
# Wrapper run
trj_wgcna <- CampaRi::generate_network(trj = trj, window = 222, method = "minkowski", minkowski_p = 4)
mst_from_trj(trj = trj_wgcna, distance_method = distance_measure, clu_radius = cluster_rad,
             birch_clu = birch_clustering, mode = "fortran", rootmax_rad = cluster_maxrad, logging = F,
             tree_height = birch_height, n_search_attempts = 100)
ret <- gen_progindex(nsnaps=nrow(trj), snap_start = 2)
# ret <- gen_progindex(adjl = adjl, snap_start = 2)
gen_annotation(ret_data = ret,snap_start = 2,local_cut_width = floor(nrow(trj)/27))
sapphire_plot(sap_file = "REPIX_000000000002.dat",timeline = T, only_timeline=F)
# -----------------------------------------------------
# plotting annotation
trj_fyc <- read.table("../../dihedral_original.fyc")
prog_index_table <- read.table("REPIX_000000000002.dat")[,c(3,4)]
xx <- seq(from=1, by=1, to=nrow(trj_fyc))
ann <- t((trj_fyc[prog_index_table[xx, 1], ] %% 360 < 120) + (trj_fyc[prog_index_table[xx, 1], ] %% 360 < 240) + 1)
title1 <- "NBU SST, 30k snapshots - WRAPPER"
sapphire_plot(sap_file = "REPIX_000000000002.dat", local_cut = T, timeline = T, ann_trace = ann[1,], title = title1)

# ----------------------------------------------------
# DELETING THE WORKING DIR where we just worked
setwd("..")
system("rm -rf output_to_delete")

# Testing the covariance construction of WGCNA
# ------------------------

library(WGCNA)
mat_wgcna <- WGCNA::adjacency(datExpr = trj, type = 'unsigned', corFnc = 'cor', power = 1)

# corOptions = "use = 'p', method = 'spearman'") # "use = 'p'" is the standard 
# NB: the spearman construction is outrageously long

homemade_mat_wgcna <- array(0, dim(mat_wgcna))

for(i in 1:nrow(mat_wgcna)){
  for(j in 1:ncol(mat_wgcna)){
    homemade_mat_wgcna[i,j] <- abs(cor(trj[,i], trj[,j]))
  }
}

# Let's check it is really as it seems! -> they are == !
sum(mat_wgcna == homemade_mat_wgcna)
42*42

# let's try the covariance matrix
cov_mat <- cov(trj)
cov_mat <- cov_mat/max(cov_mat)


