#----------------------------------------------
#         WELCOME TO CAMPARI WR TESTING
#
# Davide G.
#----------------------------------------------
# Testing SST and MST with netcdf backend (TREE DUMPING)




# Set and making working directory
#----------------------------------------------
w_dir <- getwd()
system("mkdir output_to_delete")
setwd("output_to_delete/")
w_dir <- getwd()
# remove.packages("CampaRi", lib="~/R/x86_64-pc-linux-gnu-library/3.1")
# library(devtools)
# install_github("Melkiades/CampaRi")
library(CampaRi)

# butane chain simulated data - variables init
# ------------------------------------------
library(bio3d)
trj <- read.dcd("../inst/extdata/NBU_1250fs.dcd")
options(CampaRi.data_management = "netcdf") # netcdf handling

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
# Original campari run
system('cp ../inst/extdata/NBU_1250fs.dcd .')
system('cp ../inst/extdata/nbu.key .')
system('cp ../inst/extdata/NBU_1250fs.in .')
camp_home <- "/software/campariv3/"
campari(nsnaps = nrow(trj), working_dir = w_dir, data_file = "NBU_1250fs.dcd", camp_home = camp_home, base_name = "nbu", pdb_format = 4,
        cprogindstart = 2, distance_met = distance_measure, birch_height = birch_height, cmaxrad = cluster_maxrad, cradius = cluster_rad,
        cprogindwidth = floor(nrow(trj)/27), search_attempts = 100, methodst = mst1_sst2)
sapphire_plot(sap_file = "PROGIDX_000000000002.dat",local_cut = T,timeline = T,title = "ORIGINAL CAMPARI")
# zap_ggplot(sap_file = "../base_sim/PROGIDX_000000000002.dat",local_cut = T,timeline = T,ann_trace = 2,title = "ORIGINAL CAMPARI")
# ------------------------------------------------------
# Wrapper run
if(distance_measure==1) trj <- read.table("../../dihedral_original.fyc")
mst_from_trj(trj = trj, distance_method = distance_measure, clu_radius = cluster_rad,
             birch_clu = birch_clustering, mode = "fortran", rootmax_rad = cluster_maxrad, logging = F,
             tree_height = birch_height, n_search_attempts = 100)
ret <- gen_progindex(nsnaps=nrow(trj), snap_start = 2)
# ret <- gen_progindex(adjl = adjl, snap_start = 2)
gen_annotation(ret_data = ret,snap_start = 2,local_cut_width = floor(nrow(trj)/27))
sapphire_plot(sap_file = "REPIX_000000000002.dat",local_cut = T,timeline = T,title = "WRAPPER CAMPARI")
# ----------------------------------------------------
# Adding annotation
title1 <- "NBU SST, 30k snapshots - WRAPPER"
progind2 <- read.table("REPIX_000000000002.dat")[,c(3,4)]
xx <- seq(from=1, by=1, to=nrow(trj))
ann <- t((trj[progind2[xx, 1], ] %% 360 < 120) + (trj[progind2[xx, 1], ] %% 360 < 240) + 1)
sapphire_plot(sap_file = "REPIX_000000000002.dat",local_cut = T,timeline = T,ann_trace = ann,ann_names_L = c("CCCC","HCCC","CCCH"),title = title1)
# ORIGINAL ann
title2 <- "NBU SST, 30k snapshots - ORIGINAL"
progind2 <- read.table("PROGIDX_000000000002.dat")[,c(3,4)]
xx <- seq(from=1, by=1, to=nrow(trj))
ann <- t((trj[progind2[xx, 1], ] %% 360 < 120) + (trj[progind2[xx, 1], ] %% 360 < 240) + 1)
sapphire_plot(sap_file = "PROGIDX_000000000002.dat",local_cut = T,timeline = T,ann_trace = ann,ann_names_L = c("CCCC","HCCC","CCCH"),title = title2)

# DELETING THE WORKING DIR where we just worked
setwd("..")
system("rm -rf output_to_delete")
