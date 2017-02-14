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
library(devtools)
install_github("Melkiades/CampaRi")
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

distance_measure <- 5

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
system('cp ../inst/extdata/nbu.in .')
camp_home <- "/software/campari/"
campari(nsnaps = nrow(trj), working_dir = w_dir, data_file = "NBU_1250fs.dcd", camp_home = "/software/campari/", base_name = "nbu", pdb_format = 4,
        cprogindstart = 2, distance_met = distance_measure, birch_height = birch_height, cmaxrad = cluster_maxrad, cradius = cluster_rad,
        cprogindwidth = floor(nrow(trj)/27), search_attempts = 100, methodst = mst1_sst2)
sapphire_plot(sap_file = "PROGIDX_000000000002.dat",local_cut = T,timeline = T,ann_trace = 2,title = "ORIGINAL CAMPARI")
# zap_ggplot(sap_file = "../base_sim/PROGIDX_000000000002.dat",local_cut = T,timeline = T,ann_trace = 2,title = "ORIGINAL CAMPARI")
# ------------------------------------------------------
# Wrapper run
if(distance_measure==1) trj2 <- read.table("dihedral_original.fyc")
mst_from_trj(trj = trj, distance_method = distance_measure, clu_radius = cluster_rad,
                       birch_clu = birch_clustering, mode = "fortran", rootmax_rad = cluster_maxrad, logging = F,
                       tree_height = birch_height, n_search_attempts = 100)
ret <- gen_progindex(nsnaps=nrow(trj), snap_start = 2)
# ret <- gen_progindex(adjl = adjl, snap_start = 2)
gen_annotation(ret_data = ret,snap_start = 2,local_cut_width = floor(nrow(trj)/27))
sapphire_plot(sap_file = "REPIX_000000000002.dat",local_cut = T,timeline = T, ann_trace = 2,title = "WRAPPER CAMPARI")
# ----------------------------------------------------
# Adding annotation
progind2 <- read.table("REPIX_000000000002.dat")[,c(3,4)]
xx <- seq(from=1, by=1, to=nrow(trj2))
ann <- t((trj2[progind2[xx, 1], ] %% 360 < 120) + (trj2[progind2[xx, 1], ] %% 360 < 240) + 1)
title1 <- "NBU SST, 30k snapshots - WRAPPER"
zap_ggplot(sap_file = "REPIX_000000000002.dat",local_cut = T,timeline = T,ann_trace = ann,ann_names_L = c("CCCC","HCCC","CCCH"),title = title1)
# ORIGINAL ann
title2 <- "NBU SST, 30k snapshots - ORIGINAL"
progind2 <- read.table("PROGIDX_000000000002.dat")[,c(3,4)]
xx <- seq(from=1, by=1, to=nrow(trj2))
ann <- t((trj2[progind2[xx, 1], ] %% 360 < 120) + (trj2[progind2[xx, 1], ] %% 360 < 240) + 1)
zap_ggplot(sap_file = "PROGIDX_000000000002.dat",local_cut = T,timeline = T,ann_trace = ann,ann_names_L = c("CCCC","HCCC","CCCH"),title = title2)
# ------------------------------------------------------
# Wrapper run with distance = 5
cluster_maxrad <- mean(trj)*4.0/4.0
birch_height <- 5
cluster_rad <- cluster_maxrad/birch_height
distance_measure <- 5
options(CampaRi.data_management = "netcdf")
adjl <- CampaRi::adjl_from_trj(trj = trj, distance_method = distance_measure, clu_radius = cluster_rad,
                               birch_clu = birch_clustering, mode = "fortran", rootmax_rad = cluster_maxrad, logging = F,
                               tree_height = birch_height, n_search_attempts = 7000)
ret <- gen_progindex(nsnaps=nrow(trj), snap_start = 2)
# ret <- gen_progindex(adjl = adjl, snap_start = 2)
CampaRi::gen_annotation(ret_data = ret,snap_start = 2,local_cut_width = floor(nrow(trj)/27))
zap_ggplot(sap_file = "REPIX_000000000002.dat",local_cut = T,timeline = T, ann_trace = 2,title = "WRAPPER CAMPARI")
# ------------------------------------------------------
# checking the checkable on the output tables of the two softwares
st1 <- read.table("PROGIDX_000000000002.dat")
st12 <- read.table("REPIX_000000000002.dat")
sum(st1[,3]==0);sum(st12[,3]==0);sum(st1[,4]==0);sum(st12[,4]==0) # are there any zeros?
sum(st1[,6]==0);sum(st12[,6]==0);sum(st1[,5]==0);sum(st12[,5]==0)


sum(st12[,3]!=st1[,3]);sum(st12[,6]!=st1[,6]) # neighbourlist is different? 3 connects 6
sum(st12[,4]!=st1[,4]) # is the annotation ==?
sum(st12[,4]-st1[,4])/6000 # average difference between annotation values?
mean(st12[,4])-mean(st1[,4]);mean(st1[,4]);mean(st12[,4]) # annetation differences
mean(st1[,5])==mean(st12[,5]) # is the average distance ==?
sum(st12[,5]-st1[,5]) # are the distances ==?
sum(st12[,5]!=st1[,5])

# REASSUMING THE DIFFERENCIES - MST
# the connected components are identical. The problem is in the annotation value.
# 108 are == 0, in average they differ for 2.8 and 5991 snaps are affected.
# instead the distance value is different in 6 cases for which the difference is (totally) e-16

max(abs(st12[,4]-st1[,4]))
plot(x=1:30000,y=st12[,4]-st1[,4],pch='.',cex=0.1)
points(which(st1[,3]==30000),y = 0,col="red",pch=19,cex=0.8)
points(which(st1[,3]==1),y = 0,col="blue",pch=19,cex=0.8)
#local cut
sum(st1[,10]==0);sum(st12[,10]==0);sum(st1[,12]==0);sum(st12[,12]==0)
mean(st1[,10]);mean(st12[,10]);mean(st1[,12]);mean(st12[,12])

