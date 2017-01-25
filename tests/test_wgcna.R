# testing new distance based on networks
# ------------------------------------------------------
wd <- "/home/dgarolini/projects/CampaR/trial_netcdf/"
package_dir <- "/home/dgarolini/projects/CampaR/"
setwd(wd)
# source(paste0(package_dir,"CampaRi/inst/extdata/src_netcdf/setting_up_netcdf.R"))
# setting_up_netcdf(paste0(package_dir,"CampaRi/"))
remove.packages("CampaRi", lib="~/R/x86_64-pc-linux-gnu-library/3.1")
install.packages("../CampaRi/",repos = NULL,type = 'source',
                 configure.args = c(RNetCDF = "--with-netcdf-include=/usr/include/"))
library(CampaRi)
# ------------------------------------------------------
# some graph playing
library(igraph)
g1<-erdos.renyi.game(100,0.5)
igraph::plot.igraph(g1,vertex.size = 1)



# ------------------------------------------------------
# Initialization of important variables
options(CampaRi.data_management = "netcdf")
file_dcd <- "../CampaRi/inst/extdata/NBU_250fs.dcd"
file_dcd <- "nbu_napapigiri_traj.dcd"
# file_dcd <- "../base_sim/NBU_traj.dcd"
input_trj <- CampaRi::load_trj_dcd(t_file = file_dcd)
# ------------------------------------------------------
# Torsional space projection
# library(bio3d)
# for(i in 1:24)
#   print(torsion.xyz(input_trj[1,i:(11+i)]))
# for(i in 4:14)
#   print(torsion.xyz(input_trj[1,1:(3*i)],atm.inc = i))
# 
# input_trj2 <- cbind(apply(input_trj[,1:12],MARGIN = 1,torsion.xyz),
#                     apply(input_trj[,16:27],MARGIN = 1,torsion.xyz),
#                     apply(input_trj[,24:35],MARGIN = 1,torsion.xyz))

mst <- F
cmaxrad <- 120.0
birch_hei <- 8
crad <- 60.0
b1rchclu <- T
dista <- 1
metodst <- 2
if(mst) {
  crad <- cmaxrad <-.Machine$integer.max
  metodst <- 1
  b1rchclu <- F
}
# ------------------------------------------------------
# Original campari run
campari(nsnaps = nrow(input_trj), wd = wd, data_file = file_dcd, camp_home = "/software/campari/", base_name = "nbu", pdb_format = 4,
        cprogindstart = 2,distance_met = dista, birch_height = birch_hei, cmaxrad = cmaxrad, cradius = crad,
        cprogindwidth = floor(nrow(input_trj)/27),search_attempts = 7000, methodst = metodst)
zap_ggplot(sap_file = "PROGIDX_000000000002.dat",local_cut = T,timeline = T,ann_trace = 2,title = "ORIGINAL CAMPARI")
# zap_ggplot(sap_file = "../base_sim/PROGIDX_000000000002.dat",local_cut = T,timeline = T,ann_trace = 2,title = "ORIGINAL CAMPARI")
# ------------------------------------------------------
# Wrapper run
options(CampaRi.data_management = "netcdf")
# if(dista==1) input_trj <- read.table("../base_sim/FYC.dat", header = F)[,2:4]
if(dista==1) input_trj2 <- read.table("dihedral_original.fyc")
CampaRi::adjl_from_trj(trj = input_trj2, distance_method = dista, clu_radius = crad,
                       birch_clu = b1rchclu, mode = "fortran", rootmax_rad = cmaxrad, logging = F,
                       tree_height = birch_hei, n_search_attempts = 7000)
ret <- gen_progindex(nsnaps=nrow(input_trj2), snap_start = 2)
# ret <- gen_progindex(adjl = adjl, snap_start = 2)
CampaRi::gen_annotation(ret_data = ret,snap_start = 2,local_cut_width = floor(nrow(input_trj2)/27))
zap_ggplot(sap_file = "REPIX_000000000002.dat",local_cut = T,timeline = T, ann_trace = 2,title = "WRAPPER CAMPARI")
# ----------------------------------------------------
# Adding annotation
progind2 <- read.table("REPIX_000000000002.dat")[,c(3,4)]
xx <- seq(from=1, by=1, to=nrow(input_trj2))
ann <- t((input_trj2[progind2[xx, 1], ] %% 360 < 120) + (input_trj2[progind2[xx, 1], ] %% 360 < 240) + 1)
title1 <- "NBU SST, 30k snapshots - WRAPPER"
zap_ggplot(sap_file = "REPIX_000000000002.dat",local_cut = T,timeline = T,ann_trace = ann,ann_names_L = c("CCCC","HCCC","CCCH"),title = title1)
# ORIGINAL ann
title2 <- "NBU SST, 30k snapshots - ORIGINAL"
progind2 <- read.table("PROGIDX_000000000002.dat")[,c(3,4)]
xx <- seq(from=1, by=1, to=nrow(input_trj2))
ann <- t((input_trj2[progind2[xx, 1], ] %% 360 < 120) + (input_trj2[progind2[xx, 1], ] %% 360 < 240) + 1)
zap_ggplot(sap_file = "PROGIDX_000000000002.dat",local_cut = T,timeline = T,ann_trace = ann,ann_names_L = c("CCCC","HCCC","CCCH"),title = title2)
# ------------------------------------------------------
# Wrapper run with distance = 5
cmaxrad <- mean(input_trj)*4.0/4.0
birch_hei <- 5
crad <- cmaxrad/birch_hei
dista <- 5
options(CampaRi.data_management = "netcdf")
adjl <- CampaRi::adjl_from_trj(trj = input_trj, distance_method = dista, clu_radius = crad,
                               birch_clu = b1rchclu, mode = "fortran", rootmax_rad = cmaxrad, logging = F,
                               tree_height = birch_hei, n_search_attempts = 7000)
ret <- gen_progindex(nsnaps=nrow(input_trj), snap_start = 2)
# ret <- gen_progindex(adjl = adjl, snap_start = 2)
CampaRi::gen_annotation(ret_data = ret,snap_start = 2,local_cut_width = floor(nrow(input_trj)/27))
zap_ggplot(sap_file = "REPIX_000000000002.dat",local_cut = T,timeline = T, ann_trace = 2,title = "WRAPPER CAMPARI")