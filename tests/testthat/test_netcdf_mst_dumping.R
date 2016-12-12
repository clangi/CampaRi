# testing SST and MST with netcdf backend (TREE DUMPING)
# ------------------------------------------------------
wd <- "/home/dgarolini/projects/CampaR/trial_netcdf/"
package_dir <- "/home/dgarolini/projects/CampaR/"
setwd(wd)
remove.packages("CampaRi", lib="~/R/x86_64-pc-linux-gnu-library/3.1")
install.packages("../CampaRi/",repos = NULL,type = 'source')
library(CampaRi)
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
if(dista==1) input_trj <- read.table("dihedral_original.fyc")
adjl <- CampaRi::adjl_from_trj(trj = input_trj, distance_method = dista, clu_radius = crad,
                               birch_clu = b1rchclu, mode = "fortran", rootmax_rad = cmaxrad, logging = F,
                               tree_height = birch_hei, n_search_attempts = 7000)
ret <- gen_progindex(nsnaps=nrow(input_trj), snap_start = 2)
# ret <- gen_progindex(adjl = adjl, snap_start = 2)
CampaRi::gen_annotation(ret_data = ret,snap_start = 2,local_cut_width = floor(nrow(input_trj)/27))
zap_ggplot(sap_file = "REPIX_000000000002.dat",local_cut = T,timeline = T, ann_trace = 2,title = "WRAPPER CAMPARI")

# ----------------------------------------------------
# Adding annotation
progind <- read.table("REPIX_000000000002.dat")[,c(3,4)]
xx <- seq(from=1, by=1, to=nrow(input_trj))
ann <- t((input_trj[progind[xx, 1], ] %% 360 < 120) + (input_trj[progind[xx, 1], ] %% 360 < 240) + 1)
# ann <- (ann[,as.character(xx)])
title1 <- "NBU SST, 30k snapshots - WRAPPER"
zap_ggplot(sap_file = "REPIX_000000000002.dat",local_cut = T,timeline = T,ann_trace = ann,ann_names_L = c("CCCC","HCCC","CCCH"),title = title1)
title2 <- "NBU SST, 30k snapshots - ORIGINAL"
zap_ggplot(sap_file = "PROGIDX_000000000002.dat",local_cut = T,timeline = T,ann_trace = ann,ann_names_L = c("CCCC","HCCC","CCCH"),title = title2)

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
plot(x=1:6000,y=st12[,4]-st1[,4],pch='.',cex=0.1)
points(which(st1[,3]==6000),y = 0,col="red",pch=19,cex=0.8)
points(which(st1[,3]==1),y = 0,col="blue",pch=19,cex=0.8)
#local cut
sum(st1[,10]==0);sum(st12[,10]==0);sum(st1[,12]==0);sum(st12[,12]==0)
mean(st1[,10]);mean(st12[,10]);mean(st1[,12]);mean(st12[,12])

