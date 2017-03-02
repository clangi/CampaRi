#little script to write pdb
wd <- "/home/dgarolini/projects/CampaR/pdbtrial/"
setwd(wd)
remove.packages("CampaRi", lib="~/R/x86_64-pc-linux-gnu-library/3.1")
install.packages("../CampaRi/",repos = NULL,type = 'source')
library(CampaRi)
input_trj <- matrix(stats::runif(42000),nrow = 1000,ncol = 42)
# input_trj <- matrix(rnorm(30000),nrow = 100,ncol = 300)
input_trj[500:1000,] <- input_trj[500:1000,] + 5 
# test3_bio3d <- rnorm(300)
# write.pdb.d(test3_bio3d, filename = "pdb_test3", round = T)
write.pdb.d(input_trj, base_name = "flottiglie", round = T, digit = 3)
input_trj <- round(input_trj,digits = 3)
library(bio3d)
bio3d::read.pdb("flottiglie.pdb")
cmaxrad <- mean(input_trj)*4.0/4.0
birch_hei <- 5
crad <- cmaxrad/birch_hei
campari(nsnaps = nrow(input_trj), wd = wd, data_file = "flottiglie.pdb", camp_home = "/software/campari/", base_name = "flottiglie", pdb_format = 1,
        cprogindstart = 1,distance_met = 5,birch_height = birch_hei, cmaxrad = cmaxrad, cradius = crad,
        cprogindwidth = 100,search_attempts = nrow(input_trj)/10,methodst = 2)
zap_ggplot(sap_file = "PROGIDX_000000000001.dat",local_cut = T)

adjl <- CampaRi::adjl_from_trj(trj = input_trj, distance_method = 5, clu_radius = crad,
                               birch_clu = T, mode = "fortran", rootmax_rad = cmaxrad, logging = F,
                               tree_height = birch_hei,n_search_attempts = nrow(input_trj)/10)
ret <- CampaRi::gen_progindex(adjl = adjl, snap_start = 1)
CampaRi::gen_annotation(ret_data = ret,snap_start = 1,local_cut_width = 100)
zap_ggplot(sap_file = "REPIX_000000000000.dat",local_cut = T)


# THESE TESTINGS COME FROM ANOTHER TEST (TEST_NETCDF_MST_DUMPING)
# ------------------------------------------------------
# Wrapper run with distance = 5
cluster_maxrad <- mean(trj)*4.0/4.0
birch_height <- 5
cluster_rad <- cluster_maxrad/birch_height
distance_measure <- 5
options(CampaRi.data_management = "netcdf")
adjl <- CampaRi::mst_from_trj(trj = trj, distance_method = distance_measure, clu_radius = cluster_rad,
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
