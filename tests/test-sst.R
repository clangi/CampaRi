# testing SST
wd <- "/home/dgarolini/projects/CampaR/testsst/"
package_dir <- "/home/dgarolini/projects/CampaR/"
setwd(wd)
remove.packages("CampaRi", lib="~/R/x86_64-pc-linux-gnu-library/3.1")
install.packages("../CampaRi/",repos = NULL,type = 'source')
library(CampaRi)
file_dcd <- "../CampaRi/inst/extdata/NBU_1250fs.dcd"
input_trj <- CampaRi::load_trj_dcd(t_file = file_dcd)
# CampaRi::write.pdb.d(x = input_trj, base_name = "NBU_250fs")
# library(bio3d)
# in2 <- bio3d::read.pdb("NBU250fs.pdb")
cmaxrad <- mean(input_trj)*4.0/4.0
birch_hei <- 5
crad <- cmaxrad/birch_hei
campari(nsnaps = nrow(input_trj), wd = wd, data_file = file_dcd, camp_home = "/software/campari/", base_name = "nbu", pdb_format = 4,
        cprogindstart = 2,distance_met = 5,birch_height = birch_hei, cmaxrad = cmaxrad, cradius = crad,
        cprogindwidth = floor(nrow(input_trj)/27),search_attempts = nrow(input_trj)/100,methodst = 2)
zap_ggplot(sap_file = "PROGIDX_000000000002.dat",local_cut = T,timeline = T,ann_trace = 2,title = "ORIGINAL CAMPARI")

adjl <- CampaRi::adjl_from_trj(trj = input_trj, distance_method = 1, clu_radius = crad,
                               birch_clu = T, mode = "fortran", rootmax_rad = cmaxrad, logging = F,
                               tree_height = birch_hei,n_search_attempts = nrow(input_trj)/100)
ret <- CampaRi::gen_progindex(adjl = adjl, snap_start = 2)
CampaRi::gen_annotation(ret_data = ret,snap_start = 2,local_cut_width = floor(nrow(input_trj)/27))
zap_ggplot(sap_file = "REPIX_000000000002.dat",local_cut = T,timeline = T, ann_trace = 2,title = "WRAPPER CAMPARI")
st1 <- read.table("PROGIDX_000000000002.dat")
st12 <- read.table("REPIX_000000000002.dat")

# checking the checkable on the output tables of the two softwares
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