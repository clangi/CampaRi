#THIS FILE SHOULD BE COPIED IN THE WORKING DIRECTORY. Take care about directories specifics
# wd <- "/home/dgarolini/projects2016/release1/"
wd <- "~/Projects/2016/CampaR/"
setwd(wd)
package_dir <- "CampaRi/"
install.packages(package_dir, repos = NULL, type="source")
library(CampaRi)


# REPIXING: computing sapphire plots from PROGIDX output
tree <- adjl_from_pi(fil = "CampaRi/inst/extdata/output_examples/PROGIDX_000000000001.dat")

# contrac <- contract_mst(adjl = tree,n_fold = 0) to check

r1<-gen_progindex(tree,snap_start = 7521)
r2<-gen_annotation(r1, snap_start = 7521)
sap_file <- 'REPIX_000000007521.dat'
zap_ggplot(sap_file = sap_file,
           title = paste0("Sapphire of nbu. Snapshots = ", 100000),
           timeline = F, ann_trace = F)


# SECOND TEST - butane chain2
# Load trajectory [Full 1000000 snaps]
trj2<-load_trj_dcd("CampaRi/inst/extdata/NBU.dcd")

# SUBSAMPLING. dim_reduction is the variable indicating the factor of it
# the computational complexity is O(rmsd)~O(d) having d = dimension of a snap
# *O(number of distances calculated) = O(bin(n,2)) = O(n^2) ---> O(dn^2)
object.size(trj2)/1000/1000 #Mb
dim(trj2)
dim_reduction<-10
trj<-matrix(trj2[seq(1,nrow(trj2),dim_reduction),],nrow = dim(trj2)[1]/dim_reduction, ncol = dim(trj2)[2])
dim(trj)[1]^2
dim(trj)[1]^2*dim(trj)[2]
object.size(trj)/1000 #Kb
dim(trj)


#ANALYSIS of the subsampled data-set (10 fold reduction)
install.packages("CampaRi/", repos = NULL, type="source")
library(CampaRi)

adjl<-adjl_from_trj(trj = trj, mode = "fortran")
ret<-gen_progindex(adjl,snap_start = 10)
ret2<-gen_annotation(ret,snap_start = 10,local_cut_width = 50)
sap_file <- 'REPIX_000000000010.dat'
sap_table <- read.table(sap_file)
zap_ggplot(sap_file = sap_file, 
           title = paste0("Sapphire of nbu. Snapshots = ", nrow(sap_table)),
           timeline = F, ann_trace = F)

# DIRECT COMPARISON OF THE SAME DATA-SET [original campari needed]
# less_original run
data_file <- "CampaRi/inst/extdata/NBU_1250fs.dcd"
trj<-load_trj_dcd(data_file)
adjl<-adjl_from_trj(trj = trj, mode = "fortran")
ret<-gen_progindex(adjl,snap_start = 10)
ret2<-gen_annotation(ret,snap_start = 10,local_cut_width = 50)
sap_file <- 'REPIX_000000000010.dat'
sap_table <- read.table(sap_file)
zap_ggplot(sap_file = sap_file)

# original run on the same input file
camp_home <- "/software/campari/"
campari(wd, data_file, base_name="nbu",camp_home)
sap_file2 <- "PROGIDX_000000000010.dat"
sap_table2 <- read.table(sap_file2)
zap_ggplot(sap_file2)

# checking the checkable on the output tables of the two softwares
sum(sap_table[,3]==0);sum(sap_table2[,3]==0);sum(sap_table[,4]==0);sum(sap_table2[,4]==0)
sum(sap_table[,6]==0);sum(sap_table2[,6]==0);mean(sap_table[,5]);mean(sap_table2[,5])
mean(sap_table[,5])==mean(sap_table2[,5]);mean(sap_table[,4]);mean(sap_table2[,4])
sum(sap_table2[,3]!=sap_table[,3]);sum(sap_table2[,6]!=sap_table[,6])
sum(sap_table2[,4]!=sap_table[,4]);sum(sap_table2[,5]!=sap_table[,5])
sum(sap_table[,10]==0);sum(sap_table2[,10]==0);sum(sap_table[,12]==0);sum(sap_table2[,12]==0)
mean(sap_table[,10]);mean(sap_table2[,10]);mean(sap_table[,12]);mean(sap_table2[,12])

# checking the mst output to spot the above differencies *the flag mst_print must be on INSIDE clustering.f90 (campari source)
# adjl_dis_or <- read.table("mst_original.txt",sep = "", fill=T,stringsAsFactors= F)
adjl_dis_or <- read.fwf("mst_original.txt", widths = c(15,15,15,15,15,15), fill=T) #supposing 6 column(from mst - max_degree)
adjl_dis_or[is.na(adjl_dis_or)] <- 0