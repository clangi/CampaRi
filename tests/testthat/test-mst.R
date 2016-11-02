#THIS FILE SHOULD BE COPIED IN THE WORKING DIRECTORY. Take care about directories specifics
wd <- "/home/dgarolini/projects/CampaR/"
# library(devtools)
# wd <- "~/Projects/2016/CampaR/"
setwd(wd)
# install_github("Melkiades/CampaRi")
package_dir <- "CampaRi/"
system(paste0("bash ",package_dir,"src/cleaner.sh"))
remove.packages("CampaRi", lib="~/R/x86_64-pc-linux-gnu-library/3.1")
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
           timeline = F,ann_trace = F)


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
adjl<-adjl_from_trj(trj = trj, mode = "fortran", normalize_d = FALSE, logging = T)
ret<-gen_progindex(adjl, snap_start = 10)
ret2<-gen_annotation(ret,snap_start = 10, local_cut_width = floor(10000/27))
sap_file <- 'REPIX_000000000010.dat'
sap_table <- read.table(sap_file)
zap_ggplot(sap_file = sap_file, 
           title = paste0("Sapphire of nbu. Snapshots = ", nrow(sap_table)),
           timeline = F, ann_trace = F)

# DIRECT COMPARISON OF THE SAME DATA-SET [original campari needed]
# less_original run
data_file <- "inst/extdata/NBU_1250fs.dcd"
trj<-load_trj_dcd(paste0(package_dir,data_file))
adjl<-adjl_from_trj(trj = trj, mode = "fortran", normalize_d = FALSE,birch_clu = T,logging = T,tree_height = 3)
ret<-gen_progindex(adjl, snap_start = 10)
ret2<-gen_annotation(ret, snap_start = 10, local_cut_width = floor(6000/27))
sap_file <- 'REPIX_000000000010.dat'
sap_table <- read.table(sap_file)
# sap_table[,4] <- sap_table[,4] + mean(sap_table2[,4])-mean(sap_table[,4])
zap_ggplot(sap_file = sap_file)
# zap_ggplot(sap_table = sap_table)

# original run on the same input file
system(paste0('cp ',package_dir,'inst/extdata/NBU* .'))
system(paste0('cp ',package_dir,'inst/extdata/nbu.* .'))
camp_home <- "/software/campari/"
campari(nsnaps = nrow(trj), wd, "NBU_1250fs.dcd", base_name="nbu",camp_home,
        cprogindstart = 10,pdb_format = 4,distance_met = 5,cprogindwidth=floor(6000/27))
sap_file2 <- "PROGIDX_000000000010.dat"
sap_table2 <- read.table(sap_file2)
zap_ggplot(sap_file2)
