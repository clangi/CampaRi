#THIS FILE SHOULD BE COPIED IN THE WORKING DIRECTORY. Take care about directories specifics
cur_dir <- "~/projects2016/release1/"
# cur_dir <- "~/Projects/2016/CampaR/"
setwd(cur_dir)
# setwd("..")
package_dir <- "CampaRi/"
install.packages(package_dir, repos = NULL, type="source")
library(CampaRi)


# REPIXING: computing sapphire plots from PROGIDX output
tree <- adjl_from_pi(fil = "CampaRi/data/PROGIDX_000000000001.dat")

# contrac <- contract_mst(adjl = tree,n_fold = 0) to check

r1<-gen_progindex(tree,snap_start = 7521)
r2<-gen_annotation(r1, snap_start = 7521)
sapphire_out<-read.table(file = 'REPIX_000000007521.dat')
plot(sapphire_out[,1],-log(sapphire_out[,4]/nrow(sapphire_out)), type="l", 
     ylab = "annotation function", 
     xlab = "snapshots",
     main = paste0("Sapphire of nbu; snaps = ",nrow(sapphire_out)))


# SECOND TEST - butane chain

# Load trajectory [Full 100000 snaps]
trj2<-load_trj_dcd("CampaRi/data/NBU.dcd")

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


#ANALYSIS
install.packages("CampaRi/", repos = NULL, type="source")
library(CampaRi)

trj2<-load_trj_dcd("NBU_250fs.dcd")
trj<-trj2[1:28000,]
adjl<-adjl_from_trj(trj = trj, mode = "fortran")
# adjl<-adjl_from_adjmat(adj_m) #deprecated
ret<-gen_progindex(adjl,snap_start = 10)
ret2<-gen_annotation(ret,snap_start = 10,local_cut_width = 50)
sap_file <- 'REPIX_000000000010.dat'
sap_file <- "CampaRi/data/PROGIDX_000000000001.dat"
zap_ggplot(sap_file = sap_file,timeline = T,ann_trace = 5)
