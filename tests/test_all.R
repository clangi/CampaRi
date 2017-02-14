#----------------------------------------------
#         WELCOME TO CAMPARI TESTING
#----------------------------------------------
#----------------------------------------------
# All these tests comes with netcdf integration.
# Therefore, install it.


#----------------------------------------------
# Minimum Spanning Tree test 
#----------------------------------------------

# Set working directory
#----------------------------------------------
w_dir <- getwd()
system("mkdir output_to_delete")
setwd("output_to_delete/")
w_dir <- getwd()
# remove.packages("CampaRi", lib="~/R/x86_64-pc-linux-gnu-library/3.1")
library(devtools)
install_github("Melkiades/CampaRi")
library(CampaRi)


# butane chain simulated data - load trajectory 
# ------------------------------------------
library(bio3d)
trj <- read.dcd("../inst/extdata/NBU_1250fs.dcd")

# Wrapper analysis
adjl <- mst_from_trj(trj = trj, mode = "fortran", normalize_d = FALSE, logging = F)
ret <- gen_progindex(adjl, snap_start = 10)
ret2 <- gen_annotation(ret, snap_start = 10)
sapphire_plot(sap_file = 'REPIX_000000000010.dat', title = "CAMPARI WRAPPER - MST", 
              timeline = F, ann_trace = F)

# Original campari analysis
system('cp ../inst/extdata/NBU_1250fs.dcd .')
system('cp ../inst/extdata/nbu.key .')
system('cp ../inst/extdata/nbu.in .')
camp_home <- "/software/campari/"
campari(nsnaps = nrow(trj), w_dir, "NBU_1250fs.dcd", base_name="nbu", camp_home,
        cprogindstart = 10, pdb_format = 4, distance_met = 5)
sapphire_plot(sap_file = "PROGIDX_000000000010.dat", title = "ORIGINAL CAMPARI - MST")

# REPIXING: computing sapphire plots from PROGIDX output
tree <- adjl_from_progindex(prog_index_file = "PROGIDX_000000000010.dat")


r1<-gen_progindex(tree,snap_start = 1000)
r2<-gen_annotation(r1, snap_start = 1000)
sapphire_plot(sap_file = 'REPIX_000000001000.dat',
              title = "Reconstructed timeline from progrex index matrix",
              timeline = F,ann_trace = F)

contrac <- contract_mst(adjl = tree, n_fold = 50)
r1<-gen_progindex(contrac,snap_start = 1000)
r2<-gen_annotation(r1, snap_start = 1000)
sapphire_plot(sap_file = 'REPIX_000000001000.dat',
              title = "Reconstructed timeline from progrex index matrix",
              timeline = F,ann_trace = F)

#             -----------------
# DELETING THE WORKING DIR where we just worked
setwd("..")
system("rm -rf output_to_delete")
