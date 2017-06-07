#----------------------------------------------
#         WELCOME TO CAMPARI WR TESTING
#
# Davide G.
#----------------------------------------------
# Minimum Spanning Tree test - no netcdf support needed





# Set and making working directory
#----------------------------------------------
w_dir <- getwd()
system("mkdir output_to_delete")
setwd("output_to_delete/")
w_dir <- getwd()
library(devtools)
install_github("Melkiades/CampaRi")
library(CampaRi)


# butane chain simulated data - variables init
# ------------------------------------------
library(bio3d)
trj <- read.dcd("../CampaRi/inst/extdata/NBU_1250fs.dcd")



# ------------------------------------------
#           MST - run_campari
# ------------------------------------------
system('cp ../inst/extdata/NBU_1250fs.dcd .')
system('cp ../inst/extdata/nbu.key .')
system('cp ../inst/extdata/NBU_1250fs.in .')
camp_home <- "/software/campariv3/"
campari(nsnaps = nrow(trj), w_dir, "NBU_1250fs.dcd", base_name="nbu", camp_home,
        cprogindstart = 10, pdb_format = 4, distance_met = 5)
sapphire_plot(sap_file = "PROGIDX_000000000010.dat", title = "ORIGINAL CAMPARI - MST")



# ------------------------------------------
#            MST - R backend
# ------------------------------------------
options(CampaRi.data_management = "R") # netcdf handling
adjl <- mst_from_trj(trj = trj, mode = "fortran", normalize_d = FALSE, logging = F)
ret <- gen_progindex(adjl = adjl, snap_start = 10)
ret2 <- gen_annotation(ret, snap_start = 10)
sapphire_plot(sap_file = 'REPIX_000000000010.dat', title = "CAMPARI WRAPPER - MST", 
              timeline = F, ann_trace = F)



# ------------------------------------------
#            MST - Netcdf backend
# ------------------------------------------
options(CampaRi.data_management = "netcdf") # netcdf handling
mst_from_trj(trj = trj, mode = "fortran", normalize_d = T, logging = F)
ret <- gen_progindex(nsnaps = nrow(trj), snap_start = 10)
ret2 <- gen_annotation(ret, snap_start = 10)
sapphire_plot(sap_file = 'REPIX_000000000010.dat', title = "CAMPARI WRAPPER - MST", timeline = T, ann_trace = F)










# recompute the progress index from PROGIDX/REPIX matrix output
# ------------------------------------------
tree <- adjl_from_progindex(prog_index_file = "PROGIDX_000000000010.dat")

r1<-gen_progindex(tree,snap_start = 1000)
r2<-gen_annotation(r1, snap_start = 1000)
sapphire_plot(sap_file = 'REPIX_000000001000.dat',
              title = "Reconstructed timeline from progrex index matrix",
              timeline = F,ann_trace = F)

# Contract fringe regions of the mst and re-calculate mst
# ------------------------------------------
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
