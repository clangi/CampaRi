#----------------------------------------------
#         WELCOME TO CAMPARI WR TESTING
#
# Davide G.
#----------------------------------------------
# 
#
#
#
#

# Set and making working directory
#----------------------------------------------
w_dir <- getwd()
# system("mkdir output_to_delete")
setwd("output_to_delete/")
w_dir <- getwd()
library(devtools)
install_github("Melkiades/CampaRi")
library(CampaRi)


# butane chain simulated data - variables init
# ------------------------------------------
library(bio3d)
trj <- read.dcd("../inst/extdata/NBU_1250fs.dcd")

# ------------------------------------------
# SIMULATION from tutorial 11 - run_campari
# ------------------------------------------


# reading an already formatted keyfile (it will skip the empty lines and comments)
keywords <- keywords_from_keyfile(key_file_input = "nbu.key") #list format
keywords_from_keyfile(key_file_input = "nbu.key", return_string_of_arguments = TRUE) # copy paste format


# standard run of the simulation in tutarial 11
run_campari(FMCSC_SEQFILE="NBU_1250fs.in",
            # FMCSC_BASENAME="NBU", # lets try the base_name option
            base_name = "nbu", print_status = F, # it will take 55 s in background ~
            FMCSC_SC_IPP=0.0,
            FMCSC_SC_BONDED_T=1.0,
            FMCSC_DYNAMICS=3,
            FMCSC_FRICTION=3.0,
            FMCSC_TIMESTEP=0.005,
            FMCSC_TEMP=400.0,
            FMCSC_NRSTEPS=10000000,
            FMCSC_EQUIL=0,
            FMCSC_XYZOUT=100,
            FMCSC_XYZPDB=3,
            FMCSC_TOROUT=100,
            FMCSC_COVCALC=20000000,
            FMCSC_SAVCALC=20000000,
            FMCSC_POLCALC=20000000,
            FMCSC_RHCALC=20000000,
            FMCSC_INTCALC=20000000,
            FMCSC_POLOUT=20000000,
            FMCSC_ENSOUT=20000000,
            FMCSC_ENOUT=20000000,
            FMCSC_RSTOUT=20000000
)


# ------------------------------------------
# tutorial 11 - step 2 (Simple Clustering)
# ------------------------------------------
# take the variables in a function argument format from key file exemple
keys <- keywords_from_keyfile(key_file_input = "keyfile.key", return_string_of_arguments = T, 
                              keyword_list_first = F, keyword_list = c(FMCSC_SEQFILE="NBU_1250fs.in"))
print(keys) # copy paste it in the run_campari or BETTER provide the keyfile directly

# in the following we use the key_file_input to set the majority of the keys and we override SEQFILE directly from the function
run_campari(data_file = "NBU_1250fs.dcd", key_file_input = "keyfile.key", FMCSC_SEQFILE="NBU_1250fs.in") 



# remember to use the following commad on the logfile to extract the clustering summary (or use show_clustering_summary function)
# for i in `grep -F 'CLUSTER SUMMARY' -A 10 LOGFILE | tail -n 9 | awk '{print $3}'`; do (j=`echo -e "$i + 1" | bc`; head -n $j NBU.fyc | tail -n 1 | awk '{printf("Step %10i: C1C2C3C4: %10.3f C3C2C1H11: %10.3f C2C3C4H41: %10.3f\n",$1/10,$2,$3,$4)}') done


# ------------------------------------------
#           MST - run_campari in ncminer mode
# ------------------------------------------
run_campari(trj = trj,
            FMCSC_CPROGINDMODE=1, #mst
            FMCSC_CCOLLECT=1,
            FMCSC_CMODE=4,
            FMCSC_CDISTANCE=7, #rmsd without alignment 7 - dihedral distances need a complete analysis (pdb_format dcd pdb etc...) 
            FMCSC_CPROGINDSTART=21, #starting snapshot 
            # FMCSC_CPROGINDRMAX=1000, #search att
            # FMCSC_BIRCHHEIGHT=2, #birch height
            FMCSC_CMAXRAD=6, #clustering
            FMCSC_CRADIUS=4,
            FMCSC_CCUTOFF=100,
            FMCSC_CPROGINDWIDTH=1000) #local cut
            #FMCSC_CPROGMSTFOLD 4 # b)
sapphire_plot(sap_file = "PROGIDX_000000000021.dat", title = "ORIGINAL CAMPARI - MST")



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
