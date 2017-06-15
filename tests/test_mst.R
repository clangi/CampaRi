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


# ------------------------------------------
# SIMULATION from tutorial 11 - run_campari
# ------------------------------------------


# reading an already formatted keyfile (it will skip the empty lines and comments)
keywords <- keywords_from_keyfile(key_file_input = "nbu.key") #list format
keywords_from_keyfile(key_file_input = "nbu.key", return_string_of_arguments = TRUE) # copy paste format


# standard run of the simulation in tutarial 11
run_campari(FMCSC_SEQFILE="nbu.in", # you must have it defined according to CAMPARI's rules
            # FMCSC_BASENAME="NBU", # lets try the base_name option
            base_name = "NBU", print_status = F, # it will take 55 s in background ~
            PARAMETERS="oplsaal.prm", # if this variable it is not supplied will be automatically assigned to <full path to folder>/campari/params/abs3.2_opls.prm
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
            FMCSC_RSTOUT=20000000)


# parenthesis - trials with direct reading and plotting of the variables
run_campari(data_file = "NBU.dcd", key_file_input = "key_advanced.key", base_name = "NBU_advanced") 
sapphire_plot('PROGIDX_000000000001.dat')


# ------------------------------------------
# tutorial 11 - step 2 (Simple Clustering)
# ------------------------------------------
# take the variables in a function argument format from key file exemple
keys <- keywords_from_keyfile(key_file_input = "nbu_clustering.key", return_string_of_arguments = T, 
                              keyword_list_first = F, keyword_list = c(FMCSC_SEQFILE="nbu.in"))
print(keys) # you can copy paste it in the run_campari BUT IT IS BETTER provide the keyfile directly (see next)

# in the following we use the key_file_input to set the majority of the keys and we override SEQFILE directly from the function
run_campari(data_file = "NBU.dcd", key_file_input = "nbu_clustering.key", FMCSC_SEQFILE="nbu.in",
            base_name = "nbu_clustering")

# remember to use the following commad on the logfile and fycfile to extract the clustering summary (or use show_clustering_summary function)
# for i in `grep -F 'CLUSTER SUMMARY' -A 10 LOGFILE | tail -n 9 | awk '{print $3}'`; do (j=`echo -e "$i + 1" | bc`; head -n $j FYCFILE | tail -n 1 | awk '{printf("Step %10i: C1C2C3C4: %10.3f C3C2C1H11: %10.3f C2C3C4H41: %10.3f\n",$1/10,$2,$3,$4)}') done
show_clustering_summary(log_file = "nbu_clustering.log", fyc_file = "NBU.fyc")


# ------------------------------------------
# Step 4 - Tuning Options and Choices for Clustering
# ------------------------------------------
run_campari(data_file = "NBU.dcd", key_file_input = "nbu_clustering.key", base_name = "nbu_clustering_rmsd",
            FMCSC_CDISTANCE=5, #RMSD with pairwise alignment
            FMCSC_CFILE="carbons.lst", # you MUST write it with 1,2,3,4 (see yourcalc_END.pdb for numbering scheme)
            FMCSC_CRADIUS=0.02, # changing thresholds because we changed the distance
            FMCSC_CMAXRAD=1.0,
            FMCSC_SEQFILE="nbu.in") 

# it changed the results!
show_clustering_summary(log_file = "nbu_clustering_rmsd.log", fyc_file = "NBU.fyc")

# hierarchical clustering (as there are some more keywords to copy -> keywords_from_keyfile using key_file_is_keywords)
copied_string <- "FMCSC_CMODE 3 # use hierarchical clustering
FMCSC_CDISTANCE 1 # use dihedral angles to represent each snapshot
FMCSC_CRADIUS 30.0 # threshold radius (in degrees) of clusters
FMCSC_CMAXRAD 30.0 # threshold radius (in degrees) for truncated leader algorithm used for preprocessing
FMCSC_CCUTOFF 60.0 # cutoff distance for neighbor list
FMCSC_CLINKAGE 3 # linkage criterion (3 = mean)" # KEEP THIS FORMATTING!!
keywords_from_keyfile(copied_string, 
                      return_string_of_arguments = T, # this reformat the keywords to be directly copied using the keyboard
                      key_file_is_keywords = T) # this specify that we are using a string instead of a keyfile in input. Keep the original formatting.

# now you can copy paste it in run_campari
run_campari(data_file = "NBU.dcd", key_file_input = "nbu_clustering.key", base_name = "nbu_hi_clustering",
            FMCSC_SEQFILE="nbu.in",
            FMCSC_CMODE=3, FMCSC_CDISTANCE=1, FMCSC_CRADIUS=30, FMCSC_CMAXRAD=30, FMCSC_CCUTOFF=60, FMCSC_CLINKAGE=3,
            FMCSC_CCOLLECT=10) # the algorithm is not scalable

# remember that all this runs work in append. So only the last element assigned counts

#lets see the results
show_clustering_summary(log_file = "nbu_hi_clustering.log", fyc_file = "NBU.fyc", sub_sampling_factor = 10) # this is due to CCOLLECT

# ------------------------------------------
# Step 5 - Advanced algorithms
# ------------------------------------------
# Here we apply the method of Blöchliger et al. to study the metastable states sampled in our simulation.
# We start with a new key-file from scratch. As for the case of simple clustering, the following basic keywords are required:
run_campari(data_file = "NBU.dcd", key_file_input = "nbu_advanced.key") # , FMCSC_CPROGINDRMAX=10000
            
sapphire_plot('PROGIDX_000000000001.dat')


# lets do it with pseudo energy profile
copied_string <- "
FMCSC_CPROGINDSTART -2 # consistently use largest cluster or representative thereof for reference in profiles
FMCSC_CMSMCFEP 1 # request MFPT-based cFEP
FMCSC_CRADIUS 20.0 # threshold radius (in degrees) of clusters
FMCSC_CPROGINDRMAX 1500 # number of search attempts"
keywords_from_keyfile(copied_string, return_string_of_arguments = T, key_file_is_keywords = T)
run_campari(data_file = "nbu_traj.dcd", key_file_input = "key_advanced.key", 
            FMCSC_CPROGINDSTART=-2, FMCSC_CMSMCFEP=1, FMCSC_CRADIUS=20, FMCSC_CPROGINDRMAX=1500)

# ------------------------------------------
#           MST - run_campari in ncminer mode
# ------------------------------------------

# butane chain simulated data - variables init
# ------------------------------------------
library(bio3d)
# trj <- read.dcd("nbu_traj.dcd")

# to use ncminer we need to load fyc directly (dihedral angles handling not implemented)
trj <- fread("FYC.dat", header = F, skip = 1)
head(trj)
fread("head -n 1 nbu.fyc") # head of it
trj <- as.data.frame(trj[,-1])
trj <- sapply(trj, as.numeric) # always be sure that it is numeric!
hist(trj[,2]) # this should have 3 peaks per diheadral angle

run_campari(trj = trj,
            FMCSC_CPROGINDMODE=1, #mst
            FMCSC_CCOLLECT=5,
            FMCSC_CMODE=4,
            FMCSC_CDISTANCE=1, #rmsd without alignment 7 - dihedral distances need a complete analysis (pdb_format dcd pdb etc...) 
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
options(CampaRi.data_management = "R") # R handling !!! ATTENTION MEMORY PROBLEMS ARE POSSIBLE FOR BIG DATA-SETS (USE NETCDF HANDLING)
adjl <- mst_from_trj(trj = trj, mode = "fortran", normalize_d = FALSE, logging = F)
ret <- gen_progindex(adjl = adjl, snap_start = 10)
ret2 <- gen_annotation(ret, snap_start = 10)
sapphire_plot(sap_file = 'REPIX_000000000010.dat', title = "CAMPARI WRAPPER - MST", 
              timeline = F, ann_trace = F)



# ------------------------------------------
#            MST - Netcdf backend
# ------------------------------------------
options(CampaRi.data_management = "netcdf") # netcdf handling
mst_from_trj(trj = trj[1:30000,], mode = "fortran", normalize_d = T, logging = F, distance_method = 1)
ret <- gen_progindex(nsnaps = 30000, snap_start = 10)
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
