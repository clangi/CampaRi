#----------------------------------------------
#         WELCOME TO CAMPARI WR TESTING
#
# Davide G.
#----------------------------------------------
# Testing pre-processing feature selection based on different methods:
# 1. PCA to select most valuable loadings (components weights)




# Set and making working directory
#----------------------------------------------
w_dir <- getwd()
system("mkdir output_to_delete")
setwd("output_to_delete/")
w_dir <- getwd()
library(CampaRi)

# butane chain simulated data - variables init
# ------------------------------------------
library(bio3d)
trj <- read.dcd("../inst/extdata/NBU_1250fs.dcd")
options(CampaRi.data_management = "netcdf") # netcdf handling

birch_height <- 8
cluster_maxrad <- 120.0
cluster_rad <- 60.0

distance_measure <- 5
mst1_sst2 <- 2
birch_clustering <- T

# possible preprocessing tools
# ----------------------------

# some test on pca (princomp vs PCAproj)
library(mlbench)
library(caret)
library(ggfortify)
data(PimaIndiansDiabetes)
# Check importance levels
control <- trainControl(method="repeatedcv", number=10, repeats=3)
# train the model
model <- train(diabetes~., data=PimaIndiansDiabetes, method="lvq", preProcess="scale", trControl=control)
# estimate variable importance
importance <- varImp(model, scale=FALSE)
# summarize importance
print(importance)
plot(importance)

# PCA - classic
pca <- prcomp(PimaIndiansDiabetes[,1:8], scale. = T)
princomp2 <- princomp(PimaIndiansDiabetes[,1:8])
print(princomp2$loadings)

# PCA - pcaPP
pca_pr <- pcaPP::PCAproj(PimaIndiansDiabetes[,1:8],8)
print(pca_pr$loadings)


# plotting them
autoplot(pca, data = PimaIndiansDiabetes, colour = "diabetes", size=0.05, frame=T)
autoplot(princomp2, data = PimaIndiansDiabetes, colour = "diabetes", size=0.05, frame=T)
autoplot(pca_pr, data = PimaIndiansDiabetes, colour = "diabetes", size=0.05, frame=T)

# pre-feature_selection
trj_pca <- select_features(trj, feature_selection = 'pca')

trj_pca_robust <- select_features(trj, feature_selection = 'pca', pca_method = 'robust')

# ------------------------------------------------------
# Wrapper run - no wgcna
mst_from_trj(distance_method = distance_measure, clu_radius = cluster_rad,
             birch_clu = birch_clustering, mode = "fortran", rootmax_rad = cluster_maxrad, logging = F,
             tree_height = birch_height, n_search_attempts = 100,
             # changing vars
             trj = trj_pca_robust)
ret <- gen_progindex(nsnaps=nrow(trj), snap_start = 2)
# ret <- gen_progindex(adjl = adjl, snap_start = 2)
gen_annotation(ret_data = ret,snap_start = 2,local_cut_width = floor(nrow(trj)/27))
sapphire_plot(sap_file = "REPIX_000000000002.dat",local_cut = T,timeline = T, ann_trace = 2,title = "WRAPPER CAMPARI")

# ------------------------------------------------------
# Wrapper run - wgcna
mst_from_trj(distance_method = distance_measure, clu_radius = cluster_rad,
             birch_clu = birch_clustering, mode = "fortran", rootmax_rad = cluster_maxrad, logging = F,
             tree_height = birch_height, n_search_attempts = 100,
             # changing vars
             trj = trj_pca,
             pre_process = 'wgcna', window = 1000)
ret <- gen_progindex(nsnaps=nrow(trj), snap_start = 2)
# ret <- gen_progindex(adjl = adjl, snap_start = 2)
gen_annotation(ret_data = ret,snap_start = 2,local_cut_width = floor(nrow(trj)/27))
sapphire_plot(sap_file = "REPIX_000000000002.dat",local_cut = T,timeline = T, ann_trace = 2,title = "WRAPPER CAMPARI")

# ------------------------------------------------------
# Wrapper run - wgcna with internal (post wgnca) feature selection
mst_from_trj(distance_method = distance_measure, clu_radius = cluster_rad,
             birch_clu = birch_clustering, mode = "fortran", rootmax_rad = cluster_maxrad, logging = F,
             tree_height = birch_height, n_search_attempts = 100, 
             # changing vars
             trj = trj,
             pre_process = 'wgcna', window = 1000, 
             feature_selection = 'pca', n_princ_comp = 10)
ret <- gen_progindex(nsnaps=nrow(trj), snap_start = 2)
# ret <- gen_progindex(adjl = adjl, snap_start = 2)
gen_annotation(ret_data = ret,snap_start = 2,local_cut_width = floor(nrow(trj)/27))
sapphire_plot(sap_file = "REPIX_000000000002.dat",local_cut = T,timeline = T, ann_trace = 2,title = "WRAPPER CAMPARI")

# ------------------------------------------------------
# Wrapper run - wgcna with both external (trj_pca) and internal (post wgnca) feature selection
mst_from_trj(distance_method = distance_measure, clu_radius = cluster_rad,
             birch_clu = birch_clustering, mode = "fortran", rootmax_rad = cluster_maxrad, logging = F,
             tree_height = birch_height, n_search_attempts = 100, 
             # changing vars
             trj = trj_pca,
             pre_process = 'wgcna', window = 1000, 
             feature_selection = 'pca')
ret <- gen_progindex(nsnaps=nrow(trj), snap_start = 2)
# ret <- gen_progindex(adjl = adjl, snap_start = 2)
gen_annotation(ret_data = ret,snap_start = 2,local_cut_width = floor(nrow(trj)/27))
sapphire_plot(sap_file = "REPIX_000000000002.dat",local_cut = T,timeline = T, ann_trace = 2,title = "WRAPPER CAMPARI")

# -----------------------------------------------------
# plotting annotation
trj_fyc <- read.table("../../dihedral_original.fyc")
prog_index_table <- read.table("REPIX_000000000002.dat")[,c(3,4)]
xx <- seq(from=1, by=1, to=nrow(trj_fyc))
ann <- t((trj_fyc[prog_index_table[xx, 1], ] %% 360 < 120) + (trj_fyc[prog_index_table[xx, 1], ] %% 360 < 240) + 1)
title1 <- "NBU SST, 30k snapshots - WRAPPER"
sapphire_plot(sap_file = "REPIX_000000000002.dat", local_cut = T, timeline = T, ann_trace = ann, ann_names_L = c("CCCC","HCCC","CCCH"), title = title1)

# ----------------------------------------------------
# DELETING THE WORKING DIR where we just worked
setwd("..")
system("rm -rf output_to_delete")
