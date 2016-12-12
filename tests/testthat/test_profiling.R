# testing SST and MST with netcdf backend (TREE DUMPING)
# ------------------------------------------------------
wd <- "/home/dgarolini/projects/CampaR/trial_profiling/"
if(!file.exists(wd)) dir.create(wd)
package_dir <- "/home/dgarolini/projects/CampaR/"
setwd(wd)
remove.packages("CampaRi", lib="~/R/x86_64-pc-linux-gnu-library/3.1")
install.packages("../CampaRi/",repos = NULL,type = 'source')
library(CampaRi)
# ------------------------------------------------------
# Initialization of important variables
file_dcd <- "../CampaRi/inst/extdata/nbu_napapigiri_traj.dcd"
input_trj <- CampaRi::load_trj_dcd(t_file = file_dcd)
time <- 0
for(i in seq(100,10,-5)){
  input_trj <- input_trj[seq(1,nrow(input_trj),by = i),]
  cat(i,"\n")
  mst <- F
  cmaxrad <- mean(input_trj)*4.0/4.0
  birch_hei <- 5
  crad <- cmaxrad/birch_hei
  b1rchclu <- T
  time <- c(time, system.time(adjl <- CampaRi::adjl_from_trj(trj = input_trj, distance_method = 5, clu_radius = crad,
                                                             birch_clu = b1rchclu, mode = "fortran", rootmax_rad = cmaxrad, logging = T,
                                                             tree_height = birch_hei, n_search_attempts = nrow(input_trj)/100)))
#why does it crash?
}
