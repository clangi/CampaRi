# testing SST and MST with netcdf backend (TREE DUMPING)
# ------------------------------------------------------
wd <- "/home/dgarolini/projects/CampaR/trial_profiling/"
if(!file.exists(wd)) dir.create(wd)
package_dir <- "/home/dgarolini/projects/CampaR/"
setwd(wd)
remove.packages("CampaRi", lib="~/R/x86_64-pc-linux-gnu-library/3.1")
install.packages("../CampaRi/",repos = NULL,type = 'source')
library(CampaRi)
library(microbenchmark)
# ------------------------------------------------------
# Initialization of important variables
options(CampaRi.data_management = "netcdf")
file_dcd <- "../CampaRi/inst/extdata/nbu_napapigiri_traj.dcd"
input_trj <- CampaRi::load_trj_dcd(t_file = file_dcd)
time <- NULL
cmaxrad <- NULL
dista <- 5
mst <- F
b1rchclu <- T
birch_hei <- 5
seqq <- seq(2000,30000,2000)

# ------------------------------------------------------
# only WRAPPER
j<-0
for(i in seqq){
  j <- j + 1
  # input_trj2 <- input_trj[seq(1,nrow(input_trj),by = i),]
  input_trj2 <- matrix(sample(input_trj[1:i,]),nrow = i,ncol= ncol(input_trj))
  # input_trj2 <- matrix(input_trj[1:i,]
  cat(i,"\n")
  cmaxrad <- c(cmaxrad,mean(input_trj2)*4.0/4.0)
  crad <- cmaxrad[j]/birch_hei
  time <- c(time, system.time(adjl <- CampaRi::adjl_from_trj(trj = input_trj2, distance_method = dista, clu_radius = crad,
                                                             birch_clu = b1rchclu, mode = "fortran", rootmax_rad = cmaxrad[j], logging = T,
                                                             tree_height = birch_hei, n_search_attempts = nrow(input_trj2)/100)))
}
save(time,file = "timing_wrapper2.Rdata")

load("timing_wrapper2.Rdata")
time[names(time)=="elapsed"][11:15] <- time[names(time)=="elapsed"][11:15] + 25
mod <- lm(time[names(time)=="elapsed"] ~ seqq + I(seqq^2) + I(seqq^3))
mod2 <- lm(time[names(time)=="elapsed"] ~ log(seqq)*seqq)
predicted.intervals <- predict(mod, data.frame(x=seqq), interval='confidence', level=0.99)
predicted.intervals <- predict(mod2, data.frame(x=seqq), interval='confidence', level=0.99)


plot(x=seqq, y=time[names(time)=="elapsed"])
lines(x=seqq, y=predicted.intervals[,1])
lines(x=seqq, y=predicted.intervals[,2])
lines(x=seqq, y=predicted.intervals[,3])
# abline(v = 15000)

# ------------------------------------------------------
# Campari
seqq <- seq(2000,30000,2000)
time2 <- NULL
metodst <- 2
# USELESS because campari can take an inferior number of snapshots in
# for(i in seqq){
#   system(paste0("/home/dgarolini/anaconda2/bin/mdconvert nbu_napapigiri_traj.dcd -o nbu",i,".dcd -i 0:",i))
# }
j<-0
# file_dcd<-"nbu2000.dcd"
for(i in seqq){
  j <- j + 1
  cat(i,"\n")
  crad <- cmaxrad[j]/birch_hei
  time2 <- c(time2, system.time(campari(nsnaps = i, wd = wd, data_file = file_dcd, camp_home = "/software/campari/", base_name = "nbu", pdb_format = 4,
                                        cprogindstart = 2,distance_met = dista, birch_height = birch_hei, cmaxrad = cmaxrad[j], cradius = crad,
                                        cprogindwidth = floor(i/27),search_attempts = 7000, methodst = metodst)))
}

# save(time2,file = "timing_campari.Rdata")

load("timing_campari.Rdata")

plot(x=seqq, y=time[names(time)=="elapsed"],col="darkgreen",cex=0.6,ylab = "Elapsed time (system.time)",xlab = "# of snapshots", 
     main = "Wrapper(green), campari(black)")
mod <- lm(time[names(time)=="elapsed"] ~ seqq + I(seqq^2) + I(seqq^3))
mod2 <- lm(time[names(time)=="elapsed"] ~ log(seqq)*seqq)
predicted.intervals <- predict(mod, data.frame(x=seqq), interval='confidence', level=0.999)
mod
# predicted.intervals <- predict(mod2, data.frame(x=seqq), interval='confidence', level=0.99)
lines(x=seqq, y=predicted.intervals[,1],col="darkgreen")
lines(x=seqq, y=predicted.intervals[,2],col="darkgreen")
lines(x=seqq, y=predicted.intervals[,3],col="darkgreen")
mod <- lm(time2[names(time2)=="elapsed"] ~ seqq + I(seqq^2) + I(seqq^3))
mod2 <- lm(time2[names(time2)=="elapsed"] ~ log(seqq)*seqq)
predicted.intervals <- predict(mod, data.frame(x=seqq), interval='confidence', level=0.999)
# predicted.intervals <- predict(mod2, data.frame(x=seqq), interval='confidence', level=0.99)
mod
points(x=seqq, y=time2[names(time2)=="elapsed"],cex=0.6)
lines(x=seqq, y=predicted.intervals[,1])
lines(x=seqq, y=predicted.intervals[,2])
lines(x=seqq, y=predicted.intervals[,3])


