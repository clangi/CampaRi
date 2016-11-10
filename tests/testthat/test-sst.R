# testing SST
wd <- "/home/dgarolini/projects/CampaR/sst_trial/"
setwd(wd)
remove.packages("CampaRi", lib="~/R/x86_64-pc-linux-gnu-library/3.1")
install.packages("../CampaRi/",repos = NULL,type = 'source')
library(CampaRi)
file_dcd <- "../CampaRi/inst/extdata/NBU_250fs.dcd"
input_trj <- CampaRi::load_trj_dcd(t_file = file_dcd)
# CampaRi::write.pdb.d(x = input_trj, base_name = "NBU_250fs")
# library(bio3d)
# in2 <- bio3d::read.pdb("NBU250fs.pdb")
cmaxrad <- mean(input_trj)*4.0/4.0
birch_hei <- 5
crad <- cmaxrad/birch_hei
campari(nsnaps = nrow(input_trj), wd = wd, data_file = "flottiglie.pdb", camp_home = "/software/campari/", base_name = "flottiglie", pdb_format = 1,
        cprogindstart = 1,distance_met = 5,birch_height = birch_hei, cmaxrad = cmaxrad, cradius = crad,
        cprogindwidth = 100,search_attempts = NULL,methodst = 2)
zap_ggplot(sap_file = "PROGIDX_000000000001.dat",local_cut = T)

adjl <- CampaRi::adjl_from_trj(trj = input_trj, distance_method = 5, clu_radius = crad,
                               birch_clu = T, mode = "fortran", rootmax_rad = cmaxrad, logging = F,
                               tree_height = birch_hei,n_search_attempts = nrow(input_trj)/10)
ret <- CampaRi::gen_progindex(adjl = adjl, snap_start = 1)
CampaRi::gen_annotation(ret_data = ret,snap_start = 1,local_cut_width = 100)
zap_ggplot(sap_file = "REPIX_000000000000.dat",local_cut = T)