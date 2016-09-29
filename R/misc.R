# @title .onAttach
# @description  Load required data into gloval enviroment
# @keywords internal
.onAttach<- function (libname, pkgname){
  # if (!interactive()|| stats::runif(1) > 0.1) return()
  cat(paste0(
    " ==============================================================\n",
    "    \n",
    "                            CAMPARI                           \n",
    "    \n",
    "    \n",
    " ------------------------------------------------------------\n",
    " Analysing time series.                 \n",
    " Version: ",utils::packageVersion("CampaRi"),"\n",
    " ==============================================================\n"))
}

.onLoad <- function(libname, pkgname) {
  op <- options()
  op.CampaRi <- list(
    CampaRi.path = getwd(),
    CampaRi.data_management = "R",
    CampaRi.data_filename = "dsetf.h5"
  )
  toset <- !(names(op.CampaRi) %in% names(op))
  if(any(toset)) options(op.CampaRi[toset])
  # cat("All output will be written to the current working directory: ", getwd(),"\n")
  invisible() #??
}
# .onLoad <- function() {
#   
#   cur_directory <- getwd()
#   message("Welcome to the analysis tool from CAMPARI!")
#   cat("All output will be written to the current working directory: ", cur_directory,"\n\n\n")
#   if(data_management == "R") cat("Normal memory handling selected. Without hdf5 backend file management it will be difficult for R to handle big data-sets.")
#   else if(data_management == "h5pfc") cat("Selected data support: hdf5 with mpi support\n")
#   else if(data_management == "h5fc") cat("Selected data support: hdf5 without mpi support\n")
#   else stop("Invalid data management keyword inserted. Check the available methods on the guide.")
#           
#   if(data_management != "R"){
#     warning('The dumping filename will be "dsetf.h5". If already existent it will be overwritten')
#     cat("Checking for dependencies...")
#     command_loc <- system(paste0("which ",data_management))
#     if(command_loc=="") stop("No support for hdf5. Please check installation and correct linkage of the command to your enviroment.")
#   }
#   # assign("data_man",data_management, envir = as.environment("package:CampaRi"))
#   # environment(data_management) <- as.environment("package:CampaRi")
#   CampaRi_cache<-new.env(parent = as.environment("package:CampaRi"))
#   CampaRi_cache$data_manager <- data_management
#   if(any(search()=="CampaRi_cache")) detach("CampaRi_cache")
#   attach(CampaRi_cache)
#   # As loading is no more necessary this could be done when the package is attached. *TODO
# 
# 
# 
#   # system("cp CampaRi/fortran/rerun_PIX.f90 .",intern=TRUE,ignore.stderr=TRUE,wait=TRUE,input=NULL)
# 
#   # manually copy source from tools-directory (do this by hand)
#   # then compile rerun_PIX.f90 by (on Linux) doing:
#   # system("R CMD SHLIB rerun_PIX.f90",intern=TRUE,ignore.stderr=TRUE,wait=TRUE,input=NULL)
#   # system("R CMD SHLIB CampaRi/src_fortran/f_CampaRi.f90",intern=TRUE,ignore.stderr=TRUE,wait=TRUE,input=NULL)
#   # load shared object
#   # dyn.load("rerun_PIX.so")
# }
