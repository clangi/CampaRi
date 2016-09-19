# @title .onAttach
# @description  Load required data into gloval enviroment
# @keywords internal
.onAttach<- function (libname, pkgname){
  # if (!interactive()|| stats::runif(1) > 0.1) return()
  message(paste0(
    " ==============================================================\n",
    "    \n",
    "                            CAMPARI                           \n",
    "    \n",
    "    \n",
    " ------------------------------------------------------------\n",
    " Analysing time series.                 \n",
    " Version: ",utils::packageVersion("campackage"),"\n",
    " ==============================================================\n"))
}

.onLoad <- function(libname, pkgname) {
  op <- options()
  op.campackage <- list(
    campackage.path = getwd(),
    campackage.data_management = "R",
    campackage.data_filename = "dsetf.h5"
  )
  toset <- !(names(op.campackage) %in% names(op))
  if(any(toset)) options(op.campackage[toset])
  cat("All output will be written to the current working directory: ", getwd(),"\n")
  invisible()
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
#   # assign("data_man",data_management, envir = as.environment("package:campackage"))
#   # environment(data_management) <- as.environment("package:campackage")
#   campackage_cache<-new.env(parent = as.environment("package:campackage"))
#   campackage_cache$data_manager <- data_management
#   if(any(search()=="campackage_cache")) detach("campackage_cache")
#   attach(campackage_cache)
#   # As loading is no more necessary this could be done when the package is attached. *TODO
# 
# 
# 
#   # system("cp campackage/fortran/rerun_PIX.f90 .",intern=TRUE,ignore.stderr=TRUE,wait=TRUE,input=NULL)
# 
#   # manually copy source from tools-directory (do this by hand)
#   # then compile rerun_PIX.f90 by (on Linux) doing:
#   # system("R CMD SHLIB rerun_PIX.f90",intern=TRUE,ignore.stderr=TRUE,wait=TRUE,input=NULL)
#   # system("R CMD SHLIB campackage/src_fortran/f_campackage.f90",intern=TRUE,ignore.stderr=TRUE,wait=TRUE,input=NULL)
#   # load shared object
#   # dyn.load("rerun_PIX.so")
# }
