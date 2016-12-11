# @title .onAttach
# @description  Load required data into gloval enviroment
# @keywords internal
.onAttach<- function (libname, pkgname){
  # if (!interactive()|| stats::runif(1) > 0.1) return()
  packageStartupMessage(paste0(
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
    CampaRi.data_filename = ""
  )
  toset <- !(names(op.CampaRi) %in% names(op))
  if(any(toset)) options(op.CampaRi[toset])
#   if(getOption("CampaRi.data_management")=='netcdf'){
#     nc_lib_dir <- "/usr/include/"
#     makevars_file <- paste0('MY_PKG_LIBS= -lnetcdff -I', nc_lib_dir, '
# MY_PKG_FFLAGS= -fbacktrace -fbounds-check -fcheck-array-temporaries -g
# mypackage_FFLAGS = $(FPICFLAGS) $(SHLIB_FFLAGS) $(FFLAGS)
# all: $(SHLIB)
# main_clu_adjl_mst.o: main_clu_adjl_mst.f90
#         $(FC) $(mypackage_FFLAGS) $(MY_PKG_FFLAGS) -c main_clu_adjl_mst.f90 -o main_clu_adjl_mst.o $(MY_PKG_LIBS)
# 
# utilities_netcdf.o: utilities_netcdf.f90
#         $(FC) $(mypackage_FFLAGS) $(MY_PKG_FFLAGS) -c utilities_netcdf.f90 -o utilities_netcdf.o $(MY_PKG_LIBS)
# 
# PKG_LIBS= -lnetcdff -I', nc_lib_dir, '
# ')
#     cat(makevars_file, file = paste0(".R/Makevars"))
#   }
  invisible() #no output from this function
}