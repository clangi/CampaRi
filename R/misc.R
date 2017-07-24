# @title .onAttach
# @description  Load required data into gloval enviroment
# @keywords internal
.onAttach<- function (libname, pkgname){
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

# .onLoad <- function(libname, pkgname) {
#   op <- options()
#   op.CampaRi <- list(
#     CampaRi.data_management = "R"
#   )
#   toset <- !(names(op.CampaRi) %in% names(op))
#   if(any(toset)) options(op.CampaRi[toset])
#   # .setting_up_netcdf()
# #   if(getOption("CampaRi.data_management")=='netcdf'){
# #     nc_lib_dir <- "/usr/include/"
# #     makevars_file <- paste0('MY_PKG_LIBS= -lnetcdff -I', nc_lib_dir, '
# # MY_PKG_FFLAGS= -fbacktrace -fbounds-check -fcheck-array-temporaries -g
# # mypackage_FFLAGS = $(FPICFLAGS) $(SHLIB_FFLAGS) $(FFLAGS)
# # all: $(SHLIB)
# # main_clu_adjl_mst.o: main_clu_adjl_mst.f90
# #         $(FC) $(mypackage_FFLAGS) $(MY_PKG_FFLAGS) -c main_clu_adjl_mst.f90 -o main_clu_adjl_mst.o $(MY_PKG_LIBS)
# # 
# # utilities_netcdf.o: utilities_netcdf.f90
# #         $(FC) $(mypackage_FFLAGS) $(MY_PKG_FFLAGS) -c utilities_netcdf.f90 -o utilities_netcdf.o $(MY_PKG_LIBS)
# # 
# # PKG_LIBS= -lnetcdff -I', nc_lib_dir, '
# # ')
# #     cat(makevars_file, file = paste0(".R/Makevars"))
# #   }
#   invisible() #no output from this function
# }

.lt <- function(x) return(length(x))

.check_integer <- function(x) return(is.null(x) || !is.numeric(x) || x%%1 != 0) # MUST not be null

.print_consecutio <- function(itering, total_to_iter, tot_to_print = 10, other_to_print = "", timeit = T, time_first = NULL){
  state_to_print <- floor(((itering*1.0)/total_to_iter)*tot_to_print)
  white_not_to_print <- tot_to_print - state_to_print
  if(timeit && is.null(time_first))
    stop("If you want to time me you need to give me the starting of time (time_first).")
  if(state_to_print%%1 == 0){
    if(timeit){
      time_spent <- proc.time() - time_first
      time_spent <- time_spent["elapsed"] + time_spent["user.self"] + time_spent["sys.self"]
      time_spent <- round(as.numeric(time_spent), digits = 0)
      if(time_spent>120)
        time_needed <- paste0("(needs: ", round((time_spent*1.0/itering)*total_to_iter,digits = 0)," s)")
      else
        time_needed <- ""
      if(time_spent != 0)
        time_spent <- paste0(" Time spent: ", time_spent , " s")
      else
        time_spent <- ""
      other_to_print <- paste(time_spent, time_needed, other_to_print)
    }
    string_to_print <- "\r|"
    if(state_to_print != 0){
      for(eq in 1:state_to_print)
        string_to_print <- paste0(string_to_print, "=")
    }
    if(state_to_print != tot_to_print){
      for(emp in 1:white_not_to_print)
        string_to_print <- paste0(string_to_print, " ")
    }
    string_to_print <- paste0(string_to_print,"| ", floor((state_to_print*100)/tot_to_print), "%  ", other_to_print)
    cat(string_to_print, sep = "")
  } 
}