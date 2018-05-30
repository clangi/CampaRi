# @title .onAttach
# @description  Load required data into gloval enviroment
# @keywords internal
.onAttach<- function (libname, pkgname){
  packageStartupMessage(paste0(
    " ==============================================================\n", 
    "    \n", 
    "                      CAMPARI analysis tools                   \n", 
    "    \n", 
    "    \n", 
    "           ----------------------------------------           \n", 
    " Analysing time series.                 \n", 
    " Version: ", utils::packageVersion("CampaRi"), "\n", 
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

# short hand for length
.lt <- function(x) return(length(x))

# check for single integer value
.isSingleInteger <- function(x) {
  if(!is.numeric(x) || x%%1 != 0 || (is.null(dim(x)) && length(x) != 1) || (!is.null(dim(x))))
    return(FALSE)
  else 
    return(TRUE)
}
# check for single element (e.g. character)
.isSingleElement <- function(x) {
  if((is.null(dim(x)) && length(x) != 1) || (!is.null(dim(x))))
    return(FALSE)
  else 
    return(TRUE)
}

.isSingleNumeric <- function(x) return(.isSingleElement(x) && is.numeric(x))
.isSingleChar <- function(x) return(.isSingleElement(x) && is.character(x))

# This routine is able to print a loading bar within a loop to know the work done
.print_consecutio <- function(itering, total_to_iter, tot_to_print = 10, other_to_print = "", timeit = T, time_first = NULL){
  state_to_print <- floor(((itering*1.0)/total_to_iter)*tot_to_print)
  white_not_to_print <- tot_to_print - state_to_print
  if(timeit && is.null(time_first))
    stop("If you want to time me you need to give me the starting of time (time_first).")
  if(state_to_print%%1 == 0){
    if(timeit){
      time_spent <- proc.time() - time_first
      time_spent <- time_spent["elapsed"]
      time_spent <- round(as.numeric(time_spent), digits = 0)
      if(time_spent>120)
        time_needed <- paste0("(needs: ", round((time_spent*1.0/itering)*total_to_iter, digits = 0)," s)")
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

# check for install_campari()
.get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}

# binary search for true-false. It is 100 times faster than %in%
.BiSearch <- function(table, key, start.idx = 1, end.idx = length(table),
                      tol = .Machine$double.eps ^ 0.5,
                      check = TRUE) {
  # Takes sorted (in ascending order) vectors
  if (check) stopifnot(is.vector(table), is.numeric(table))
  m <- as.integer(ceiling((end.idx + start.idx) / 2)) # Midpoint
  if (table[m] > key + tol) {
    if (start.idx == end.idx) return(FALSE)
    Recall(table, key, start.idx = start.idx, end.idx = m - 1L, tol = tol, check = FALSE)
  } else if (table[m] < key - tol) {
    if (start.idx == end.idx) return(FALSE)
    Recall(table, key, start.idx = m + 1L, end.idx = end.idx, tol = tol, check = FALSE)
  } else return(TRUE)
}

# normalize
.normalize <- function(x, xmax = NULL, xmin = NULL) {
  # if(.isSingleElement(x)) return(1)
  if(is.null(xmax)) xmax <- max(x)
  if(is.null(xmin)) xmin <- min(x)
  if(xmin == xmax) return(x/xmax)
  else return((x*1.0 - xmin)/(xmax - xmin))
}

# cutting the 0s
.cut0 <- function(x) { x[x<0] <- 0; x }

# fitting a quadratic model
.denaturate <- function(yyy, xxx, polydeg = 7, plotit = FALSE){
  # yyy <- kin.pl
  # xxx <- seq(lpi)
  q.mod <- lm(yyy ~ poly(xxx, polydeg, raw=TRUE))
  new_yyy <- .normalize(yyy - stats::predict(q.mod))
  if(plotit){
    plot(xxx, yyy, type = 'l', ylim = c(0,1))
    abline(lm(yyy ~ xxx), col = 'darkblue') # linear
    lines(xxx, stats::predict(q.mod), col = "darkgreen", lwd = 2)
    lines(xxx, new_yyy, col = "red", lwd = 3)
  }
  return(new_yyy)
}

# checking function for output command (intern or not, i.e. talkative or not)
.check_cmd.out <- function(cmd.out, intern = FALSE){
  if(intern) {
    return(!is.null(attr(cmd.out, which = 'status')))
  }else{
    if(!.isSingleInteger(cmd.out)) stop('Something went wrong. The command output is not an integer even if intern is FALSE.')
    return(cmd.out > 0)
  } 
}

# simple double separator creator; something like the following.
# > a <- .get_adjusted_separator(' a ', 11)
# > cat(a$sep1, '\n', a$sep2, sep = '')
# ~~~~ a ~~~~
# ~~~~~~~~~~~
.get_adjusted_separators <- function(tit, length_it = 68){
  if(length_it%%2) length_it <- length_it + 1 # simple check on the lengths
  if(.lt(tit) > length_it) stop('Please insert correctly the final length. It was too short in comparison to the title length.')
  if(nchar(tit)%%2 != 0) nchar_sep <- length_it - 1
  else nchar_sep <- length_it
  half_it <- (nchar_sep - nchar(tit))/2
  stp1 <- ""
  for(u in 1:half_it) stp1 <- paste0(stp1, "~")
  stp1 <- paste0(stp1, tit, stp1)
  stp2 <- ""
  for(u in 1:nchar_sep) stp2 <- paste0(stp2, "~")
  stp2 <- paste0(stp2)
  return(list('sep1' = stp1, 'sep2' = stp2))
}

# create annotation vector from barrier vector # not used at the moment
.vec_from_barriers <- function(bar.vec, label.vec = NULL, end.point = NULL){
  
  stopifnot(all(end.point > bar.vec))
  if(is.null(label.vec)) label.vec <- 1:length(bar.vec)
  
  if(sum(bar.vec) != end.point) {
    stopifnot(!is.null(end.point))
    bar.vec <- sort(c(bar.vec, end.point))
  }
  bar.vec <- c(bar.vec[1], diff(bar.vec))
  
  stopifnot(sum(bar.vec) == end.point)
  stopifnot(.lt(label.vec) == .lt(bar.vec))
  
  return(c(unlist(sapply(label.vec, function(x) rep(x, bar.vec[x])))))
}
