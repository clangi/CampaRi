#' @title Compute the progress index along with parent and distance-to-parent information
#' @description
#'      This function is able to create the progrex index of a time series from its minimum spanning tree of the distances.
#'      In order to understand how the progrex index is generated please refer to \url{http://www.sciencedirect.com/science/article/pii/S0010465513002038}.
#'
#'
#' @param adjl A list of three elements: degree list, connectivity matrix and weights. This input must be a minimum spanning tree, otherwise
#' the function will not be able to resolve the progress index in a unique solution. If backend netcdf data handling is active, this input can be left
#' empty. In this case please set the general option \code{CampaRi.data_management="netcdf"}.
#' @param nsnaps Number of snashots in the trajectory. If \code{adjl} is inserted this specification is not necessary. Instead if the netcdf backend data handling
#' is active this value will be necessary and precise in order to avoid unpredictable behaviours.
#' @param snap_start The snapshot from which the progrex index will start.
#' @param read_from_netcdf If \code{FALSE} the netcdf support will be used. The minimum spanning tree will be readed from file (see \code{\link{mst_from_trj}}).
#' @param mute If \code{TRUE} will mute the code messages
#' @return The output of this function MUST be used as input for the \code{\link{gen_annotation}} function.
#'
#' @details For details, please refer to the main documentation of the original campari software \url{http://campari.sourceforge.net/documentation.html}.
#'
#' @examples
#' adjl <- mst_from_trj(trj = matrix(rnorm(1000), nrow = 100, ncol = 10))
#' ret <- gen_progindex(adjl = adjl)
#'
#' @export gen_progindex
#' @useDynLib CampaRi

gen_progindex <- function(adjl=NULL, nsnaps = NULL, snap_start = 1, read_from_netcdf = FALSE, mute = FALSE){

  if(!is.logical(read_from_netcdf))
    stop('read_from_netcdf must be a logical.')
  if(!is.logical(mute))
    stop('mute must be a logical.')

  if(read_from_netcdf)
    if(!mute) message('The netcdf mode is active. A file called MST_DUMPLING.nc will be searched in the current directory (output from adjl_from_trj).')

  if(is.null(adjl))
    if(!read_from_netcdf) stop("The input must follow the standard used in 'R' data management protocol (output from adjl_from_trj).")
  if(!is.null(adjl)) {
    if(!is.list(adjl)||length(adjl[[1]])!=nrow(adjl[[2]])||length(adjl[[1]])!=nrow(adjl[[3]])||
       !is.numeric(adjl[[1]])||!is.numeric(adjl[[2]])||!is.numeric(adjl[[3]]))
      stop("This function accepts only a list of 3 elements: a degree list, an adjacency list and an associated list of distances.")
  }
  if(!is.null(adjl)&&is.null(nsnaps)) nsnaps <- length(adjl[[1]])
  if(is.null(adjl)){
    if(is.null(nsnaps)||!is.numeric(nsnaps))
      stop("It is possible to avoid using the direct insertion of the adjacency list only if there is netcdf backend file handling.
           For this you need to specify the number of snapshots anyway. If the number inserted is not corrected a crash is highly probable.")
  }
  if(!is.numeric(snap_start) || snap_start <= 0 || snap_start > nsnaps) stop("Starting snapshot is out of snapshot bounds.")
  #    this generates the new index along with parent and distance-to-parent information
  o_invvec <- array(as.integer(0),c(nsnaps+2)) # mistery
  o_iv2 <- array(as.integer(0),c(nsnaps)) # mistery
  o_progind <- array(as.integer(0),c(nsnaps))
  o_distv <- array(as.single(0.0),c(nsnaps))
  if(!read_from_netcdf){
    # Working only if the graph it is complete -crash with not connected components
    maxnb <- max(as.integer(adjl[[1]]))
    adjl[[2]] <- matrix(as.integer(adjl[[2]]), nrow = nsnaps, ncol = maxnb)
    attr(adjl[[3]],"Csingle") <- TRUE
    ret_data <- .Fortran("gen_progind_from_adjlst",
                         n_snaps=as.integer(nsnaps),
                         starter=as.integer(snap_start),
                         mnb=as.integer(maxnb),
                         alnbs=as.integer(adjl[[1]]),
                         alst=adjl[[2]],
                         aldis=adjl[[3]],
                         progind=as.integer(o_progind),
                         distv=as.single(o_distv),
                         invvec=as.integer(o_invvec),
                         iv2=as.integer(o_iv2),
                         mute_in=as.logical(mute))
  }else{
    if(!mute) cat("Going netcdf...\n")
    if(is.null(nsnaps) || !is.numeric(nsnaps)) stop("For netcdf data management you must write the number of snapshots in your trj")
    ret_data <- tryCatch(.Fortran("gen_progind_from_adjlst_r",
                                  n_snaps_in=as.integer(nsnaps),
                                  starter=as.integer(snap_start),
                                  progind=as.integer(o_progind),
                                  distv=as.single(o_distv),
                                  invvec=as.integer(o_invvec),
                                  iv2=as.integer(o_iv2),
                                  mute_in=as.logical(mute)),
                         error = function(e) stop('ERROR: the read_from_netcdf mode is active but the CampaRi package was built without netcdf support'))
  }
  return(ret_data)
}
