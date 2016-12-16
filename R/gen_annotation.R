#' @title Compute the annotation of the progrex index
#' @description
#'      \code{gen_annotation} is able to annotate the progrex index using the mean first passage time between two basins on the progrex index.
#'      This function will produce a file (e.g. \code{"REPIX_000000000001.dat"}) that contains all the information about progress index and annotations
#'      and it can be used to plot the sapphire plot using \code{sapphire_plot}.
#'
#' @param ret_data This input must be the output from the function \code{gen_progindex}.
#' @param local_cut_width This variable defines the width for a 3 states annotation of the progrex index.
#' @param snap_start Starting snapshot. It must be the same as in gen_progindex. Best to leave not assigned.
#' @details For details, please refer to the main documentation of the original campari software \url{http://campari.sourceforge.net/documentation.html}.
#'
#' @seealso 
#' \code{\link{gen_progindex}}, \code{\link{sapphire_plot}}.
#'
#' @examples 
#' adjl <- mst_from_trj(trj = matrix(rnorm(1000), nrow = 100, ncol = 10))
#' adjl2 <- contract_mst(adjl = adjl)
#' ret<-gen_progindex(adjl = adjl2)
#' gen_annotation(ret_data = ret, local_cut_width = 10)
#' 
#' 
#' @export gen_annotation
#' @useDynLib CampaRi

gen_annotation <- function(ret_data, local_cut_width = NULL, snap_start = NULL){
  warning("if you use an input different in format to the output of the other CampaRi functions there is an high probability of crashing")
  # local_cut_widtg is a (potentially) new setting for the width for the local cut
  n_breaks <- 0 #not a clue ?
  brklst <- array(as.integer(0),c(1))
  nsnaps <- ret_data$n_snaps
  if((is.null(local_cut_width))||
     (is.numeric(local_cut_width)&&(local_cut_width>nsnaps||local_cut_width<1))) {
    local_cut_width <- as.integer(floor(nsnaps-1)/2)
    warning(paste0(
      "Local cut width has not been assigned or it has been wrogly assigned.
       It will be set at: ",local_cut_width))
  }
  if(is.null(snap_start)||
     (is.numeric(snap_start)&&snap_start>nsnaps)) {
    snap_start <- as.integer(ret_data$starter)
    warning(paste0(
      "Starting snapshot has not been assigned or it has been wrogly assigned.
      It will be set at: ",snap_start))
  }
  o_progind <- ret_data$progind
  o_distv <- ret_data$distv
  o_invvec <- ret_data$invvec
  o_iv2 <- ret_data$iv2
  
  ret_data2 <- .Fortran("gen_manycuts",
                        n_snaps=as.integer(nsnaps),
                        start=as.integer(snap_start), #starting snapshot?
                        ntbrks2=as.integer(n_breaks), #?
                        pwidth=as.integer(local_cut_width),
                        setis=as.integer(o_progind),
                        distv=as.single(o_distv),
                        invvec=as.integer(o_invvec),
                        ivec2=as.integer(o_iv2),
                        trbrkslst=as.integer(brklst)) #?
  invisible(ret_data2)
}