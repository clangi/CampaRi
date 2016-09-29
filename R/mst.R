  #' @title Compute the contraction of the extremes of the mst
#' @description
#'      This is a wonderful description(X)
#'
#' @param adjl A list of 3 adj lists: degree list, connectivity matrix and weights
#' @param n_fold Number of links contractions (folds)
#' Default: \code{0}
#'
#'
#' @return adjl modified
#'
#' @export contract_mst
#' @useDynLib CampaRi

# because the function is built in this way we need also treennb. Then we will keep the 3 output from
# adjl_from_piOut in the adjl variable
contract_mst <- function(adjl,n_fold=0){

  # n_fold is number of leaf-folding operations (FMCSC_CPROGMSTFOLD)
  nsnaps <- nrow(adjl[[1]]) # Working only if the graph it is complete -crash with not connected components
  maxnb <- max(as.numeric(adjl[[1]])) # If we avoid 0s in the loading function this could be not working
  # output initialization
  o_istats <- array(as.integer(0),c(2))
  if (!is.loaded('contract_mst')) {
    warning("The fortran function contract_mst was not loaded. Please reconsider the package loading.")
    dyn.load("src/rerun_PIX.so")
  }

  ret_data <- .Fortran("contract_mst", PACKAGE="CampaRi",
                       n_snaps=as.integer(nsnaps),
                       mnb=as.integer(maxnb),
                       alnbs=as.integer(adjl[[1]]), #why should I need this thing
                       alst=as.matrix(as.integer(adjl[[2]]),nsnaps,maxnb),
                       aldis=as.single(adjl[[3]]),
                       nrnds=as.integer(n_fold),
                       istats=as.integer(o_istats))

  o_istats <- ret_data$istats
  treedis2 <- ret_data$aldis
  cat(all(treedis2==adjl[[3]])) #true if cprogmstfold = 0, False if != 0, reaching the upper limit at 60% ==
  cat(sum(treedis2==adjl[[3]])*100.0/(nsnaps*maxnb)); cat("% are the same")
  rm(treedis2)


  return(ret_data)
}


#' @title Compute the new index along with parent and distance-to-parent information
#' @description
#'      This is a wonderful description(X)
#'
#' @param adjl A list of 3 adj lists: degree list, connectivity matrix and weights
#' @param snap_start Starting snapshot
#' Default: \code{0}
#'
#'
#' @return ret_data
#'
#' @export gen_progindex
#' @useDynLib CampaRi

gen_progindex<-function(adjl, snap_start = 0){
  warning("if you use an input different in format to the output of the other CampaRi functions there is an high probability of crashing")
  nsnaps <- length(adjl[[1]]) # Working only if the graph it is complete -crash with not connected components
  maxnb <- max(as.numeric(adjl[[1]]))
  
  if(!is.numeric(snap_start) || snap_start<0 || snap_start>nsnaps) stop("Starting snapshot is out of snapshot bounds.")

  #    this generates the new index along with parent and distance-to-parent information
  o_invvec <- array(as.integer(0),c(nsnaps+2)) # mistery
  o_iv2 <- array(as.integer(0),c(nsnaps)) # mistery
  o_progind <- array(as.integer(0),c(nsnaps))
  o_distv <- array(as.single(0.0),c(nsnaps))
  
  adjl[[2]] <- matrix(as.integer(adjl[[2]]),nsnaps,maxnb)
  if(is.null(attributes(adjl[[3]])$Csingle)) attr(adjl[[3]],"Csingle") <- TRUE
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
                       iv2=as.integer(o_iv2))

  # writeLines(sprintf("Check your directory, something happened! %g %f",o_progind[1],o_distv[12]))
  #    this generates the kinetic annotation function both globally and locally
  #    it also writes (to disk) the new output files
  return(ret_data)
}


#' @title Compute the new index along with parent and distance-to-parent information
#' @description
#'      This is a wonderful description(X)
#'
#' @param ret_data asdrubali che corrono (X)
#' @param local_cut_width local cut(X)
#' Default: \code{15000}
#' @param snap_start Starting snapshot
#' Default: \code{0}
#'
#' @return ret_data
#'
#' @export gen_annotation
#' @useDynLib CampaRi

gen_annotation<-function(ret_data, local_cut_width=NULL, snap_start = NULL){
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
    snap_start <- as.integer(0)
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
  return(ret_data2)
}