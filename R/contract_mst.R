#' @title Compute the contraction of the leaves in the minimum spanning tree extreme branches
#' @description
#'      This function aims to reduce the fringe regions by collapsing the tree leaves 
#'      on their branches. In this way we have less grouping of fringe regions
#'
#' @param adjl A list of three elements: degree list, connectivity matrix and weights
#' @param n_fold Number of links contractions (folds)
#' @details For details, please refer to the main documentation of the original campari software \url{http://campari.sourceforge.net/documentation.html}
#'
#'
#' @return This function will return an adjacency list modified
#' @examples 
#' adjl <- mst_from_trj(trj = matrix(rnorm(1000), nrow = 100, ncol = 10))
#' adjl2 <- contract_mst(adjl = adjl)
#'
#' @export contract_mst
#' @useDynLib CampaRi
 
contract_mst <- function(adjl, n_fold = 1){
  # n_fold is number of leaf-folding operations (FMCSC_CPROGMSTFOLD)
  nsnaps <- nrow(adjl[[2]]) # Working only if the graph it is complete -crash with not connected components
  maxnb <- max(as.numeric(adjl[[1]])) # If we avoid 0s in the loading function this could be not working
  # output initialization
  output_fin <- list()
  o_istats <- array(as.integer(0),c(2))
  attr(adjl[[3]],"Csingle") <- TRUE
  output <- .Fortran("contract_mst", PACKAGE="CampaRi",
                       n_snaps=as.integer(nsnaps),
                       mnb=as.integer(maxnb),
                       alnbs=as.integer(adjl[[1]]), #why should I need this thing
                       alst=matrix(as.integer(adjl[[2]]), nrow = nsnaps, ncol = maxnb),
                       aldis=adjl[[3]],
                       nrnds=as.integer(n_fold),
                       istats=as.integer(o_istats))
  cat(sum(output$aldis==adjl[[3]])*100.0/(nsnaps*maxnb)); cat("% are the same")
  output_fin[[1]] <- output$alnbs
  output_fin[[2]] <- output$alst[,1:output$mnb]
  output_fin[[3]] <- output$aldis[,1:output$mnb]
  return(output_fin)
}
