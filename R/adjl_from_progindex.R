#' @title Build the network from the already processed Progress Index file
#' @description
#'      \code{adjl_from_progindex} is able to use the output file from original campari software
#'      (or the output file from \code{\link{gen_annotation}}) in order to generate again the minimum spanning tree.
#'
#' @param prog_index_file Progress index file location. This should be of the kind \code{"PROGIDX_000000000001.dat"} (original campari)
#' or \code{"REPIX_000000000001.dat"} (CampaRi).
#'
#'
#' @return \code{adjl_from_progindex} will return a minimum spanning tree: degree list, connectivity matrix and weights
#' @details For details, please refer to the main documentation of the original campari software \url{http://campari.sourceforge.net/documentation.html}
#'
#' @seealso
#' \code{\link{gen_progindex}}, \code{\link{gen_annotation}}
#'
#' @examples
#' \dontrun{
#' adjl <- adjl_from_progindex(fil = "PROGIDX_000000000001.dat")
#' }
#'
#' @importFrom data.table fread
#' @export adjl_from_progindex

adjl_from_progindex <- function(prog_index_file){
  # extract the SST or MST from the output of the analysis already made with campari.
  # Here we will reconstruct a bit of the tree in order to be able to find again the MST/SST
  piOut <- data.frame(fread(prog_index_file))

  # number of snapshots
  nsnaps <- nrow(piOut)
  nbl<-piOut[,6] #number list
  itl<-piOut[,3] #index teo? list
  ditl<-piOut[,5] #distance of 6 from 3
  #max number of connections
  maxnb <- max(hist(breaks=seq(from=0.5,by=1,to=nsnaps+0.5),x=nbl,plot=FALSE)$counts) + 1
  # empty vector of future number of connections for each node in column 6
  treennb <- array(as.integer(0),c(nsnaps))
  # adjlist matrix with obvious limit in maxnb of connections
  treenbl <- array(as.integer(0),c(nsnaps,maxnb))
  # distance of each connection
  treedis <- array(as.single(0.0),c(nsnaps,maxnb))
  # slow
  for (i in 2:nsnaps) {
    treennb[nbl[i]] <- treennb[nbl[i]] + 1
    treenbl[nbl[i],treennb[nbl[i]]] <- itl[i]
    treedis[nbl[i],treennb[nbl[i]]] <- ditl[i]
    treennb[itl[i]] <- treennb[itl[i]] + 1
    treenbl[itl[i],treennb[itl[i]]] <- nbl[i]
    treedis[itl[i],treennb[itl[i]]] <- ditl[i]
    # list of breaks (must be passed at least of size 1, but n_breaks can be zero) eh?
  }
  return(list(treennb,treenbl,treedis))
}
