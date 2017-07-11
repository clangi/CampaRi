#' @title From adjacency matrix to adjacency list
#' @description
#'      This function is able to transform a matrix of distances into an adjacency list that can be used for the analysis pipeline.
#'      Please remember that \code{gen_progindex} accepts only minimum spanning trees.
#'
#' @param adj_m Input matrix (adjacency matrix).
#'
#' @return A list of three elements: degree list, connectivity matrix and weights.
#'
#' @seealso
#' \code{\link{gen_progindex}}, \code{\link{mst_from_trj}}.
#' @export adjl_from_adjmat
adjl_from_adjmat<-function(adj_m){ #deprecated
  # extract the SST or MST from the output of the analysis already made with campari.
  # Here we will reconstruct a bit of the tree in order to be able to find again the MST/SST
  adjl_nmbrs<-c()
  adjl<-list()
  adjl_dis<-list()
  for (i in 1:nrow(adj_m)){
    tmp <- sort(adj_m[i,adj_m[i,]!=0], index.return = TRUE)
    adjl_nmbrs[i] <- length(tmp$ix)
    adjl[[i]] <- tmp$ix
    adjl_dis[[i]] <- tmp$x
  }
  adjl <- array(unlist(adjl),dim = c(length(adjl_nmbrs),max(adjl_nmbrs)))
  adjl_dis <- array(unlist(adjl_dis),dim = c(length(adjl_nmbrs),max(adjl_nmbrs)))
  return(list(array(adjl_nmbrs), adjl, adjl_dis))
}
