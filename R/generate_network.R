#' @title Create a network for each snapshots 
#' @description
#'      \code{generate_network} creates a sequence of network (one for each snapshot) using a sliding window of trajectory (see \code{window}). 
#'      It is also possible to reduce the dimensionality using an overlapping reduction (see \code{overlapping_reduction})
#'
#' @param trj Input trajectory (variables on the columns and equal-time spaced snpashots on the row). It must be a \code{matrix} or a \code{data.frame} of numeric.
#' @param window
#' @param overlapping_reduction
#' 
#' 
#' @details #For more details, please refer to the main documentation of the original campari software \url{http://campari.sourceforge.net/documentation.html}.
#'
#' @return If no netcdf support is available the function will return a trj with vectorized variables
#' @seealso
#' \code{\link{adjl_from_progindex}}, \code{\link{gen_progindex}}, \code{\link{gen_annotation}}.
# @examples
#' 
#'
#' @export generate_network
#' @import wgcna

generate_network <- function(trj, window = NULL, overlapping_reduction = NULL){
  # initialising the variables
  trj_out <- matrix(NA, nrow = nrow(trj), ncol = ncol(trj)*ncol(trj)) 
  mode2 <- "R"
  if(is.character(mode2)&&mode2 == "R"){
    if(dim(trj)[1] > 1000) warning('the computation could be incredibly long. We advise to use mode = "fortran" option')
    
    # Initiate cluster
    # cl <- makeCluster(n_cores)
    # clusterExport(cl, "trj")
    # cl <- makeCluster(mc <- getOption("cl.cores", 4))
    # clusterExport(cl=cl, varlist=c("text.var", "ntv", "gc.rate", "pos"))
    # clusterEvalQ(cl, library(rms))

    for(i in 1:dim(trj)[1]){
      # if(n_cores>1) clusterExport(cl, "i")
      if(i - window < 0)
        trj_out[i,] <- as.vector(WGCNA::adjacency(trj[1:(i+window),]))
      else if(window + i > nrow(trj))
        trj_out[i,] <- as.vector(WGCNA::adjacency(trj[(i-window):dim(trj)[1],]))
      else
        trj_out[i,] <- as.vector(WGCNA::adjacency(trj[(i-window):(i+window),]))
    }
    
    # stopCluster(cl)
    return(trj_out)
    # warning('No MST made. Please use igraph package and in particular "mst" function in order to continue the analysis')
  
  }
  
}