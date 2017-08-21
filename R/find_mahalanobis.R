#' @title Find the mahalanobis distance for a set of triplets
#' @description
#'      \code{find_mahalanobis} uses a group of vectors belonging to a cluster and a second group of vector belonging to another cluster and learns
#'      a distance matrix which maximise the distance between different clusters and minimize the distance between elements of the same cluster
#'
#' @param clu1 Elements of the first clustering. Needed shape: (features, number of vectors). n_vectors must be > of features.
#' @param clu1e1 Elements of the first clustering. Needed shape: (features, number of vectors). n_vectors must be > of features.
#' @param clu2 As above. The number of elements should be identical for clu1e1, clu1e2 and clu2
#' @param must_be_positive If FALSE, it will consider also negative values for the Mahalanobis distance (the convergence can be much longer).
# @param ... Various variables. Possible values are \code{c('wgcna_type', 'wgcna_power', 'wgcna_corOp')}.
#' 
#' @details For more details on the SAPPHIRE plot, please refer to the main documentation of the original 
#' campari software \url{http://campari.sourceforge.net/documentation.html}.
#'
#' @return It will return a matrix n_features*n_features.
#' @seealso
#' \code{\link{adjl_from_progindex}}, \code{\link{gen_progindex}}, \code{\link{gen_annotation}}.
# @examples
#' 
#' @export find_mahalanobis
#' @useDynLib CampaRi, .registration = TRUE

find_mahalanobis <- function(clu1e1, clu1e2, clu2, must_be_positive = TRUE){
  
  
  # -----------------------------
  # input checks
  if(!is.numeric(clu1e1))
    stop('clu1e1 must be a numeric.')
  if(!is.numeric(clu1e2))
    stop('clu1e2 must be a numeric.')
  if(!is.numeric(clu2))
    stop('clu2 must be a numeric.')
  if(is.null(dim(clu1e1)))
    stop('clu1e1 must be a matrix with n_features rows and n_elements columns')
  if(is.null(dim(clu1e2)))
    stop('clu1e2 must be a matrix with n_features rows and n_elements columns')
  if(is.null(dim(clu2)))
    stop('clu2 must be a matrix with n_features rows and n_elements columns')
  if(nrow(clu1e1) != nrow(clu1e2) || nrow(clu1e1) != nrow(clu2))
    stop('The input matrices must have same number of rows (features).')
  if(ncol(clu1e1) != ncol(clu1e2) || ncol(clu1e1) != ncol(clu2))
    stop('The input matrices must have same number of cols (elements).')
  if(!is.logical(must_be_positive))
    stop('must_be_positive must be a logical.')
  
  
  # Long calculation warning
  # if(dim(trj)[1] > 20000) warning('The network generation can be really long. Please consider multi-threads options of the WGCNA package.')
  # if(dim(trj)[2] > 50) warning('The network generation can create an exagerated number of variables')

  # pre-defining the variables for Fortran talker
  n_features <- nrow(clu1e1)
  n_elements <- ncol(clu1e1)
  clu1e1 <- matrix(as.single(clu1e1), nrow = n_features, ncol = n_elements)
  clu1e2 <- matrix(as.single(clu1e2), nrow = n_features, ncol = n_elements)
  clu2 <- matrix(as.single(clu2), nrow = n_features, ncol = n_elements)
  X_MA <- matrix(as.single(rep(0.0,n_features*n_features)), nrow = n_features, ncol = n_features)
  attr(clu1e1, "Csingle") <- TRUE
  attr(clu1e2, "Csingle") <- TRUE
  attr(clu2, "Csingle") <- TRUE
  attr(X_MA, "Csingle") <- TRUE
  
  #main fortran talker
  output <- .Fortran("find_mahalanobis_F", PACKAGE="CampaRi",
                     #input
                     clu1e1=clu1e1,
                     clu1e2=clu1e2,
                     clu2=clu2,
                     n_feature=as.integer(n_features),
                     n_elements=as.integer(n_elements),
                     must_be_positive=as.logical(must_be_positive),
                     #output
                     X_MA=X_MA)
  return(output$X_MA)
}
