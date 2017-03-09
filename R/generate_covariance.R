#' @title Create a covariance matrix for each snapshots 
#' @description
#'      \code{generate_covariance} creates a sequence of covariance matrices (one for each snapshot) using a sliding window of trajectory (see \code{window}). 
#'      It is also possible to reduce the dimensionality using an overlapping reduction (see \code{overlapping_reduction})
#'
#' @param trj Input trajectory (variables on the columns and equal-time spaced snpashots on the row). It must be a \code{matrix} or a \code{data.frame} of numeric.
#' @param window Number of snapshots taken before and after the considered snapshot. 
#' @param overlapping_reduction Not yet supported (It will if a great number of snapshots will be considered in this analysis - snapshot selection)
#' 
#' 
#' @details 
#' For more details on the SAPPHIRE plot, please refer to the main documentation of the original 
#' campari software \url{http://campari.sourceforge.net/documentation.html}.
#'
#' @return It will return a modified trajectory matrix.
#' @seealso
#' \code{\link{adjl_from_progindex}}, \code{\link{gen_progindex}}, \code{\link{gen_annotation}}.
# @examples
#' 
#'
#' @export generate_covariance
#' @importFrom stats cov

generate_covariance <- function(trj, window = NULL, overlapping_reduction = NULL){
  
  # Checking input variables (again - it is also a stand alone function)
  if(!is.null(window) && (length(window) != 1 || !is.numeric(window) || window <= 3 || window > nrow(trj)/2))
    stop('The used window (distance 12) is too small or too big (must be less than half to have sense) or it is simply an erroneus insertion.')
  if((!is.null(overlapping_reduction) && (length(overlapping_reduction) != 1 ||!is.numeric(overlapping_reduction) ||
                                          overlapping_reduction <= 0 || overlapping_reduction > 1)))
    stop('The used overlapping_reduction is not correctly defined. It must be a number between 0 and 1.')
  
  # Long calculation warning
  if(dim(trj)[1] > 20000) warning('The network generation can be really long. Please consider multi-threads options of the WGCNA package.')
  if(dim(trj)[2] > 50) warning('The network generation can create an exagerated number of variables')
  
  # setting standard window size
  if(is.null(window)) window <- nrow(trj)/100
  
  # Check for NA
  if(any(is.na(trj))) stop('There are NA values in the input trajectory')
  
  # setting window left and window right (if it is not divisible by 2)
  if(((window-1)/2)%%1 == 0) {
    window_r <- window_l <- (window-1)/2
  }else{
    window_r <- floor((window-1)/2)
    window_l <- ceiling((window-1)/2)
  } 
  
  # initialising the variables
  trj_out <- matrix(NA, nrow = nrow(trj), ncol = (ncol(trj)*ncol(trj))) 
  
  # Main transformation
  for(i in 1:dim(trj)[1]){
    if((i - window_l) <= 0)
      cov_mat <- cov(trj[1:(i+window_r),])
    else if((window_r + i) > nrow(trj))
      cov_mat <- cov(trj[(i-window_l):dim(trj)[1],])
    else
      cov_mat <- cov(trj[(i-window_l):(i+window_r),])
    
    trj_out[i,] <- array(cov_mat/max(abs(cov_mat)))
  }
  return(trj_out)
}