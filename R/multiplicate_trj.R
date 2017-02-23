#' @title Multiplicate the trajectory variables
#' @description
#'      \code{multiplicate_trj} copy the future and past variables in the snapshot considered. This leads to a multiplication of the variables 
#'      for each snapshot, following the variable \code{window}.
#'
#' @param trj Input trajectory (variables on the columns and equal-time spaced snpashots on the row). It must be a \code{matrix} or a \code{data.frame} of numeric.
#' @param window  Number of snapshots taken before and after the considered snapshot.
#' @param overlapping_reduction
#' 
#' 
#' @details #For more details, please refer to the main documentation of the original campari software \url{http://campari.sourceforge.net/documentation.html}.
#'
#' @return It will return a modified trajectory matrix.
#' @seealso
#' \code{\link{adjl_from_progindex}}, \code{\link{gen_progindex}}, \code{\link{gen_annotation}}.
# @examples
#' 
#'
#' @export multiplicate_trj
#' 

multiplicate_trj <- function(trj, window = NULL, overlapping_reduction = NULL){
  
  # Checking input variables (again - it is also a stand alone function)
  if(!is.null(window) && (length(window) != 1 || !is.numeric(window) || window <= 0 || window > nrow(trj)/2))
    stop('The used window (distance 12) is too small or too big (must be less than half to have sense) or it is simply an erroneus insertion.')
  if((!is.null(overlapping_reduction) && (length(overlapping_reduction) != 1 ||!is.numeric(overlapping_reduction) ||
                                          overlapping_reduction <= 0 || overlapping_reduction > 1)))
    stop('The used overlapping_reduction is not correctly defined. It must be a number between 0 and 1.')
  
  # Long calculation warning
  if(dim(trj)[1] > 20000) warning('The network generation can be really long. Please consider multi-threads options of the WGCNA package.')
  
  # setting standard window size
  if(is.null(window)) window <- 3
  
  # initialising the variables
  trj_out <- matrix(NA, nrow = nrow(trj), ncol = ncol(trj)*window) 
  
  # setting window left and window right (if it is not divisible by 2)
  if(((window-1)/2)%%1 == 0) {
    window_r <- window_l <- (window-1)/2
  }else{
    window_r <- floor((window-1)/2)
    window_l <- ceiling((window-1)/2)
  } 
  
  # Control check on the unfeasibility of certain windows
  if(window >= ncol(trj)) 
    stop('The constructed dimensionality (n-plication of the dimentions) has encountered a value that is really high (>dim^2). It could be slow.')
  
  # Main transformation
  for(i in 1:dim(trj)[1]){
    if((i - window_l) <= 0)
      trj_out[i,] <- as.vector(t(trj[1:(window),]))
    else if((window_r + i) > nrow(trj))
      trj_out[i,] <- as.vector(t(trj[(dim(trj)[1]-window+1):dim(trj)[1],]))
    else
      trj_out[i,] <- as.vector(t(trj[(i-window_l):(i+window_r),]))
  }
  return(trj_out)
}