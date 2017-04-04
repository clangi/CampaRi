#' @title Create a network for each snapshots 
#' @description
#'      \code{generate_network} creates a sequence of network (one for each snapshot) using a sliding window of trajectory (see \code{window}). 
#'      It is also possible to reduce the dimensionality using an overlapping reduction (see \code{overlapping_reduction})
#'
#' @param trj Input trajectory (variables on the columns and equal-time spaced snpashots on the row). It must be a \code{matrix} or a \code{data.frame} of numeric.
#' @param window Number of snapshots taken before and after the considered snapshot. 
#' @param overlapping_reduction Not yet supported (It will if a great number of snapshots will be considered in this analysis - snapshot selection)
#' @param transpose_trj Defaults to F. If set T the junk of trj (i.e. with a specific window) is transposed so to infer a network with dimensions window*(window-1)/2
#' @param ... Various variables. Possible values are \code{c('wgcna_type', 'wgcna_power', 'wgcna_corOp')}.
#' 
#' @details From WGCNA::adjacency: Correlation and distance are transformed as follows: for type = "unsigned", 
#' adjacency = |cor|^power; for type = "signed", adjacency = (0.5 * (1+cor) )^power; for type = "signed hybrid",
#' adjacency = cor^power if cor>0 and 0 otherwise; and for type = "distance", adjacency = (1-(dist/max(dist))^2)^power.
#' For more details on the SAPPHIRE plot, please refer to the main documentation of the original 
#' campari software \url{http://campari.sourceforge.net/documentation.html}.
#'
#' @return It will return a modified trajectory matrix.
#' @seealso
#' \code{\link{adjl_from_progindex}}, \code{\link{gen_progindex}}, \code{\link{gen_annotation}}.
# @examples
#' 
#'
#' @export generate_network
#' @importFrom WGCNA adjacency

generate_network <- function(trj, window = NULL, overlapping_reduction = NULL, transpose_trj = FALSE,... ){
 
  # checking additional variable
  input_args <- list(...)
  avail_extra_argoments <- c('wgcna_type', 'wgcna_power', 'wgcna_corOp')
  if(any(!(names(input_args) %in% avail_extra_argoments))) 
    warning('There is a probable mispelling in one of the inserted variables. Please check the available extra input arguments.')
  
  # Default handling
  if(!('wgcna_type' %in% names(input_args))) wgcna_type <- 'unsigned' else wgcna_type <- input_args[['wgcna_type']]
  if(!('wgcna_power' %in% names(input_args))) wgcna_power <- 1 else wgcna_power <- input_args[['wgcna_power']]
  if(!('wgcna_corOp' %in% names(input_args))) wgcna_corOp <- "use = 'p'" else wgcna_corOp <- input_args[['wgcna_corOp']]
  
  # Additional input (...) checks
  if(!is.character(wgcna_type) || !is.character(wgcna_corOp) || !is.numeric(wgcna_power)) stop('Inserted values for wgcna specifics not correct.')
  
  if(wgcna_corOp == 'pearson') wgcna_corOp <- "use = 'p'"
  else if(wgcna_corOp == 'spearman') wgcna_corOp <- "use = 'p', method = 'spearman'"
  
  # Checking input variables (again - it is also a stand alone function)
  if(!is.null(window) && (length(window) != 1 || !is.numeric(window) || window <= 3 || window > nrow(trj)/2))
    stop('The used window (distance 12) is too small or too big (must be less than half to have sense) or it is simply an erroneus insertion.')
  if((!is.null(overlapping_reduction) && (length(overlapping_reduction) != 1 ||!is.numeric(overlapping_reduction) ||
                                          overlapping_reduction <= 0 || overlapping_reduction > 1)))
    stop('The used overlapping_reduction is not correctly defined. It must be a number between 0 and 1.')
  if(!is.logical(transpose_trj))
    stop('transpose_trj must be a logical value.')
  
  # Long calculation warning
  if(dim(trj)[1] > 20000) warning('The network generation can be really long. Please consider multi-threads options of the WGCNA package.')
  if(dim(trj)[2] > 50) warning('The network generation can create an exagerated number of variables')
  
  # setting standard window size
  if(is.null(window)) window <- nrow(trj)/100
  
  # setting window left and window right (if it is not divisible by 2)
  if(((window-1)/2)%%1 == 0) {
    window_r <- window_l <- (window-1)/2
  }else{
    window_r <- floor((window-1)/2)
    window_l <- ceiling((window-1)/2)
  } 
  
  # Check for NA
  if(any(is.na(trj))) stop('There are NA values in the input trajectory')
    
  
  # initialising the variables
  trj_out <- matrix(NA, nrow = nrow(trj), ncol = ((ncol(trj)-1)*ncol(trj)/2)) 
  
  # Main transformation
  message('Network construction started.')
  for(i in 1:dim(trj)[1]){
    if((i - window_l) <= 0)
      tmp_trj <- trj[1:(i+window_r),]

    else if((window_r + i) > nrow(trj))
      tmp_trj <- trj[(i-window_l):dim(trj)[1],]
  
    else
      tmp_trj <- trj[(i-window_l):(i+window_r),]
    
    if(transpose_trj)
      tmp_trj <- transpose(tmp_trj)
    
    built_net <- WGCNA::adjacency(tmp_trj, type = wgcna_type, corFnc = 'cor', power = wgcna_power, 
                                    corOptions = wgcna_corOp)
    # Taking only the upper.tri  
    trj_out[i,] <- built_net[upper.tri(built_net)]
  }
  
  if(any(is.na(trj_out))){
    n_na_trj <- sum(is.na(trj_out))
    warning('Attention: NA generated. Probably it is due to too short window of time used. Fraction of NA: ', 
            (n_na_trj*100/(nrow(trj_out)*ncol(trj_out))), ' %')
    if((n_na_trj*100/(nrow(trj_out)*ncol(trj_out))) > 5) 
      warning('ATTENTION!!! The generated NA overcame the threshold of 5%. Please consider alternative methods and vars.')
    else message('Assigning 0 to the generated NA.')
    trj_out[is.na(trj_out)] <- 0 
  }
  
  message('Network construction completed.')
  return(trj_out)
}