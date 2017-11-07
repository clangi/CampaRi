#' @title Create a network for each snapshots 
#' @description
#'      \code{generate_network} creates a sequence of network (one for each snapshot) using a sliding window of trajectory (see \code{window}). 
#'      It is also possible to reduce the dimensionality using an overlapping reduction (see \code{overlapping_reduction})
#'
#' @param trj Input trajectory (variables on the columns and equal-time spaced snpashots on the row). It must be a \code{matrix} or a \code{data.frame} of numeric.
#' @param window Number of snapshots taken before and after the considered snapshot. 
#' @param overlapping_reduction Not yet supported (It will if a great number of snapshots will be considered in this analysis - snapshot selection)
#' @param post_processing_method 'path_euclidean', 'path_minkowski', 'path_manhattan', 'path_maximum', 'path_canberra', 'path_binary' will find the distance path to make
#' an hamiltonian cycle. 'svd' if a singular vector decomposition will be used to have the diagional of the d matrix as output. With 'SymmetricUncertainty' you will have
#' AdjacencySymmetricUncertainty=2*I(X;Y)/(H(X)+H(Y)) [method='MI_*' must be active]. 'tsne' will make a tsne transformation.
#' @param method Supported pairwise similarity measures are 'wgcna', 'binary', 'euclidean', 'maximum', 'canberra', 'minkowski', 'covariance', 
#' 'MI_MM'. If set 'none' the windows will be used only for the post_processing_method (works only with path_* or tsne/svd).
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
#' @export generate_network
#' @importFrom WGCNA adjacency
#' @importFrom WGCNA mutualInfoAdjacency
#' @importFrom data.table transpose
#' @importFrom PairViz find_path
#' @importFrom Rtsne Rtsne
#' @importFrom TSA periodogram

generate_network <- function(trj, window = NULL, method = 'wgcna', post_processing_method = NULL, overlapping_reduction = NULL, transpose_trj = FALSE, ...){
 
  # checking additional variable
  input_args <- list(...)
  avail_extra_argoments <- c('wgcna_type', 'wgcna_power', 'wgcna_corOp', 'minkowski_p', 'cov_method', 'tsne_dimensions', 'tsne_perplexity', 'tsne_maxIter')
  avail_methods <- c('wgcna', 'none', 
                     'binary', 'euclidean', 'manhattan', 'maximum', 'canberra', 'minkowski', 'mahalanobis',
                     'covariance',
                     'MI', 'MI_MM', 'MI_ML', 'MI_shrink', 'MI_SG',
                     'fft')
  
  # Default handling of extra inputs
  if(any(!(names(input_args) %in% avail_extra_argoments))) 
    stop('There is a probable mispelling in one of the inserted variables. Please check the available extra input arguments.')
  if(!('wgcna_type' %in% names(input_args))) wgcna_type <- 'unsigned' else wgcna_type <- input_args[['wgcna_type']]
  if(!('wgcna_power' %in% names(input_args))) wgcna_power <- 1 else wgcna_power <- input_args[['wgcna_power']]
  if(!('wgcna_corOp' %in% names(input_args))) wgcna_corOp <- "use = 'p'" else wgcna_corOp <- input_args[['wgcna_corOp']]
  if(!('minkowski_p' %in% names(input_args))) minkowski_p <- 3 else minkowski_p <- input_args[['minkowski_p']]
  if(!('cov_method' %in% names(input_args))) cov_method <- 'pearson' else cov_method <- input_args[['cov_method']]
  if(!('tsne_dimensions' %in% names(input_args))) tsne_dimensions <- 2 else tsne_dimensions <- input_args[['tsne_dimensions']]
  if(!('tsne_perplexity' %in% names(input_args))) tsne_perplexity <- 30 else tsne_perplexity <- input_args[['tsne_perplexity']]
  if(!('tsne_maxIter' %in% names(input_args))) tsne_maxIter <- 500 else tsne_maxIter <- input_args[['tsne_maxIter']]
  
  # Method checks
  if(!is.null(method) && (!is.character(method) || length(method) != 1))
    stop('method must be a single character.')
  if(!(method %in% avail_methods))
    stop('method not considered.')
  
  # Checking MI method
  if(method == 'MI')
    method <- 'MI_MM'
  
  # checks for available post_processing_method
  avail_post_proc_met <- c('path_euclidean', 'path_minkowski', 'path_manhattan', 'path_maximum', 'path_canberra', 'path_binary',
                           'svd', 'tsne',
                           'SymmetricUncertainty',
                           'amplitude', 'amplitude_maxfreq', 'maxfreq')
  if(!is.null(post_processing_method) && (!is.character(post_processing_method) || length(post_processing_method) != 1))
    stop('post_processing_method must be a single character.')
  if(!is.null(post_processing_method) && is.character(post_processing_method) &&
     !(post_processing_method %in% avail_post_proc_met))
    stop("post_processing_method must be choosen between these: 
c('path_euclidean', 'path_minkowski', 'path_manhattan', 'path_maximum', 'path_canberra', 'path_binary',
                           'svd', 'tsne',
                           'SymmetricUncertainty',
                           'amplitude', 'amplitude_maxfreq', 'maxfreq')")
  if(!is.null(post_processing_method) && grepl(post_processing_method, pattern = 'SymmetricUncertainty', fixed = T) && 
     !grepl(method, pattern = 'MI', fixed = T)){
    warning('You did not insert any Mutual Information option but you are trying to use SymmetricUncertainty. The method will be set to "MI".')
    method <- 'MI'
  }
  if(!is.null(post_processing_method) && post_processing_method %in% c('amplitude', 'amplitude_maxfreq', 'maxfreq')){
    if(method != 'fft'){
      warning('You inserted post_processing_methods which depends on the fourier transformation of the dataset. "fft" method set automatically.')
      method <- 'fft'
      }
  }
  
  if(method == 'none' && is.null(post_processing_method))
    stop('It is needed a post_processing_method if you want to use directly the window (without any network construction).')
  
  # Checking input variables 
  if(!is.null(window) && (length(window) != 1 || !is.numeric(window) || window <= 3 || window > nrow(trj)/2))
    stop('The used window (distance 12) is too small or too big (must be less than half to have sense) or it is simply an erroneus insertion.')
  
  # Additional input (...) checks
  # ---------------------------
  # wgcna specific inputs
  if(!.isSingleInteger(wgcna_power) || 
     !.isSingleElement(wgcna_type) || 
     !.isSingleElement(wgcna_corOp) || 
     !is.character(wgcna_type) || !is.character(wgcna_corOp) || !is.numeric(wgcna_power)) 
    stop('Inserted values for wgcna specifics inserted.')
  if(wgcna_corOp == 'pearson') wgcna_corOp <- "use = 'p'"
  else if(wgcna_corOp == 'spearman') wgcna_corOp <- "use = 'p', method = 'spearman'"
  
  # minkowski specific checks
  if(!is.numeric(minkowski_p) || minkowski_p %% 1 != 0) stop('minkowski_p must be an integer')
  else minkowski_p <- paste0(", p = '",minkowski_p,"'") 
  
  # specifiic checkings: the cov_method
  if(!cov_method %in% c('pearson', 'spearman', 'kendall'))
    stop('cov_method not supported.')

  # specific checks for tsne
  if(!is.null(post_processing_method) && post_processing_method == 'tsne'){
    if(!.isSingleInteger(tsne_dimensions))
      stop('tsne_dimensions must be a single integer')
    if(!.isSingleInteger(tsne_perplexity))
      stop('tsne_perplexity must be a single integer')
    if(!.isSingleInteger(tsne_maxIter))
      stop('tsne_maxIter must be a single integer')
    if(tsne_perplexity > (window - 1)/3){
      warning('tsne_perplexity too big. It will be set to tsne_perplexity <- (window/2 - 1)/3')
      tsne_perplexity <- (window/2 - 1)/3
    }  
  }
  
  if(!is.null(overlapping_reduction)) warning('overlapping_reduction functionality is not implemented still. Not use it.')
  # if((!is.null(overlapping_reduction) && (length(overlapping_reduction) != 1 ||!is.numeric(overlapping_reduction) ||
  #                                         overlapping_reduction <= 0 || overlapping_reduction > 1)))
  #   warning('The used overlapping_reduction is not correctly defined. It must be a number between 0 and 1.')
  
  # -----------------------------
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
  # ============================================================================================================
  if(transpose_trj){
    trj_out <- matrix(NA, nrow = nrow(trj), ncol = ((window*(window-1))/2)) 
  }else if(method == 'fft'){
    if(!is.null(post_processing_method) && post_processing_method %in% c('amplitude', 'amplitude_maxfreq', 'maxfreq')){
      if(post_processing_method == 'amplitude')
        trj_out <- matrix(NA, nrow = nrow(trj), ncol = ncol(trj))
      else if(post_processing_method == 'maxfreq')
        trj_out <- matrix(NA, nrow = nrow(trj), ncol = ncol(trj))
      else if(post_processing_method == 'amplitude_maxfreq')
        trj_out <- matrix(NA, nrow = nrow(trj), ncol = 2*ncol(trj))
    }else{
      tester <- TSA::periodogram(trj[1:window, 1], plot = F)
      tester <- length(tester$spec)
      trj_out <- matrix(NA, nrow = nrow(trj), ncol = tester*ncol(trj))   
    }
  }else if(is.null(post_processing_method) || (grepl(method, pattern = 'MI', fixed = T) && 
                                               (!grepl(post_processing_method, pattern = 'path_', fixed = T) || 
                                                !grepl(post_processing_method, pattern = 'svd', fixed = T)))){
    trj_out <- matrix(NA, nrow = nrow(trj), ncol = ((ncol(trj)-1)*ncol(trj)/2)) 
  }else if(!is.null(post_processing_method) && grepl(post_processing_method, pattern = 'path_', fixed = T)){
    net_test <- WGCNA::adjacency(trj[1:100,])
    standard_path_length <- length(c(PairViz::find_path(net_test, path = function(xu) dist(xu, method = strsplit(post_processing_method, split = '_')[[1]][2]))))
    trj_out <- matrix(NA, nrow = nrow(trj), ncol = standard_path_length) 
  }else if(!is.null(post_processing_method) && grepl(post_processing_method, pattern = 'svd', fixed = T)){
    trj_out <- matrix(NA, nrow = nrow(trj), ncol = ncol(trj))
  }else if(!is.null(post_processing_method) && grepl(post_processing_method, pattern = 'tsne', fixed = T)){
    trj_out <- matrix(NA, nrow = nrow(trj), ncol = tsne_dimensions*window)
  }else{
    stop('impossible initilization of the variables. Please check the modes.')
  }
  # ============================================================================================================
  
  
  # time keeping variables
  if(nrow(trj)*ncol(trj) > 50000)
    timing_it <- TRUE
  else
    timing_it <- FALSE
  
  # Main transformation
  # ============================================================================================================
  message('Network construction started (selected window: ', window, ').')
  for(i in 1:dim(trj)[1]){
    # adding time if necessary
    if(i==1&&timing_it)
      time1 <- proc.time()
    .print_consecutio(if(i==1) 0 else i, dim(trj)[1], 70, timeit=timing_it, time_first = time1)
    
    if((i - window_l) <= 0) {
      tmp_trj <- trj[1:(i + window_r),]
      if(transpose_trj || method == 'fft') tmp_trj <- rbind(tmp_trj, tmp_trj[1:(window - nrow(tmp_trj)),])
    }else if((window_r + i) > nrow(trj)){
      tmp_trj <- trj[(i - window_l):dim(trj)[1],]
      if(transpose_trj || method == 'fft') tmp_trj <- rbind(tmp_trj, tmp_trj[1:(window - nrow(tmp_trj)),])
    }else{
      tmp_trj <- trj[(i - window_l):(i + window_r),]
    } 
    
    if(transpose_trj)
      tmp_trj <- suppressWarnings(transpose(data.frame(tmp_trj))) # I suppress the warnings for lost names
    
    # main adjacency constructor.
    # --------------------------
    if(method != 'none'){
      
      # WGCNA
      if(method == 'wgcna'){
        built_net <- WGCNA::adjacency(tmp_trj, type = wgcna_type, corFnc = 'cor', power = wgcna_power, corOptions = wgcna_corOp)
      
        # covariance
      }else if(method == 'covariance'){
        built_net <- cov(tmp_trj, method = cov_method)
      
        # mahalanobis (stat implementation)
      }else if(method == 'covariance'){
        built_net <- cov(tmp_trj, method = cov_method)
      
        # MI based
      }else if(grepl(method, pattern = 'MI', fixed = T)){
        while(dim(tmp_trj)[1] < 4) tmp_trj <- rbind(tmp_trj, tmp_trj)
        built_net <- WGCNA::mutualInfoAdjacency(tmp_trj, entropyEstimationMethod = strsplit(method, split = '_')[[1]][2])
      
        # construction based on Fourier Transform
      }else if( method == 'fft'){
        vec_freq <- c()
        for(fr in 1:ncol(trj)){
          if(is.null(post_processing_method) || (!is.null(post_processing_method) && post_processing_method != 'amplitude')) # just to avoid further calc
            freq_spectr <- periodogram(tmp_trj[,fr], plot = F)
          if(!is.null(post_processing_method) && post_processing_method %in% c('amplitude', 'amplitude_maxfreq', 'maxfreq')){
            if(post_processing_method == 'amplitude')
              vec_freq <- c(vec_freq, diff(range(tmp_trj[,fr])))
            else if(post_processing_method == 'maxfreq')
              vec_freq <- c(vec_freq, freq_spectr$freq[which(freq_spectr$spec == max(freq_spectr$spec))])
            else if(post_processing_method == 'amplitude_maxfreq')
              vec_freq <- c(vec_freq, diff(range(tmp_trj[,fr])), freq_spectr$freq[which(freq_spectr$spec == max(freq_spectr$spec))[1]])
          }else{
            vec_freq <- c(vec_freq, freq_spectr$spec)
          }
        }
      
        # distance inverse (e.g. Minkowskij similarity)
      }else if(!(method %in% c('wgcna', 'covariance', 'MI_MM', 'MI_ML', 'MI_shrink', 'MI_SG'))){
        distOptions_manual <- paste0("method = '",method,"'")
        if(method == 'minkowski') distOptions_manual <- paste0(distOptions_manual, minkowski_p)
        built_net <- WGCNA::adjacency(tmp_trj, type = "distance", distOptions = distOptions_manual)
      
        # eventual misunderstanding
      }else{
        stop('Something in the method construction went wrong. Please refer to the developers.')
      }
    }else{
      built_net <- tmp_trj
    }
          
    # Post-processing
    # -------------------------
    # special case: MI
    if(grepl(method, pattern = 'MI', fixed = T)){
      if(!is.null(post_processing_method) && grepl(post_processing_method, pattern = 'SymmetricUncertainty', fixed = T))
        built_net_fin <- built_net$AdjacencySymmetricUncertainty
      else
        built_net_fin <- built_net$MutualInformation
    }else if(method != 'fft'){
      built_net_fin <- built_net
    }
    # other cases    
    if(method == 'fft'){
      trj_out[i,] <- vec_freq
    }else if(!is.null(post_processing_method) && grepl(post_processing_method, pattern = 'path_', fixed = T)){
      tmp_path <- c(PairViz::find_path(built_net_fin, path = function(xu) dist(xu, method = strsplit(post_processing_method, split = '_')[[1]][2])))
      if(standard_path_length != length(tmp_path))
        stop('the path length changed during the analysis. This is not possible.')
      trj_out[i,] <- tmp_path
    }else if(!is.null(post_processing_method) && grepl(post_processing_method, pattern = 'svd', fixed = T)){
      trj_out[i,] <- svd(built_net_fin, nu = 0, nv = 0)$d
    }else if(!is.null(post_processing_method) && grepl(post_processing_method, pattern = 'tsne', fixed = T)){
      tmp_to_expand <- c(Rtsne::Rtsne(X = built_net_fin, dims = tsne_dimensions, perplexity = tsne_perplexity, verbose = FALSE, max_iter = tsne_maxIter)$Y)
      if(length(tmp_to_expand) < ncol(trj_out))
        tmp_to_expand <- c(rep(tmp_to_expand, 3))[ncol(trj_out)]
      trj_out[i,] <- tmp_to_expand
    }else{
      # Taking only the upper.tri
      trj_out[i,] <- built_net_fin[lower.tri(x = built_net_fin, diag = FALSE)]
    }
  }
  # ============================================================================================================
  #
  
  # Checking the creation of NAs in the inference
  # -------------------------
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
