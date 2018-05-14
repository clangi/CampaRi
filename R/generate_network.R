#' @title Create a network for each snapshots 
#' @description
#'      \code{generate_network} creates a sequence of network (one for each snapshot) using a sliding window of trajectory (see \code{window}). 
#'      It is also possible to reduce the dimensionality using an overlapping reduction (see \code{overlapping_reduction})
#'
#' @param trj Input trajectory (variables on the columns and equal-time spaced snpashots on the row). It must be a \code{matrix} or a \code{data.frame} of numeric.
#' @param window Number of snapshots taken before and after the considered snapshot. 
#' @param overlapping_reduction Not yet supported (It will if a great number of snapshots will be considered in this analysis - snapshot selection)
#' @param pp_method 'path_euclidean', 'path_minkowski', 'path_manhattan', 'path_maximum', 'path_canberra', 'path_binary' will find the distance path to make
#' an hamiltonian cycle. 'svd' if a singular vector decomposition will be used to have the diagional of the d matrix as output. With 'SymmetricUncertainty' you will have
#' AdjacencySymmetricUncertainty=2*I(X;Y)/(H(X)+H(Y)) [method='MI_*' must be active]. 'tsne' will make a tsne transformation.
#' @param method Supported pairwise similarity measures are 'wgcna', 'binary', 'euclidean', 'maximum', 'canberra', 'minkowski', 'covariance', 
#' 'MI_MM'. If set 'none' the windows will be used only for the pp_method (works only with path_* or tsne/svd).
#' @param transpose_trj Defaults to F. If set T the junk of trj (i.e. with a specific window) is transposed so to infer a 
#' network with dimensions window*(window-1)/2
#' @param silent A logical value indicating whether the function has to remain silent or not. Default value is \code{FALSE}.
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
#' @importFrom WGCNA adjacency
#' @importFrom WGCNA mutualInfoAdjacency
#' @importFrom WGCNA enableWGCNAThreads
#' @importFrom WGCNA disableWGCNAThreads 
#' @importFrom data.table transpose
#' @importFrom PairViz find_path
#' @importFrom Rtsne Rtsne
#' @importFrom TSA periodogram
#' @importFrom minerva mine
#' @importFrom parallel detectCores
#' 
#' 
#' @export generate_network

generate_network <- function(trj, window = NULL, method = 'wgcna', pp_method = NULL, 
                             overlapping_reduction = NULL, transpose_trj = FALSE, silent = FALSE, ...){
 
  # checking additional variable
  input_args <- list(...)
  
  avail_extra_argoments <- c('wgcna_type', 'wgcna_power', 'wgcna_corOp',                 # correlation options (WGCNA)
                             'minkowski_p',                                              # power for minkowski
                             'cov_method',                                               # cov met
                             'tsne_dimensions', 'tsne_perplexity', 'tsne_maxIter',       # tsne opts
                             'dbg_gn', 'do_multithreads')      
  
  avail_methods <- c('none',                                                                                   # only post-proc
                     'wgcna',                                                                                  # cor
                     'binary', 'euclidean', 'manhattan', 'maximum', 'canberra', 'minkowski', #'mahalanobis',   # distances
                     'covariance',                                                                             # cov
                     'MI', 'MI_MM', 'MI_ML', 'MI_shrink', 'MI_SG',                                             # MI based
                     'MIC', 'MAS', 'MEV', 'MCN', 'MICR2', 'GMIC', 'TIC',                                       # mine based measures (minerva)
                     'fft')                                                                                    # fft 
  
  which_are_distances <- c('binary', 'euclidean', 'manhattan', 'maximum', 'canberra', 'minkowski')
  MIC_mets <- c('MIC', 'MAS', 'MEV', 'MCN', 'MICR2', 'GMIC', 'TIC')
  wgcna_dep_mets <- c('wgcna', 'binary', 'euclidean', 'manhattan', 'maximum', 'canberra', 'minkowski', 'MI', 'MI_MM', 'MI_ML', 'MI_shrink', 'MI_SG')
  
  # default handling of extra inputs
  if(any(!(names(input_args) %in% avail_extra_argoments))) 
    stop('There is a probable mispelling in one of the inserted variables. Please check the available extra input arguments.')
  if(!('wgcna_type' %in% names(input_args))) wgcna_type <- 'unsigned' else wgcna_type <- input_args[['wgcna_type']]
  if(!('wgcna_power' %in% names(input_args))) wgcna_power <- 2 else wgcna_power <- input_args[['wgcna_power']]
  if(!('wgcna_corOp' %in% names(input_args))) wgcna_corOp <- "pearson" else wgcna_corOp <- input_args[['wgcna_corOp']]
  if(!('minkowski_p' %in% names(input_args))) minkowski_p <- 3 else minkowski_p <- input_args[['minkowski_p']]
  if(!('cov_method' %in% names(input_args))) cov_method <- 'pearson' else cov_method <- input_args[['cov_method']]
  if(!('tsne_dimensions' %in% names(input_args))) tsne_dimensions <- 2 else tsne_dimensions <- input_args[['tsne_dimensions']]
  if(!('tsne_perplexity' %in% names(input_args))) tsne_perplexity <- 30 else tsne_perplexity <- input_args[['tsne_perplexity']]
  if(!('tsne_maxIter' %in% names(input_args))) tsne_maxIter <- 500 else tsne_maxIter <- input_args[['tsne_maxIter']]
  if(!('dbg_gn' %in% names(input_args))) dbg_gn <- F else dbg_gn <- input_args[['dbg_gn']]
  if(!('do_multithreads' %in% names(input_args))) do_multithreads <- F else do_multithreads <- input_args[['do_multithreads']]
  
  # checks for fundamental variables: trj, window, method, (overlapping reduction)
  # trj checks ------------------------------------------------------------------------------------------------------------------------------
  
  # input trj checks
  if(!is.data.frame(trj)){
    if(!is.character(trj)) stop('The trj must be the tsv of the trj or a data.frame.')
    if(!silent) cat("Reading trj file...\n")
    trj <- as.matrix(data.table::fread(file = trj, data.table = F)) 
  }
  
  # long calculation warning
  if(dim(trj)[1] > 20000) warning('The network generation can be really long. Please consider multi-threads options of the WGCNA package.')
  if(dim(trj)[2] > 50) warning('The network generation can create an exagerated number of variables')
  
  # check for NA in trj
  if(any(is.na(trj))) stop('There are NA values in the input trajectory')
    
  # window checks if inserted
  if(!is.null(window) && .isSingleInteger(window) && (window <= 3 || window > nrow(trj)/2))
    stop('The used window (distance 12) is too small or too big (must be less than half to have sense) or it is simply an erroneus insertion.')
  
  # setting standard window size
  if(is.null(window)) window <- round(nrow(trj)/100)

  # setting window left and window right (if it is not divisible by 2)
  if(((window-1)/2)%%1 == 0) {
    window_r <- window_l <- (window-1)/2
  }else{
    window_r <- floor((window-1)/2)
    window_l <- ceiling((window-1)/2)
  } 
  
  # method checks
  if(!is.null(method) && (!is.character(method) || length(method) != 1)) stop('method must be a single character.')
  if(!(method %in% avail_methods)) stop(paste0('The method must be chosen between', paste0(avail_methods, collapse = ' '),'.'))
  
  # checking MI method
  if(method == 'MI') method <- 'MI_MM'
  
  # enabling the multi-threading if using wgcna
  if(method %in% wgcna_dep_mets && do_multithreads) WGCNA::enableWGCNAThreads(nThreads = parallel::detectCores() - 1)
  if(method %in% MIC_mets && do_multithreads) ncores <- parallel::detectCores() - 1 else ncores <- 1
  
  # checks for logicals: transpose_trj, silent
  if(!is.logical(transpose_trj)) stop('transpose_trj must be a logical value.')
  if(!is.logical(silent)) stop('silent must be a logical value.')
  
  # overlapping reduction checks - not available
  if(!is.null(overlapping_reduction)) stop('overlapping_reduction functionality is not implemented still. Not use it.')
  # if((!is.null(overlapping_reduction) && (length(overlapping_reduction) != 1 ||!is.numeric(overlapping_reduction) ||
  #                                         overlapping_reduction <= 0 || overlapping_reduction > 1)))
  #   warning('The used overlapping_reduction is not correctly defined. It must be a number between 0 and 1.')
  
  # checks for available pp_method
  # pp_method checks ------------------------------------------------------------------------------------------------------------------------------
  avail_post_proc_met <- c('path_euclidean', 'path_minkowski', 'path_manhattan', 'path_maximum', 'path_canberra', 'path_binary', # find_path
                           'svd', 
                           'tsne',
                           'SymmetricUncertainty',
                           'amplitude', 'amplitude_maxfreq', 'maxfreq')
  
  # general post processing checks
  if(!is.null(pp_method) && !.isSingleChar(pp_method)) stop('pp_method must be a single character.')
  if(!is.null(pp_method) && !(pp_method %in% avail_post_proc_met)) 
    stop(paste0( "pp_method must be choosen between ", paste0(avail_post_proc_met, collapse = ' '))) 
  
  # post processing checks for MI
  if(!is.null(pp_method) && grepl(pp_method, pattern = 'SymmetricUncertainty', fixed = T) && !grepl(method, pattern = 'MI_', fixed = T)){
    warning('You did not insert any Mutual Information method (MI) but you are trying to use SymmetricUncertainty. The method will be set to "MI".')
    method <- 'MI_MM'
  }
  
  # post processing checks for fft
  if(!is.null(pp_method) && pp_method %in% c('amplitude', 'amplitude_maxfreq', 'maxfreq')){
    if(method != 'fft'){
      warning('You inserted pp_methods which depends on the fourier transformation of the dataset. "fft" method set automatically.')
      method <- 'fft'
      }
  }
  
  # checks if you inserted the pp_method when method is none (you do only post processing on windows)
  if(method == 'none' && is.null(pp_method)) stop('It is needed a pp_method if you want to use directly the window (without any network construction).')
  
  # checks for extra inputs 
  # extra input (...) checks ------------------------------------------------------------------------------------------------------------

  # wgcna specific inputs
  if(!.isSingleInteger(wgcna_power)) stop('wgcna_power must be a single integer.')
  if(!.isSingleChar(wgcna_type)) stop('wgcna_type must be a single string.')
  if(!.isSingleChar(wgcna_corOp)) stop('wgcna_corOp must be a single string.')
  if(wgcna_corOp == 'pearson') wgcna_corOp <- "use = 'p'"
  else if(wgcna_corOp == 'spearman') wgcna_corOp <- "use = 'p', method = 'spearman'"
  else stop('Inserted wgcna_corOp not valid. choose between pearson and sperman.')
  
  # minkowski specific checks
  if(!.isSingleInteger(minkowski_p)) stop('minkowski_p must be a single integer.')
  else minkowski_p <- paste0(", p = '", minkowski_p, "'") 
  
  # specifiic checkings: the cov_method
  if(!.isSingleChar(cov_method)) stop('cov_method must be a single string.')
  if(!cov_method %in% c('pearson', 'spearman', 'kendall')) stop('cov_method not supported.')

  # specific checks for tsne
  if(!is.null(pp_method) && pp_method == 'tsne'){
    if(!.isSingleInteger(tsne_dimensions)) stop('tsne_dimensions must be a single integer')
    if(!.isSingleInteger(tsne_perplexity)) stop('tsne_perplexity must be a single integer')
    if(!.isSingleInteger(tsne_maxIter)) stop('tsne_maxIter must be a single integer')
    if(tsne_perplexity > (window - 1)/3){
      warning('tsne_perplexity too big. It will be set to tsne_perplexity <- (window/2 - 1)/3')
      tsne_perplexity <- (window/2 - 1)/3
    }  
  }
  
  # initialising the variables
  # __init ============================================================================================================
  
  # transpose option 
  if(transpose_trj){
    trj_out <- matrix(NA, nrow = nrow(trj), ncol = ((window*(window-1))/2)) 
    
  # fft methods
  }else if(method == 'fft'){
    if(!is.null(pp_method) && pp_method %in% c('amplitude', 'amplitude_maxfreq', 'maxfreq')){
      if(pp_method == 'amplitude')
        trj_out <- matrix(NA, nrow = nrow(trj), ncol = ncol(trj))
      else if(pp_method == 'maxfreq')
        trj_out <- matrix(NA, nrow = nrow(trj), ncol = ncol(trj))
      else if(pp_method == 'amplitude_maxfreq')
        trj_out <- matrix(NA, nrow = nrow(trj), ncol = 2*ncol(trj))
    }else{
      tester <- TSA::periodogram(trj[1:window, 1], plot = F)
      tester <- length(tester$spec)
      trj_out <- matrix(NA, nrow = nrow(trj), ncol = tester*ncol(trj))   
    }
    
  # standard run methods (no post processing)
  }else if(is.null(pp_method) || 
           (grepl(method, pattern = 'MI', fixed = T) && (!grepl(pp_method, pattern = 'path_', fixed = T) || !grepl(pp_method, pattern = 'svd', fixed = T)))){
    trj_out <- matrix(NA, nrow = nrow(trj), ncol = ((ncol(trj)-1)*ncol(trj)/2)) 
    
  # path construction for final matrix
  }else if(!is.null(pp_method) && grepl(pp_method, pattern = 'path_', fixed = T)){
    net_test <- WGCNA::adjacency(trj[1:100,])
    standard_path_length <- length(c(PairViz::find_path(net_test, path = function(xu) dist(xu, method = strsplit(pp_method, split = '_')[[1]][2]))))
    trj_out <- matrix(NA, nrow = nrow(trj), ncol = standard_path_length) 
  
  # final svd
  }else if(!is.null(pp_method) && grepl(pp_method, pattern = 'svd', fixed = T)){
    trj_out <- matrix(NA, nrow = nrow(trj), ncol = ncol(trj))
    
  # final tse
  }else if(!is.null(pp_method) && grepl(pp_method, pattern = 'tsne', fixed = T)){
    trj_out <- matrix(NA, nrow = nrow(trj), ncol = tsne_dimensions*window)
  
  # other options are not available
  }else{
    stop('impossible initilization of the variables. Please check the modes.')
  }

  # time keeping variables
  if(nrow(trj)*ncol(trj) > 50000) timing_it <- TRUE
  else timing_it <- FALSE
  
  if(dbg_gn) browser()
  
  # main transformation
  # __ loop ============================================================================================================
  
  # general loop
  if(!silent) cat('Network construction started (selected window: ', window, ').\n')
  for(i in 1:dim(trj)[1]){
    
    # adding time if necessary
    if(i==1 && timing_it) time1 <- proc.time()
    
    # print the bar if necessary
    if(!silent) .print_consecutio(if(i==1) 0 else i, dim(trj)[1], 70, timeit=timing_it, time_first = time1)
    
    # taking the piece of trj to transform
    if((i - window_l) <= 0) {                  # left border collision
      tmp_trj <- trj[1:(i + window_r),]
      if(transpose_trj || method == 'fft') tmp_trj <- rbind(tmp_trj, tmp_trj[1:(window - nrow(tmp_trj)),])
    }else if((window_r + i) > nrow(trj)){      # right border collision
      tmp_trj <- trj[(i - window_l):dim(trj)[1],]
      if(transpose_trj || method == 'fft') tmp_trj <- rbind(tmp_trj, tmp_trj[1:(window - nrow(tmp_trj)),])
    }else{                                     # usual center piece
      tmp_trj <- trj[(i - window_l):(i + window_r),]
    } 
    
    # transpose case
    if(transpose_trj) tmp_trj <- suppressWarnings(transpose(data.frame(tmp_trj))) # suppress the warnings for lost names
    
    # main adjacency constructor.
    # main --------------------------
    if(method != 'none'){
      
      # WGCNA
      if(method == 'wgcna'){
        built_net <- WGCNA::adjacency(tmp_trj, type = wgcna_type, corFnc = 'cor', power = wgcna_power, corOptions = wgcna_corOp)
      
      # covariance
      }else if(method == 'covariance'){
        built_net <- stats::cov(tmp_trj, method = cov_method)
      
      # mahalanobis (stat implementation) # todo
      # }else if(method == 'covariance'){
      #   built_net <- stats::cov(tmp_trj, method = cov_method)
      
      # MI based
      }else if(grepl(method, pattern = 'MI_', fixed = T)){
        while(dim(tmp_trj)[1] < 4) tmp_trj <- rbind(tmp_trj, tmp_trj)
        built_net <- WGCNA::mutualInfoAdjacency(tmp_trj, entropyEstimationMethod = strsplit(method, split = '_')[[1]][2])
      
      # construction based on Fourier Transform
      }else if(method == 'fft'){
        vec_freq <- c()
        for(fr in 1:ncol(trj)){
          if(is.null(pp_method) || (!is.null(pp_method) && pp_method != 'amplitude')) # just to avoid further calc
            freq_spectr <- periodogram(tmp_trj[,fr], plot = F)
          if(!is.null(pp_method) && pp_method %in% c('amplitude', 'amplitude_maxfreq', 'maxfreq')){
            if(pp_method == 'amplitude')
              vec_freq <- c(vec_freq, diff(range(tmp_trj[,fr])))
            else if(pp_method == 'maxfreq')
              vec_freq <- c(vec_freq, freq_spectr$freq[which(freq_spectr$spec == max(freq_spectr$spec))])
            else if(pp_method == 'amplitude_maxfreq')
              vec_freq <- c(vec_freq, diff(range(tmp_trj[,fr])), freq_spectr$freq[which(freq_spectr$spec == max(freq_spectr$spec))[1]])
          }else{
            vec_freq <- c(vec_freq, freq_spectr$spec)
          }
        }
      # distance inverse (e.g. Minkowskij similarity)
      }else if(method %in% which_are_distances){
        distOptions_manual <- paste0("method = '", method, "'")
        if(method == 'minkowski') distOptions_manual <- paste0(distOptions_manual, minkowski_p)
        built_net <- WGCNA::adjacency(tmp_trj, type = "distance", distOptions = distOptions_manual)
      
      # eventual misunderstanding
      }else if(method %in% MIC_mets){
        built_net <- minerva::mine(x = tmp_trj, n.cores = ncores, alpha = 0.4, C = 2)
        if(method == 'MIC') built_net <- built_net$MIC
        else if(method == 'MAS') built_net <- built_net$MAS
        else if(method == 'MEV') built_net <- built_net$MEV
        else if(method == 'MCN') built_net <- built_net$MCN
        else if(method == 'MICR2') built_net <- built_net$MICR2
        else if(method == 'GMIC') built_net <- built_net$GMIC
        else if(method == 'TIC') built_net <- built_net$TIC
      }else{
        stop('Something in the method construction went wrong. Please refer to the developers.')
      }
    # none case -> pass on the post processing without modifications
    }else{
      built_net <- tmp_trj
    }
          
    # Post-processing
    # pp_method -------------------------
    # special case: MI
    if(grepl(method, pattern = 'MI_', fixed = T)){
      if(!is.null(pp_method) && grepl(pp_method, pattern = 'SymmetricUncertainty', fixed = T))
        built_net_fin <- built_net$AdjacencySymmetricUncertainty
      else
        built_net_fin <- built_net$MutualInformation
    }else if(method != 'fft'){
      built_net_fin <- built_net
    }
    # other cases    
    if(method == 'fft'){
      trj_out[i,] <- vec_freq
    }else if(!is.null(pp_method) && grepl(pp_method, pattern = 'path_', fixed = T)){
      tmp_path <- c(PairViz::find_path(built_net_fin, path = function(xu) dist(xu, method = strsplit(pp_method, split = '_')[[1]][2])))
      if(standard_path_length != length(tmp_path)) stop('the path length changed during the analysis. This is not possible.')
      trj_out[i,] <- tmp_path
    }else if(!is.null(pp_method) && grepl(pp_method, pattern = 'svd', fixed = T)){
      trj_out[i,] <- svd(built_net_fin, nu = 0, nv = 0)$d
    }else if(!is.null(pp_method) && grepl(pp_method, pattern = 'tsne', fixed = T)){
      tmp_to_expand <- c(Rtsne::Rtsne(X = built_net_fin, dims = tsne_dimensions, perplexity = tsne_perplexity, verbose = FALSE, max_iter = tsne_maxIter)$Y)
      if(length(tmp_to_expand) < ncol(trj_out))
        tmp_to_expand <- c(rep(tmp_to_expand, 3))[ncol(trj_out)]
      trj_out[i,] <- tmp_to_expand
    }else{
      # Taking only the upper.tri
      trj_out[i,] <- built_net_fin[lower.tri(x = built_net_fin, diag = FALSE)]
    }
  }
  # __ end loop ============================================================================================================
  #
  if(method %in% wgcna_dep_mets && do_multithreads) WGCNA::disableWGCNAThreads()
  
  # Checking the creation of NAs in the inference
  # final output -------------------------
  n_na_trj <- sum(is.na(trj_out))
  if(n_na_trj > 0){
    warning('Attention: NA generated. Maybe too short window of time used? Fraction of NA: ', (n_na_trj*100/(nrow(trj_out)*ncol(trj_out))), ' %')
    if((n_na_trj*100/(nrow(trj_out)*ncol(trj_out))) > 5) 
      stop('ATTENTION!!! The generated NA overcame the threshold of 5%. Please consider alternative methods and vars.')
    if(!silent) cat('Assigning 0 to the generated NA.\n')
    trj_out[is.na(trj_out)] <- 0 
  }
  if(!silent) cat('Network construction completed.\n')
  invisible(list('trj_out' = trj_out, 'n_na_trj' = n_na_trj))
}
