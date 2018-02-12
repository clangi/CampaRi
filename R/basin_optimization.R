#' @title Basin optimization of the SBR
#' @description
#'     \code{basin_optimization} uses the information provided by the SAPPHIRE basin recognition in order to define the optimal partition of the data (uses in 
#'     particulare the number of bins as a tuning variable).
#'     
#' @param the_sap Name of the PROGIDX_<...>.dat or REPIX_<...>.dat file or Data frame with three columns containing, in order, the (sorted) progress index,
#'  the relative time indices and the cut function. 
#' @param basin_optimization_method Method name for the selection of precise basin. This is a costly procedure which is based on an heuristical obj function w
#' with the following subgroups.
#' \itemize{
#'    \item "\code{uniformity}" not available
#'    \item "\code{minimal_sd_entropy}" not available yet
#' }
#' @param how_fine_search This variable define the length of the vector of bins in which is splitted the 7-nbins_x interval. As the algorithm undergo binary 
#' search this is a fundamental variable to define the sensitivity of the search.
#' @param plot_basin_identification A logical value indicating whether to display the SAPPHIRE plot with the resulting partitions or not. 
#' Black partitions are the matched ones, blue ones derive only from the dynamic analysis and orange ones only from the kinetic analysis. 
#' The green curve is the kinetic annotation (black curve) where the parabolic shape has been subtracted, i.e. the actual curve used for 
#' the peaks identification. Default value is \code{FALSE}.
#' @param nbins_x Number of bins on  x-axis of the 2-D histogram. Default to sqrt(nrow(the_sap)).
#' @param nbins_y Number of bins on the y-axis of the 2-D histogram. Default to sqrt(nrow(the_sap)).
#' @param merge_cluster Logical that allow clusters to be merged automatically if consecutives
#' @param number_of_clusters if basin_optimization_method is active accordingly this must be set to integer.
#' @param force_matching Please refer to \code{\link{basin_recognition}} for further details about the match option.
#' @param silent A logical value indicating whether the function has to remain silent or not. Default value is \code{FALSE}.
#'      
#' @return The resulting score (0-1)
#' @examples
#' 
#' adjl <- mst_from_trj(trj = matrix(rnorm(1000), nrow = 100, ncol = 10))
#' ret<-gen_progindex(adjl = adjl)
#' gen_annotation(ret_data = ret, local_cut_width = 10)
#' \dontrun{
#' score_sapphire(the_sap = "PROGIDX_000000000001.dat", ann = rnorm(100))
#' CampaRi::basin_optimization(the_sap = "PROGIDX_000000000001.dat",  how_fine_search = 10, number_of_clusters = 2, force_matching = T, silent = F)
#' }
#' 
#' @details For details regarding the SAPPHIRE plot, please refer to the relative publications \url{http://www.nature.com/articles/srep06264}. 
#' Main documentation of the original campari software \url{http://campari.sourceforge.net/documentation.html}.
#' 
#' @importFrom data.table fread
#' @export basin_optimization

basin_optimization <- function(the_sap, basin_optimization_method = NULL, how_fine_search = 100,
                               plot_basin_identification = FALSE, nbins_x = NULL, nbins_y = nbins_x,
                               number_of_clusters = NULL, force_matching = FALSE, silent = FALSE){
  # general input checking
  if(!is.character(the_sap) && !is.data.frame(the_sap)) stop("the_sap must be a string or a data frame")
  if(is.character(the_sap) && (!all(grepl("PROGIDX", the_sap)) && !all(grepl("REPIX", the_sap)))) stop("Please provide a data name starting with 'PROGIDX' or 'REPIX'" )
  if(!is.null(basin_optimization_method) && !is.character(basin_optimization_method)) stop('basin_optimization_method must be a character')
  if(!is.null(nbins_x) && !.isSingleInteger(nbins_x)) stop('nbins_x must be a single integer')
  if(!is.null(nbins_y) && !.isSingleInteger(nbins_x)) stop('nbins_y must be a single integer')
  if(!.isSingleInteger(how_fine_search)) stop('how_fine_search must be a single integer')
  if(!is.null(number_of_clusters) && !.isSingleInteger(number_of_clusters)) stop('number_of_clusters must be a single integer')
  if(!is.logical(silent)) stop('silent must be a logical')
  if(!is.logical(plot_basin_identification)) stop('plot_basin_identification must be a logical')
  if(!is.logical(force_matching)) stop('force_matching must be a logical')
  
  # methods input check
  basin_optimization.opt <- c("uniformity", "minimal_sd_entropy")
  if(!is.null(basin_optimization_method) && !(basin_optimization_method[1] %in% basin_optimization.opt)) stop("basin_optimization_method method option not valid")
  
  # sapphire table loading
  if(!is.data.frame(the_sap))
    st <- as.data.frame(data.table::fread(the_sap, data.table = F))
  else
    st <- the_sap
  
  # checking the number of bins inserted
  if(is.null(nbins_x)) nbins_x <- round(sqrt(nrow(st)*10))
  if(is.null(nbins_y)) nbins_y <- round(sqrt(nrow(st)*10))
  if(nbins_x != nbins_y) {
    if(!silent) cat('For simplicity we do not allow yet to use different number of nx and ny (bins for x and y).')
    nbins_x <- nbins_y <- max(nbins_x, nbins_y)
  }
  
  
  # -------------------------------------------------------------------------------------------------------- number of clusters
  # If the number of clusters is inserted this method consist only into finding the right number of clusters 
  # without splitting the possible subdivision (of the same number of clustering)
  # and afterwards it is possible to score them per their relevance (e.g. uniformity or minimal entropy).
  #
  if(!is.null(number_of_clusters)){
    
    # checks
    if(!silent) cat('Automatic optimization of basins based on number of cluster inserted... \n')
    if(number_of_clusters < 2 && number_of_clusters >= nrow(st)) stop('Number of cluster too high or too low.')
    if(number_of_clusters > nbinsxy) stop('Please set the number of clusters lower than the number of initial nbins (nbinsxy).')
    
    # init
    lin_scale <- unique(round(seq(7, nbins_x, length.out = how_fine_search)))
    if(!silent) cat('Selected a convergence step of', how_fine_search, 'subdivisions, ranging from', nbins_x, 'to 7. \n')
    bisbr.out <- .BiSBR(st = st, ncl_found = 1, ncl_teo = number_of_clusters, start.idx = 1,
                        end.idx = length(lin_scale), lin_scale =  lin_scale, force_matching = force_matching, silent = silent)
    
    # if not found some handling - stop for now
    if(!bisbr.out$found) stop('We could not find optimal separation for the number of selected clusters. Try to put force_matching = F or more how_fine_search.')
    
    # final plot and saving results
    bas <- CampaRi::basins_recognition(st, nx = bisbr.out$nbins, new.dev = F, out.file = F, match = force_matching, plot = plot_basin_identification)
    if(!silent) cat('Number of (automatically) selected bins for the basin recognition step is', bisbr.out$nbins, '\n')
    fin_nbins <- bisbr.out$nbins
  }else{
    # initial linear scale for the basin_optimization_method
    lin_scale <- unique(round(seq(7, nbins_x, length.out = how_fine_search)))
  }
  
  # -------------------------------------------------------------------------------------------------------- optimization of basins
  if(!is.null(basin_optimization_method)){
    if(!is.null(number_of_clusters)){
      cat('Looking for plateu around', fin_nbins, 'nbins. In this way we can optimize even for identical number of clusters.\n')
      #!!
    }
    
    ################################ uniformity
    if(basin_optimization_method[1] == 'uniformity') stop('uniformity basin_optimization_method option not yet ready.')
    ################################ minimal_sd_entropy
    if(basin_optimization_method[1] == 'minimal_sd_entropy') stop('minimal_entropy basin_optimization_method option not yet ready.')
  }else{
    if(is.null(number_of_clusters)) stop('One between number_of_clusters and basin_optimization_method must be active to optimize the number of bins.')
  }
  
  # -------------------------------------------------------------------------------------------------------- final return
  return(list('Optimal_nbins' = bas$nbins, 'bas' = bas, 'ncl' = bisbr.out$ncl, 'idx' = bisbr.out$idx))
}


# abba <- unique(round(seq(2, nbinsxy, length.out = 100)))
# .BiSBR(st = st, ncl_found = 1, ncl_teo = 3, start.idx = 1, end.idx = length(abba), lin_scale =  abba, force_matching = F, silent = F)
# .BiSBR(st = st, ncl_found = 1, ncl_teo = 3, start.idx = 1, end.idx = length(abba), lin_scale =  abba, force_matching = T, silent = F)
# # testing for minimum number of bins without crashing
# CampaRi::basins_recognition(st[1:1000,], nx = 4, new.dev = F, out.file = F, plot = T)
# CampaRi::basins_recognition(st, nx = 12, new.dev = F, out.file = F, plot = T)
# CampaRi::basins_recognition(st, nx = 88, new.dev = F, out.file = F, plot = T, match = T)
# CampaRi::basins_recognition(st, nx = 97, new.dev = F, out.file = F, plot = T, match = T)
# CampaRi::basins_recognition(st, nx = 343, new.dev = F, out.file = F, plot = T, match = T)
# alu <- CampaRi::basins_recognition(st, nx = 905, new.dev = F, out.file = F, plot = T, match = T)



# help function for binary search on sapphire basin recognition
.BiSBR <- function(st, ncl_found, ncl_teo, start.idx = 1, end.idx = NULL, lin_scale = NULL, 
                   force_matching = FALSE, tol = .Machine$double.eps ^ 0.5, check = TRUE, silent = TRUE) {
  # Takes sorted (in ascending order) vectors
  # if (check) stopifnot(is.vector(table), is.numeric(table))
  if(check) stopifnot(!is.null(end.idx))
  m <- as.integer(ceiling((end.idx + start.idx) / 2)) # Midpoint
  if(!is.null(lin_scale)) nbins <- lin_scale[m] else nbins <- m
  if(!silent) cat('Looking for right divisions using', nbins, 'nbins...')
  bas <- CampaRi::basins_recognition(st, nx = nbins, plot = F, match = force_matching, out.file = F, new.dev = F, silent = T)
  ncl_found <- nrow(bas$tab.st)
  if(!silent) cat(' found', ncl_found,'clusters. \n')
  if (ncl_found > ncl_teo + tol) {
    if (start.idx == end.idx) return(list('found' = FALSE, 'ncl' = ncl_found, 'start.idx' = start.idx, 'nbins' = lin_scale[start.idx]))
    Recall(st, ncl_found, ncl_teo, start.idx = start.idx, end.idx = m - 1L, lin_scale = lin_scale,
           force_matching = force_matching, tol = tol, check = FALSE, silent = silent)
  } else if (ncl_found < ncl_teo - tol) {
    if (start.idx == end.idx) return(list('found' = FALSE, 'ncl' = ncl_found, 'start.idx' = start.idx, 'nbins' = lin_scale[start.idx]))
    Recall(st, ncl_found, ncl_teo, start.idx = m + 1L, end.idx = end.idx, lin_scale = lin_scale,
           force_matching = force_matching, tol = tol, check = FALSE, silent = silent)
  } else return(list('found' = TRUE, 'ncl' = ncl_found, 'idx' = m, 'nbins' = lin_scale[m]))
}
