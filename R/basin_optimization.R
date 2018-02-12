#' @title Basin optimization of the SBR
#' @description
#'     \code{score_sapphire} uses the information provided by the SAPPHIRE basin recognition in order to score the identified free energy basins
#'     against the original partition (annotation).
#'     
#' @param the_sap Name of the PROGIDX_<...>.dat or REPIX_<...>.dat file or Data frame with three columns containing, in order, the (sorted) progress index,
#'  the relative time indices and the cut function. 
#' @param ann Annotation must be a single line with the true cluster labels.
#' @param scoring_method Precise scoring method (final comparison between cluster labels)
#'       \itemize{
#'            \item "\code{mni}" 
#'            \item "\code{adjusted_rand_index}" 
#'            \item "\code{jaccard_index}" 
#'            \item "\code{purity}" 
#'       }
#' @param silent A logical value indicating whether the function has to remain silent or not. Default value is \code{FALSE}.
#' @param plot_basin_identification A logical value indicating whether to display the SAPPHIRE plot with the resulting partitions or not. 
#' Black partitions are the matched ones, blue ones derive only from the dynamic analysis and orange ones only from the kinetic analysis. 
#' The green curve is the kinetic annotation (black curve) where the parabolic shape has been subtracted, i.e. the actual curve used for 
#' the peaks identification. Default value is \code{FALSE}.
#' @param nbinsxy Number of bins on the y and x axis of the 2-D histogram. Default to sqrt(nrow(the_sap)).
#' @param merge_cluster Logical that allow clusters to be merged automatically if consecutives
#' @param basin_optimization Method name for the selection of precise basin. This is a costly procedure which is based on an heuristical obj function w
#' with the following subgroups.
#' \itemize{
#'    \item "\code{uniformity}" not available
#'    \item "\code{number_of_clusters}" optimize on the basis of the number of clusters
#'    \item "\code{minimal_entropy}" not available yet
#' }
#' @param number_of_clusters if basin_optimization is active accordingly this must be set to integer.
#' @param force_matching Please refer to \code{\link{basin_recognition}} for further details about the match option.
#'      
#' @return The resulting score (0-1)
#' @examples
#' 
#' adjl <- mst_from_trj(trj = matrix(rnorm(1000), nrow = 100, ncol = 10))
#' ret<-gen_progindex(adjl = adjl)
#' gen_annotation(ret_data = ret, local_cut_width = 10)
#' \dontrun{
#' score_sapphire(the_sap = "PROGIDX_000000000001.dat", ann = rnorm(100))
#' }
#' 
#' @details For details regarding the SAPPHIRE plot, please refer to the relative publications \url{http://www.nature.com/articles/srep06264}. 
#' Main documentation of the original campari software \url{http://campari.sourceforge.net/documentation.html}.
#' 
#' @importFrom ClusterR external_validation
#' @importFrom data.table fread
#' @export score_sapphire

basin_optimization <- function(the_sap, basin_optimization_method = NULL, how_fine_search = 100,
                               plot_basin_identification = FALSE, nbins_x = NULL, nbins_y = nbins_x,
                               number_of_clusters = NULL, force_matching = FALSE, silent = FALSE){
  # general input checking
  if(!is.character(the_sap) && !is.data.frame(the_sap)) stop("the_sap must be a string or a data frame")
  if(is.character(the_sap) && (!all(grepl("PROGIDX", the_sap)) && !all(grepl("REPIX", the_sap)))) stop("Please provide a data name starting with 'PROGIDX' or 'REPIX'" )
  if(!is.null(basin_optimization) && !is.character(basin_optimization)) stop('basin_optimization must be a character')
  if(!is.null(nbins_x) && !.isSingleInteger(nbins_x)) stop('nbins_x must be a single integer')
  if(!is.null(nbins_y) && !.isSingleInteger(nbins_x)) stop('nbins_y must be a single integer')
  if(!is.null(number_of_clusters) && !.isSingleInteger(number_of_clusters)) stop('number_of_clusters must be a single integer')
  if(!is.logical(silent)) stop('silent must be a logical')
  if(!is.logical(plot_basin_identification)) stop('plot_basin_identification must be a logical')
  
  # methods input check
  basin_optimization.opt <- c("uniformity", "number_of_clusters", "minimal_entropy")
  if(!is.null(basin_optimization) && !(basin_optimization[1] %in% basin_optimization.opt)) stop("basin_optimization method option not valid")
  
  # sapphire table loading
  if(!is.data.frame(the_sap))
    st <- as.data.frame(data.table::fread(the_sap, data.table = F))
  else
    st <- the_sap
  
  if(is.null(nbins_x)) nbins_x <- round(sqrt(nrow(st)*10))
  if(is.null(nbins_y)) nbins_y <- round(sqrt(nrow(st)*10))
  if(nbins_x != nbins_y) {
    if(!silent) cat('For simplicity we do not allow yet to use different number of nx and ny (bins for x and y).')
    nbins_x <- nbins_y <- max(nbins_x, nbins_y)
  }
  
  if(!silent) cat('Number of (automatically) selected bins for the basin recognition step is', nbins_x, '\n')
  if(!is.null(basin_optimization)){
    if(basin_optimization[1] == 'uniformity') stop('uniformity basin_optimization option not yet ready.')
    if(basin_optimization[1] == 'minimal_entropy') stop('minimal_entropy basin_optimization option not yet ready.')
    
    # ------------------------------------------------------- basin_opt - number of clusters
    # This method consist only into finding the right number of clusters 
    # without splitting the possible subdivision (of the same number of clustering)
    # and score them per their relevance (e.g. uniformity or minimal entropy).
    #
    if(basin_optimization_method[1] == 'number_of_cluster'){
      
      # checks
      if(!silent) cat('Automatic optimization of basins based on number of cluster selected... \n')
      if(is.null(number_of_clusters)) stop('To use basin_optimization with number_of_clusters the correspective variable must be set to a single integer.')
      if(number_of_clusters < 2 && number_of_clusters >= nrow(st)) stop('Number of cluster too high or too low.')
      if(number_of_clusters > nbinsxy) stop('Please set the number of clusters lower than the number of initial nbins (nbinsxy).')
      
      # init
      how_fine <- 100
      lin_scale <- unique(round(seq(7, nbins_x, length.out = how_fine)))
      if(!silent) cat('Selected a convergence step of', how_fine, 'subdivisions, ranging from', nbinsxy, 'to 2. \n')
      .BiSBR(st = st, ncl_found = 1, ncl_teo = 3, start.idx = 1, end.idx = length(abba), lin_scale =  abba, force_matching = T, silent = F)
      
    }
  }
}
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
    
    require(microbenchmark); require(ggplot2); v <- 1:10e7
    a <- microbenchmark("bs" = .BiSearch(v, 432), "classic" = (432 %in% v), times = 10)
    autoplot(a)
    
    nbinsxy <- sqrt(nrow(st)*10)
    abba <- unique(round(seq(2, nbinsxy, length.out = 10)))
    .BiSBR(st = st, ncl_found = 1, ncl_teo = 3, start.idx = 1, end.idx = length(abba), lin_scale =  abba, force_matching = F, silent = F)
    .BiSBR(st = st, ncl_found = 1, ncl_teo = 3, start.idx = 1, end.idx = length(abba), lin_scale =  abba, force_matching = T, silent = F)
    # testing for minimum number of bins without crashing
    CampaRi::basins_recognition(st[1:1000,], nx = 4, new.dev = F, out.file = F, plot = T)
    CampaRi::basins_recognition(st, nx = 12, new.dev = F, out.file = F, plot = T)
    CampaRi::basins_recognition(st, nx = 88, new.dev = F, out.file = F, plot = T, match = T)
    CampaRi::basins_recognition(st, nx = 97, new.dev = F, out.file = F, plot = T, match = T)
    CampaRi::basins_recognition(st, nx = 343, new.dev = F, out.file = F, plot = T, match = T)
    alu <- CampaRi::basins_recognition(st, nx = 905, new.dev = F, out.file = F, plot = T, match = T)