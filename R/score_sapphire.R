#' @title Scoring the SAPPHIRE-plot basin division
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

score_sapphire <- function(the_sap, ann, scoring_method = 'nmi', silent = FALSE,
                           plot_basin_identification = FALSE, nbinsxy = NULL, merge_clusters = FALSE,
                           basin_optimization = NULL, number_of_clusters = NULL, force_matching = FALSE){
  
  # general input checking
  if(!is.character(the_sap) && !is.data.frame(the_sap)) stop("the_sap must be a string or a data frame")
  if(is.character(the_sap) && (!all(grepl("PROGIDX", the_sap)) && !all(grepl("REPIX", the_sap)))) stop("Please provide a data name starting with 'PROGIDX' or 'REPIX'" )
  if(!is.null(basin_optimization) && !is.character(basin_optimization)) stop('basin_optimization must be a character')
  if(!is.numeric(ann) && (!is.null(dim(ann)))) stop('Please provide an integer vector for ann')
  if(!is.null(nbinsxy) && !.isSingleInteger(nbinsxy)) stop('nbinsxy must be a single integer')
  if(!is.null(number_of_clusters) && !.isSingleInteger(number_of_clusters)) stop('number_of_clusters must be a single integer')
  if(!is.logical(silent)) stop('silent must be a logical')
  if(!is.logical(plot_basin_identification)) stop('plot_basin_identification must be a logical')
  if(!is.logical(merge_clusters)) stop('merge_clusters must be a logical')
  
  # methods input check
  scoring_method.opt <- c("adjusted_rand_index", "jaccard_index", "purity", "nmi")
  basin_optimization.opt <- c("uniformity", "number_of_clusters", "minimal_entropy")
  if(!(scoring_method[1] %in% scoring_method.opt)) stop("Scoring method option not valid")
  if(!is.null(basin_optimization) && !(basin_optimization[1] %in% basin_optimization.opt)) stop("basin_optimization method option not valid")
  
  # sapphire table loading
  if(!is.data.frame(the_sap))
    st <- as.data.frame(data.table::fread(the_sap, data.table = F))
  else
    st <- the_sap
  
  # hist(st[,5], breaks = 1000) # hist of the distances
  if(is.null(nbinsxy)) nbins_x <- nbins_y <- nbinsxy <- round(sqrt(nrow(st)*10))
  else nbins_x <- nbins_y <- nbinsxy

  if(!silent) cat('Number of (automatically) selected bins for the basin recognition step is', nbins_x, '\n')
  if(!is.null(basin_optimization)){
    if(basin_optimization[1] == 'uniformity') stop('uniformity basin_optimization option not yet ready.')
    if(basin_optimization[1] == 'minimal_entropy') stop('minimal_entropy basin_optimization option not yet ready.')
    
    # ------------------------------------------------------- basin_opt - number of clusters
    # This method consist only into finding the right number of clusters 
    # without splitting the possible subdivision (of the same number of clustering)
    # and score them per their relevance (e.g. uniformity or minimal entropy).
    #
    
    # if(basin_optimization[1] == 'number_of_clusters') { 
      
      # CampaRi::basin_optimization(...)
      #####################################
    #   
    #   # checks
    #   if(!silent) cat('Automatic optimization of basins based on number of cluster selected... \n')
    #   if(is.null(number_of_clusters)) stop('To use basin_optimization with number_of_clusters the correspective variable must be set to a single integer.')
    #   if(number_of_clusters < 2 && number_of_clusters >= nrow(st)) stop('Number of cluster too high or too low.')
    #   if(number_of_clusters > nbinsxy) stop('Please set the number of clusters lower than the number of initial nbins (nbinsxy).')
    #   
    #   # init
    #   how_fine <- 10
    #   lin_scale <- unique(round(seq(2, nbinsxy, length.out = how_fine)))
    #   if(!silent) cat('Selected a convergence step of', how_fine, 'subdivisions, ranging from', nbinsxy, 'to 2. \n')
    #   n_cl <- nbins_x
    #   whch <- how_fine
    #   
    #   # loop over the first coarse search of best partitioning
    #   while(n_cl != number_of_clusters){
    #     
    #     nbins_x <- lin_scale[whch] # linear scale from 2 to nbinsxy
    #     if(!silent) cat('Looking for right divisions using', nbins_x, ' nbins...\n')
    #     bas <- CampaRi::basins_recognition(st, nx = nbins_x, plot = F, match = force_matching, out.file = F, new.dev = F, silent = T)
    #     n_cl <- nrow(bas$tab.st) # take the number of clusters found
    #     
    #     # normal finer search - descent
    #     if(n_cl > number_of_clusters){
    #       old_whch <- whch
    #       whch <- round(whch/2)
    #       if(old_whch == whch) whch <- old_whch - 1
    #       if(whch < 1){
    #         n_fin_cl <- n_cl
    #         n_cl <- number_of_clusters
    #       } 
    #     # if fine partitioning went wrong, we want to try a finer search
    #     }else if(n_cl < number_of_clusters){
    #       how_fine2 <- 10
    #       if(whch+1 > how_fine) stop('Initial guess of bins brought you directly to have less barriers found than needed. Please consider an increment of nbinsxy parameter.')
    #       # case in which it is simply to select the upper part!
    #       
    #       if(!silent) cat('Found that our split was too coarse. Trying to split it again in', how_fine2,
    #                       'parts in the found ', nbins_x, '-', lin_scale[whch+1] , 'range.\n')
    #       lin_scale <- unique(round(seq(nbins_x, lin_scale[whch+1], length.out = how_fine2)))
    #       if(.lt(lin_scale) == 1){
    #         if(!silent) cat('Perfect number of divisions not found. In particular, we finally found', n_cl, 'partitions. Probably a finer search could work.', 
    #                         'Otherwise, it is necessary to suppose another number of clusters. In particular, we stopped before starting internal loop because',
    #                         'the span between the last two splits was not sufficiently large.\n')
    #         n_fin_cl <- n_cl
    #         n_cl <- number_of_clusters
    #       }
    #       whch <- how_fine2
    #       while(n_cl != number_of_clusters){
    #         nbins_x <- lin_scale[whch]
    #         bas <- CampaRi::basins_recognition(st, nx = nbins_x, plot = F, match = force_matching, out.file = F, new.dev = F, silent = T)
    #         n_cl <- max(bas$tab.st[,1]) 
    #         
    #         if(n_cl > number_of_clusters){
    #           whch <- whch - 1
    #         } else if(n_cl < number_of_clusters){
    #           if(!silent) cat('Perfect number of divisions not found. In particular, we finally found', n_cl, 'partitions. Probably a finer search could work.', 
    #                           'Otherwise, it is necessary to suppose another number of clusters. \n')
    #           n_fin_cl <- n_cl
    #           n_cl <- number_of_clusters
    #         }else{
    #           if(!silent) cat('We found a perfect binning using', nbins_x, 'number of bins.\n')
    #           n_fin_cl <- n_cl
    #         }
    #       }
    #     }else{ # end of the finer part, i.e. n_cl == number of clusters without having to go inside the inner loop
    #       if(!silent) cat('We found a perfect binning using', nbins_x, 'number of bins.\n')
    #       n_fin_cl <- n_cl
    #     }
    #   }
    # }
    # # final call if you want to plot!
    # bas <- CampaRi::basins_recognition(st, nx = nbins_x, dyn.check = 1, 
    #                                    plot = plot_basin_identification, match = force_matching, out.file = F, new.dev = F, silent = silent)
    # 
    # 
    # 
    cat('bog')
      
    ##########################################
  }else{
    bas <- CampaRi::basins_recognition(st, nx = nbins_x, dyn.check = 1, 
                                       plot = plot_basin_identification, match = force_matching, out.file = F, new.dev = F, silent = silent)
    n_fin_cl <- nrow(bas$tab.st) # here we suppose that our split is the one we wanted and we define the number of resulting clusters
  }
  
  
  # ann is reordered following progrex index and it is checked for the absence of 0s
  pin <- ann[c(st[,3])]
  uni_ann <- unique(pin)
  if(anyNA(pin)) stop('We can not handle NAs in the annotation. Check it!')
  # min of ann correction
  if(min(uni_ann) < 1){
    if(!silent) cat('Found negative or 0 values in the annotation. Correcting by sum of the abs of the min +1 per each annotation value. \n')
    pin <- pin + abs(min(uni_ann)) + 1
  }else if(min(uni_ann) > 1){
    if(!silent) cat('Found minimum values > 1 in the annotation. Correcting it automatically. \n')
    pin <- pin - min(uni_ann) + 1
  }
  if(any(seq(1:max(uni_ann)) != uni_ann)) stop('Please provide an annotation withouth gaps.')
  
  # init 
  n_labels <- .lt(uni_ann) # At this point I know how many annotation points we found.
  if(n_labels < n_fin_cl){
    if(!silent) cat('ATTENTION: found', n_fin_cl, 'clusters using basin recognition while the inserted annotation has only', n_labels, 'number of labels.',
                    'It would be the case to reduce the number of barriers. To do so, please consider reducing the number of bins (nbinsxy).\n')
    
    # This is the problematic case. 
  }
  if(n_labels > n_fin_cl){
    # no problem we can continue...
    if(!silent) cat('ATTENTION: found', n_fin_cl, 'clusters while in ann we have', n_labels, 'number of labels.',
                    'This is not problematic but it will lead to a possible overestimation of the error. Please consider more bins (nbinsxy).\n')
  }
  # plot(pin, pch='.')
  # sort(table(pin), decreasing=TRUE) # sorted
  
  # we choose the representative label
  major_freq <- list()
  for(jk in 1:n_fin_cl){
    maj_selected <- sort(table(pin[bas$tab.st[jk,2]:bas$tab.st[jk,3]]), decreasing=TRUE) # count and sort for each cluster found.
    maj_sel_filtered <- maj_selected[1:n_labels]  # select only first 3!
    major_freq[[jk]] <- (maj_sel_filtered*1.0)/bas$tab.st[jk,4] # calculating the density of major label
    major_freq[[jk]] <- rbind('d' = major_freq[[jk]], 'n' = maj_sel_filtered)
    major_freq[[jk]][is.na(major_freq[[jk]])] <- 0 
    # if a group does not contain some of the labels a <NA> will appear
    if(anyNA(colnames(major_freq[[1]]))){
      # no scream will be done. This part should be merged with the other one. Lets say that there is a better entropy or something else
      
    }
  }
  major_freq
  
  for(jk in 1:length(major_freq)){
    it <- c()
    for(kj in 1:ncol(major_freq[[jk]])) {
      ulm <- major_freq[[jk]][1, kj]
      if(!is.na(names(ulm))) kholn <- - ulm * log(ulm) else kholn <- 0
      it <- c(it,  kholn)
    }
    major_freq[[jk]] <- cbind(major_freq[[jk]], sh_en = sum(it))
  }
  
  lab <- list()
  size <- list()
  sh_en <- list()
  for(jk in 1:length(major_freq)){
    lab[[jk]] <- as.integer(colnames(major_freq[[jk]])[1])
    size[[jk]] <- as.integer(major_freq[[jk]][2, 1])
    sh_en[[jk]] <- major_freq[[jk]][1,ncol(major_freq[[jk]])]
    
  }  
  
  fs <- data.frame(cbind(label = unlist(lab), size = unlist(size), sh_en = unlist(sh_en)))
  # fs
  # sum(fs[,2])
  # sum(unlist(size))
  # fs[order(fs$sh_en),]
  
  # principle: uniformity
  res_bound <- list()
  res_label <- list()
  diff_vec <- diff(unlist(lab))
  lt_sec <- bas$tab.st[,4]
  ul <- 1 
  a <- lt_sec[1]
  res_label[1] <- lab[1]
  for(ini in 1:length(diff_vec)){
    if(diff_vec[ini] == 0) {
      a <- a + lt_sec[ini + 1]
    }
    if(diff_vec[ini] != 0){
      res_label[ul] <- lab[ini]
      res_bound[ul] <- a
      a <- lt_sec[ini + 1]
      ul <- ul + 1
    }
    if(ini == length(diff_vec)){
      res_label[ul] <- lab[ini]
      res_bound[ul] <- a
    }
  }
  res_bound <- unlist(res_bound)
  # res_bound
  # sum(res_bound)
  # lt_sec
  res_label <- unlist(res_label)
  # res_label
  
  # now creating the output vector
  predicted_div <- list()
  for(i in 1:length(res_label)){
    predicted_div[[i]] <- rep(res_label[i], res_bound[i])
  }
  predicted_div <- unlist(predicted_div)
  plot(predicted_div)
  
  # calculating the entropy of the new selection
  major_freq <- list()
  res_b <- 0
  for(jk in 1:length(res_label)){
    major_freq[[jk]] <- sort(table(ann[c(st[,3])][(res_b+1):(res_b + res_bound[i])]), decreasing=TRUE)[1:4] * 1.0 / res_bound[i] 
    # (res_bound[i]-res_b)
    major_freq[[jk]] <- rbind(major_freq[[jk]], sort(table(ann[c(st[,3])][(res_b+1):(res_b + res_bound[i])]), decreasing=TRUE)[1:4])
    major_freq[[jk]][is.na(major_freq[[jk]])] <- 0 
    res_b <- res_b + res_bound[i]
  }
  
  # major_freq
  for(jk in 1:length(major_freq)){
    it <- c()
    for(kj in 1:ncol(major_freq[[jk]])) {
      ulm <- major_freq[[jk]][1, kj]
      if(!is.na(names(ulm))) kholn <- - ulm * log(ulm) else kholn <- 0
      it <- c(it,  kholn)
    }
    major_freq[[jk]] <- cbind(major_freq[[jk]], sh_en = sum(it))
  }
  lab <- list()
  size <- list()
  sh_en <- list()
  for(jk in 1:length(major_freq)){
    lab[[jk]] <- as.integer(colnames(major_freq[[jk]])[1])
    size[[jk]] <- as.integer(major_freq[[jk]][2, 1])
    sh_en[[jk]] <- major_freq[[jk]][1,ncol(major_freq[[jk]])]
  }  
  
  fs <- data.frame(cbind(label = unlist(lab), size = unlist(size), sh_en = unlist(sh_en)))
  # fs
  # fs[order(fs$sh_en), ]
  
  # scoring it with accuracy?
  # external_validation(true_labels = ann[st[,3]], clusters = predicted_div, method = "adjusted_rand_index", summary_stats = FALSE)
  # external_validation(true_labels = ann[st[,3]], clusters = predicted_div, method = "jaccard_index", summary_stats = FALSE)
  # external_validation(true_labels = ann[st[,3]], clusters = predicted_div, method = "purity", summary_stats = FALSE)
  return(external_validation(true_labels = ann[st[,3]], clusters = predicted_div, method = "nmi", summary_stats = FALSE))
}

