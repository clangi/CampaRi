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
#' @param basin_optimization If \code{TRUE} it will use \code{\link{basin_optimization}} to optimize the clusters. Please consider adding its relevant variables.
#' @param plot_basin_identification A logical value indicating whether to display the SAPPHIRE plot with the resulting partitions or not. 
#' Black partitions are the matched ones, blue ones derive only from the dynamic analysis and orange ones only from the kinetic analysis. 
#' The green curve is the kinetic annotation (black curve) where the parabolic shape has been subtracted, i.e. the actual curve used for 
#' the peaks identification. Default value is \code{FALSE}.
#' @param plot_pred_true_resume Defaults tp \code{FALSE}. If set to true it plots the predicted and true labels one along the other.
#' @param nbins_x Number of bins on  x-axis of the 2-D histogram. Default to sqrt(nrow(the_sap)).
#' @param nbins_y Number of bins on the y-axis of the 2-D histogram. Default to sqrt(nrow(the_sap)).
#' @param merge_cluster Logical that allow clusters to be merged automatically if consecutives
#' @param force_matching Please refer to \code{\link{basin_recognition}} for further details about the match option.
#' @param silent A logical value indicating whether the function has to remain silent or not. Default value is \code{FALSE}.
#' @param ... This variables will be sent to \code{\link{basin_optimization}} without any checking
#' 
#' @return A list containing
#'       \itemize{
#'         \item "\code{score.out}" Resulting score between 0 and 1 (good).
#'         \item "\code{max_freq_table}" Matrix with a description of the labels selected and their cluster characteristics
#'         \item "\code{label_freq_list}" Total representation of the clusters
#'       }   
#'       
#'       
#' @examples
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
#' @import ggplot2
#' @export score_sapphire

score_sapphire <- function(the_sap, ann, scoring_method = 'nmi', merge_clusters = FALSE, 
                           basin_optimization = FALSE, plot_basin_identification = FALSE, plot_pred_true_resume = FALSE,
                           nbins_x = NULL, nbins_y = nbins_x, force_matching = FALSE, silent = FALSE, ...){
  
  # general input checking
  if(!is.character(the_sap) && !is.data.frame(the_sap)) stop("the_sap must be a string or a data frame")
  if(is.character(the_sap) && (!all(grepl("PROGIDX", the_sap)) && !all(grepl("REPIX", the_sap)))) stop("Please provide a data name starting with 'PROGIDX' or 'REPIX'" )
  if(!is.numeric(ann) && (!is.null(dim(ann)))) stop('Please provide an integer vector for ann')
  if(!is.null(nbins_x) && !.isSingleInteger(nbins_x)) stop('nbins_x must be a single integer')
  if(!is.null(nbins_y) && !.isSingleInteger(nbins_x)) stop('nbins_y must be a single integer')
  # if(!is.null(number_of_clusters) && !.isSingleInteger(number_of_clusters)) stop('number_of_clusters must be a single integer')
  if(!is.logical(silent)) stop('silent must be a logical')
  if(!is.logical(basin_optimization)) stop('basin_optimization must be a logical')
  if(!is.logical(plot_basin_identification)) stop('plot_basin_identification must be a logical')
  if(!is.logical(plot_pred_true_resume)) stop('plot_pred_true_resume must be a logical')
  if(!is.logical(merge_clusters)) stop('merge_clusters must be a logical')
  
  # methods input check
  scoring_method.opt <- c("adjusted_rand_index", "jaccard_index", "purity", "nmi")
  if(!(scoring_method[1] %in% scoring_method.opt)) stop("Scoring method option not valid")
  
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
  
  if(!silent) cat('Number of (automatically) selected bins for the basin recognition step is', nbins_x, '\n')
  if(basin_optimization){
    optim_bas <- CampaRi::basin_optimization(the_sap = the_sap, plot_basin_identification = plot_basin_identification, 
                                             nbins_x = nbins_x, nbins_y = nbins_x, force_matching = force_matching, silent = silent, ...)
    bas <- optim_bas$bas
    n_fin_cl <- nrow(bas$tab.st)
  }else{
    bas <- CampaRi::basins_recognition(st, nx = nbins_x, ny = nbins_x, dyn.check = 1, 
                                       plot = plot_basin_identification, match = force_matching, out.file = F, new.dev = F, silent = silent)
    n_fin_cl <- nrow(bas$tab.st) # here we suppose that our split is the one we wanted and we define the number of resulting clusters
  }
  
  
  # ann is reordered following progrex index and it is checked for the absence of 0s
  pin <- ann[c(st[,3])]
  lpin <- .lt(pin)
  uni_ann <- unique(pin)
  if(anyNA(pin)) stop('We can not handle NAs in the annotation. Check it!')
  
  # min of ann correction
  if(min(uni_ann) < 1){
    if(!silent) cat('Found negative or 0 values in the annotation. Correcting by sum of the abs of the min +1 per each annotation value. \n')
    pin <- pin + abs(min(uni_ann)) + 1
    uni_ann <- unique(pin)
  }else if(min(uni_ann) > 1){
    if(!silent) cat('Found minimum values > 1 in the annotation. Correcting it automatically. \n')
    pin <- pin - min(uni_ann) + 1
    uni_ann <- unique(pin)
  }
  if(any(seq(1:max(uni_ann)) != sort(uni_ann))) stop('Please provide an annotation withouth gaps.')
  
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
  label_freq_list <- list() # each element of this list is a cluster
  for(jk in 1:n_fin_cl){
    maj_selected <- sort(table(pin[bas$tab.st[jk,2]:bas$tab.st[jk,3]]), decreasing=TRUE) # count and sort for each cluster found.
    maj_sel_filtered <- maj_selected[1:n_labels]  # select only first 3!
    label_freq_list[[jk]] <- (maj_sel_filtered*1.0)/bas$tab.st[jk,4] # calculating the density of major label
    label_freq_list[[jk]] <- rbind('d' = label_freq_list[[jk]], 'n' = as.integer(maj_sel_filtered))
    label_freq_list[[jk]][is.na(label_freq_list[[jk]])] <- 0 
    # if a group does not contain some of the labels a <NA> will appear
    if(anyNA(colnames(label_freq_list[[jk]]))){
      # no scream will be done. This part should be merged with the other one. Lets say that there is a better entropy or something else
      clnms <- colnames(label_freq_list[[jk]])
      clnms[is.na(clnms)] <- setdiff(c(1:n_labels), c(colnames(label_freq_list[[jk]])))
      colnames(label_freq_list[[jk]]) <- clnms
    }
    label_freq_list[[jk]] <- rbind('lab' = as.integer(colnames(label_freq_list[[jk]])), 'clu' = rep(jk, n_labels), label_freq_list[[jk]])
    
  }
  # label_freq_list
  
  # calculating the entropy
  for(jk in 1:length(label_freq_list)){
    it <- c()
    for(kj in 1:ncol(label_freq_list[[jk]])) {
      ulm <- label_freq_list[[jk]]['d', kj]
      if(ulm != 0) kholn <- - ulm * log(ulm) else kholn <- 0
      it <- c(it,  kholn)
    }
    label_freq_list[[jk]] <- rbind(label_freq_list[[jk]], 'sh_en' = it)
  }
  # label_freq_list
  
  # collecting the various elected labels - selection policy: unique labels!
  # selection_policy <- 'unique' # this is like saying do not merge clusters
  if(!merge_clusters){ # TO DO BETTER - spawning policy
    # now I construct the grouping for afterwards selection
    # group_lab_b_cl <- rep(list(NULL), n_fin_cl)
    # for(jk in 1:n_fin_cl){
    #   for(ul in 1:n_fin_cl){
    #     whc <- (label_freq_list[[jk]][1,]*1.0 == ul*1.0)
    #     group_lab_b_cl[[ul]] <- cbind(group_lab_b_cl[[ul]], label_freq_list[[jk]][, whc])
    #   }
    # }
    # group_lab_b_cl <- lapply(X = group_lab_b_cl, FUN = function(x) rbind('cl' = 1:n_fin_cl, x, 'weight' = (1-x[2,])*x[4,]))
    # group_lab_b_cl
    lab <- list()
    slct <- array(1, dim = .lt(label_freq_list))
    for(jk in 1:.lt(label_freq_list)){
      lab[[jk]] <- as.integer(colnames(label_freq_list[[jk]])[slct[jk]])
      if(jk != 1){
        while(lab[[jk]] %in% c(lab[1:(jk-1)])){
          slct[jk] <- slct[jk] + 1
          # label_freq_list[r]; unlist(lab[r]); slct[r]
          if(slct[jk] > n_labels || label_freq_list[[jk]]['n', slct[jk]] == 0){
            if(!silent) cat('Unfortunately we found that all the labels in cluster', jk,
                        'are colliding with others (we did not count empty buckets). Kept the duplication.\n')
            slct[jk] <- 1
            print_new_lab <- FALSE
            break
          }else{
            print_new_lab <- TRUE
            lab[[jk]] <- as.integer(colnames(label_freq_list[[jk]])[slct[jk]])
          }
          if(!silent && print_new_lab) cat('Label', lab[[jk]], 'has been found already in another cluster. We will automatically',
                                           'select the second most present value in this cluster.\n')
        }
      }
    }
    # df_lab_freq <- data.frame(label_freq_list)
    # df_lab_freq[2,] <- 1 - df_lab_freq[2,]
    # sort policy: density of elements and entropy
    # df_lab_freq[,order(df_lab_freq[4,], decreasing = T)]
    # df_lab_freq
  } else {
    slct <- array(1, dim = .lt(label_freq_list))
  } 
  
  # main constructor loop for the selected (slct) labels
  lab <- list()
  size <- list()
  sh_en <- list()
  for(jk in 1:.lt(label_freq_list)){
    lab[[jk]] <- as.integer(colnames(label_freq_list[[jk]])[slct[jk]])
    size[[jk]] <- as.integer(label_freq_list[[jk]]['n', slct[jk]])
    sh_en[[jk]] <- label_freq_list[[jk]]['sh_en', slct[jk]]
  }  
  
  max_freq_table <- data.frame(cbind(label = unlist(lab), size = unlist(size), sh_en = unlist(sh_en)))
  if(!silent) cat('We found the following splits, accounting for', sum(unlist(size)),'of the elemnts. This means there are roughly', 
                  round((lpin - sum(unlist(size))) * 100 / lpin, digits = 1), '% missassignment.\n\n')
  # now attaching it to the bas output
  max_freq_table <- cbind(max_freq_table, bas$tab.st)
  if(!silent) print(max_freq_table); cat('\n')
  # max_freq_table[order(max_freq_table$sh_en),] # ordering based on internal entropy
  
  # merging policy - inputs: (label, size, sh_en) from major_freq_table
  if(merge_clusters && !silent) cat('Merging policy applied. For the moment only consecutive identically labeled clusters are merged. \n')
  
  # init
  res_bound <- list()
  res_label <- list()
  diff_vec <- diff(unlist(max_freq_table$label)) # barriers
  lt_sec <- bas$tab.st[,'.lt'] # length of the cl
  ul <- 1 # hlp loop
  res_bound[1] <- lt_cl <- lt_sec[1]
  res_label[1] <- max_freq_table$label[1]
  
  # loop on the barriers
  for(ini in 1:length(diff_vec)){
    
    # clusters have same lebels merging then
    if(diff_vec[ini] == 0){
      if(merge_clusters) {
        lt_cl <- lt_cl + lt_sec[ini + 1]
      } else{ # standard run diff = 0
        res_label[ul] <- max_freq_table$label[ini]
        res_bound[ul] <- lt_cl
        lt_cl <- lt_sec[ini + 1]
        ul <- ul + 1
      }
    }
    
    # standard run
    if(diff_vec[ini] != 0){
      res_label[ul] <- max_freq_table$label[ini]
      res_bound[ul] <- lt_cl
      lt_cl <- lt_sec[ini + 1]
      ul <- ul + 1
    }
    
    # final split
    if(ini == length(diff_vec)){
      res_label[ul] <- max_freq_table$label[ini + 1]
      res_bound[ul] <- lt_cl
    }
  }
  res_bound <- unlist(res_bound)
  res_label <- unlist(res_label)
  # res_bound
  # unlist(res_bound); max_freq_table$.lt; diff_vec; unlist(lab)
  # sum(res_bound)
  # lt_sec
  # res_label
  
  # creating the predicted vector
  predicted_div <- list()
  for(i in 1:length(res_label)){
    predicted_div[[i]] <- rep(res_label[i], res_bound[i])
  }
  predicted_div <- unlist(predicted_div)
  if(plot_pred_true_resume){
    plot_df <- data.frame(predicted = as.factor(predicted_div), true = ann[st[,3]]-0.05)
    gg <- ggplot(data = plot_df) + 
          geom_point(aes(y = predicted, x = 1:lpin, colour = 'red'), size = 0.3) + 
          geom_point(aes(y = true, x = 1:lpin, colour = 'blue'), size = 0.3) + 
          theme_minimal() + xlab('Progress Index') + ylab('Cluster') + 
          scale_color_manual(name = "Labels", labels = c("Predicted", "True"), values = c("red", "blue")) +
          guides(color = guide_legend(override.aes = list(size=5)))
    print(gg)
  }
  # calculating the entropy of the new selection
  # label_freq_list <- list()
  # res_b <- 0
  # for(jk in 1:length(res_label)){
  #   label_freq_list[[jk]] <- sort(table(ann[c(st[,3])][(res_b+1):(res_b + res_bound[i])]), decreasing=TRUE)[1:4] * 1.0 / res_bound[i] 
  #   # (res_bound[i]-res_b)
  #   label_freq_list[[jk]] <- rbind(label_freq_list[[jk]], sort(table(ann[c(st[,3])][(res_b+1):(res_b + res_bound[i])]), decreasing=TRUE)[1:4])
  #   label_freq_list[[jk]][is.na(label_freq_list[[jk]])] <- 0 
  #   res_b <- res_b + res_bound[i]
  # }
  # 
  # # label_freq_list
  # for(jk in 1:length(label_freq_list)){
  #   it <- c()
  #   for(kj in 1:ncol(label_freq_list[[jk]])) {
  #     ulm <- label_freq_list[[jk]][2, kj]
  #     if(!is.na(names(ulm))) kholn <- - ulm * log(ulm) else kholn <- 0
  #     it <- c(it,  kholn)
  #   }
  #   label_freq_list[[jk]] <- cbind(label_freq_list[[jk]], sh_en = sum(it))
  # }
  # lab <- list()
  # size <- list()
  # sh_en <- list()
  # for(jk in 1:length(label_freq_list)){ # check label_freq
  #   lab[[jk]] <- as.integer(colnames(label_freq_list[[jk]])[1])
  #   size[[jk]] <- as.integer(label_freq_list[[jk]][2, 1])
  #   sh_en[[jk]] <- label_freq_list[[jk]][1,ncol(label_freq_list[[jk]])]
  # }
  # 
  # max_freq_table <- data.frame(cbind(label = unlist(lab), size = unlist(size), sh_en = unlist(sh_en)))
  # max_freq_table
  # max_freq_table[order(max_freq_table$sh_en), ]
  
  # scoring it with accuracy?
  # external_validation(true_labels = ann[st[,3]], clusters = predicted_div, method = "adjusted_rand_index", summary_stats = FALSE)
  # external_validation(true_labels = ann[st[,3]], clusters = predicted_div, method = "jaccard_index", summary_stats = FALSE)
  # external_validation(true_labels = ann[st[,3]], clusters = predicted_div, method = "purity", summary_stats = FALSE)
  score.out <- ClusterR::external_validation(true_labels = ann[st[,3]], clusters = predicted_div, method = scoring_method, summary_stats = FALSE)
  if(!silent) cat('Using', scoring_method,'we obtained a final score of', score.out, '\n')
  invisible(list('score.out' = score.out, 'max_freq_table' = max_freq_table, 'label_freq_list' = label_freq_list))
}



#####################################
# ------------------------------------------------------- basin_opt - number of clusters
# This method consist only into finding the right number of clusters 
# without splitting the possible subdivision (of the same number of clustering)
# and score them per their relevance (e.g. uniformity or minimal entropy).
#

# if(basin_optimization[1] == 'number_of_clusters') { 
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
##########################################


