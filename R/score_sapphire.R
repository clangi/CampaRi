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
#' @param basin_obj Output of \code{\link{basins_recognition}}. If you used \code{\link{basin_optimization}} to optimize the clusters please
#'  insert only the resulting bas object
#' @param manual_barriers If an integer vector is inserted, it is used as the barrier locations.
#' @param plot_pred_true_resume Defaults tp \code{FALSE}. If set to true it plots the predicted and true labels one along the other.
#' @param merge_clusters Logical that allow clusters to be merged automatically if consecutives
#' @param silent A logical value indicating whether the function has to remain silent or not. Default value is \code{FALSE}.
#' @param merge_clusters A logical value for indicating if merging consecutive clusters with same labels.
#' @param ... Not yet in use.
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

score_sapphire <- function(the_sap, ann, manual_barriers = NULL, basin_obj = NULL,                     # fundamental inputs
                           scoring_method = 'nmi', merge_clusters = FALSE,                                      # scoring details
                           plot_pred_true_resume = FALSE, silent = FALSE,                                       # it refers to this function
                           ...){                                                                                # to add
  
  # ----------------------------------------------------------------------------------------------- general input checking
  have_the_bas <- FALSE
  manual_mode <- FALSE
  
  # Fundamental inputs: the_sap and optimizations 
  # - basin_obj = bas: uses bas object from basins_recognition() or basin_optimization()
  # - manual_barriers = c(1,2,3,4): insert manually the barriers [it creates a fake bas obj]
  if(!is.character(the_sap) && !is.data.frame(the_sap)) stop("the_sap must be a string or a data frame")
  if(is.character(the_sap) && (!grepl("PROGIDX", the_sap) && !grepl("REPIX", the_sap))) stop("Please provide a the_sap name starting with 'PROGIDX' or 'REPIX'" )
  
  if(!is.null(basin_obj) && !is.logical(basin_obj)){
    if(is.null(names(basin_obj)) || names(basin_obj)[1] != 'tab.st') 
      stop('Use the basin_optimization (optimal_bas) or the basins_recognition output for basin_optimization var.
           It must be of the shape tab.st nbins seq.st etc.. (or put it TRUE/FALSE).')
    if(!is.null(manual_barriers)) stop('Use either manual mode or bas mode')
    have_the_bas <- TRUE
  }
  if(!is.null(manual_barriers)){
    if(!all(sapply(manual_barriers, .isSingleInteger))) stop('All the values in manual_barriers must be integers')
    if(any(manual_barriers < 1)) stop('the manual barriers should be above 0')
    manual_mode <- TRUE
  }
  # - basin_optimization = TRUE: runs the optimization of the basins
  # if(!is.null(basin_optimization) && is.logical(basin_optimization)) {
  #   if(!silent) cat('You are using basin_optimization as a logical. We see it. \n')
  #   do_optimization <- TRUE
  #   if(!is.null(manual_barriers) && basin_optimization) stop('use manual mode or basin optimization mode')
  # }
  # 
  
  # Other input checks
  if(!is.numeric(ann) && (!is.null(dim(ann)))) stop('Please provide an integer vector for ann')
  if(!is.logical(merge_clusters)) stop('merge_clusters must be a logical')
  if(!is.logical(plot_pred_true_resume)) stop('plot_pred_true_resume must be a logical')
  if(!is.logical(silent)) stop('silent must be a logical')
  
  # methods input check
  scoring_method.opt <- c("adjusted_rand_index", "jaccard_index", "purity", "nmi")
  if(!(scoring_method[1] %in% scoring_method.opt)) stop("Scoring method option not valid")
  
  # sapphire table loading
  if(!is.null(the_sap) && !is.data.frame(the_sap)){
    st <- as.data.frame(data.table::fread(the_sap, data.table = F))
  }else{
    st <- as.data.frame(the_sap)
  }

  # - basin_optimization = TRUE: runs the optimization of the basins
  # if(do_optimization){
  #   if(!silent) cat('It is advisable to use the basin_optimization function EXTERNALLY to this one and feed the output to basin_optimization.\n')
  #   optim_bas <- CampaRi::basin_optimization(the_sap = the_sap, silent = silent, ...)
  #   bas <- optim_bas$bas
  #   n_fin_cl <- nrow(bas$tab.st)
    
  # Main switcher for the final output  
  # - basin_obj = bas: uses bas object from basins_recognition() or basin_optimization()
  # - manual_barriers = c(1,2,3,4): insert manually the barriers [it creates a fake bas obj]
  if(have_the_bas){
    if(!silent) cat('You inserted the bas object directly. No specific check of your wrong doing is applied. Therefore USE THE RIGHT ONE.\n')
    bas <- basin_obj
    n_fin_cl <- nrow(bas$tab.st)
    pifilename <- bas$filename
  }else if(manual_mode){
    if(!silent) cat('You inserted manual barriers for testing.\n')
    n_fin_cl <- .lt(manual_barriers) + 1
    if(is.null(the_sap)) stop('with manual_barriers it is needed to specify the SAPPHIRE table in the_sap')
  }else{
    stop('We did not understood what to do.')
  }
  
  # ---------------------------------------------------------------------------------------------- Annotation analysis
  # Using the prog. idx for reordering and lpi
  if(!silent) cat('Having inserted the_sap we reorder the ann using it.\n')
  if(.lt(ann) != nrow(st)) stop('Annotation and progress index must have same length. ')
  piann <- ann[c(st[,3])]
  lpiann <- .lt(piann)
  uni_ann <- unique(piann)
  if(anyNA(piann)) stop('We can not handle NAs in the annotation. Check it!')
  
  # min of ann correction
  if(min(uni_ann) < 1){
    if(!silent) cat('Found negative or 0 values in the annotation. Correcting by sum of the abs of the min +1 per each annotation value. \n')
    piann <- piann + abs(min(uni_ann)) + 1
    uni_ann <- unique(piann)
  }else if(min(uni_ann) > 1){
    if(!silent) cat('Found minimum values > 1 in the annotation. Correcting it automatically. \n')
    piann <- piann - min(uni_ann) + 1
    uni_ann <- unique(piann)
  }
  
  # check for gaps in the annotation
  if(any(seq(1:max(uni_ann)) != sort(uni_ann))) stop('Please provide an annotation withouth gaps.')
  
  # Checks for the number of labels VS number of clusters
  n_labels <- .lt(uni_ann)
  
  if(n_labels < n_fin_cl){
    if(!silent) cat('ATTENTION: found', n_fin_cl, 'clusters using basin recognition while the inserted annotation has only', n_labels, 'number of labels.',
                    'It would be the case to reduce the number of barriers. To do so, please consider reducing the number of bins (nbinsxy).\n')
  } else {
    if(!silent) cat('ATTENTION: found', n_fin_cl, 'clusters while in ann we have', n_labels, 'number of labels.',
                    'This is not problematic but it will lead to a possible overestimation of the error. Please consider more bins (nbinsxy).\n')
  }

  # creating a fake bas to fit the following analysis
  if(manual_mode){
    n_cl.b <- 1:n_fin_cl
    if(any(manual_barriers >= lpiann)) stop('one manual_barriers or more are higher than the number of snapshots!.')
    start.b <- c(1, sort(manual_barriers))
    end.b <- c(sort(manual_barriers), lpiann)
    lt.b <- diff(c(1, sort(manual_barriers), lpiann))
    lt.b[1] <- lt.b[1] + 1
    if(sum(lt.b) != lpiann) stop('The inserted barriers dont sum up to the length of the annotation.')
    bastbl <- cbind('n_cl' = n_cl.b, 'start' = start.b, 'end' = end.b, '.lt' = lt.b)
    bas <- list('tab.st' = bastbl)
  }
  
  
  # ---------------------------------------------------------------------------------------------- Creation of the Entropy levels
  # we choose the representative label
  label_freq_list <- list() # each element of this list is a cluster
  for(jk in 1:n_fin_cl){
    maj_selected <- sort(table(piann[bas$tab.st[jk,2]:bas$tab.st[jk,3]]), decreasing=TRUE) # count and sort for each cluster found.
    maj_sel_filtered <- maj_selected[1:n_labels]  # select only first n_labels!
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
  
  # ---------------------------------------------------------------------------------------------- Creation of Predicted vector
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
                  round((lpiann - sum(unlist(size))) * 100 / lpiann, digits = 1), '% missassignment.\n\n')
  # now attaching it to the bas output
  max_freq_table <- cbind(max_freq_table, bas$tab.st)
  if(!silent) {print(max_freq_table); cat('\n')}
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
  if(nrow(max_freq_table) > 1){
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
  } else{
    if(merge_clusters && !silent) cat('Only one cluster found. Nothing to be merged really. \n')
    res_bound <- lt_cl[1]
    res_label <- max_freq_table$label[1]
  }
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
    plot_df <- data.frame('predicted' = as.factor(predicted_div), 'true' = ann[st[,3]]-0.05)
    gg <- ggplot(data = plot_df) + 
          geom_point(aes(y = 'predicted', x = 1:lpiann, colour = 'red'), size = 0.3) + 
          geom_point(aes(y = 'true', x = 1:lpiann, colour = 'blue'), size = 0.3) + 
          theme_minimal() + xlab('Progress Index') + ylab('Cluster') + 
          scale_color_manual(name = "Labels", labels = c("True", "Predicted"), values = c("red", "blue")) +
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
  
  # ---------------------------------------------------------------------------------------------- Final scoring
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


