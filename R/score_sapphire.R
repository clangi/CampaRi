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
#'            \item "\code{nmi}" 
#'            \item "\code{adjusted_rand_index}" 
#'            \item "\code{jaccard_index}" 
#'            \item "\code{purity}" 
#'       }
#' @param basin_obj Output of \code{\link{basins_recognition}}. If you used \code{\link{basin_optimization}} to optimize the clusters please
#'  insert only the resulting bas object
#' @param manual_barriers If an integer vector is inserted, it is used as the barrier locations.
#' @param plot_pred_true_resume Defaults tp \code{FALSE}. If set to true it plots the predicted and true labels one along the other.
#' @param return_predicted If true return the predicted vector.
#' @param multi_cluster_policy This string decides how the clusters must be handled in the case of more clusters found in comparison to the annotation 
#' inserted. The default is \code{'popup'}. The available values are:
#'      \itemize{
#'         \item "\code{popup}" Creates new clusters (ordered by shannon weight).
#'         \item "\code{keep}" Keeps the duplications and so on. This is a bias feature.
#'         \item "\code{merge_previous}" Takes in account the last label inserted in the shannon weigthed table.
#'       }   
#'       
#' @param silent A logical value indicating whether the function has to remain silent or not. Default value is \code{FALSE}.
#' @param ... \code{max_number_of_elements} for plotting can be supplied. It defaults to 20k.
#' 
#' @return A list containing
#'       \itemize{
#'         \item "\code{score.out}" Resulting score between 0 and 1 (good).
#'         \item "\code{label_freq_list}" Total representation of the clusters along with the single shannon entropies.
#'         \item "\code{main_desc}" Description of the main clusters found.
#'         \item "\code{miss_perc}" Percentage of not correctly defined divisions.
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

score_sapphire <- function(the_sap, ann, manual_barriers = NULL, basin_obj = NULL,                              # fundamental inputs
                           scoring_method = 'adjusted_rand_index', multi_cluster_policy = 'popup',              # scoring  and merging details
                           plot_pred_true_resume = FALSE, silent = FALSE, return_plot = FALSE,                  # it refers to this function
                           return_predicted = FALSE,                                                              # return the prediction and true vector
                           ...){                                                                                # to add
  
  # ----------------------------------------------------------------------------------------------- general input checking
  input.args <- list(...)
  avail.extra.arg <- c('dbg_score_sapphire', 'max_number_of_elements')
  avail.policies <- c('popup', 'keep', 'merge_previous')
  
  if(!is.null(names(input.args)) && any(!(names(input.args) %in% avail.extra.arg))){
    if(!silent) warning('There is a probable mispelling in one of the inserted variables. Please check the available extra input arguments.')
    if(!silent) cat('!!!!!!!!! We found the following variables without a father (not between our extra input arguments) !!!!!!!!!!\n')
    if(!silent) cat(names(input.args)[!(names(input.args) %in% avail.extra.arg)], '\n')
  }
  if(!('dbg_score_sapphire' %in% names(input.args))) dbg_score_sapphire <- FALSE else dbg_score_sapphire <- input.args[['dbg_score_sapphire']]
  if(!('max_number_of_elements' %in% names(input.args))) max_number_of_elements <- NULL else max_number_of_elements <- input.args[['max_number_of_elements']]
  if(!is.logical(dbg_score_sapphire)) stop('dbg_score_sapphire must be a logical')
  if(is.null(max_number_of_elements)) max_number_of_elements <- 20000
  if(!.isSingleInteger(max_number_of_elements)) stop('max_number_of_elements must be a single integer')
  # methods input check
  scoring_method.opt <- c("adjusted_rand_index", "jaccard_index", "purity", "nmi")
  if(!(scoring_method[1] %in% scoring_method.opt)) stop("Scoring method option not valid")

  ###
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
  if(!is.logical(plot_pred_true_resume)) stop('plot_pred_true_resume must be a logical')
  if(!is.logical(return_plot)) stop('return_plot must be a logical')
  if(!is.logical(silent)) stop('silent must be a logical')
  if(!is.logical(return_predicted)) stop("return_predicted must be a logical value")
  
  # check over the merging or popping policy
  if(!.isSingleChar(multi_cluster_policy)) stop('multi_cluster_policy must be a logical')
  if(!(multi_cluster_policy %in% avail.policies)) stop(paste0('multi_cluster_policy must be between the following: ', paste(avail.policies, collapse = " ")))
  
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
  piann_true_v <- ann[c(st[,3])]
  lpiann <- .lt(piann_true_v)
  uni_ann <- unique(piann_true_v)
  if(anyNA(piann_true_v)) stop('We can not handle NAs in the annotation. Check it!')
  
  # min of ann correction
  if(min(uni_ann) < 1){
    if(!silent) cat('Found negative or 0 values in the annotation. Correcting by sum of the abs of the min +1 per each annotation value. \n')
    piann_true_v <- piann_true_v + abs(min(uni_ann)) + 1
    uni_ann <- unique(piann_true_v)
  }else if(min(uni_ann) > 1){
    if(!silent) cat('Found minimum values > 1 in the annotation. Correcting it automatically. \n')
    piann_true_v <- piann_true_v - min(uni_ann) + 1
    uni_ann <- unique(piann_true_v)
  }
  
  # check for gaps in the annotation
  if(any(seq(1:max(uni_ann)) != sort(uni_ann))) stop('Please provide an annotation withouth gaps.')
  
  # Checks for the number of labels VS number of clusters
  n_labels <- .lt(uni_ann)
  
  do_cl_spawn <- FALSE
  if(n_labels < n_fin_cl){
    if(!silent) cat('ATTENTION: found', n_fin_cl, 'clusters using basin recognition while the inserted annotation has only', n_labels, 'number of labels.',
                    '\nIt would be the case to reduce the number of barriers. To do so, please consider reducing the number of bins (nbinsxy).\n')
    if(!silent) cat('Spawning procedure (if merging is not TRUE) will be applied. This will create new cluster labels for the less populated or more hetereogenic clusters.\n')
    do_cl_spawn <- TRUE
  } else if(n_labels > n_fin_cl){
    if(!silent) cat('ATTENTION: found', n_fin_cl, 'clusters while in ann we have', n_labels, 'number of labels.',
                    '\nThis is not problematic but it will lead to a possible overestimation of the error. Please consider more bins (nbinsxy).\n')
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
    if(jk != 1) tarts <- bas$tab.st[jk,2] + 1
    else tarts <- bas$tab.st[jk,2]
    maj_selected <- sort(table(piann_true_v[tarts:bas$tab.st[jk,3]]), decreasing=TRUE) # count and sort for each cluster found.
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
  for(jk in 1:.lt(label_freq_list)){
    it <- c()
    for(kj in 1:ncol(label_freq_list[[jk]])) {
      ulm <- label_freq_list[[jk]]['d', kj]
      if(ulm != 0) kholn <- - ulm * log(ulm) else kholn <- 0
      it <- c(it,  kholn)
    }
    label_freq_list[[jk]] <- rbind(label_freq_list[[jk]], 'sh_en' = it)
  }
  # label_freq_list
  
  if(F) {
    ltest <- 500
    n_tcl <- 4
    tvec <- sample(1:n_tcl, size = ltest, replace = T)
    tvec2 <- rep(1,10)
    entr <- .myShEn(tvec); entr
    entr <- .myShEn(tvec2); entr
    
    .myShEn <- function(tvec){
      x2 <- sapply(1:max(tvec), FUN = function(x) sum(tvec == x))
      x2 <- x2 / .lt(tvec)
      return(-sum(x2*log(x2)))
    }
  }
  
  
  # ---------------------------------------------------------------------------------------------- Creation of the Entropy levels
  # collecting the various elected labels - selection policy: unique labels!
    
  # looking for collisions and candidates
  top_label <- sapply(label_freq_list, function(x) x['lab', 1])
  top_lab_n <- sapply(label_freq_list, function(x) x['n', 1])
  top_lab_sh <- sapply(label_freq_list, function(x) x['sh_en', 1])
  n_labels_in_cl <- sapply(label_freq_list, function(x) sum(x['n',] != 0))
  tot_dim <- sapply(label_freq_list, function(x) sum(x['n',]))
  tot_shen <- sapply(label_freq_list, function(x) sum(x['sh_en',]))
  n_cl.b <- 1:n_fin_cl
  main_desc <- data.frame('pos' = n_cl.b, 
                          'top_lab' = top_label, 'tl_sh' = top_lab_sh, 'tl_n' = top_lab_n, 
                          'freq_pos' = rep(1, n_fin_cl), 
                          'res_lab' = rep(NA, n_fin_cl), #'res_sh' = top_lab_sh, 'res_n' = top_lab_n,
                          'n_lab_in_cl' = n_labels_in_cl, 
                          'tl_shen_o_tot_n' = (top_lab_sh + 1) / tot_dim, 
                          'tl_shen_o_tl_n' = (top_lab_sh + 1) / top_lab_n, 
                          'tot_d' =  tot_dim, 'tot_shen' = tot_shen,
                          'shen_frac_d' = (tot_shen + 1) / tot_dim, 'd_frac_shen' = tot_dim / (tot_shen + 1))
  
  # ordering policy - weighted SH EN 
  # ordering_principle <- 'tl_shen_o_tot_n'
  # ordering_principle <- 'ori_ord'
  ordering_principle <- 'tl_shen_o_tl_n'
  # ordering_principle <- 'shen_frac_d'
  main_desc <- main_desc[order(main_desc[, ordering_principle]),]
  main_desc <- cbind(main_desc, 'ori_ord' = n_cl.b)

  # main_desc <- cbind(main_desc, 'problematic' = rep(FALSE, n_fin_cl))
  # main_desc[main_desc$res_label %in% c(6, 19, 22, 23, 25, 3), 'problematic'] <- TRUE
  
  # loop for choosing the label - using the number of clusters defined
  ncl.i <- 0
  extreme_search <- F
  
  while(ncl.i < n_fin_cl){
    
    if(ncl.i == 0) n_cl_popped <- n_labels
    ncl.i <- ncl.i + 1
    
    # select first candidate
    candida <- main_desc[ncl.i, 'top_lab']
    cl_pos <- main_desc$pos[ncl.i]
    h_diply <- main_desc$freq_pos[ncl.i]
    
    # repeat untill res_label is set
    while (TRUE) {
      h_diply <- h_diply + 1 # search level (it starts from two because first is in the init)
      
      # candidate not yet picked
      if(!(candida %in% main_desc$res_lab)){
        main_desc$res_lab[ncl.i] <- candida
        main_desc$freq_pos[ncl.i] <- h_diply - 1
        break
      
      # candidate already present
      } else {
        
        # redefining a candidate if the selection was not optimal of the first label
        if(h_diply <=  main_desc$n_lab_in_cl[ncl.i]){
          candida <- label_freq_list[[cl_pos]]['lab', h_diply]
          
          if(extreme_search){
            tls <- label_freq_list[[cl_pos]]['sh_en', h_diply]
            tln <- label_freq_list[[cl_pos]]['n', h_diply]
            main_desc[ncl.i, 'top_lab'] <- candida
            main_desc[ncl.i, 'tl_sh'] <- tls
            main_desc[ncl.i, 'tl_n'] <- tln
            main_desc[ncl.i, 'tl_shen_o_tot_n'] <- (tls + 1) / main_desc[ncl.i, 'tot_d']
            main_desc[ncl.i, 'tl_shen_o_tl_n'] <- (tls + 1) / tln
            main_desc$freq_pos[ncl.i] <- h_diply 
            main_desc <- main_desc[order(main_desc[, 'tl_shen_o_tot_n']),]
            main_desc$res_lab <- rep(NA, n_fin_cl)
            ncl.i <- 0
            break
          }
          
        # we finished the available labels and we have to define new clusters!
        }else{
          
          n_cl_popped <- n_cl_popped + 1 # it is valid also for the merge (just for the check)
          
          # define new clusters 
          if(multi_cluster_policy == 'popup'){
            main_desc$res_lab[ncl.i] <- n_cl_popped
            main_desc$freq_pos[ncl.i] <- n_cl_popped
            # cat(n_cl_popped, '\n')
            
          # keep the most freq - bias
          }else if(multi_cluster_policy == 'keep'){
            main_desc$res_lab[ncl.i] <- main_desc[ncl.i, 'top_lab'] # candida is different here
            # main_desc$freq_pos[ncl.i] <- 1
            
          # merge previous on the ordered main_desc
          }else if(multi_cluster_policy == 'merge_previous'){
            main_desc$res_lab[ncl.i] <- main_desc[ncl.i-1, 'top_lab'] 
            main_desc$freq_pos[ncl.i] <- main_desc$freq_pos[ncl.i-1] 
          }
          break # we must exit somehow
        }
      }
    }
  }
  # stopifnot(n_cl_popped == n_fin_cl) # crashing if n_labels > n_fin_cl
  main_desc <- main_desc[order(main_desc$pos),]
  
  res_bound <- bas$tab.st[,4]
  res_label <- main_desc$res_lab
  
  # DEPRECATED -------------------------------------------------------------------------------------------------------
  # ftab_for_coll <- data.frame('labs' = 1:n_labels, 'freq' = sapply(1:n_labels, function(x) sum(x == possible_coll)))
  # 
  # not_to_consider <- ftab_for_coll$labs[which(ftab_for_coll$freq == 1)]
  # ftab_for_coll <- ftab_for_coll[ftab_for_coll$labs[-not_to_consider],]
  
  # if(F){
  #   if(multi_cluster_policy == 'select_non_conflict' || multi_cluster_policy == 'select_non_conflict_spawn'){
  #     
  #     # init - broken
  #     lab <- list()
  #     slct <- array(1, dim = .lt(label_freq_list))
  #     
  #     for(jk in 1:.lt(label_freq_list)){
  #       lab[[jk]] <- as.integer(colnames(label_freq_list[[jk]])[slct[jk]])
  #       if(jk != 1){
  #         while(lab[[jk]] %in% c(lab[1:(jk-1)])){
  #           slct[jk] <- slct[jk] + 1
  #           # label_freq_list[r]; unlist(lab[r]); slct[r]
  #           if(slct[jk] > n_labels || label_freq_list[[jk]]['n', slct[jk]] == 0){
  #             if(!silent) cat('Unfortunately we found that all the labels in cluster', jk,
  #                         'are colliding with others (we did not count empty buckets). Kept the duplication.\n')
  #             if(multi_cluster_policy == 'select_non_conflict') slct[jk] <- 1
  #             else if(multi_cluster_policy == 'select_non_conflict_spawn') slct[jk] <- 'pop'
  #             print_new_lab <- FALSE
  #             break
  #           }else{
  #             print_new_lab <- TRUE
  #             lab[[jk]] <- as.integer(colnames(label_freq_list[[jk]])[slct[jk]])
  #           }
  #           if(!silent && print_new_lab) cat('Label', lab[[jk]], 'has been found already in another cluster. We will automatically',
  #                                            'select the second most present value in this cluster.\n')
  #         }
  #       }
  #     }
  #   } else if(multi_cluster_policy == 'merge_consecutive'){
  #     if(!silent) cat('Merging policy applied. For the moment only consecutive identically labeled clusters are merged. \n')
  #     slct <- array(1, dim = .lt(label_freq_list))
  #   } 
  #   
  #   # main constructor loop for the selected (slct) labels
  #   lab <- list()
  #   size <- list()
  #   sh_en <- list()
  #   for(jk in 1:.lt(label_freq_list)){
  #     lab[[jk]] <- as.integer(colnames(label_freq_list[[jk]])[slct[jk]])
  #     size[[jk]] <- as.integer(label_freq_list[[jk]]['n', slct[jk]])
  #     sh_en[[jk]] <- label_freq_list[[jk]]['sh_en', slct[jk]]
  #   }  
  #   # fin def
  #   res_bound <- bas$tab.st[,4]
  #   res_label <- unlist(lab)
  # }
  
  # ---------------------------------------------------------------------------------------------- Creation of Predicted vector
  # creating the predicted vector
  predicted_div <- unlist(sapply(1:.lt(res_label), function(x) rep(res_label[x], res_bound[x])))
  # pred_test <- .vec_from_barriers(bar.vec = res_bound, label.vec = res_label)
  # piann_true_v <- ann[st[,3]] # defined before
  missass <- as.numeric(predicted_div != piann_true_v)
  miss_perc <- round(sum(missass)*100/lpiann, 2)
  
  # printing number of missass
  if(!silent) cat('We found roughly', miss_perc, '% missassignment.\n')
  
  # plotting
  if(plot_pred_true_resume || return_plot){
    if(lpiann > max_number_of_elements) {
      s_pi <- round(seq(1, lpiann, length.out = max_number_of_elements))
      fac <- round(lpiann/max_number_of_elements)
      if(!silent) cat('Reducing plotting size by a factor of', fac, '\n')
    } else {
      s_pi <- 1:lpiann
    }
    plot_df <- data.frame('pi' = s_pi, 'predicted' = as.factor(predicted_div[s_pi]), 'true' = as.factor(piann_true_v[s_pi]), 
                          'misass' = as.factor(missass[s_pi]))
    
    if(dbg_score_sapphire) browser()
    gg <- ggplot(data = plot_df) + 
          geom_point(aes_string(y = 'true', x = 'pi', col = shQuote("grey")), size = 8, shape = 108) + 
          geom_point(aes_string(y = 'predicted', x = 'pi', col = shQuote("black")), size = 2, shape = 108) + 
          theme_minimal() + xlab('Progress Index') + ylab('State') + scale_y_discrete(limits = sort(unique(res_label)))
    # scale_color_manual(name = "", labels = c("Predicted", "True"), values = c("black", "lightblue")) +
    #   guides(color = guide_legend(override.aes = list(size=5))) + 
    #   theme(panel.grid = element_blank())
    
    for(gg.i in unique(res_label)) gg <- gg + geom_line(aes_string(y = gg.i, x = 'pi'), size = 0.1, alpha = 0.5)
    gg <- gg + geom_segment(aes_string(y = 0, yend = 0.5, x = 'pi', xend = 'pi', col = 'misass')) + 
          scale_color_manual(name = "", labels = c("Correct", "Miss", "Predicted", "True"), values = c('green4', 'red3', "black", "grey")) +
          guides(color = guide_legend(override.aes = list(size=7))) + 
          theme(panel.grid = element_blank()) 
    gg <- gg + annotate('text', x = lpiann/7.5, y = 0.25, label = paste0('Misses: ', miss_perc, '%'), col = 'white')
    gg <- gg + geom_vline(xintercept = bas$tab.st[,2][-1], size = 0.1, linetype = 'dashed')  
    # gg + geom_ribbon(aes_string(ymin = -0.1, ymax = 0.1, x = 'pi', fill = 'misass')) + 
    #   scale_fill_manual(name = "", labels = c("Correct", "Miss"), values = c("darkgreen", "darkred"))
    
    # cp <- cowplot::plot_grid(gg_popup + theme(legend.position="none") + ggtitle('popup'), 
    #                          gg_merge + theme(legend.position="none") + ggtitle('merge'), gg_keep + ggtitle('keep'), nrow = 1)
    # ggsave('test_diff_policies.png', plot = cp, width = 25, height = 9)
  }
  if(.lt(piann_true_v) != .lt(predicted_div)) 
    stop('Something went wrong. The pred and true vec have differenct lengths and respectively ', .lt(piann_true_v),' and ', .lt(predicted_div))
  
  # ---------------------------------------------------------------------------------------------- Final scoring
  # scoring it with accuracy?
  score.out <- ClusterR::external_validation(true_labels = piann_true_v, clusters = predicted_div, method = scoring_method, summary_stats = FALSE)
  if(!silent) cat('Using', scoring_method,'we obtained a final score of', score.out, '\n')
  
  
  # final output
  return_list <- list('score.out' = score.out, 'label_freq_list' = label_freq_list, 'main_desc' = main_desc, 'perc_miss' = miss_perc)
  if(plot_pred_true_resume) print(gg)
  if(return_predicted) return_list[['predicted_div']] <- predicted_div
  if(return_plot) return_list[['plot']] <- gg
  invisible(return_list)
}

# deprecated
# ----------------------------------------------------------------------------------------
# calculating the entropy of the new selection
# label_freq_list <- list()
# res_b <- 0
# for(jk in 1:.lt(res_label)){
#   label_freq_list[[jk]] <- sort(table(ann[c(st[,3])][(res_b+1):(res_b + res_bound[i])]), decreasing=TRUE)[1:4] * 1.0 / res_bound[i] 
#   # (res_bound[i]-res_b)
#   label_freq_list[[jk]] <- rbind(label_freq_list[[jk]], sort(table(ann[c(st[,3])][(res_b+1):(res_b + res_bound[i])]), decreasing=TRUE)[1:4])
#   label_freq_list[[jk]][is.na(label_freq_list[[jk]])] <- 0 
#   res_b <- res_b + res_bound[i]
# }
# 
# # label_freq_list
# for(jk in 1:.lt(label_freq_list)){
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
# for(jk in 1:.lt(label_freq_list)){ # check label_freq
#   lab[[jk]] <- as.integer(colnames(label_freq_list[[jk]])[1])
#   size[[jk]] <- as.integer(label_freq_list[[jk]][2, 1])
#   sh_en[[jk]] <- label_freq_list[[jk]][1,ncol(label_freq_list[[jk]])]
# }
# 
# max_freq_table <- data.frame(cbind(label = unlist(lab), size = unlist(size), sh_en = unlist(sh_en)))
# max_freq_table
# max_freq_table[order(max_freq_table$sh_en), ]


# external_validation(true_labels = ann[st[,3]], clusters = predicted_div, method = "adjusted_rand_index", summary_stats = FALSE)
# external_validation(true_labels = ann[st[,3]], clusters = predicted_div, method = "jaccard_index", summary_stats = FALSE)
# external_validation(true_labels = ann[st[,3]], clusters = predicted_div, method = "purity", summary_stats = FALSE)



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


