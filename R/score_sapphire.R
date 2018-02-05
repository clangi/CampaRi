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
#' 
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

score_sapphire <- function(the_sap, ann, scoring_method = 'nmi', silent = FALSE, plot_basin_identification = FALSE, nbinsxy = NULL){
  
  # general input checking
  if(!is.character(the_sap) && !is.data.frame(the_sap)) stop("the_sap must be a string or a data frame")
  if(is.character(the_sap) && (!all(grepl("PROGIDX", the_sap)) && !all(grepl("REPIX", the_sap)))) stop("Please provide a data name starting with 'PROGIDX' or 'REPIX'" )
  if(!is.numeric(ann) && (!is.null(dim(ann)))) stop('Please provide an integer vector for ann')
  if(!is.null(nbinsxy) && !.isSingleInteger(nbinsxy)) stop('nbinsxy must be a single integer')
  if(!is.logical(silent)) stop('silent must be a logical')
  if(!is.logical(plot_basin_identification)) stop('plot_basin_identification must be a logical')
  
  # methods input check
  scoring_method.opt <- c("adjusted_rand_index", "jaccard_index", "purity", "nmi")
  if(!(scoring_method[1] %in% scoring_method.opt)) stop("Scoring method option not valid")
  
  # sapphire table loading
  if(!is.data.frame(the_sap))
    st <- as.data.frame(data.table::fread(the_sap, data.table = F))
  else
    st <- the_sap
  
  # hist(st[,5], breaks = 1000) # hist of the distances
  if(is.null(nbinsxy)) nbins_x <- nbins_y <- round(sqrt(nrow(st)*10))
  else nbins_x <- nbins_y <- nbinsxy
  if(!silent) cat('Number of (automatically) selected bins for the basin recognition step is', nbins_x, '\n')
  bas <- CampaRi::basins_recognition(st, nx = nbins_x, ny = nbins_y, dyn.check = 1, plot = plot_basin_identification, out.file = F, silent = silent)
  
  # plot(ann[c(st[,3])], pch='.')
  # sort(table(ann[c(st[,3])]), decreasing=TRUE)[1:4]
  
  major_freq <- list()
  for(jk in 1:nrow(bas$tab.st)){
    major_freq[[jk]] <- sort(table(ann[c(st[,3])][bas$tab.st[jk,2]:bas$tab.st[jk,3]]), decreasing=TRUE)[1:4]/bas$tab.st[jk,4]
    major_freq[[jk]] <- rbind(major_freq[[jk]], sort(table(ann[c(st[,3])][bas$tab.st[jk,2]:bas$tab.st[jk,3]]), decreasing=TRUE)[1:4])
    major_freq[[jk]][is.na(major_freq[[jk]])] <- 0 
  }
  
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

