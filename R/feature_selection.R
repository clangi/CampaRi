#' @title Multiple pre-processing methods for feature selection
#' @description
#'      \code{select_features} is able to select input variables on the basis of the trajectory input. For the moment only PCA-based feature selection
#'      is supported. Moreover, this tool is meant to be used with the total trajectory input.
#'
#' @param trj Input trajectory (variables on the columns and equal-time spaced snpashots on the row). It must be a \code{matrix} or a \code{data.frame} of numeric.
#' @param feature_selection Available method is 'pca'
#' @param n_princ_comp number of principal components to use
#' @param pca_method If set 'R' (default) it will use \code{\link{princomp}}. The other (slower) option is 'robust' which is using \code{\link{PCAproj}}.
#' @param plotit Plot the PCA components if two are selected.
#' @param frameit Add a frame (shaded clustering) of the whole performed PCA.
#' @param return_plot This option is usually used to add layers to the ggplot (made using autoplot).
#' @param cluster_vector This option can be used to set the clusters you want and show them with different colors (and shades if \code{frameit = TRUE}). 
#' Please set this option with the same dimensionality of the trj (n_snapshots) and use integer numbers (to define the clusters).
#' @param points_size It must be a number and it defines the size of the points.
#'
#' @details This function is based primarly on the basic R function \code{pricomp} and on \code{PCAproj} from the package pcaPP. Insead, for more details on the SAPPHIRE anlysis, please refer to the main documentation
#' of the original campari software \url{http://campari.sourceforge.net/documentation.html}.
#'
#' @return It will return a modified trajectory matrix and print the principal components vector.
#' @seealso
#' \code{\link{princomp}}, \code{\link{PCAproj}}, \code{\link{adjl_from_progindex}}, \code{\link{gen_progindex}}, \code{\link{gen_annotation}}.
# @examples
#'
#'
#' @export select_features
#' @importFrom data.table fwrite
#' @importFrom pcaPP PCAproj
#' @importFrom stats princomp
#' @importFrom stats loadings
#' @importFrom plotly ggplotly
#' @import ggplot2
#' @import ggfortify

select_features <- function(trj, feature_selection = 'pca', n_princ_comp = floor(ncol(trj)/10), pca_method = 'R',
                            plotit = FALSE, frameit = FALSE, return_plot = FALSE, cluster_vector = NULL, plotly_it = FALSE,
                            points_size = 1, specific_palette = NULL, plot_legend = FALSE, legend_title = NULL, legend_labels = NULL){

  # Checking input variables (again - it is also a stand alone function)
  if(is.null(n_princ_comp) || (length(n_princ_comp) != 1 || !is.numeric(n_princ_comp) || n_princ_comp < 0 || n_princ_comp > ncol(trj)))
    stop('The number of principal components to use has not been inserted correctly.')
  if(!is.integer(n_princ_comp)) n_princ_comp <- as.integer(n_princ_comp)
  if(n_princ_comp < 2) n_princ_comp <- 2
  if(is.null(feature_selection) || (length(feature_selection) != 1 ||!is.character(feature_selection) || !(feature_selection %in% c('pca'))))
    stop('The used feature selection method is not correctly defined.')
  if(is.null(pca_method) || (length(pca_method) != 1 ||!is.character(pca_method) || !(pca_method %in% c('robust','R'))))
    stop('The inserted pca method is not correctly defined.')
  
  # Checking the size of the points
  if(!is.numeric(points_size) || length(points_size) != 1)
    stop('The points_size must be a single numeric.')
  
  # checking the palette
  if(!is.null(specific_palette) &&
     (!all(is.character(specific_palette)) || !all(nchar(specific_palette)==7) || !all(sapply(specific_palette, function(x){substring(x,1,1)}) == "#")))
    stop('specific_palette must be of the type "#b47b00", "#D9E000" (7 characters starting with #')
  
  # Checking the logicals
  if(!is.logical(plotit))
    stop('plotit must be a logical value.')
  if(!is.logical(frameit))
    stop('frameit must be a logical value.')
  if(!is.logical(return_plot))
    stop('return_plot must be a logical value.')
  if(!is.logical(plotly_it))
    stop('plotly_it must be a logical value.')
  
  # adding the cross checks
  if(plotit){
    # checking the cluster vector
    if(!is.null(cluster_vector)){
      if(length(cluster_vector) != nrow(trj))
        stop('The cluster vector length is different from the number of snapshots.')
      if(!is.numeric(cluster_vector))
        stop('The cluster vector must be numeric.')
      if(any(cluster_vector %% 1 != 0))
        stop('Only integers can be inserted.')
      col <- 'Clusters'
    }else{
      col <- NULL
      cluster_vector <- rep(1, nrow(trj))
    }
    if(!plot_legend && (!is.null(legend_title) || !is.null(legend_labels))){
      warning('Inserted legend_title or legend_labels or annotation_type variables WITHOUT activating plot_legend. This option will be turned on automatically.')
      plot_legend <- TRUE
    }
    
    # checking the legend
    if(plot_legend){
      # checking the legend title
      if(!is.null(legend_title)){
        if(!is.character(legend_title))
          stop('legend_title must be a character.')
        if(length(legend_title) != 1)
          stop('legend_title must be a SINGLE character.')
        
        leg_tit <- legend_title
      }else{
        leg_tit <- 'Clusters'
      }
      col <- leg_tit # to make them overlap (shade and legend)
      
      # checking the labels
      if(!is.null(legend_labels)){
        if(!is.character(legend_labels))
          stop('legend_labels must be a vector of characters.')
        if(!is.null(dim(legend_labels)))
          stop('legend_labels must be a vector of characters (found more than 1 dimension).')
        leg_lab <- c(legend_labels)
      }else{
        warning('No label inserted for the legend. they will be assigned to integers.')
        leg_lab <- NULL
      }
      # checking the legend vs the palette etc
      cluster_vector <- factor(cluster_vector)
      if(!is.null(specific_palette)){
        if(length(specific_palette) != length(levels(cluster_vector)))
          stop('When using a specific palette please set it to the same number of factors in the cluster_vector.')
      }else{
        specific_palette <- levels(cluster_vector)
      }
      if(is.null(leg_lab)) leg_lab <- seq(1, length(specific_palette))
      if(!is.null(leg_lab) && length(leg_lab) != length(specific_palette))
        stop('Inserted different number of labels than colors. In this mode this cannot be done.')
      if(!is.null(leg_lab) && length(leg_lab) != length(levels(cluster_vector)))
        stop('Inserted different number of labels than factors in the cluster_vectors (e.g. different numbers 1-2-3).') 
      if(length(levels(cluster_vector)) > 10)
        warning('There are a lot of levels (clusters > 10).')
      if(length(levels(cluster_vector)) > 20){
        stop('There are too many levels (clusters > 20).')
      }
    }
    if(plotly_it && frameit)
      stop('plotly_it and frameit creates weird legends if this last is plotted. Avoid to use these options toghether')
  }else{
    if(!is.null(cluster_vector))
      stop('Please put plotit variables TRUE if you want to define the clusters to plot.')
    if(frameit)
      stop('Please put plotit variables TRUE if you want to frame the clusters to plot.')
    if(plotly_it)
      stop('Please put plotit variables TRUE if you want to plot the pca using plotly.')
    if(return_plot)
      stop('Please put plotit variables TRUE if you want to return the plot object.')
    if(plot_legend || !is.null(legend_title) || !is.null(specific_palette) || !is.null(legend_labels))
      stop('You inserted something between plot_legend, legend_title, specific_palette, legend_labels but you did not activate the plotting.')
  }
  
  # General messages
  cat('Feature selection mode active.', feature_selection, 'dimensionality reduction will be performed.\n')

  # Long calculation warning
  if(dim(trj)[1] > 20000) warning('The dimensionality reduction can be really long.')

  # Check for NA
  if(any(is.na(trj))) stop('There are NA values in the input trajectory')
    
  # Starting the PCA
  if(feature_selection == 'pca'){
    if(pca_method == 'robust')
      pcahah <- pcaPP::PCAproj(trj, k = n_princ_comp, method = "mad", CalcMethod = "eachobs",
                               nmax = 1000, update = TRUE, scores = TRUE, maxit = 5, maxhalf = 5, scale = NULL, zero.tol = 1e-16)
    else pcahah <- princomp(x = trj)
    if(pca_method == 'robust') cat('Robust method selected. The pca-selection will be using robust PCA by projection pursuit.\n')
    selected_components <- array(0, n_princ_comp)
    for(i in seq(n_princ_comp)){
      selected_components[i] <- which.max(abs(loadings(pcahah)[,i]))
    }
    # For printing the first 10 most important variables (most indipendent)
    if(n_princ_comp > 10) more_selected <- selected_components[1:10]
    else more_selected <- selected_components
    # impossible check
    if(any(selected_components == 0)) stop('An internal error occured. Please refer it to the mantainers.')
  
    # plotting it eventually
    if(plotit || return_plot){
      data_to_plot <- as.data.frame(cbind(trj, cls = cluster_vector))
      data_to_plot[, ncol(data_to_plot)] <- as.factor(data_to_plot[, ncol(data_to_plot)])
      names(data_to_plot)[length(names(data_to_plot))] <- col 
      appca <- autoplot(pcahah, data = data_to_plot, colour = col, size=0.1*points_size, frame=frameit, frame.type = 'norm', frame.colour = col) + theme_minimal()
      # plot the legend
      if(plot_legend){
        appca <- appca + scale_color_manual(name = leg_tit, values = specific_palette, labels = leg_lab) + 
          guides(color = guide_legend(override.aes = list(size=points_size)))
      }else{
        if(!is.null(specific_palette))
          appca <- appca + scale_color_gradientn(colours = specific_palette, guide = FALSE)
        else
          appca <- appca + theme(legend.position='none')
      }
    } 
    if(plotly_it) {
      if(nrow(trj) > 15000)
        warning('nrow(trj) is more than 15000. Ggplotly would be stupidly slow and crash. Consider reducing the dimensionality. It will plotted normally.')
      else
        appca <- ggplotly(appca)
    }
    # Printing simple output and to file
    cat('PCA performed successfully.\n')
    cat('The most indipendent values selected are:', more_selected)
    message('The details of the pca output will be oscurated but the selected values will be written in "selected_pca.out" ASCII file.')
    fwrite(x = as.list(selected_components), file = "selected_pca.out", sep = '\t')
  }
  
  if(return_plot){
    cat('Returning the plot object. No selected features will be returned with this mode active.\n')  
    invisible(appca)
  }else{
    if(plotit) print(appca)
    cat('Returning the selected features...\n')
    if(!plotit)
      return(trj[,selected_components])
    else
      invisible(trj[,selected_components])
  }
}
