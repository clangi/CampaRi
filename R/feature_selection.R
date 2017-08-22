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
#' 

select_features <- function(trj, feature_selection = 'pca', n_princ_comp = floor(ncol(trj)/10), pca_method = 'R',
                            plotit = FALSE, frameit = FALSE, return_plot = FALSE, cluster_vector = NULL, plotly_it = FALSE,
                            points_size = 1){

  # Checking input variables (again - it is also a stand alone function)
  if(is.null(n_princ_comp) || (length(n_princ_comp) != 1 || !is.numeric(n_princ_comp) || n_princ_comp < 0 || n_princ_comp > ncol(trj)))
    stop('The number of principal components to use has not been inserted correctly.')
  if(!is.integer(n_princ_comp)) n_princ_comp <- as.integer(n_princ_comp)
  if(n_princ_comp < 2) n_princ_comp <- 2
  if(is.null(feature_selection) || (length(feature_selection) != 1 ||!is.character(feature_selection) || !(feature_selection %in% c('pca'))))
    stop('The used feature selection method is not correctly defined.')
  if(is.null(pca_method) || (length(pca_method) != 1 ||!is.character(pca_method) || !(pca_method %in% c('robust','R'))))
    stop('The inserted pca method is not correctly defined.')
  
  # checking the cluster vector
  if(!is.null(cluster_vector)){
    if(length(cluster_vector) != nrow(trj))
      stop('The cluster vector length is different from the number of snapshots.')
    if(!is.numeric(cluster_vector))
      stop('The cluster vector must be numeric.')
    if(any(cluster_vector %% 1 != 0))
      stop('Only integers can be inserted.')
    col <- 'clusters'
  }else{
    col <- NULL
  }
  
  # Checking the size of the points
  if(!is.numeric(points_size) || length(points_size) != 1)
    stop('The points_size must be a single numeric.')
  
  # Checking the logicals
  if(!is.logical(plotit))
    stop('plotit must be a logical value.')
  if(!is.logical(frameit))
    stop('frameit must be a logical value.')
  if(!is.logical(return_plot))
    stop('return_plot must be a logical value.')
  if(!is.logical(plotly_it))
    stop('plotly_it must be a logical value.')
  
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
    if(plotit) appca <- autoplot(pcahah, data = cbind(trj, clusters = cluster_vector), colour = col, size=0.1*points_size, frame=frameit) + theme_minimal()
    if(return_plot) appca <- autoplot(pcahah, data = cbind(trj, clusters = cluster_vector), colour = col, size=0.1*points_size, frame=frameit) + theme_minimal()
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
  if(plotit) plot(appca)
  if(return_plot){
    cat('Returning the plot object. No selected features will be returned with this mode active.\n')  
    invisible(appca)
  }else{
    cat('Returning the selected features...\n')  
    return(trj[,selected_components])
  }
}
