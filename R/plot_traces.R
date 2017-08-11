#' @title Trace plotting (neural mainly)
#' @description The following functions will plot a trace of valors from a matrix with the index on the y and the snpashots on x (principally used in neurons).
#' @param voltage_traces_mat A logical indicating whether to write the plot to file.
#' @param return_plot if \code{TRUE} it returns the plot object (ggplot).
#' @param sub_sampling_factor factor of subsampling.
#' @param size_line size of the plotted line.
#' @param title Title of the plot (default "")
#' @param xlab x-axis label of the plot (default "")
#' @param ylab y-axis label of the plot (default "")
#' @param highlight_rows which rows should be highlighted in dark red.
#' @param nrows_to_plot The number of vertical neurons (rows) to plot.
#' @details For details, please refer to the main documentation of the original campari software \url{http://campari.sourceforge.net/documentation.html}.
#'
#' @return If \code{return_plot} is active it will return the plot.
#' @seealso
#' \code{\link{mst_from_trj}}, \code{\link{gen_progindex}}, \code{\link{gen_annotation}}.
#'
#' 
#' @export plot_traces
#' @import ggplot2

plot_traces <- function(voltage_traces_mat, return_plot = F, sub_sampling_factor = NULL, title = "", xlab = "", ylab = "", highlight_rows = NULL,
                          size_line = 1, nrows_to_plot = NULL){
  
  # Analysis of extra args
  # input_args <- list(...)
  # avail_extra_argoments <- c('only_timeline',
  #                            'annotate_snap_dist',
  #                            'size_points_on_timeline')
  
  
  # if(any(!(names(input_args) %in% avail_extra_argoments))) 
  #   warning('There is a probable mispelling in one of the inserted variables. Please check the available extra input arguments.')
  # 
  # # Default handling
  # if(!('only_timeline' %in% names(input_args))) only_timeline <- FALSE else only_timeline <- input_args[['only_timeline']]
  # if(!('annotate_snap_dist' %in% names(input_args))) annotate_snap_dist <- FALSE  else annotate_snap_dist <- input_args[['annotate_snap_dist']]
  # if(!('size_points_on_timeline' %in% names(input_args))) size_points_on_timeline <- 0.01 else size_points_on_timeline <- input_args[['size_points_on_timeline']]
  # 
  # ============================
  #          CHECKS
  # ============================
  #
  # ---------------------
  # check on title and axis labels
  if(!is.character(title)) stop("title var must be a string")
  if(!is.character(xlab)) stop("xlab var must be a string")
  if(!is.character(ylab)) stop("ylab var must be a string")
  # ---------------------
  # checking the logicals
  if(!is.logical(return_plot))
    stop('return_plot must be a logical.')
  # ---------------------  
  # size_line = 1 must be a numeric
  if(!is.numeric(size_line) || length(size_line) != 1)
    stop('size_line must be a single numeric.')  
  # ---------------------
  # checking single numbers 
  if(!is.null(sub_sampling_factor) && (!is.numeric(sub_sampling_factor) || length(sub_sampling_factor) != 1 || (n_snaps%%sub_sampling_factor != 0))){
    warning('Wrong sub_sampling_factor insertion. It should be a number divisible by the number of snapshots. Otherwise, it will be used anyway with a truncated ending.')
    # cat('Checking a different value for sub_sampling factor (the first divisible number).')
    # divisible_sub_sampling values
    # for(i in 2:n_snaps){
    #   {if(n_snaps%%i==0) print(sum(trial_len)/i)}
    # }
  }else if(!is.null(sub_sampling_factor)){
    do_subsam_ann <- TRUE
  }else{
    do_subsan_ann <- FALSE
    sub_sampling_factor <- 1
  }
  # ---------------------
  # checking the main input (the traces)
  if(!is.matrix(voltage_traces_mat)){
    if(!is.data.frame(voltage_traces_mat)) stop('voltage_traces_mat input must be a matrix or a data.frame')
    voltage_traces_mat <- as.matrix(voltage_traces_mat)
  }
  if(!is.numeric(voltage_traces_mat)) stop('voltage_traces_mat must be numeric.')
  n_traces <- nrow(voltage_traces_mat)
  n_snaps <- ncol(voltage_traces_mat)
  xx <- seq(from=1, by=1, to=n_snaps)[seq(1, n_snaps, sub_sampling_factor)]
  # ---------------------
  # nrows_to_plot 
  if(!is.null(nrows_to_plot) && is.numeric(nrows_to_plot) && length(nrows_to_plot) == 1L){
    if(nrows_to_plot < 0 || nrows_to_plot > n_traces){
      warning("Inserted background height too small or too big.")
    }else{
      n_traces <- nrows_to_plot
    }
  }
  # ---------------------  
  # checking the row to highlight
  if(!is.null(highlight_rows)){
    if(!is.numeric(highlight_rows) || !all(highlight_rows %in% c(1:n_traces)))
      stop('highlight_rows must be a numeric between 1 and the number of rows to plot.')
  } 
  # ---------------------
  # initial creation of the plot
  xlabel <- ylabel <- ""
  if(xlab != "") xlabel <- xlab
  if(ylab != "") ylabel <- ylab
  gg <- ggplot() +
    xlab(xlabel) + ylab(ylabel) + theme_minimal() 
  # theme(panel.grid.minor = element_line(colour="gray80"))
  
  # ---------------------
  # main vectorization and plot
  
  # var init
  height_one_band <- 1./n_traces
  y_multilines <- NULL
  x_multilines <- NULL
  col_grps <- NULL
  grps <- NULL
  
  # sequential add of the height of the horizontal band
  for(i in 1:n_traces){
    if(any(voltage_traces_mat[i,][seq(1, n_snaps, sub_sampling_factor)] < 0))
      neural_v <- voltage_traces_mat[i,][seq(1, n_snaps, sub_sampling_factor)] + abs(min(voltage_traces_mat[i,][seq(1, n_snaps, sub_sampling_factor)]))
    else
      neural_v <- voltage_traces_mat[i,][seq(1, n_snaps, sub_sampling_factor)]
    max_v <- max(neural_v)
    y_multilines <- c(y_multilines, (neural_v*height_one_band)/max_v + height_one_band*(i-1))
    x_multilines <- c(x_multilines, xx)
    grps <- c(grps, rep(i,n_snaps)[seq(1, n_snaps, sub_sampling_factor)])
    if(!is.null(highlight_rows) && (i %in% highlight_rows))
      col_grps <- c(col_grps, rep('B', n_snaps)[seq(1, n_snaps, sub_sampling_factor)])
    else
      col_grps <- c(col_grps, rep('A', n_snaps)[seq(1, n_snaps, sub_sampling_factor)])
  }
  gg <- gg + geom_line(aes(x = x_multilines,
                              y = y_multilines, colour=col_grps, group = grps),
                          size = 0.1*size_line) 
  gg <- gg + scale_colour_manual(values=c(A="#000000", B="chocolate3")) + theme(legend.position="none")
  gg <- gg + scale_y_continuous(breaks = seq(0, (1 - 1./n_traces), length.out = n_traces), labels = 1:n_traces)
  if(title != "") gg <- gg + ggtitle(title)
  if(return_plot) return(gg)
  else plot(gg)
}
