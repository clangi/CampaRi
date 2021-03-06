% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_traces.R
\name{plot_traces}
\alias{plot_traces}
\title{Trace plotting (neural mainly)}
\usage{
plot_traces(voltage_traces_mat, return_plot = F, sub_sampling_factor = NULL,
  title = "", xlab = "", ylab = "", highlight_rows = NULL,
  size_line = 1, nrows_to_plot = NULL, ...)
}
\arguments{
\item{voltage_traces_mat}{A logical indicating whether to write the plot to file.}

\item{return_plot}{if \code{TRUE} it returns the plot object (ggplot).}

\item{sub_sampling_factor}{factor of subsampling.}

\item{title}{Title of the plot (default "")}

\item{xlab}{x-axis label of the plot (default "")}

\item{ylab}{y-axis label of the plot (default "")}

\item{highlight_rows}{which rows should be highlighted in dark red.}

\item{size_line}{size of the plotted line.}

\item{nrows_to_plot}{The number of vertical neurons (rows) to plot.}

\item{...}{'add_to_this_plot'}
}
\value{
If \code{return_plot} is active it will return the plot.
}
\description{
The following functions will plot a trace of valors from a matrix with the index on the y and the snpashots on x (principally used in neurons).
}
\details{
For details, please refer to the main documentation of the original campari software \url{http://campari.sourceforge.net/documentation.html}.
}
\seealso{
\code{\link{mst_from_trj}}, \code{\link{gen_progindex}}, \code{\link{gen_annotation}}.
}
