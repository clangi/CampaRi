#' @title Shows the clustering summary from log file
#' @description
#'      Using sed and awk to print a formatted cluster summary. log file and FYC file are needed.
#' @param log_file log output of a campari run in which a clustering was performed
#' @param fyc_file the fyc file which was output of the same campari run
#' @param verbose if verbose is \code{TRUE} the clustering summary will be printed directly in the console (first 10 elements).
#' @param number_of_clusters It defaults to 10 and defines the number of centers to consider in the clustering summary.
#' @param sub_sampling_factor It defaults to 1 and defines the subsampling that was performed on the trajectory (FMCSC_CCOLLECT). Use it only if 
#' the analysis was running with that variable different from 1 (use everysnapshot). Otherwise the results will be inconsistent or even erroneous.
#' @param return_centers it returns the cluster centers.
#' @param return_angles it returns the cluster centers' angles.
#' @return A named list or a table of keywords readed from the keyfile in input
#' @seealso
#' \code{\link{run_campari}}, \code{\link{sapphire_plot}}.
#' @examples
#' \dontrun{
#'  keywords_from_keyfile("keyfile.key")
#' }
#'
#' 
#' @export show_clustering_summary
#' 


show_clustering_summary <- function(log_file, fyc_file, verbose = TRUE, number_of_clusters = 10, sub_sampling_factor = 1,
                                    return_centers = FALSE, return_angles = FALSE){
  
  # some check
  if(!is.logical(verbose))
    stop('verbose must be a logical.')
  if(!is.logical(return_centers))
    stop('return_centers must be a logical.')
  if(!is.logical(return_angles))
    stop('return_angles must be a logical.')
  if(!is.character(log_file) || !file.exists(log_file))
    stop('Log file not existant.')
  if(!is.character(fyc_file) || !file.exists(fyc_file))
    stop('FYC file not existant.')
  if(!is.numeric(number_of_clusters) || length(number_of_clusters) != 1 || number_of_clusters%%1 != 0)
    stop('Wrong insertion of the argument number_of_clusters. It must be a single integer.')
  if(!is.numeric(sub_sampling_factor) || length(sub_sampling_factor) != 1 || sub_sampling_factor%%1 != 0)
    stop('Wrong insertion of the argument sub_sampling_factor. It must be a single integer.')
  if(verbose)
    system(paste("grep -F 'CLUSTER SUMMARY' -A", number_of_clusters+1, log_file))
  # checking the number of snashot is correct:
  incorrect_n_snaps <- suppressWarnings(system(paste0("cat ", log_file, " | grep Warning | grep snapshots | grep incorrect"), intern = TRUE))
  if(length(incorrect_n_snaps) != 0)
    stop('It seems that the number of snapshots used for the analysis it is different from the simulated dimensions (check log file for snapshots incorrect).')
    
  
  # selecting the centers:
  clu_centers <- as.numeric(system(paste0("grep -F 'CLUSTER SUMMARY' -A ", number_of_clusters + 1, " ", log_file, 
                                          " | tail -n ", number_of_clusters, " | awk '{print $3}'"), intern = TRUE))*sub_sampling_factor + 1
  # selecting the line from fyc_file:
  angles <- system(paste0('head -n', clu_centers[1],' ', fyc_file, ' | tail -n 1'), intern = TRUE) 
  
  for(i in clu_centers[-1])
    angles <- c(angles, system(paste0('head -n', i,' ', fyc_file, ' | tail -n 1'), intern = TRUE))
  
  angles_num <- lapply(strsplit(angles, split = " "), FUN = function(x) as.numeric(x[!is.na(as.numeric(x))]))
  cat('\n')
  for(j in 1:number_of_clusters)  
    cat(paste0("Center ", clu_centers[j]-1, ":    C1C2C3C4: ", angles_num[[j]][2],
               "    C3C2C1H11: ", angles_num[[j]][3], "    C2C3C4H41: ", angles_num[[j]][4]), "\n")
    # cat(paste0("Step ", as.integer(angles_num[[j]][1]), ":    C1C2C3C4: ", angles_num[[j]][2],
  
  if(return_angles && return_centers)
    stop('Only one between return_angles and return_centers must be provided.')
  if(return_centers)
    invisible(clu_centers)
  if(return_angles)
    invisible(angles_num)
}
