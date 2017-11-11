#' @title Shows the clustering summary from log file
#' @description
#'      Using sed and awk to print a formatted cluster summary. log file and FYC file are needed.
#' @param log_file log output of a campari run in which a clustering was performed
#' @param fyc_file the fyc file which was output of the same campari run. This is set only to work with the nbutane example (tutorial 11)
#' @param verbose if verbose is \code{TRUE} the clustering summary will be printed directly in the console (first 10 elements).
#' @param which_first_clusters It defaults to 10 and defines the number of centers to consider in the clustering summary.
#' @param which_cluster you can select a specific cluster element with this option
# @param sub_sampling_factor It defaults to 1 and defines the subsampling that was performed on the trajectory (FMCSC_CCOLLECT). Use it only if
# the analysis was running with that variable different from 1 (use everysnapshot). Otherwise the results will be inconsistent or even erroneous. (fyc option)
#' @param return_centers it returns the cluster centers. (fyc option)
#' @param return_angles it returns the cluster centers' angles. (fyc option)
#' @return see return_*
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


show_clustering_summary <- function(log_file, fyc_file = NULL, verbose = TRUE, which_first_clusters = NULL, which_cluster = NULL, sub_sampling_factor = 1,
                                    return_centers = FALSE, return_angles = FALSE){
  
  if(is.null(which_first_clusters) && is.null(which_cluster)){
    warning('WARNING: neither a number of first cluster to show (which_first_clusters), neither a specific cluster (which_cluster) options used.
            we will continue with the first 10 clusters (which_first_clusters = 10).')
    which_first_clusters <- 10
  }
  if(!is.null(which_first_clusters) && !is.null(which_cluster)){
    warning('WARNING: inserted both a number of first cluster to show (which_first_clusters) and a specific cluster (which_cluster).
            We will continue with the first 10 clusters (which_first_clusters = 10).')
    which_first_clusters <- 10
  }
  
  # some check
  if(!is.logical(verbose))
    stop('verbose must be a logical.')
  if(!is.logical(return_centers))
    stop('return_centers must be a logical.')
  if(!is.logical(return_angles))
    stop('return_angles must be a logical.')
  if(!is.character(log_file) || !file.exists(log_file))
    stop('Log file not existant.')
  if(!is.null(fyc_file) && (!is.character(fyc_file) || !file.exists(fyc_file)))
    stop('FYC file not existant.')
  if(!is.null(which_first_clusters) && !.isSingleInteger(which_first_clusters))
    stop('Wrong insertion of the argument which_first_clusters. It must be a single integer.')
  if(!is.null(which_cluster) && !.isSingleInteger(which_cluster))
    stop('Wrong insertion of the argument which_cluster It must be a single integer.')
  # if(!.isSingleInteger(sub_sampling_factor))
  #   stop('Wrong insertion of the argument sub_sampling_factor. It must be a single integer.')
  
  
  # MAIN PRINTER
  if(verbose || is.null(fyc_file)){
    # Overview of clusters: printing first part of the summary (levels)
    cat('\nPrinting the general overview of the tree... \n')
    cat(.print.cluster.summary(log_file), sep = '\n')
    
    # specific cluster printing
    if(!is.null(which_cluster)){
      cat('\nPrinting the details of one cluster... \n')
      # checking its existance
      visual_string <- suppressWarnings(system(paste("grep -F 'CLUSTER SUMMARY' -A", which_cluster + 1, log_file), intern = TRUE))
      if(any(visual_string == ' ----------------------------------------------------')){
        clu_showable <- which(visual_string == " ----------------------------------------------------") - 3 # number of stuff before the starting
        warning('WARNING: the log file found only ', clu_showable, ' clusters in the log file while you asked for ', which_cluster, '.
                Only the last one we can show will be shown.')
        which_cluster <- clu_showable
      }
      # printing it
      cat(.print.cluster.details(log_file, which_cluster), sep = '\n')
      
    # first tot cluster printing
    }else if(!is.null(which_first_clusters)){
      cat('\nPrinting first', which_first_clusters, 'number of clusters: \n')
      # checking their existance
      visual_string <- suppressWarnings(system(paste("grep -F 'CLUSTER SUMMARY' -A", which_first_clusters + 1, log_file), intern = TRUE))
      # problem: found the end of it
      if(any(visual_string == ' ----------------------------------------------------')){
        clu_showable <- which(visual_string == " ----------------------------------------------------") - 3 # number of stuff before the starting
        warning('WARNING: the log file found only ', clu_showable, ' clusters in the log file while you asked for ', which_first_clusters, '.
                Only the ones we can show will be shown.')
        which_first_clusters <- clu_showable
        # printing them
        cat(visual_string[1:(clu_showable + 3)], sep = '\n')
      }else{
        # printing them
        cat(visual_string, sep = '\n')
      }      
    }
  }
  
  # finding stuff in the FYC file
  if(!is.null(fyc_file)){
    if(is.null(which_first_clusters)) which_first_clusters <- which_cluster
    # checking the number of snashot is correct:
    incorrect_n_snaps <- suppressWarnings(system(paste0("cat ", log_file, " | grep Warning | grep snapshots | grep incorrect"), intern = TRUE))
    if(length(incorrect_n_snaps) != 0)
      stop('It seems that the number of snapshots used for the analysis it is different from the simulated dimensions (check log file for snapshots incorrect).')
    
    # see if there are enough lines
    n_lines <- strsplit(suppressWarnings(system(paste0('wc -l ', fyc_file), intern = TRUE)), split = ' ')[[1]]
    n_lines <- as.numeric(n_lines[n_lines != ''][1])
    
    # selecting the centers:
    clu_centers <- as.numeric(suppressWarnings(system(paste0("grep -F 'CLUSTER SUMMARY' -A ", which_first_clusters + 1, " ", log_file, 
                                                             " | tail -n ", which_first_clusters, " | awk '{print $3}'"), intern = TRUE)))
    
    # handling for the header (selecting and collapsing all the dihedral angles with its residue)
    strings_angles <- system(paste0("printf \"%s\n\" `head -n 1 ", fyc_file, " | sed -e 's/||/ /g' | sed -e 's/|//g' | sed -e 's/I 0/I0/g' | sed -e 's/C 0/C0/g' | \
                                    sed -e 's/\\([A-Z]\\) OME/\\1_OME/g' | \
                                    sed -e 's/\\([A-Z]\\) PHI/\\1_PHI/g' | \
                                    sed -e 's/\\([A-Z]\\) PSI/\\1_PSI/g' | \
                                    sed -e 's/\\([A-Z]\\) CHI/\\1_CHI/g' | \
                                    sed -e 's/\\([A-Z]\\) NUC/\\1_NUC/g'` | \
                              awk -v rs=0 '{if ((NR > 1) && (length($1) > 5)) {rs = rs + 1; split($1, resname, \"_\"); printf(\"%5d : %10s %5d\\n\", NR, $1, rs)};\
                              if ((NR > 1) && (length($1) <= 5)) {printf(\"%5d : %10s %5d\\n\", NR, resname[1]\"_\"$1, rs);}}'"),
                   intern = TRUE)
    # selecting the 4 columns
    names_angs <- unlist(lapply(strsplit(strings_angles, split = "[ ]+"), function(x) x[4]))
    n_angles <- length(names_angs)
    
    # selecting the line from fyc_file:
    if((clu_centers[1] + 1) > n_lines) stop('It seems that you inserted a fyc_file with less clusters than the one you are pointing at from cluster centers.
                                            Maybe wrong fyc_file log_file coupling?')
    angles <- system(paste0('head -n', clu_centers[1] + 1,' ', fyc_file, ' | tail -n 1'), intern = TRUE) 
    for(i in clu_centers[-1]){
      if((i + 1) > n_lines) stop('It seems that you inserted a fyc_file with less clusters than the one you are pointing at from cluster centers.
                                            Maybe wrong fyc_file log_file coupling?')
      angles <- c(angles, system(paste0('head -n', i + 1,' ', fyc_file, ' | tail -n 1'), intern = TRUE))
    }
    
    angles_num <- lapply(strsplit(angles, split = " "), FUN = function(x) as.numeric(x[!is.na(as.numeric(x))]))
    
    cat('\n\t\t')
    for(i in 1:n_angles) cat(names_angs[i], " ")
    cat('\n')
    for(j in 1:which_first_clusters){
      cat(paste0("Center ", clu_centers[j], ":   "))
      for(k in 1:n_angles)
        cat(angles_num[[j]][k + 1], " ")
      cat('\n')
    } 
      
    # cat(paste0("Step ", as.integer(angles_num[[j]][1]), ":    C1C2C3C4: ", angles_num[[j]][2],
    
    if(return_angles && return_centers)
      stop('Only one between return_angles and return_centers must be provided.')
    if(return_centers)
      invisible(clu_centers)
    if(return_angles)
      invisible(angles_num)
  }
}

# Overview of clusters: printing first part of the summary (levels)
.print.cluster.summary <- function(file.log) {
  line <- as.numeric(system(paste("grep -n MAXIMAL", file.log, "| awk -F : '{print $1}'"), intern=TRUE))  ## Line at which clustering information appears
  line.stop <- as.numeric(system(paste("sed -n ", line,"p ", file.log, " | awk '{print $1}'", sep=""), intern=TRUE))  ## On the same line it's written how many levels (!= tree height)
  system(paste0("sed -n ", line-1, ",", line+line.stop-1, "p ", file.log), intern=TRUE)
}

## Print the #clust cluster detail from the log file
.print.cluster.details <- function(file.log, clust) {
  tree.height <- as.numeric(system(paste0("grep \"Number of levels in tree\" ", file.log, " | awk -F : '{print $2}'"), intern=TRUE)) 
  if(clust > tree.height) return(print("Cluster information not available, increase BIRCHMULTI"))
  line <- as.numeric(system(paste0("grep -n \"CLUSTER SUMMARY\" ", file.log, " | awk -F : '{print $1}' | sed -n ", clust,"p"), intern=TRUE))
  return(print(system(paste0("sed -n ", line, ",/Total/p ", file.log), intern=TRUE)))
}