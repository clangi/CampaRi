#' @title Create the minimum spanning tree from time series
#' @description
#'      \code{mst_from_trj} creates a minimum spanning tree from a time series (e.g. a trajectory in molecular dynamics) using different distance metrics
#'      between pairwise snapshots.
#'
#' @param trj Input trajectory (variables on the columns and equal-time spaced snpashots on the row). It must be a \code{matrix} or a \code{data.frame} of numeric.
#' @param distance_method Distance metric between snapshots. This value can be set 1 (dihedral angles) or 5 (root mean square deviation) or 11 (balistic distance).
#' @param distance_weights Vector of weights to be applied in order to compute averaged weighted distance using multiple \code{distance_method}.
#' Each value must be between 0 and 1. This option works only if \code{birch_clu=F} and is highly NOT reccomended.
#' @param clu_radius This numeric argument is used in the clustering step in order to make clusters of the same radius at the base level.
#' @param clu_hardcut This option is used only with \code{birch_clu=F} and defines the inter-clusters distance threshold.
#' @param normalize_d A logical that indicates whether the distances must be normalized or not. Usually used with averaging.
#' @param birch_clu A logical that indicates whether the algorithm will use a birch-tree like clustering step (short spanning tree - fast) or it will be generated
#' using a simple leader clustering algorithm (minimum spanning tree).
#' @param min_span_tree This option is used only with \code{birch_clu=F} and defines if the returning adjacency list must be a minimum spanning tree.
#' @param mode It takes a string in input and can be either "fortran" (highly advised and default) or "R".
#' @param rootmax_rad If \code{birch_clu=T} this option defines the maximum radius at the root level of the tree in the advanced clustering algorithm.
#' @param tree_height If \code{birch_clu=T} this option defines the height of the tree in the advanced clustering algorithm.
#' @param n_search_attempts If \code{birch_clu=T} a number of search attempts must be provided for the minimum spanning tree search.
#' @param cores If \code{mode="R"} a complete adjacency matrix can be created in parallel using multiple cores (anyhow slower than "fortran" mode).
#' @param logging If \code{logging=T} the function will print to file the fortran messages ("campari.log").
#' @param ... Various variables. Possible values are \code{c('pre_process', 'window', 'overlapping_reduction','wgcna_type', 'wgcna_power', 'wgcna_corOp','feature_selection', 'n_princ_comp')}
#'
#' @details For more details, please refer to the main documentation of the original campari software \url{http://campari.sourceforge.net/documentation.html}.
#'
#' @return If no netcdf support is available the function will return a list with 3 arguments: node degrees, adjacency list and associated distances.
#' If netcdf support is activated the function will dump the mst in the file "DUMPLING.nc".
#' @seealso
#' \code{\link{adjl_from_progindex}}, \code{\link{gen_progindex}}, \code{\link{gen_annotation}}.
#' @examples
#' adjl <- mst_from_trj(trj = matrix(rnorm(1000), nrow = 100, ncol = 10))
#'
#' \dontrun{
#' adjl <- mst_from_trj(trj = matrix(rnorm(1000),ncol=10,nrow=100),
#' distance_method = 5, clu_radius = 100, clu_hardcut = 100,
#' birch_clu = FALSE, mode = "fortran", logging = FALSE)
#'
#' adjl <- adjl_from_trj(trj = matrix(rnorm(1000),ncol=10,nrow=100),
#' distance_method = 5, clu_radius = 0.1,
#' birch_clu = TRUE, mode = "fortran", rootmax_rad = 1.3, logging = FALSE,
#' tree_height = 5, n_search_attempts = 50)
#' }
#'
#'
#' @importFrom bio3d rmsd
#' @export mst_from_trj
#' @import parallel

mst_from_trj<-function(trj, distance_method = 5, distance_weights = NULL, clu_radius = NULL, clu_hardcut = NULL, #inputs
                       normalize_d = TRUE, birch_clu = FALSE, min_span_tree = TRUE, mode = "fortran", #algo modes
                       rootmax_rad = NULL, tree_height = NULL, n_search_attempts = NULL, #sst default
                       cores = NULL, logging = FALSE, ...){ #misc
  # Checking additional inputs
  input_args <- list(...)
  avail_extra_argoments <- c('pre_process', 'window', 'overlapping_reduction',
                             # specific wgcna vars
                             'wgcna_type', 'wgcna_power', 'wgcna_corOp',
                             # feature selection
                             'feature_selection', 'n_princ_comp')
  if(any(!(names(input_args) %in% avail_extra_argoments)))
    warning('There is a probable mispelling in one of the inserted variables. Please check the available extra input arguments.')

  # wgcna and multiplication pre-processing
  if(!('pre_process' %in% names(input_args))) pre_process <- NULL else pre_process <- input_args[['pre_process']]
  if(!('window' %in% names(input_args))) window <- NULL else window <- input_args[['window']]
  if(!('overlapping_reduction' %in% names(input_args))) overlapping_reduction <- NULL else overlapping_reduction <- input_args[['overlapping_reduction']]

  # feature selection
  if(!('feature_selection' %in% names(input_args))) feature_selection <- NULL else feature_selection <- input_args[['feature_selection']]
  if(!('n_princ_comp' %in% names(input_args))) n_princ_comp <- NULL else n_princ_comp <- input_args[['n_princ_comp']]

  # checking trajectory input
  if(!is.matrix(trj)){
    if(!is.data.frame(trj)) stop('trj input must be a matrix or a data.frame')
    trj <- as.matrix(trj)
  }
  # Memory handling (data management checks)
  data_management <- getOption("CampaRi.data_management")
  if(data_management == "R") cat("Normal memory handling selected. Without hdf5/netcdf backend file management it will be difficult for R to handle big data-sets.\n")
  else if(data_management == "netcdf") cat("Selected data support: netcdf data management.\n")
  else stop("Invalid data management keyword inserted. Check the available methods on the guide.\n")

  if(data_management == "R") cat("To set new data_management method: options(list(CampaRi.data_management = 'R'))\n")

  if(data_management != "R"){
    # warning(paste0('The dumping filename will be ',getOption("CampaRi.data_filename"),'. If already existent it will be overwritten'))
    # TODO in installing routine the check!!
    #     cat("Checking for netcdf support...")
    #     command_loc <- system(paste0("which ",data_management))
    #     if(command_loc=="") stop("No support for hdf5. Please check installation and correct linkage of the command to your enviroment.")
  }

  #Input setting
  n_snaps <- nrow(trj)
  n_xyz <- ncol(trj)



  # -----------------------------------------------------------------------
  # Preprocessing
  #
  # Available options: 'wgcna' 'multiplication'
  #
  #
  #

  preprocessing_opts <- c('wgcna', 'multiplication')
  if(!is.null(pre_process) && (length(pre_process)!=1 || !is.character(pre_process) || !(pre_process %in% preprocessing_opts)))
    stop('Inserted preprocessing method (string) not valid.')
  if(!is.null(pre_process)) cat('Preprocessing mode activated. ')
  if(!is.null(pre_process) && pre_process == 'wgcna'){
    # checking network construction varoables
    if(!is.null(window) && (length(window) != 1 || !is.numeric(window) || window <= 3 || window > n_snaps/2))
      stop('The used window (distance 12) is too small or too big (must be less than half to have sense) or it is simply an erroneus insertion.')
    if((!is.null(overlapping_reduction) && (length(overlapping_reduction) != 1 ||!is.numeric(overlapping_reduction) ||
                                            overlapping_reduction <= 0 || overlapping_reduction > 1)))
      stop('The used overlapping_reduction is not correctly defined. It must be a number between 0 and 1.')

    # setting standard window size
    if(is.null(window)) window <- nrow(trj)/100
    cat('A network will be generated using the WGCNA correlation algorithm and using a sliding window of', window, 'snapshots.\n')
    # Calling generate_network
    trj <- generate_network(trj, window, overlapping_reduction, ...)
    n_xyz <- ncol(trj)
  }
  if(!is.null(pre_process) && pre_process == 'multiplication'){
    cat('A multiplication will be generated copy-pasting dimensionalities from a sliding window of', window, 'snapshots.\n')
    # Calling multiplicate_trj (the )
    trj <- multiplicate_trj(trj, window, overlapping_reduction)
    n_xyz <- ncol(trj)
  }


  # -----------------------------------------------------------------------
  # Feature selection
  #
  # This method has been implemented after the common pre-processing analysis
  # function because this step of feature selection coul be performed only
  # internally while the a priori feature selection was passible also outside
  # this function.
  #
  # Supported algorithms: 'pca'
  #
  if(!is.null(feature_selection)){
    if(is.null(n_princ_comp)) n_princ_comp <- floor(ncol(trj)/10)
    trj <- select_features(trj, feature_selection = feature_selection, n_princ_comp = n_princ_comp)
    n_xyz <- ncol(trj)
  }


  # -----------------------------------------------------------------------
  # -----------------------------------------------------------------------
  # Normal mode (R).
  #
  # This feature is using only R in parallel to calculate all the paired
  # distances. It is incredible LONGY.
  #
  #


  if(is.character(mode)&&mode == "R"){
    n_cores <- detectCores() - 1
    if(!is.null(cores)&&(cores%%1==0)) n_cores=cores
    else warning("No or wrong entry for number of cores: using all of them -1")

    dim<-attributes(trj)$dim
    if(dim>1000) warning('the computation could be incredibly long. We advise to use mode = "fortran" option')

    # Initiate cluster
    cl <- makeCluster(n_cores)
    # clusterExport(cl, "trj")
    # cl <- makeCluster(mc <- getOption("cl.cores", 4))
    # clusterExport(cl=cl, varlist=c("text.var", "ntv", "gc.rate", "pos"))
    # clusterEvalQ(cl, library(rms))

    adjl<-c()
    for(i in 1:dim[1]){
      # if(n_cores>1) clusterExport(cl, "i")
      adjl<-c(adjl,parLapply(cl = cl,X = trj[(i+1):dim[1],],fun = function(x){
        bio3d::rmsd(trj[i,],x)
      }))
    }

    stopCluster(cl)
    warning('No MST made. Please use igraph package and in particular "mst" function in order to continue the analysis')


  # -----------------------------------------------------------------------
  # -----------------------------------------------------------------------
  # Fortran mode.
  #
  # This is the main routine which is using extracted code from campari.
  # Here all the variables specific to the fortran interface are defined.
  #
  #

  }else if(is.character(mode)&&(mode=="fortran")){
    # --------------
    # Default vars
    # --------------

    # Distance value
    max_supported_dist <- 11
    sup_dist <- c(1, 5, 11)
    tmp_dis <- distance_method
    if(!birch_clu){
      # The multiple distance insertion (and weights) are available only with MST
      if(!is.numeric(tmp_dis) ||
         (!all(tmp_dis %in% sup_dist)) || length(tmp_dis)>max_supported_dist)
        stop("The distance values that have been inserted are not supported. The supported values can be checked in the documentation.")
      if(length(tmp_dis) > 1) cat("More than one distance selected. They will be averaged (this feature is available only for MST).")
    }else{
      if(length(distance_method) != 1) stop('When using the birch clustering algorithm (SST) only one distance is available per time.')
      if(!is.numeric(distance_method) || !(distance_method %in% sup_dist)) stop("The distance inserted is not valid. Check the documentation for precise values.")
    }
    distance_method <- rep(0, max_supported_dist)
    distance_method[1:length(tmp_dis)] <- tmp_dis

    # Distance weights
    if(!is.null(distance_weights)){
      tmp_dis_w <- distance_weights
      distance_weights <- rep(1,max_supported_dist)
      if(!is.numeric(tmp_dis_w) || length(tmp_dis_w)>max_supported_dist || tmp_dis_w < 0 || tmp_dis_w > 1)
        warning("Distances are not num or they are not in [0,1]. The option will be turned off")
      else distance_weights[1:length(tmp_dis_w)] <- tmp_dis_w
    }else{
      distance_weights <- rep(1,max_supported_dist)
    }

    # Thresholds for radius and inter radius values. This is MST leader clustering
    if(is.null(clu_radius) || clu_radius <= 0){
      clu_radius <- 214748364
      warning(paste("clu_radius variable (a priori fixed clustering radius) has not been selected.
              A standard value of", clu_radius, "will be used."))
    }
    if(is.null(clu_hardcut) || clu_hardcut <= 0){
      clu_hardcut <- 214748364
      warning(paste("clu_hardcut variable (a priori fixed distance threshold between different cluster members) has not been selected.
              A standard value of", clu_hardcut, "will be used."))
    }
    if(!is.numeric(clu_radius) || length(clu_radius)!=1 || !is.numeric(clu_hardcut) || length(clu_hardcut)!=1)
      stop("clu_radius and clu_hardcut must be a real number.")

    # Logical inputs check
    if(!is.logical(normalize_d))
      stop("Normalization mode must be activated using T/F inputs.")
    if(!is.logical(min_span_tree))
      stop("MST must be enabled using T/F inputs. Using the SST (birch_clu) it is not needed.")
    if(!is.logical(birch_clu))
      stop("SST(birch_clu) mode must be enabled using T/F inputs.")
    if(!is.logical(logging))
      stop("logging mode must be a T/F input.")
    if(birch_clu&&min_span_tree)
      message("MST option is automatically used when birch_clu is activated.")

    # sst checks
    if(birch_clu){
      if(is.null(rootmax_rad))
        rootmax_rad <- max(trj)*2
      else if(!is.numeric(rootmax_rad) || length(rootmax_rad)!=1)
        stop('rootmax_rad must be a numeric of length 1.')

      if(is.null(tree_height) || (is.numeric(tree_height) && length(tree_height)==1 && tree_height<2))
        tree_height <- 5
      else if(!is.numeric(tree_height) || length(tree_height)!=1)
        stop('tree_heigth must be a numeric of length 1.')

      if(is.null(n_search_attempts))
        n_search_attempts <- ceiling(nrow(trj)/10)
      else if(!is.numeric(n_search_attempts) || length(n_search_attempts)!=1)
        stop('n_search_attempts must be a numeric of length 1.')

      if(is.null(clu_radius) || clu_radius==214748364)
        clu_radius <- rootmax_rad/tree_height
    }else{
      rootmax_rad <- 0
      tree_height <- 0
      n_search_attempts <- 0
    }

    # ------------------------------------------------------------------------------
    # Main functions for internal calling of Fortran code
    #
    #

    # -------------
    #    R - old
    # -------------
    if(data_management == "R"){
      #input-output initialization
      if(n_snaps > 25000)
        stop("Using more than 25000 snapshots with no memory handling will generate a memory overflow (tested with 16gb).
             Please set the option data.management to netcdf handler.")
      output_fin <- list()
      max_d <- 0
      adj_deg <- as.integer(rep(0,n_snaps))
      adj_ix <- matrix(as.integer(rep(0,n_snaps*n_snaps)),n_snaps,n_snaps)
      adj_dis <- matrix(as.single(rep(0.0,n_snaps*n_snaps)),n_snaps,n_snaps)
      #double Cstyle deginitions
      attr(trj,"Csingle") <- TRUE
      attr(adj_dis,"Csingle") <- TRUE
      attr(distance_weights,"Csingle") <- TRUE
      #main fortran talker
      output<-.Fortran("generate_neighbour_list", PACKAGE="CampaRi",
                              #input
                              trj_data=trj,
                              n_xyz_in=as.integer(n_xyz),
                              n_snaps_in=as.integer(n_snaps),
                              clu_radius_in=as.single(clu_radius),
                              clu_hardcut_in=as.single(clu_hardcut),
                              #output
                              adjl_deg=adj_deg,
                              adjl_ix=adj_ix,
                              adjl_dis=adj_dis,
                              max_degr=as.integer(max_d),
                              #algorithm details
                              dis_method_in=as.integer(distance_method),
                              dis_weight_in=distance_weights,
                              birch_in=as.logical(birch_clu),
                              mst_in=as.logical(min_span_tree),
                              #sst details
                              rootmax_rad_in=as.single(rootmax_rad),
                              tree_height_in=as.integer(tree_height),
                              n_search_attempts_in=as.integer(n_search_attempts),
                              #modes
                              normalize_dis_in=as.logical(normalize_d),
                              log_print_in=as.logical(logging),
                              verbose_in=as.logical(TRUE))

      #output adjustment
      output_fin[[1]] <- output$adjl_deg
      output_fin[[2]] <- output$adjl_ix[,1:output$max_degr]
      output_fin[[3]] <- output$adjl_dis[,1:output$max_degr]
      return(output_fin)

    # -------------
    #    NetCDF
    # -------------
    }else if(data_management=="netcdf"){
      output_fin <- list()
      #double Cstyle deginitions
      attr(trj,"Csingle") <- TRUE
      attr(distance_weights,"Csingle") <- TRUE
      #main fortran talker
      invisible(.Fortran("generate_neighbour_list_w", PACKAGE="CampaRi",
               #input
               trj_data=trj,
               n_xyz_in=as.integer(n_xyz),
               n_snaps_in=as.integer(n_snaps),
               clu_radius_in=as.single(clu_radius),
               clu_hardcut_in=as.single(clu_hardcut),
               #output
               #none -> it is printed to file
               #algorithm details
               dis_method_in=as.integer(distance_method),
               dis_weight_in=distance_weights,
               birch_in=as.logical(birch_clu),
               # mst_in=as.logical(min_span_tree), #no more specificable
               #sst details
               rootmax_rad_in=as.single(rootmax_rad),
               tree_height_in=as.integer(tree_height),
               n_search_attempts_in=as.integer(n_search_attempts),
               #modes
               normalize_dis_in=as.logical(normalize_d),
               log_print_in=as.logical(logging),
               verbose_in=as.logical(TRUE)))
    }else if(data_management=="h5fc"){
      stop("still to-do")
    }else{
      stop("Data_management assigned to an unknown value. Please complain to the developers.")
    }
  }else{
      stop("Mode entry not correct.")
    }
}
