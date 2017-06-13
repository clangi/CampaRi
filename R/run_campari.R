#' @title Frontend for the original CAMPARI library
#' @description
#'      \code{run_campari} will call the original campari software, specified by its specific system-wide commands (e.g. campari, campari_threads).
#'      It is possible to use this function following the instructions in \url{http://campari.sourceforge.net/tutorial11.html}.
#'      For that purpose the keyfile is generated with a standard file name (\code{keyfile.key}) and overwritten automatically in the workind directory.
#'      Please consider that this function (therefore the CAMPARI library) will generate all the files in the working directory.
#'      As policy, the standard arguments of this functions will be overwritten in the specific FMSCS_* if they are supplied (e.g. data_file always
#'      overrides FMCSC_*FILE).
#'
#' @param trj Trajectory to be analysed by campari. It must be in the format (snapshots per variables). It will be written as an .tsv in the current working directory.
#' @param base_name This string can be used for the input/output files, such as (\code{base_name.key}, \code{base_name.in}) and others.
#' @param data_file Input file (e.g. \code{trajectory.dcd}) location. This or \code{trj} must be set.
#' @param nsnaps Number of snapshots in the trajectory file. If the data_file format is not ASCII this variable must be set to the number of snapshots in the trajectory.
#' @param multi_threading Default to FALSE. It will run campari_threads with openmp directives if installed.
#' @param mpi Default to FALSE. It will run campari_mpi with installed mpi-compiler. Please consider also campari_threads_mpi (also multi_threading=TRUE).
#' @param key_file_input If you provide an already formatted keyfile to this argument, every variable defined in the provided keyfile will be overriden by the ones in this function and CAMPARI will 
#' be run after the generation of a new keyfile. Please note that in the case of multivariable keywords the algorithm will keep only the first value. Consider using \code{\link{keywords_from_keyfile}} function.
#' @param ... Analysis variables (similarly to \code{\link{mst_from_trj}}). You can check all of these in the original documentation (\url{http://campari.sourceforge.net/documentation.html}). 
#'
#' @details For details, please refer to the main documentation of the original campari software \url{http://campari.sourceforge.net/documentation.html}.
#' @seealso
#' \code{\link{mst_from_trj}}, \code{\link{gen_progindex}}, \code{\link{gen_annotation}}.
#' @examples
#' \dontrun{
#'  run_campari(...)
#' }
#'
#' @importFrom data.table fread
#' @importFrom tools file_ext
#' 
#' @export run_campari
#' 

run_campari <- function(trj=NULL, base_name='base_name', data_file=NULL, nsnaps=NULL,
                        multi_threading=FALSE, mpi=FALSE, print_status=TRUE, key_file_input=NULL, ...){
  
  # -----------------------
  #        CHECKS  
  # -----------------------
  args_list <- list(...)
  ascii_mode <- FALSE
  netcdf_mode <- FALSE
  analysis_mode <- FALSE
  simulation_mode <- FALSE
  particular_analy <- FALSE
  nvars <- NULL
  
  # set upper case for args names
  names(args_list) <- toupper(toupper(names(args_list)))
  args_names <- names(args_list)
  if(!is.null(key_file_input)){
    cat('Keyfile manually inserted for enhanced keywords feeding. Attention: some standard checks will be overriden and others could be misbehaving.\n')
    args_list <- keywords_from_keyfile(key_file_input = key_file_input, keyword_list = args_list, keyword_list_first = FALSE)
    args_names <- names(args_list)
  }
  
  # -----------------------
  # base input checks
  if(!is.logical(multi_threading))
    stop('multi_threading must be a logical.')
  if(!is.logical(mpi))
    stop('mpi must be a logical.')
  if(!is.logical(print_status))
    stop('print_status must be a logical.')
  if(!is.null(nsnaps) && (!is.numeric(nsnaps) || nsnaps%%1 != 0))
    stop('nsnaps must be an integer.')
  if(!is.null(data_file) && (!is.character(data_file) || length(nchar(data_file)) > 1))
    stop('data_file must be a single character string.')
  if(!is.null(trj) && (!is.numeric(trj) || length(dim(trj)) != 2))
    stop("trj must be a numeric dataframe/matrix (2 rank).")
  if(!is.null(trj) && nrow(trj) < ncol(trj))
    warning('The inserted trajectory has more variables (columns) than snapshots (rows). Please consider checking again the input data.')
  
  # -----------------------
  # Base name checks and needed file setting
  if(!is.character(base_name) || length(nchar(base_name)) > 1)
    stop('base_name must be a single character string.')
  message('Selected base name for every output/input files (WARNING: they will be overwritten): ', base_name)
  if(!"FMCSC_BASENAME" %in% args_names) args_list <- c(args_list, FMCSC_BASENAME=base_name)
  
  # seq_in <- paste0(base_name, '.in') # should it be inserted or not??
  log_f <- paste0(base_name, '.log')
  key_f <- paste0(base_name, '.key')
  
  # checking existance of files
  # if(file.exists(seq_in)) warning(seq_in, " is present in the working directory. It will be overwritten.")
  if(file.exists(log_f)) warning(log_f, " is present in the working directory. It will be overwritten.")
  if(file.exists(key_f)) warning(key_f, " is present in the working directory. It will be overwritten.")
  
  
  # -----------------------
  # time series input check
  if(!is.null(trj)){
    if(!is.null(data_file))
      stop('Only one between trj and data_file must be provided.')
    message('The trajectory inserted will be written for the analysis in: ', base_name, '.tsv')
    data_file <- paste0(base_name, '.tsv')
    if(file.exists(data_file)) warning(data_file, " is present in the working directory. It will be overwritten.")
    data.table::fwrite(data.frame(trj), data_file, sep = '\t', row.names = FALSE, col.names = FALSE) 
    if(is.null(nsnaps)) nsnaps <- nrow(trj)
    if(is.null(nvars)) nvars <- ncol(trj) 
    message('Input trajectory dimensions: (', nsnaps, ', ', nvars, ')')
    analysis_mode <- TRUE
  }else if(!is.null(data_file)){
    if(!file.exists(data_file))
      stop('Inserted file for the analysis is not present in the directory.')
    if(is.null(nsnaps) && file_ext(data_file) %in% c('tsv', 'dat')) nsnaps <- as.numeric(strsplit(suppressWarnings(system(paste0('wc -l ', base_name)), " ")[[1]][1]))
    if(is.null(nvars) && file_ext(data_file) %in% c('tsv', 'dat')) nvars <- as.numeric(strsplit(suppressWarnings(system(paste0('wc -n1 ', base_name, '| wc -l')), " ")[[1]][1]))
    if(is.null(nsnaps) && "FMCSC_NCDM_NRFRMS" %in% args_names) nsnaps <- args_list[["FMCSC_NCDM_NRFRMS"]]
    if(is.null(nsnaps) && !"FMCSC_NRSTEPS" %in% args_names)
      stop('FMCSC_NCDM_NRFRMS or nsnaps must be provided for non-ASCII file data analysis.')
    analysis_mode <- TRUE
  }else{
    if(any(c('FMCSC_PDB_FORMAT') %in% args_names)){
      warning('We found analysis keywords while no trj nor data_file was supplied. All the checks for the analysis inputs will be disabled. No NCMINER mode active.')
      cat('ANALYSIS MODE (for manual insertion of FMCSC_PDB_FORMAT keyword)\n')
      particular_analy <- TRUE
    }
    if(any(c('FMCSC_NCDM_ASFILE', 'FMCSC_NCDM_NCFILE') %in% args_names)){
      warning('We found analysis keywords while no trj nor data_file was supplied. All the checks for the analysis inputs will be disabled. NCMINER mode active.')
      cat('ANALYSIS MODE (for manual insertion of FMCSC_NCDM_* keyword)\n')
      particular_analy <- TRUE
      if('FMCSC_NCDM_ASFILE' %in% args_names) ascii_mode <- TRUE
      if('FMCSC_NCDM_NCFILE' %in% args_names) netcdf_mode <- TRUE
    }
    simulation_mode <- TRUE
  }
  
  # printing the running mode
  if(simulation_mode && !particular_analy)
    cat('SIMULATION RUN: time series (trj) or trajectory file (data_file) not inserted.\n')
  if(analysis_mode)
    cat('ANALYSIS MODE\n')
  
  # -----------------------
  # checking analysis mode file specifics
  if(analysis_mode){
    if(is.null(trj) && is.null(data_file))
      stop('Either a trajectory table or a data_file character must be supplied if in analysis mode.')
    
    # checking the file extension
    if(!file_ext(data_file) %in% c('tsv', 'dat', 'nc', 'xtc', 'dcd', 'pdb'))
      stop("Please provide a data_file with a clear formatted extension (one between 'tsv', 'dat', 'nc', 'xtc', 'dcd', 'pdb').")
    if(file_ext(data_file) %in% c('tsv', 'dat')){
      message('File mode: ASCII (e.g. tab separated values)')
      ascii_mode <- TRUE
      args_list <- c(args_list,
                     FMCSC_NCDM_ASFILE=data_file,
                     FMCSC_NCDM_NRFRMS=nsnaps, 
                     FMCSC_NCDM_NRFEATS=nvars)
    }else{
      args_list <- c(args_list, FMCSC_NRSTEPS=nsnaps)
      if(file_ext(data_file) %in% c('pdb')){
        message('File mode: PDB (consider changing pdb convention with FMCSC_PDB_R_CONV)')
        args_list <- c(args_list, 
                       FMCSC_PDBFILE=data_file,
                       FMCSC_PDB_FORMAT=1) # option two not ready                
      }else if(file_ext(data_file) %in% c('xtc')){
        message('File mode: XTC (gromacs convention)')
        args_list <- c(args_list, 
                       FMCSC_XTCFILE=data_file,
                       FMCSC_PDB_FORMAT=3)      
      }else if(file_ext(data_file) %in% c('dcd')){
        message('File mode: DCD (CHARMM/NAMD style)')
        args_list <- c(args_list, 
                       FMCSC_DCDFILE=data_file,
                       FMCSC_PDB_FORMAT=4) 
      }else if(file_ext(data_file) %in% c('nc')){
        message('File mode: NetCDF (AMBER style)')
        args_list <- c(args_list, 
                       FMCSC_NCDM_NCFILE=data_file) 
        #nsnaps?
        #nvars?
        #ncdm_ananas?
        netcdf_mode <- TRUE
      }
    }
  }else if(simulation_mode && !particular_analy){
    warning('All the check for the simulation mode are shallow and the real error messages are delivered directly by the campari library.')
  }
  
  # -----------------------
  # checking and setting the executable
  if(netcdf_mode || ascii_mode){
    base_exe <- "camp_ncminer"
    if(analysis_mode) args_list <- c(args_list,
                                     FMCSC_NCDM_ANONAS=1) # This simple logical (1 is true) allows the user to activate the analysis (vs file conversion)
  }else{
    base_exe <- "campari"
  }
  
  # eventual add of bashrc PATH exports
  cat('Looking for additional campari bin locations (exports) in ~./bashrc file... ')
  camp_bin_path <- suppressWarnings(system('cat ~/.bashrc | grep PATH | grep camp', intern = T))
  if(length(camp_bin_path) != 0){
    camp_bin_path <- strsplit(paste(camp_bin_path, collapse = 'PATH'), split = "PATH|:", fixed = FALSE)[[1]]
    if(length(camp_bin_path) > 1){
      cat('found.\n')
      camp_bin_path <- paste(camp_bin_path[grep(x = camp_bin_path, pattern = '/')], collapse = ":")
    }else{
      cat('not correct format.\n')
    }
  }else{
    cat('not found.\n')
    camp_bin_path <- ""
  }
  # eventual add of additional aliases
  cat('Looking for additional campari bin locations (aliases) in ~./bashrc file... ')
  camp_bin_alias <- suppressWarnings(system('cat ~/.bashrc | grep alias | grep camp', intern = T))
  if(length(camp_bin_alias) != 0){
    camp_bin_alias <- strsplit(paste(camp_bin_alias, collapse = '='), split = "=|'|\"", fixed = FALSE)[[1]]
    if(length(camp_bin_alias) > 1){
      cat('found.\n')
      camp_bin_alias <- paste(dirname(camp_bin_alias[grep(x = camp_bin_alias, pattern = '/')]), collapse = ":")
    }else{
      cat('not correct format.\n')
    }
  }else{
    cat('not found.\n')
    camp_bin_alias <- ""
  }
  
  # adding the paths to the std PATH variable
  Sys.setenv(PATH=paste(Sys.getenv("PATH"), camp_bin_path, camp_bin_alias, sep=":"))
  
  
  # standard exe
  campari_exe <- suppressWarnings(system(paste0('which ',base_exe), intern = T))
  if(!multi_threading && !mpi && length(campari_exe) == 0) 
    stop('No campari (or camp_ncminer) executable found between the standard commands. Please install the library and alias the bin executables.')
  else if(!multi_threading && !mpi)
    campari_main_exe <- campari_exe
  
  # multi_threading exe
  if(multi_threading){
    campari_threads_exe <- suppressWarnings(system(paste0('which ', base_exe, '_threads'), intern = T))
    if(!mpi && length(campari_threads_exe) == 0) 
      stop('No campari_threads (or camp_ncminer_threads) executable found between the standard commands. Please install the library and alias the bin executables.')
    else if(!mpi)
      campari_main_exe <- campari_threads_exe
  }
  
  # mpi exe
  if(mpi){
    campari_mpi_exe <- suppressWarnings(system(paste0('which ', base_exe, '_mpi'), intern = T))
    if(!multi_threading && length(campari_mpi_exe) == 0) 
      stop('No campari_mpi (or camp_ncminer_mpi) executable found between the standard commands. Please install the library and alias the bin executables.')
    else if(!multi_threading)
      campari_main_exe <- campari_mpi_exe
  }
  
  # mpi and multi_threading exe
  if(mpi && multi_threading){
    campari_mpi_threads_exe <- suppressWarnings(system(paste0('which ', base_exe, '_mpi_threads'), intern = T))
    if(length(campari_mpi_threads_exe) == 0) 
      stop('No campari_mpi_threads (or camp_ncminer_mpi_threads) executable found between the standard commands. Please install the library and alias the bin executables.')
    else
      campari_main_exe <- campari_mpi_threads_exe
  }
  # printing the exe which will be used
  cat(campari_main_exe, 'will be used for this CAMPARI run. \n')
  
  # -----------------------
  # set the default directory
  camp_home <- strsplit(campari_main_exe, split = "bin")[[1]][1]
  
  # -----------------------
  # must exist checks - PARAMETERS
  if("PARAMETERS" %in% args_names)
    paramiters <- args_list[["PARAMETERS"]]
  else
    paramiters <- paste0(camp_home,"/params/abs3.2_opls.prm")  # file defining system energies. Irrelevant fuer blosse Analyse.
  
  if(!file.exists(paramiters)) stop('Parameter file not found in: ', paramiters)
  args_list <- c(args_list, PARAMETERS=paramiters)
  
  # -----------------------
  # Certain variables CANNOT be repeated (or at least must exist specifically)
  # therefore we coerce to have only the last one existing
  must_exist_unique_keys <- c('PARAMETERS', 'FMCSC_SEQFILE')
  for(i in must_exist_unique_keys){
    if(i %in% args_names){ # bug was for SEQFILE when it was not needed
      tmp_k <- args_list[which(names(args_list)==i)] # keep all the elements on the side
      tmp_k <- tmp_k[length(tmp_k)] # take the last element
      args_list <- args_list[-which(names(args_list)==i)] # delete all the duplicate
      args_list <- c(args_list, tmp_k)
    }
  }
  
  
  # -----------------------
  # keyfile writer
  cat('writing of keyfile ', key_f, '\n\n')
  args_names <- names(args_list)
  cat(args_names[1],args_list[[1]],'\n', file = key_f)
  for(i in 2:length(args_list))
    cat(args_names[i],args_list[[i]],'\n', file = key_f, append = TRUE)
  
  
  # -----------------------
  # final run of campari
  cat('Starting campari... \n')
  cat('-------------------------------------------------------\n')
  cat('                         CAMPARI                       \n')
  cat('-------------------------------------------------------\n')
  if(print_status){
    suppressWarnings(system(paste0(campari_main_exe, " -k ", key_f, " | tee ", log_f)))
  }else{
    cat('Direct console printing disabled (it will run in background). Please check', log_f, ' file for real time logging.')
    suppressWarnings(system(paste0(campari_main_exe, " -k ", key_f, " >& ", log_f, "&")))
    
  }
}
