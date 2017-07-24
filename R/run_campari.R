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
#' @param seq_in This input file specifies the sequence of the system to be simulated or analyzed. Every residue (the base organizational unit) occupies one line and uses a base three-letter representation.
#' This file is needed only in case of a simulation run or a non-ncminer(also ascii then) run.
#' @param multi_threading Default to \code{FALSE}. It will run campari_threads with openmp directives if installed.
#' @param mpi Default to \code{FALSE}. It will run campari_mpi with installed mpi-compiler. Please consider also campari_threads_mpi (also multi_threading=TRUE).
#' @param print_status Default to \code{TRUE}. It will print CAMPARI output on the R console or not otherwise.
#' @param run_in_background Default to \code{FALSE}. It will run CAMPARI in background. This option will disable print_status automatically.
#' @param return_log Default to \code{FALSE}. It will return CAMPARI log as a string (also).
#' @param key_file_input If you provide an already formatted keyfile to this argument, every variable defined in the provided keyfile will be overriden by the ones in this function and CAMPARI will 
#' be run after the generation of a new keyfile. Please note that in the case of multivariable keywords the algorithm will keep only the first value. Consider using \code{\link{keywords_from_keyfile}} function.
#' @param silent Defaults to \code{FALSE}. It will silent all outputs (not the warnings).
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
#' @importFrom parallel detectCores
#' @export run_campari
#' 

run_campari <- function(trj=NULL, base_name='base_name', data_file=NULL, nsnaps=NULL, seq_in=NULL,
                        multi_threading=FALSE, mpi=FALSE, print_status=TRUE, run_in_background=FALSE,
                        key_file_input=NULL, return_log = FALSE, silent = FALSE, ...){
  
  # -----------------------
  #        CHECKS  
  # -----------------------
  args_list <- list(...)
  ascii_mode <- FALSE
  netcdf_mode <- FALSE
  netcdf_mode_amber <- FALSE
  analysis_mode <- FALSE
  simulation_mode <- FALSE
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
  if(!is.logical(run_in_background))
    stop('run_in_background must be a logical.')
  if(!is.logical(return_log))
    stop('return_log must be a logical.')
  if(!is.logical(silent))
    stop('silent must be a logical.')
  if(run_in_background){
    warning('run_in_background option manually forced print_status to FALSE.\n')
    print_status <- FALSE
  }
  if(return_log && (run_in_background)){
    warning('If return_log is active it is not possible to run_in_background. run_in_background is set to FALSE.\n')
    run_in_background <- FALSE
  }
    
  # checks on nsnaps and eventual other inputs.
  if(!is.null(nsnaps) && (!is.numeric(nsnaps) || nsnaps%%1 != 0))
    stop('nsnaps must be an integer.')
  if(!is.null(nsnaps) && ("FMCSC_NCDM_NRFRMS" %in% args_names || "FMCSC_NRSTEPS" %in% args_names))
    stop('Please define nsnaps without using also FMCSC_NCDM_NRFRMS or FMCSC_NRSTEPS. With nsnaps the number of steps will be selected accordingly to the type of input dataset.')
  if("FMCSC_NCDM_NRFRMS" %in% args_names && "FMCSC_NRSTEPS" %in% args_names)
    stop('Only one between FMCSC_NCDM_NRFRMS and FMCSC_NRSTEPS vars must be provided.')  
  
  if(!is.null(data_file) && (!is.character(data_file) || length(nchar(data_file)) > 1))
    stop('data_file must be a single character string.')
  if(!is.null(trj) && (!is.numeric(trj) || length(dim(trj)) != 2))
    stop("trj must be a numeric dataframe/matrix (2 rank).")
  if(!is.null(trj) && nrow(trj) < ncol(trj))
    warning('The inserted trajectory has more variables (columns) than snapshots (rows). Please consider checking again the input data.\n')
  
  # -----------------------
  # Base name checks and needed file setting
  if(!is.character(base_name) || length(nchar(base_name)) > 1)
    stop('base_name must be a single character string.')
  if(!silent) cat('Selected base name for every output/input files (WARNING: they will be overwritten): ', base_name, '\n')
  if(!"FMCSC_BASENAME" %in% args_names) args_list <- c(args_list, FMCSC_BASENAME=base_name)
  
  # seq_in <- paste0(base_name, '.in') # should it be inserted or not??
  log_f <- paste0(base_name, '.log')
  key_f <- paste0(base_name, '.key')
  
  # checking existance of files
  # if(file.exists(seq_in)) warning(seq_in, " is present in the working directory. It will be overwritten.")
  if(file.exists(log_f)) if(!silent) cat(log_f, " is present in the working directory. It will be overwritten.\n")
  if(file.exists(key_f)) if(!silent) cat(key_f, " is present in the working directory. It will be overwritten.\n")
  
  
  # -----------------------
  # time series input check - if not specified trj and data_file this chunck looks for it in the variables. 
  if(!is.null(trj)){
    analysis_mode <- TRUE
    if(!is.null(data_file))
      stop('Only one between trj and data_file must be provided.')
    if(!silent) cat('The trajectory inserted will be written for the analysis in: ', base_name, '.tsv\n')
    data_file <- paste0(base_name, '.tsv')
    if(file.exists(data_file)) warning(data_file, " is present in the working directory. It will be overwritten. \n")
    data.table::fwrite(data.frame(trj), data_file, sep = '\t', row.names = FALSE, col.names = FALSE) 
    if(is.null(nsnaps)) nsnaps <- nrow(trj)
    if(is.null(nvars)) nvars <- ncol(trj) 
    if(!silent) cat('Input trajectory dimensions: (', nsnaps, ', ', nvars, ')\n')
    if(!silent) cat('ANALYSIS MODE\n')
  }else if(is.null(data_file) && 
           any(c('FMCSC_XTCFILE', 'FMCSC_DCDFILE', 'FMCSC_PDBFILE', 'FMCSC_NETCDFFILE', 'FMCSC_NCDM_ASFILE', 'FMCSC_NCDM_NCFILE') %in% args_names)){
    if(any(c('FMCSC_XTCFILE', 'FMCSC_DCDFILE', 'FMCSC_PDBFILE', 'FMCSC_NETCDFFILE') %in% args_names)){
      analysis_mode <- TRUE
      warning('We found analysis keywords while no trj nor data_file was supplied. All the checks for the analysis inputs will be disabled. No NCMINER mode active.\n')
      if(!silent) cat('ANALYSIS MODE (for manual insertion of FMCSC_*FILE keywords)\n')
      if(any(c('FMCSC_NCDM_ASFILE', 'FMCSC_NCDM_NCFILE') %in% args_names)) stop('Once inserted the FMCSC_PDB_FORMAT no NCMINER mode keywords are usable.')
      if(sum(c('FMCSC_XTCFILE', 'FMCSC_DCDFILE', 'FMCSC_PDBFILE', 'FMCSC_NETCDFFILE') %in% args_names) < 1)
        stop('Please supply data_file to this function (or some FMCSC_*FILE following FMCSC_PDB_FORMAT).')
      if(sum(c('FMCSC_XTCFILE', 'FMCSC_DCDFILE', 'FMCSC_PDBFILE', 'FMCSC_NETCDFFILE') %in% args_names) > 1)
        stop('Please use only one file mode (following FMCSC_*FILE described in FMCSC_PDB_FORMAT).')
      if('FMCSC_PDBFILE' %in% args_names){
        data_file <- args_list[['FMCSC_PDBFILE']]
        if(!args_list[['FMCSC_PDB_FORMAT']] %in% c(1, 2)){
          warning('Even if FMCSC_PDB_FORMAT was used for FMCSC_PDBFILE no standard (1 or 2) was selected. It will be defaulted to 1.\n')
          args_list <- c(args_list, FMCSC_PDB_FORMAT=1)
        }
      }
      if('FMCSC_XTCFILE' %in% args_names){
        data_file <- args_list[['FMCSC_XTCFILE']]
        if(!args_list[['FMCSC_PDB_FORMAT']] != 3){
          warning('Even if FMCSC_PDB_FORMAT was used for FMCSC_XTCFILE no standard (3) was selected. It will be defaulted to 3.\n')
          args_list <- c(args_list, FMCSC_PDB_FORMAT=3)
        }
      }
      if('FMCSC_DCDFILE' %in% args_names){
        data_file <- args_list[['FMCSC_DCDFILE']]
        if(!args_list[['FMCSC_PDB_FORMAT']] != 4){
          warning('Even if FMCSC_PDB_FORMAT was used for FMCSC_DCDFILE no standard (4) was selected. It will be defaulted to 4.\n')
          args_list <- c(args_list, FMCSC_PDB_FORMAT=4)
        }
      }
      if('FMCSC_NETCDFFILE' %in% args_names){
        data_file <- args_list[['FMCSC_NETCDFFILE']]
        netcdf_mode_amber <- TRUE
        if(!args_list[['FMCSC_PDB_FORMAT']] != 5){
          warning('Even if FMCSC_PDB_FORMAT was used for FMCSC_NETCDFFILE no standard (5) was selected. It will be defaulted to 5.\n')
          args_list <- c(args_list, FMCSC_PDB_FORMAT=5)
        }
      }
    }else if(any(c('FMCSC_NCDM_ASFILE', 'FMCSC_NCDM_NCFILE') %in% args_names)){
      analysis_mode <- TRUE
      warning('We found analysis keywords for the analysis while no trj R-object was supplied. These NCMINER mode active.\n')
      if(!silent) cat('ANALYSIS MODE (for manual insertion of FMCSC_NCDM_* keyword)\n')
      if('FMCSC_NCDM_ASFILE' %in% args_names && 'FMCSC_NCDM_NCFILE' %in% args_names)
        stop('Use only one between FMCSC_NCDM_ASFILE and FMCSC_NCDM_NCFILE (or use data_file input).')
      if('FMCSC_NCDM_ASFILE' %in% args_names) data_file <- args_list[['FMCSC_NCDM_ASFILE']]
      if('FMCSC_NCDM_NCFILE' %in% args_names) data_file <- args_list[['FMCSC_NCDM_NCFILE']]
      if(is.null(nsnaps) && "FMCSC_NCDM_NRFRMS" %in% args_names) nsnaps <- args_list[["FMCSC_NCDM_NRFRMS"]]
      if('FMCSC_NCDM_ASFILE' %in% args_names) ascii_mode <- TRUE
      if('FMCSC_NCDM_NCFILE' %in% args_names) netcdf_mode <- TRUE
    }
  }else if(!is.null(data_file)){
    analysis_mode <- TRUE
    if(!silent) cat('ANALYSIS MODE\n')
    if(any(c('FMCSC_XTCFILE', 'FMCSC_DCDFILE', 'FMCSC_PDBFILE', 'FMCSC_NETCDFFILE') %in% args_names))
      stop("We found both data_file and one of 'FMCSC_XTCFILE', 'FMCSC_DCDFILE', 'FMCSC_PDBFILE', 'FMCSC_NETCDFFILE' keywords. Please select only one of those options.")
  }else{
    simulation_mode <- TRUE
    if(!silent) cat('SIMULATION RUN: time series (trj) or trajectory file (data_file) not inserted (not even some FMCSC_*FILE).\n')
  }
  
  if(!is.null(data_file)){
    if(!file.exists(data_file))
      stop('Inserted file for the analysis is not present in the directory.')
    if(is.null(nsnaps) && file_ext(data_file) %in% c('tsv', 'dat')) 
      nsnaps <- as.numeric(strsplit(x = suppressWarnings(system(paste0('wc -l ', base_name), intern = TRUE)), split = " ")[[1]][1])
    if(is.null(nvars) && file_ext(data_file) %in% c('tsv', 'dat')) 
      nvars <- as.numeric(strsplit(x = suppressWarnings(system(paste0('wc -n1 ', base_name, '| wc -l '), intern = TRUE)), split =  " ")[[1]][1])
    if(is.null(nsnaps) && "FMCSC_NCDM_NRFRMS" %in% args_names) nsnaps <- args_list[["FMCSC_NCDM_NRFRMS"]]
    if(is.null(nsnaps) && "FMCSC_NRSTEPS" %in% args_names) nsnaps <- args_list[["FMCSC_NRSTEPS"]]
    if(is.null(nsnaps))
      stop('FMCSC_NCDM_NRFRMS or FMCSC_NRSTEPS or nsnaps must be provided for data analysis.')
  }


  # -----------------------
  # checking analysis mode file specifics # reduntant
  if(analysis_mode){
    if(is.null(trj) && is.null(data_file))
      stop('Either a trajectory table or a data_file character must be supplied if in analysis mode.')
    
    # checking the file extension
    if(!file_ext(data_file) %in% c('tsv', 'dat', 'nc', 'xtc', 'dcd', 'pdb'))
      stop("Please provide a data_file with a clear formatted extension (one between 'tsv', 'dat', 'nc', 'xtc', 'dcd', 'pdb').")
    # ascii
    if(file_ext(data_file) %in% c('tsv', 'dat')){
      if(!silent) cat('File mode: ASCII (e.g. tab separated values)\n')
      ascii_mode <- TRUE
      args_list <- c(args_list, FMCSC_NCDM_ASFILE=data_file, FMCSC_NCDM_NRFRMS=nsnaps, FMCSC_NCDM_NRFEATS=nvars)
    # xtc, dcd, pdb
    }else if(file_ext(data_file) %in% c('xtc', 'dcd', 'pdb')){
      args_list <- c(args_list, FMCSC_NRSTEPS=nsnaps)
      # pdb
      if(file_ext(data_file) %in% c('pdb')){
        if(!silent) cat('File mode: PDB (consider changing pdb convention with FMCSC_PDB_R_CONV)\n')
        args_list <- c(args_list, FMCSC_PDBFILE=data_file)
        if(!'FMCSC_PDB_FORMAT' %in% args_names){
          warning('Selected option 1 for FMCSC_PDB_FORMAT because it was not provided.\n')
          args_list <- c(args_list, FMCSC_PDB_FORMAT=1)  
        }
      # xtc
      }else if(file_ext(data_file) %in% c('xtc')){
        if(!silent) cat('File mode: XTC (gromacs convention)\n')
        args_list <- c(args_list, FMCSC_XTCFILE=data_file, FMCSC_PDB_FORMAT=3)      
      # dcd
      }else if(file_ext(data_file) %in% c('dcd')){
        if(!silent) cat('File mode: DCD (CHARMM/NAMD style)\n')
        args_list <- c(args_list, FMCSC_DCDFILE=data_file, FMCSC_PDB_FORMAT=4)
      # nc
      }else if(file_ext(data_file) %in% c('nc')){
        if(netcdf_mode_amber){
          if(!silent) cat('File mode: NetCDF (AMBER style)\n')
          args_list <- c(args_list, FMCSC_NETCDFFILE=data_file, FMCSC_NRSTEPS=nsnaps, FMCSC_PDB_FORMAT=5)
          netcdf_mode_amber <- TRUE
        }else if(netcdf_mode){
          if(!silent) cat('File mode: NetCDF (NCMINER style)\n')
          args_list <- c(args_list, FMCSC_NCDM_NCFILE=data_file, FMCSC_NCDM_NRFRMS=nsnaps)
          netcdf_mode <- TRUE
          # ps: no need of NRFEATS in netcdf
        }
      }
    }
  }else if(simulation_mode){
    warning('All the check for the simulation mode are shallow and the real error messages are delivered directly by the campari library.\n')
  }else{
    stop('mode not valid')
  }
  args_names <- names(args_list) # update of the names
  
  # -----------------------
  # checking and setting the executable
  if((netcdf_mode || ascii_mode) && !netcdf_mode_amber){
    base_exe <- "camp_ncminer"
    if(analysis_mode){
      if('FMCSC_NCDM_ANONAS' %in% args_names && args_list[['FMCSC_NCDM_ANONAS']] != 1)
        stop('ATTENTION: the conversion feature of FMCSC_NCDM_ANONAS (0) is not fully supported and it could have some minor problems (e.g. due to the automatic setting of analysis mode).')
      if(!'FMCSC_NCDM_ANONAS' %in% args_names)
        args_list <- c(FMCSC_NCDM_ANONAS=1, args_list) # This simple logical (1 is true) allows the user to activate the analysis (vs file conversion)
    }
  }else{
    base_exe <- "campari"
    if(analysis_mode){
      # default things for analysis non-ncminer run
      args_list <- c(FMCSC_PDBANALYZE=1,
                     FMCSC_COVCALC=20000000,
                     FMCSC_SAVCALC=20000000,
                     FMCSC_POLCALC=20000000,
                     FMCSC_RHCALC=20000000,
                     FMCSC_INTCALC=20000000,
                     FMCSC_POLOUT=20000000,
                     FMCSC_ENSOUT=20000000,
                     FMCSC_ENOUT=20000000,
                     FMCSC_RSTOUT=20000000,
                     FMCSC_EQUIL=0,
                     FMCSC_SC_IPP=0.0,
                     args_list)
    }
  }
  # eventual add of bashrc PATH exports
  if(!silent) cat('Looking for additional campari bin locations (exports) in ~./bashrc file... ')
  camp_bin_path <- suppressWarnings(system('cat ~/.bashrc | grep PATH | grep camp', intern = T))
  if(length(camp_bin_path) != 0){
    camp_bin_path <- strsplit(paste(camp_bin_path, collapse = 'PATH'), split = "PATH|:", fixed = FALSE)[[1]]
    if(length(camp_bin_path) > 1){
      if(!silent) cat('found.\n')
      camp_bin_path <- paste(camp_bin_path[grep(x = camp_bin_path, pattern = '/')], collapse = ":")
    }else{
      if(!silent) cat('not correct format.\n')
    }
  }else{
    if(!silent) cat('not found.\n')
    camp_bin_path <- ""
  }
  # eventual add of additional aliases
  if(!silent) cat('Looking for additional campari bin locations (aliases) in ~./bashrc file... ')
  camp_bin_alias <- suppressWarnings(system('cat ~/.bashrc | grep alias | grep camp', intern = T))
  if(length(camp_bin_alias) != 0){
    camp_bin_alias <- strsplit(paste(camp_bin_alias, collapse = '='), split = "=|'|\"", fixed = FALSE)[[1]]
    if(length(camp_bin_alias) > 1){
      if(!silent) cat('found.\n')
      camp_bin_alias <- paste(dirname(camp_bin_alias[grep(x = camp_bin_alias, pattern = '/')]), collapse = ":")
    }else{
      if(!silent) cat('not correct format.\n')
    }
  }else{
    if(!silent) cat('not found.\n')
    camp_bin_alias <- ""
  }
  
  # adding the paths to the std PATH variable
  Sys.setenv(PATH=paste(Sys.getenv("PATH"), camp_bin_path, camp_bin_alias, sep=":"))
  
  
  # standard exe
  campari_exe <- suppressWarnings(system(paste0('which ', base_exe), intern = T))
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
  if(!silent) cat(campari_main_exe, 'will be used for this CAMPARI run. \n')
  
  # -----------------------
  # set the default directory
  camp_home <- strsplit(campari_main_exe, split = "bin")[[1]][1]
  # -----------------------
  # defining number of cores
  if(multi_threading){
    n_cores <- parallel::detectCores()
    if("FMCSC_NRTHREADS" %in% args_names){
      if(!silent) cat('The number of threads defined for the openMP run has been set manually to:', args_list[["FMCSC_NRTHREADS"]], '\n')
      if(args_list[["FMCSC_NRTHREADS"]] > n_cores || args_list[["FMCSC_NRTHREADS"]] < 2 || !is.numeric(args_list[["FMCSC_NRTHREADS"]])){
        warning('Selected more/less/notcorrectly cores than the available number of units. It will be set to n_cores - 1')
        args_list <- c(args_list, FMCSC_NRTHREADS=n_cores-1)
      }
    }else{
      if(!silent) cat('Not finding the variable FMCSC_NRTHREADS it will be assigned to n_cores - 1:', n_cores - 1, '\n')
      args_list <- c(args_list, FMCSC_NRTHREADS=n_cores-1)
    }
  }
  
  # -----------------------
  # must exist checks - PARAMETERS
  if("PARAMETERS" %in% args_names){
    paramiters <- args_list[["PARAMETERS"]]
    if(grepl(pattern = "/", paramiters)){
      if(!silent) cat('Found PARAMETERS variable with full path to the parameter file. Please use simply the filename to look directly into the exe directory.\n')
    }else{
      if(!silent) cat('Inserted only filename in the PARAMETERS variable. This file will be searched in the exe directory. To use current directory please add "./" in front of the filename.\n')      
      paramiters <- paste0(camp_home, "params/", paramiters)
    }
  }else{
    paramiters <- paste0(camp_home,"params/abs3.2_opls.prm")  # file defining system energies. Irrelevant fuer blosse Analyse.
    warning('PARAMETERS variable not found. It MUST be supplied, therefore we automatically assign it to: ', paramiters, '.\n')
    if(!silent) cat('PARAMITERS variable AUTOMATICALLY assigned to', paramiters, '\nATTENTION! This file should be correctly assigned to avoid spurious behaviours. \n')
  }
  if(!silent) cat('Using the following specific PARAMETERS:', paramiters, '\n')
  
  if(!file.exists(paramiters)) stop('Parameter file not found in: ', paramiters)
  args_list <- c(args_list, PARAMETERS=paramiters)
  
  # -----------------------
  # must exist checks - SEQFILE
  if("FMCSC_SEQFILE" %in% args_names && !is.null(seq_in))
    stop('Either "FMCSC_SEQFILE" or seq_in file can be supplied each run.')
  if(!is.null(seq_in)){
    args_list <- c(args_list, FMCSC_SEQFILE=seq_in)
    args_names <- names(args_list)
  }
  if(simulation_mode || (analysis_mode && (!any(c(ascii_mode, netcdf_mode))))){
    if("FMCSC_SEQFILE" %in% args_names){
      if(!file.exists(args_list[['FMCSC_SEQFILE']]))
        stop('In this mode an existing sequence file must be specified. It was not found (use se_in).')
      if(!silent) cat('Using the following specific sequence file:', args_list[['FMCSC_SEQFILE']], '\n')
    }else{
      stop('A seq_in file must be provided when in simulation mode or in analysis (not ncminer - ascii nc (not amber)) mode.')
    }
  }else if((ascii_mode || netcdf_mode)){
    if("FMCSC_SEQFILE" %in% args_names)
      warning('Even if in ncminer mode a sequence file was provided. It will not be used because superfluous for this modes.\n')
  }else{
    stop('Sequence file must be supplied when not using the ncminer mode (analysis using *NC* variables - see documentation).')
  }
  
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
  if(!silent) cat('writing of keyfile', key_f, 'completed. \n\n')
  args_names <- names(args_list)
  cat(args_names[1],args_list[[1]],'\n', file = key_f)
  for(i in 2:length(args_list))
    cat(args_names[i],args_list[[i]],'\n', file = key_f, append = TRUE)
  
  
  # -----------------------
  # final run of campari
  if(!silent) cat('Starting campari... \n')
  if(!silent) cat(' --------------------------------------------------------------------------\n')
  if(!silent) cat('                                    CAMPARI                                \n')
  if(!silent) cat(' --------------------------------------------------------------------------\n')
  if(!silent) cat('If not in bakground mode, an error in CAMPARI will be reflected in R.\n')
  if(print_status){
    suppressWarnings(system(paste0(campari_main_exe, " -k ", key_f, " | tee ", log_f)))
    if(any(grepl(x = suppressWarnings(system(paste0("tail -n8 ", log_f), intern = TRUE)), pattern = "CAMPARI CRASHED")))
      stop('======== detected CAMPARI CRASHED in log. Please check it for details ========')
  }else if(run_in_background){
    if(!silent) cat('Direct console printing disabled (it will run in background). Please check', log_f, ' file for real time logging.\n')
    suppressWarnings(system(paste0(campari_main_exe, " -k ", key_f, " > ", log_f, "&")))
  }else{
    if(!silent) cat('Direct console printing disabled. Please check', log_f, ' file for real time logging (at the end it will be tailed).\n')
    suppressWarnings(system(paste0(campari_main_exe, " -k ", key_f, " > ", log_f)))
    if(!silent) cat('\n')
    suppressWarnings(system(paste0("tail -n7 ", log_f)))
    if(any(grepl(x = suppressWarnings(system(paste0("tail -n8 ", log_f), intern = TRUE)), pattern = "CAMPARI CRASHED")))
      stop('======== detected CAMPARI CRASHED in log. Please check it for details ========')
  }
  if(return_log) invisible(readChar(log_f, file.info(log_f)$size))
}
