#' @title Installing the original Fortran code
#' @description
#'      \code{install_campari} is able to install the original campari Fortran code in the default directory (package/inst) or in a specified one. 
#'      Please remember that make and configure are two fundamental steps in this process.
#'
#' @param installation_location It defaults to the package installation path (package/inst/campari) but this can be copied in another directory and installed from there.
#' NB: it will contain a source directory with its configure file and therefore must be executable.
#' @param install_ncminer If true the executable for netcdf and ascii (tsv, csv) file handling will be installed on the top of normal installation.
#' @param install_threads If this option is true and you have some multithreading fortran compiler (e.g. openmp) the campari_threads will be installed.
#' @param install_mpi If this option is true and you have some MPI fortran compiler the campari_mpi will be installed.
#' @param silent_built The configuration and make step will not print anything to console.
#' @param no_optimization The configuration and make step will not print anything to console.
#' @details For details, please refer to the main documentation of the original campari software \url{http://campari.sourceforge.net/documentation.html}.
#' 
#' @export install_campari


install_campari <- function(installation_location = NULL, install_ncminer = FALSE, install_threads = FALSE, install_mpi = FALSE, silent_built = FALSE, no_optimization = FALSE){
  
  # insertion checks
  campari_source <- paste0(system.file('extdata/', "for_campari/", package = "CampaRi"), '/')
  if(!is.null(installation_location)){
    if(!dir.exists(installation_location))
      stop('Inserted directory does not exist. Please take care about the possibility to execute a makefile.')
  }
  if(!is.logical(install_ncminer)) stop('install_ncminer must be a logical value.')
  if(!is.logical(install_threads)) stop('install_threads must be a logical value.')
  if(!is.logical(install_mpi)) stop('install_mpi must be a logical value.')
  if(!is.logical(silent_built)) stop('silent_built must be a logical value.')
  if(!is.logical(no_optimization)) stop('silent_built must be a logical value.')
  if(install_ncminer && install_mpi) stop('ncminer option cannot be installed using mpi.')
  
  
  # checking for install-sh config.sub config.guess and add them eventually
  if(.get_os() == 'linux'){
    automake_path <- '/usr/share/'
    automake_name <- dir(automake_path, pattern = 'automake-*')[1]
  }else if(.get_os() == 'osx'){
    automake_path <- '/usr/local/Cellar/automake/'
    automake_name <- dir(automake_path, pattern = '[0-9].[0.9]*')[1]
  }
  if(!silent_built) cat('Checking for the presence of install-sh, config.guess and config.sub ....')
  if(!file.exists(paste0(campari_source, 'source/install-sh'))){
    if(!file.exists(paste0(automake_path, automake_name, '/install-sh')))
      stop('Please install automake on your computer (e.g. apt install automake) and be sure that there is this file: /usr/share/automake-*/install-sh .')
    file.copy(from = paste0(automake_path, automake_name, '/install-sh'), to = paste0(campari_source, 'source/'), overwrite = TRUE)
  }
  if(!file.exists(paste0(campari_source, 'source/config.sub'))){
    if(!file.exists(paste0(automake_path, automake_name, '/config.sub')))
      stop('Please install automake on your computer (e.g. apt install automake) and be sure that there is this file: /usr/share/automake-*/config.sub .')
    file.copy(from = paste0(automake_path, automake_name, '/config.sub'), to = paste0(campari_source, 'source/'), overwrite = TRUE)
  }
  if(!file.exists(paste0(campari_source, 'source/config.guess'))){
    if(!file.exists(paste0(automake_path, automake_name, '/config.guess')))
      stop('Please install automake on your computer (e.g. apt install automake) and be sure that there is this file: /usr/share/automake-*/config.guess .')
    file.copy(from = paste0(automake_path, automake_name, '/config.guess'), to = paste0(campari_source, 'source/'), overwrite = TRUE)
  }
  if(!silent_built) cat('done\n')
  
  # specific installation handling
  if(!is.null(installation_location)){
    if(!silent_built) cat('Copying source files in installation directory...')
    installation_location <- suppressWarnings(system(paste0('pwd ', installation_location), intern = T)) # check to unfold '..' which would break make
    suppressWarnings(system(paste0('cp -r ', campari_source, ' ', installation_location), ignore.stdout = silent_built))
    installing_place <- paste0(installation_location, '/for_campari/')
    if(!silent_built) cat('done\n')
  }else{
    installing_place <- campari_source
  }
  
  cat(paste0('Specifing the directory: ', installing_place, ' ...\n'))
  confit <- paste0(installing_place, '/source/configure --with-campari-home=', installing_place)
  if(no_optimization) confit <- paste(confit, '--enable-fast-compilation')
  
  # stanard installations
  .installer_conf.make(make_name = 'campari', configure_base = confit, extra_directives = '', silent_built = silent_built)
  if(install_threads)
    .installer_conf.make(make_name = 'campari_threads', configure_base = confit, extra_directives = '--enable-threads', silent_built = silent_built)
  if(install_mpi)
    .installer_conf.make(make_name = 'campari_mpi', configure_base = confit, extra_directives = '--enable-mpi', silent_built = silent_built)
  if(install_threads && install_mpi)
    .installer_conf.make(make_name = 'campari_mpi_threads', configure_base = confit, extra_directives = '--enable-threads --enable-mpi', silent_built = silent_built)
  
  # ncminer installations
  if(install_ncminer){
    .installer_conf.make(make_name = 'camp_ncminer', configure_base = confit, extra_directives = '', silent_built = silent_built)
    if(install_threads)
      .installer_conf.make(make_name = 'camp_ncminer_threads', configure_base = confit, extra_directives = '--enable-threads', silent_built = silent_built)
    if(install_mpi)
      .installer_conf.make(make_name = 'camp_ncminer_mpi', configure_base = confit, extra_directives = '--enable-mpi', silent_built = silent_built)
    if(install_threads && install_mpi)
      .installer_conf.make(make_name = 'camp_ncminer_mpi_threads', configure_base = confit, extra_directives = '--enable-threads --enable-mpi', silent_built = silent_built)
  }
  if(!silent_built) cat('Installation completed.\n')
}

.installer_conf.make <- function(make_name, configure_base, extra_directives = '', silent_built = FALSE){
  
  # simple printing routine
  seps.l <- .get_adjusted_separators(tit = paste0(' installing ', make_name,' executable '), length_it = 68)
  
  # definition and show of the commands
  if(!silent_built) cat('\n\n', seps.l$sep1, '\n', sep = '')
  conf_cmd <- paste(configure_base, extra_directives)
  make_cmd <- paste('make', make_name)
  if(!silent_built) cat('Running the following commands in sequence.\n')
  if(!silent_built) cat('CONFIGURE:     ', conf_cmd, '\n')
  if(!silent_built) cat('MAKE:          ', make_cmd, '\n')
  if(!silent_built) cat(seps.l$sep2, '\n')   
  
  # configure
  cmd.out <- suppressWarnings(system(conf_cmd, ignore.stdout = silent_built, ignore.stderr = silent_built, intern = silent_built)) 
  if(.check_cmd.out(cmd.out, intern = silent_built)) stop('Error during ./configure. Please check the logging. Following command failed:\n', conf_cmd)
  if(!silent_built) cat('\n\n')
  
  # make
  cmd.out <- suppressWarnings(system(make_cmd, ignore.stdout = silent_built, ignore.stderr = silent_built, intern = silent_built))
  if(.check_cmd.out(cmd.out, intern = silent_built)) stop('Error during make. Please check the logging. Following command failed:\n', make_cmd)
  
  # end
  if(!silent_built) cat('\nInstallation ended successfully.\n')
  if(!silent_built) cat(seps.l$sep2, '\n')   
}

