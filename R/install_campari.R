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
#' @details For details, please refer to the main documentation of the original campari software \url{http://campari.sourceforge.net/documentation.html}.
#' 
#' @export install_campari


install_campari <- function(installation_location = NULL, install_ncminer = FALSE, install_threads = FALSE, install_mpi = FALSE, silent_built = FALSE){
  
  # insertion checks
  campari_source <- system.file('extdata/', "for_campari/", package = "CampaRi")
  cat(campari_source, '\n')
  if(!is.null(installation_location)){
    if(!dir.exists(installation_location))
      stop('Inserted directory does not exist. Please take care about the possibility to execute a makefile.')
  }
  if(!is.logical(install_ncminer))
    stop('install_ncminer must be a logical value.')
  if(!is.logical(install_threads))
    stop('install_threads must be a logical value.')
  if(!is.logical(install_mpi))
    stop('install_mpi must be a logical value.')
  if(!is.logical(silent_built))
    stop('silent_built must be a logical value.')
  if(install_ncminer && install_mpi)
    stop('ncminer option cannot be installed using mpi.')
  
  
  # checking for install-sh config.sub config.guess and add them eventually
  if(.get_os() == 'linux'){
    automake_path <- '/usr/share/'
    automake_name <- dir(automake_path, pattern = 'automake-*')[1]
  }else if(.get_os() == 'osx'){
    automake_path <- '/usr/local/Cellar/automake/'
    automake_name <- dir(automake_path, pattern = '[0-9].[0.9]*')[1]
  }
  cat('Checking for the presence of install-sh, config.guess and config.sub ....')
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
  cat('done\n')
  
  ori_wd <- getwd()
  
  if(!is.null(installation_location)){
    cat('Copying source files in installation directory...')
    suppressWarnings(system(paste0('cp -r ', campari_source, ' ', installation_location), ignore.stdout = silent_built))
    installing_place <- paste0(installation_location, '/for_campari/')
    cat('done\n')
  }else{
    installing_place <- campari_source
  }
  
  cat('Changing directory...\n')
  setwd(paste0(installing_place, 'source/'))
  cat('\n\n############ installing classic campari executable ############\n')
  # some notes on the following: cat is sending to console out the intern = T (needed for error catching). ignore.stderr = T is necessary for
  # not showing the error in console.
  tryCatch(cat(suppressWarnings(system('./configure', ignore.stdout = silent_built, ignore.stderr = T, intern = T)), sep = '\n'), 
           error = function(err) stop('Error during ./configure. Please check the logging.'))
  cat('\n#\n')
  tryCatch(cat(suppressWarnings(system('make campari', ignore.stdout = silent_built, ignore.stderr = T, intern = T)), sep = '\n'), 
           error = function(err) stop('Error during make campari. Please check the logging.'))
  cat('\n############ INSTALLED ############\n\n')
  if(install_threads){
    cat('\n\n############ installing campari_threads executable ############\n')
    tryCatch(cat(suppressWarnings(system('./configure --enable-threads', ignore.stdout = silent_built, ignore.stderr = T, intern = T)), sep = '\n'), 
             error = function(err) stop('Error during ./configure --enable-threads. Please check the logging.'))
    cat('\n#\n')
    tryCatch(cat(suppressWarnings(system('make campari_threads', ignore.stdout = silent_built, ignore.stderr = T, intern = T)), sep = '\n'), 
             error = function(err) stop('Error during make campari_threads. Please check the logging.'))
    cat('\n############ INSTALLED ############\n\n')
  }
  if(install_mpi){
    cat('\n\n############ installing campari_mpi executable ############\n')
    tryCatch(cat(suppressWarnings(system('./configure --enable-mpi', ignore.stdout = silent_built, ignore.stderr = T, intern = T)), sep = '\n'), 
             error = function(err) stop('Error during ./configure --enable-mpi. Please check the logging.'))
    cat('\n#\n')
    tryCatch(cat(suppressWarnings(system('make campari_mpi', ignore.stdout = silent_built, ignore.stderr = T, intern = T)), sep = '\n'), 
             error = function(err) stop('Error during make campari_mpi. Please check the logging.'))
    cat('\n############ INSTALLED ############\n\n')
  }
  if(install_threads && install_mpi){
    cat('\n\n############ installing campari_mpi_threads executable ############\n')
    tryCatch(cat(suppressWarnings(system('./configure --enable-threads --enable-mpi', ignore.stdout = silent_built, ignore.stderr = T, intern = T)), sep = '\n'), 
             error = function(err) stop('Error during ./configure --enable-threads --enable-mpi. Please check the logging.'))
    cat('\n#\n')
    tryCatch(cat(suppressWarnings(system('make campari_mpi_threads', ignore.stdout = silent_built, ignore.stderr = T, intern = T)), sep = '\n'), 
             error = function(err) stop('Error during make campari_mpi_threads. Please check the logging.'))
    cat('\n############ INSTALLED ############\n\n')
  }
  
  # nc_minare
  if(install_ncminer){
    cat('\n\n############ installing classic camp_ncminer executable ############\n')
    tryCatch(cat(suppressWarnings(system('./configure', ignore.stdout = silent_built, ignore.stderr = T, intern = T)), sep = '\n'), 
             error = function(err) stop('Error during ./configure. Please check the logging.'))
    cat('\n#\n')
    tryCatch(cat(suppressWarnings(system('make camp_ncminer', ignore.stdout = silent_built, ignore.stderr = T, intern = T)), sep = '\n'), 
             error = function(err) stop('Error during make camp_ncminer. Please check the logging.'))
    cat('\n############ INSTALLED ############\n\n')
    if(install_threads){
      cat('\n\n############ installing camp_ncminer_threads executable ############\n')
      tryCatch(cat(suppressWarnings(system('./configure --enable-threads', ignore.stdout = silent_built, ignore.stderr = T, intern = T)), sep = '\n'), 
               error = function(err) stop('Error during ./configure --enable-threads. Please check the logging.'))
      cat('\n#\n')
      tryCatch(cat(suppressWarnings(system('make camp_ncminer_threads', ignore.stdout = silent_built, ignore.stderr = T, intern = T)), sep = '\n'), 
               error = function(err) stop('Error during make camp_ncminer_threads. Please check the logging.'))
      cat('\n############ INSTALLED ############\n\n')
    }
    if(install_mpi){
      stop('nc_miner with mpi is not ready.')
      cat('\n\n############ installing camp_ncminer_mpi executable ############\n')
      tryCatch(cat(suppressWarnings(system('./configure --enable-mpi', ignore.stdout = silent_built, ignore.stderr = T, intern = T)), sep = '\n'), 
               error = function(err) stop('Error during ./configure --enable-mpi. Please check the logging.'))
      cat('\n#\n')
      tryCatch(cat(suppressWarnings(system('make camp_ncminer_mpi', ignore.stdout = silent_built, ignore.stderr = T, intern = T)), sep = '\n'), 
               error = function(err) stop('Error during make camp_ncminer_mpi. Please check the logging.'))
      cat('\n############ INSTALLED ############\n\n')
    }
    if(install_threads && install_mpi){
      stop('nc_miner with mpi is not ready.')
      cat('\n\n############ installing camp_ncminer_mpi_threads executable ############\n')
      tryCatch(cat(suppressWarnings(system('./configure --enable-threads --enable-mpi', ignore.stdout = silent_built, ignore.stderr = T, intern = T)), sep = '\n'), 
               error = function(err) stop('Error during ./configure --enable-threads --enable-mpi. Please check the logging.'))
      cat('\n#\n')
      tryCatch(cat(suppressWarnings(system('make camp_ncminer_mpi_threads', ignore.stdout = silent_built, ignore.stderr = T, intern = T)), sep = '\n'), 
               error = function(err) stop('Error during make camp_ncminer_mpi_threads. Please check the logging.'))
      cat('\n############ INSTALLED ############\n\n')
    }
  }
  cat('Installation completed.\n')
  cat('Going to the original directory...\n')
  setwd(ori_wd)
}