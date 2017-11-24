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
#' @details For details, please refer to the main documentation of the original campari software \url{http://campari.sourceforge.net/documentation.html}.
#' 
#' @export install_campari


install_campari <- function(installation_location = NULL, install_ncminer = FALSE, install_threads = FALSE, install_mpi = FALSE, silent = FALSE){
  
  # insertion checks
  campari_source <- system.file("for_campari/", package = "CampaRi")
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

  ori_wd <- getwd()
  
  if(!is.null(installation_location)){
    cat('Copying source files in installation directory...')
    suppressWarnings(system(paste0('cp -r ', campari_source, ' ', installation_location)))
    installing_place <- installation_location
  }else{
    installing_place <- campari_source
  }
  
  cat('Changing directory...\n')
  setwd(paste0(installing_place, 'source/'))
  cat('\n\n############ installing classic campari executable ############\n\n')
  suppressWarnings(system(command = './configure'))
  suppressWarnings(system(command = 'make campari'))
  if(install_threads && !install_mpi){
    cat('\n\n############ installing campari_threads executable ############\n\n')
    suppressWarnings(system(command = './configure --enable-threads'))
    suppressWarnings(system(command = 'make campari_threads'))
  }
  if(install_mpi && !install_threads){
    cat('\n\n############ installing campari_mpi executable ############\n\n')
    suppressWarnings(system(command = './configure --enable-mpi'))
    suppressWarnings(system(command = 'make campari_mpi'))
  }
  if(install_threads && install_mpi){
    cat('\n\n############ installing campari_mpi_threads executable ############\n\n')
    suppressWarnings(system(command = './configure --enable-threads --enable-mpi'))
    suppressWarnings(system(command = 'make campari_mpi_threads'))
  }
  
  # nc_minare
  if(install_ncminer){
    cat('\n\n############ installing classic campari executable ############\n\n')
    suppressWarnings(system(command = './configure'))
    suppressWarnings(system(command = 'make camp_ncminer'))
    if(install_threads && !install_mpi){
      cat('\n\n############ installing campari_threads executable ############\n\n')
      suppressWarnings(system(command = './configure --enable-threads'))
      suppressWarnings(system(command = 'make camp_ncminer_threads'))
    }
    if(install_mpi && !install_threads){
      stop('nc_miner with mpi is not ready.')
      cat('\n\n############ installing campari_mpi executable ############\n\n')
      suppressWarnings(system(command = './configure --enable-mpi'))
      suppressWarnings(system(command = 'make camp_ncminer_mpi'))
    }
    if(install_threads && install_mpi){
      stop('nc_miner with mpi is not ready.')
      cat('\n\n############ installing campari_mpi_threads executable ############\n\n')
      suppressWarnings(system(command = './configure --enable-threads --enable-mpi'))
      suppressWarnings(system(command = 'make camp_ncminer_mpi_threads'))
    }
  }
  cat('Installation completed.\n')
  cat('Going to the original directory...\n')
  setwd(ori_wd)
}