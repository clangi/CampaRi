setting_up_netcdf<-function(package_dir, autolocating = T, no_ncdf = F){
  if(!autolocating){
    cat("Press enter to continue without netcdf support or insert the directory 
        location of the netcdf.mod file (you can use 'locate netcdf.mod'): ")
    dir_netcdf_mod_in <- readLines(con=stdin(),1)
    if(!file.exists(dir_netcdf_mod_in)) stop("This directory does not exist.")
    dir_netcdf_mod <- dirname(suppressWarnings(system(paste0("find . -name=netcdf.mod"), intern = T)))
  }else{
    warning("The automatic location of netcdf.mod file uses the locate function in terminal. It must be updated.")
    dir_netcdf_mod <- dirname(suppressWarnings(system(paste0("locate netcdf.mod"), intern = T)))
  }
  lines<-readLines(paste0(package_dir,"src/Makevars"), file.info(paste0(package_dir,"src/Makevars"))$size)
  if(length(dir_netcdf_mod) > 0){
    lines[1]<-paste0("MY_PKG_LIBS= -I",dir_netcdf_mod," -L",dirname(dir_netcdf_mod),"/lib -lnetcdff")
    file.copy(paste0(package_dir,"inst/extdata/src_netcdf/4m_mst_dumping.f90"),paste0(package_dir,"src/"))
    file.copy(paste0(package_dir,"inst/extdata/src_netcdf/main_clu_adjl_mst_dumping.f90"), paste0(package_dir,"src/"))
    file.copy(paste0(package_dir,"inst/extdata/src_netcdf/utilities_netcdf.f90"), paste0(package_dir,"src/"))
    cat("Netcdf support available and activated (lib directory: ", dir_netcdf_mod, ")\n")
  }
  if(!(length(dir_netcdf_mod) > 0)||no_ncdf){
    cat("Warning: netcdf.mod not found or decided not to use netcdf backend data handling even if it is highly advised.\n
        The code has been cleaned up.")
    lines[1]<-paste0("MY_PKG_LIBS= ")
    if(exists(paste0(package_dir,"src/4m_mst_dumping.f90"))) file.remove(paste0(package_dir,"src/4m_mst_dumping.f90"))
    if(exists(paste0(package_dir,"src/main_clu_adjl_mst_dumping.f90"))) file.remove(paste0(package_dir,"src/main_clu_adjl_mst_dumping.f90"))
    if(exists(paste0(package_dir,"src/utilities_netcdf.f90"))) file.remove(paste0(package_dir,"src/utilities_netcdf.f90"))
  }
  writeLines(lines,paste0(package_dir,"src/Makevars"))
}