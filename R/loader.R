  #' @title Load trajectory
#' @description
#'      This is a wonderful description(X)
#'
#' @param fl file location
#' @param n_fold Number of links contractions (folds)
#' Default: \code{0}
#'
#' @examples
#' tree <- adjl_from_pi(fil = "data/PROGIDX_000000000001.dat")
#'
#' @return tree: degree list, connectivity matrix and weights
#'
#' @export load_trj_dcd
#' @import bio3d
load_trj_dcd<-function(t_file){
 return(read.dcd(trjfile = t_file)) 
}

#' @title Build the network from the trajectory file
#' @description
#'      This is a wonderful description(X)
#'
#' @param fl file location
#' @param n_fold Number of links contractions (folds)
#' Default: \code{0}
#'
#' @examples
#' tree <- adjl_from_trj(fil = "data/PROGIDX_000000000001.dat")
#'
#' @return tree: degree list, connectivity matrix and weights
#'
#' @export adjl_from_trj
#' @import parallel bio3d

adjl_from_trj<-function(trj, clu_radius=NULL, clu_hardcut=NULL, cores=NULL, mode=NULL, min_span_tree=TRUE){ # to add functional of distance
  if(!is.matrix(trj)){
    if(!is.data.frame(trj)) stop('trj input must be a matrix or a data.frame')
    trj <- as.matrix(trj)
  }
  #memory handling
  data_management <- getOption("campackage.data_management")
  if(data_management == "R") cat("Normal memory handling selected. Without hdf5 backend file management it will be difficult for R to handle big data-sets.\n")
  else if(data_management == "h5pfc") cat("Selected data support: hdf5 with mpi support\n")
  else if(data_management == "h5fc") cat("Selected data support: hdf5 without mpi support\n")
  else stop("Invalid data management keyword inserted. Check the available methods on the guide.\n")
  
  cat("To set new data_management method: options(list(campackage.data_management = 'R'))\n")

  if(data_management != "R"){
    warning(paste0('The dumping filename will be ',getOption("campackage.data_filename"),'. If already existent it will be overwritten'))
    cat("Checking for hdf5 support...")
    command_loc <- system(paste0("which ",data_management))
    if(command_loc=="") stop("No support for hdf5. Please check installation and correct linkage of the command to your enviroment.")
  }
  
  #Input setting
  r1 <- nrow(trj)
  c1 <- ncol(trj)
  
  #normal -LONG- mode (parallel)
  if(is.null(mode)){
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
    rm(trj)
    
    
  }else if(is.character(mode)&&(mode=="fortran")){
    #Fortran mode
    #default vars
    distance_method <- 5
    if(is.null(clu_radius)){
      clu_radius <- 50
      warning(paste("clu_radius variable (a priori fixed clustering radius) has not been selected. 
              A standard value of",clu_radius,"will be used."))
    }
    if(is.null(clu_hardcut)){
      clu_hardcut <- 150
      warning(paste("clu_hardcut variable (a priori fixed distance threshold between different cluster members) has not been selected. 
              A standard value of",clu_hardcut,"will be used."))
    }
    if(!is.numeric(clu_radius)||!is.numeric(clu_hardcut)) stop("clu_radius and clu_hardcut must be numeric")
    
    if(data_management == "R"){
      #input-output initialization
      output_fin <- list()
      max_d <- 0
      adj_deg <- as.integer(rep(0,r1))
      adj_ix <- matrix(as.integer(rep(0,r1*r1)),r1,r1)
      adj_dis <- matrix(as.single(rep(0.0,r1*r1)),r1,r1)
      trj <- matrix(as.single(trj),ncol = c1,nrow = r1)
      #double Cstyle deginitions
      attr(trj,"Csingle") <- TRUE
      attr(adj_dis,"Csingle") <- TRUE
      #main fortran talker
      output<-.Fortran("generate_neighbour_list", PACKAGE="campackage",
                      trj_data=trj,
                      n_xyz_in=as.integer(c1),
                      n_snaps_in=as.integer(r1),
                      clu_radius_in=as.single(clu_radius),
                      clu_hardcut_in=as.single(clu_hardcut), 
                      adjl_deg=adj_deg, 
                      adjl_ix=adj_ix, 
                      adjl_dis=adj_dis, 
                      max_degr=as.integer(max_d),
                      dis_method_in=as.integer(distance_method),
                      mst_log_in=as.logical(min_span_tree),
                      data_meth_in=as.integer(1),
                      verbose_in=as.logical(TRUE))
      #output adjustment
      output_fin[[1]] <- output$adjl_deg
      output_fin[[2]] <- output$adjl_ix[,1:output$max_degr]
      output_fin[[3]] <- output$adjl_dis[,1:output$max_degr]
      return(output_fin)
    }else if(data_management=="h5pfc"){
      stop("still to-do")
    }else if(data_management=="h5fc"){
      stop("still to-do")
    }else{
      stop("Data_management assigned to an unknown value. Please complain to the developers.")
    }
  }else{
      stop("Mode entry not correct.")
    }
}

#' @title Build the network from the already processed Progress Index file
#' @description
#'      This is a wonderful description(X)
#'
#' @param fl file location
#' @param n_fold Number of links contractions (folds)
#' Default: \code{0}
#'
#' @examples
#' tree <- adjl_from_pi(fil = "data/PROGIDX_000000000001.dat")
#'
#' @return tree: degree list, connectivity matrix and weights
#'
#' @export adjl_from_pi
#' @useDynLib campackage

adjl_from_pi<-function(fil){
  # extract the SST or MST from the output of the analysis already made with campari.
  # Here we will reconstruct a bit of the tree in order to be able to find again the MST/SST
  piOut<-read.table(file = fil)

  # number of snapshots
  nsnaps <- nrow(piOut)
  nbl<-piOut[,6] #number list
  itl<-piOut[,3] #index teo? list
  ditl<-piOut[,5] #distance of 6 from 3
  #max number of connections
  maxnb <- max(hist(breaks=seq(from=0.5,by=1,to=nsnaps+0.5),x=nbl,plot=FALSE)$counts) + 1
  # empty vector of future number of connections for each node in column 6
  treennb <- array(as.integer(0),c(nsnaps))
  # adjlist matrix with obvious limit in maxnb of connections
  treenbl <- array(as.integer(0),c(nsnaps,maxnb))
  # distance of each connection
  treedis <- array(as.single(0.0),c(nsnaps,maxnb))
  # slow
  for (i in 2:nsnaps) {
    treennb[nbl[i]] <- treennb[nbl[i]] + 1
    treenbl[nbl[i],treennb[nbl[i]]] <- itl[i]
    treedis[nbl[i],treennb[nbl[i]]] <- ditl[i]
    treennb[itl[i]] <- treennb[itl[i]] + 1
    treenbl[itl[i],treennb[itl[i]]] <- nbl[i]
    treedis[itl[i],treennb[itl[i]]] <- ditl[i]
    # list of breaks (must be passed at least of size 1, but n_breaks can be zero) eh?
  }
  return(list(treennb,treenbl,treedis))
}

#' @title Transform from an adjmatrix to an adjlist of variables
#' @description
#'      This is a wonderful description(X)
#'
#' @param adj_mat Matrix 
#' @param n_fold Number of links contractions (folds)
#' Default: \code{0}
#'
#' @examples
#' tree <- adjl_from_adjmat(adj_m)
#'
#' @return tree: degree list, connectivity matrix and weights
#'
#' @export adjl_from_adjmat
#' @useDynLib campackage

adjl_from_adjmat<-function(adj_m){
  # extract the SST or MST from the output of the analysis already made with campari.
  # Here we will reconstruct a bit of the tree in order to be able to find again the MST/SST
  adjl_nmbrs<-c()
  adjl<-list()
  adjl_dis<-list()
  for (i in 1:nrow(adj_m)){
    tmp <- sort(adj_m[i,adj_m[i,]!=0], index.return = TRUE)
    adjl_nmbrs[i] <- length(tmp$ix)
    adjl[[i]] <- tmp$ix
    adjl_dis[[i]] <- tmp$x 
  }
  adjl <- array(unlist(adjl),dim = c(length(adjl_nmbrs),max(adjl_nmbrs)))
  adjl_dis <- array(unlist(adjl_dis),dim = c(length(adjl_nmbrs),max(adjl_nmbrs)))
  return(list(array(adjl_nmbrs), adjl, adjl_dis))
}