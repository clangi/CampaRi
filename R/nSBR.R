#' @title New identifying basins in the SAPPHIRE plot
#' @description
#'     Todo.
#'     
#' @param data Name of the PROGIDX_<...>.dat file or Data frame with three columns containing, in order, the (sorted) progress index, the relative time indices and the cut function. 
#' @param nx Number of bins on the x-axis of the 2-D histogram.
#' @param ny Number of bins on the y-axis of the 2-D histogram. Default to nx.
#' @param ny.aut Logical value indicating whether a suitable number of bins on the y-axis has to be identified automatically. Default value is \code{FALSE}.
#' @param local.cut Logical value indicating whether the localized cut function has to be used (see references). Default value is \code{FALSE}.
#' @param plot A logical value indicating whether to display the SAPPHIRE plot with the resulting partitions or not. Black partitions 
#' are the matched ones, blue ones derive only from the dynamic analysis and orange ones only from the kinetic analysis. The green curve is the kinetic annotation (black curve) 
#' where the parabolic shape has been subtracted, i.e. the actual curve used for the peaks identification. Default value is \code{FALSE}.
#' @param silent A logical value indicating whether the function has to remain silent or not. Default value is \code{FALSE}.
#' @param ...
#'      \itemize{
#'          \item "\code{time.series}" File name. If specified, it substitutes the time series of the PROGIDX_<..> file with the one provided by the file.
#'          \item "\code{cl.stat}" If \code{TRUE} it will activate the cluster analysis with Mutual Information, Hellinger distance and Shannon Entropy
#'          as defaults.
#'          \item "\code{cl.stat.weight.barriers}" If \code{TRUE} (standard once cl.stat is active) it will calculate two different Mutual Information: 
#'          one based on the clusters founds and one based on uniform subsampling (done using \code{cl.stat.nBreaks} or the number of cluster founds using 
#'          SBR multiplied for 10). The second one will be inverted to find the discontinuities and the barriers will be weighted simply multiplying the two values 
#'          found on the specific progress index data point.  
#'          \item "\code{cl.stat.nBreaks}" Integer. This variable defines the number of splits for the analysis. If set to 0 it will use the breaks found in the SBR.
#'          \item "\code{cl.stat.denat}" This value can be set to \code{"process_subtraction"} or \code{"poly_interpolation"} and it is defining the removal of 
#'          parabolic artifacts in the kinetic trace. The polynomial fit is by default 7 and 12 in degree for the kinetic annotation and uniform MI curves. The process 
#'          option istead is referring to the simulated process which is behind the kin ann. This option can be set to \code{TRUE} for using the process way.
#'          \item "\code{cl.stat.denat.MI}" Integer. If NULL its default is 7 for poly_interpolation of the kinetic ann but nothing is done on the MI score. 
#'          \item "\code{cl.stat.MI_comb}" This value is the representative of how the various calculation in barrier weighting are combined. In particular the options 
#'          are the following:
#'          \itemize{
#'            \item "\code{mean}" Mean between the MI based on the SBR and the MI based on uniform subdivisions
#'            \item "\code{multip}" Multiplication between the MI based on the SBR and the MI based on uniform subdivisions
#'            \item "\code{kin_MI}" Mean between the kinetic annotaton and the MI based on the SBR and the MI based on uniform subdivisions
#'            \item "\code{kin}" Only kinetic annotation value on the barrier
#'            \item "\code{MI}" Only the MI based on the SBR and the MI based on uniform subdivisions
#'          }
#'          \item "\code{cl.stat.nUni}" Integer. This variable defines the number of splits for the uniform sampling. If not set the algorithm will use 40 divisions which 
#'          could not be optimal for clogged data. Generally speaking, this value is useful for the MI - ratio weighting of the barriers. It is also possible to insert a 
#'          vector of integer and the resulting value is an average of these expanded results. We advice to use \code{c(5,10,15,20,25,30,40,50)} generally.
#'          \item "\code{plot.cl.stat}" Logical for plotting the statistics of the SBR basins using various colors for the methods along with 
#'          the temporal annotation (points in a progress index / time plot). Time has been normalized between 0 and 1 as all the statistics.
#'          \item "\code{cl.stat.entropy}" Logical. Calculate Shannon entropy.
#'          \item "\code{cl.stat.stft}" Logical. Calculate short time fourier transform.
#'          \item "\code{cl.stat.TE}" Logical. Calculate symmetric Transfer Entropy.
#'          \item "\code{cl.stat.KL}" Logical. Calculate symmetric and non-symmetric Kullback-Leibler divergence.
#'          \item "\code{cl.stat.wMI}" Logical. Use the number of breaks (\code{cl.stat.nBreaks}) to find the Mutual Information for 10 slided divisions.
#'          This will result in 10*nBreaks values.
#'          \item "\code{data.out.it}" Should I output the data.frame in input?
#'      }
#'      
#' @return A list containing
#'       \itemize{
#'         \item "\code{tab.st}" Data frame containing the boundaries of each state, their lengths and the type of the right boundary. If the cl.stat option
#'         is active this will have other columns with the barrier statistics (e.g. Hellinger distance).
#'         Type=1 means that it is a matched partition, type=2 is only dynamic, type=3 is only kinetic. The type of the last state is always equal to 1.
#'         \item "\code{nbins}" 2-D vector containing number of bins on x-axis and y-axis.
#'         \item "\code{seq.st}" The time-ordered discretized trajectory.
#'         \item "\code{statistics}" If \code{cl.stat} is \code{TRUE} this element will contain all the cluster statistics (unbound to the found barriers). 
#'         Otherwise it is \code{NULL}. If no barrier has been found this value defaults to \code{FALSE}.
#'         \item "\code{call}" The matched call. 
#'         \item "\code{data}" The sapphire file name if it is provided, otherwise the data obj. It is used by \code{\link{score_sapphire}}.
#'       }   
#' @examples
#' 
#' adjl <- mst_from_trj(trj = matrix(rnorm(1000), nrow = 100, ncol = 10))
#' ret<-gen_progindex(adjl = adjl)
#' gen_annotation(ret_data = ret, local_cut_width = 10)
#' \dontrun{
#' basins_recognition("PROGIDX_000000000001.dat", nx=500)
#' }
#' 
#' @details For details regarding the SAPPHIRE plot, please refer to the relative publications \url{http://www.nature.com/articles/srep06264}. 
#' Main documentation of the original campari software \url{http://campari.sourceforge.net/documentation.html}.
#' 
#' @importFrom data.table fread fwrite
#' @importFrom RColorBrewer brewer.pal
#' @importFrom TransferEntropy computeTE
#' @importFrom infotheo mutinformation
#' @importFrom e1071 stft
# @importFrom RcppRoll roll_mean
#' 
#' @export nSBR

nSBR <- function(data, nx, ny=nx, ny.aut=FALSE, local.cut=FALSE, plot=FALSE, silent=FALSE, ...) {
  call <- match.call()
  
  if(!is.character(data) && !is.data.frame(data)) stop("data must be a string or a data frame")
  if(is.character(data) && (!all(grepl("PROGIDX", data)) && !all(grepl("REPIX", data)))) stop("Please provide a data name starting with 'PROGIDX' or 'REPIX'" )
  if((nx %% 1) != 0) stop("nx must be an integer")
  if((ny %% 1) != 0) stop("ny must be an integer")
  if((nx < 7 || ny < 7)) stop('nx and ny must be > 7')
  if(!is.logical(local.cut)) stop("local.cut must be a logical value")
  if(!is.logical(plot)) stop("plot must be a logical value")
  
  # Extra arguments checks
  input.args <- list(...)
  avail.extra.arg <- c("pol.degree", "time.series",
                       "cl.stat", "plot.cl.stat", 'cl.stat.entropy', 'cl.stat.weight.barriers',
                       'cl.stat.stft', 'cl.stat.TE', 'cl.stat.KL', 'cl.stat.wMI',
                       'cl.stat.nUni', 'cl.stat.nBreaks', 'cl.stat.MI_comb', 'cl.stat.denat', 'cl.stat.denat.MI', 
                       'dbg_basins_recognition', 'data.out.it')
  
  if(!is.null(names(input.args)) && any(!(names(input.args) %in% avail.extra.arg))){
    warning('There is a probable mispelling in one of the inserted variables. Please check the available extra input arguments.')
    if(!silent) cat('!!!!!!!!! We found the following variables without a father (not between our extra input arguments) !!!!!!!!!!\n')
    if(!silent) cat(names(input.args)[!(names(input.args) %in% avail.extra.arg)], '\n')
  }
  
  # dbg_basins_recognition for stats - dgarol
  if("dbg_basins_recognition" %in% names(input.args)) { # dgarol
    dbg_basins_recognition <- input.args$dbg_basins_recognition
    stopifnot(is.logical(dbg_basins_recognition))
  } else dbg_basins_recognition <- FALSE
  
  
  if("data.out.it" %in% names(input.args)) { # dgarol
    data.out.it <- input.args$data.out.it
    stopifnot(is.logical(data.out.it))
  } else data.out.it <- FALSE
  
  # cluster statistics check - dgarol
  
  if("cl.stat" %in% names(input.args)) { # dgarol
    cl.stat <- input.args$cl.stat
    cl.stat.weight.barriers <- TRUE
    stopifnot(is.logical(cl.stat))
  } else cl.stat <- FALSE
  
  # Checking the weight of the barriers
  if("cl.stat.weight.barriers" %in% names(input.args)) { # dgarol
    cl.stat.weight.barriers <- input.args$cl.stat.weight.barriers
    stopifnot(is.logical(cl.stat.weight.barriers))
    if(!cl.stat && cl.stat.weight.barriers) cl.stat <- TRUE
  } else if(!cl.stat) cl.stat.weight.barriers <- FALSE
  
  # Checking the number of breaks for the statistics (0 is the standard barriers)
  if("cl.stat.nBreaks" %in% names(input.args)) { # dgarol
    cl.stat.nBreaks <- input.args$cl.stat.nBreaks
    stopifnot(.isSingleInteger(cl.stat.nBreaks))
    if(!cl.stat) cl.stat <- TRUE
  } else cl.stat.nBreaks <- 0
  
  # checking the uniform sampling nsplits
  if("cl.stat.denat" %in% names(input.args)) { # dgarol
    cl.stat.denat <- input.args$cl.stat.denat
    denat_opt.ava <- c("process_subtraction", "poly_interpolation")
    if(is.logical(cl.stat.denat) && cl.stat.denat){
      cl.stat.denat <- denat_opt.ava[1]
    }else if(is.logical(cl.stat.denat) && !cl.stat.denat) cl.stat.denat <- NULL 
    if(!is.null(cl.stat.denat)){
      stopifnot(is.character(cl.stat.denat))
      if(!(cl.stat.denat[1] %in% denat_opt.ava)) 
        stop("cl.stat.denat method option not valid.")
      if(!cl.stat) {
        cl.stat <- TRUE
        cl.stat.weight.barriers <- TRUE
      }
    }
  } else cl.stat.denat <- NULL
  if("cl.stat.denat.MI" %in% names(input.args)) { # dgarol
    cl.stat.denat.MI <- input.args$cl.stat.denat.MI
    if(!is.null(cl.stat.denat.MI)){
      if(is.logical(cl.stat.denat.MI) && cl.stat.denat.MI){
        cl.stat.denat.MI <- 7
      }else if(is.logical(cl.stat.denat.MI) && !cl.stat.denat.MI) cl.stat.denat.MI <- NULL 
      stopifnot(.isSingleInteger(cl.stat.denat.MI))
      if(cl.stat.denat.MI < -1 || cl.stat.denat.MI == 0 || cl.stat.denat.MI > 50) stop('cl.stat.denat.MI must be between 1 and 50.')
      if(!cl.stat) cl.stat <- TRUE
      if(!cl.stat.weight.barriers) cl.stat.weight.barriers <- TRUE
      if(is.null(cl.stat.denat)) cl.stat.denat <- 'poly_interpolation'
    }
  } else cl.stat.denat.MI <- NULL
  
  if("cl.stat.MI_comb" %in% names(input.args)) { # dgarol
    cl.stat.MI_comb <- input.args$cl.stat.MI_comb
    stopifnot(is.character(cl.stat.MI_comb))
    stopifnot(cl.stat.MI_comb %in% c('MI', 'kin', 'mean', 'kin_MI', 'multip'))
    if(!cl.stat) {
      cl.stat <- TRUE
      cl.stat.weight.barriers <- TRUE
    }
  } else cl.stat.MI_comb <- 'kin_MI'
  
  if("cl.stat.nUni" %in% names(input.args)) { # dgarol
    cl.stat.nUni <- input.args$cl.stat.nUni
    if(!is.null(cl.stat.nUni)){
      stopifnot(all(sapply(cl.stat.nUni, function(x) x%%1) == 0))
      if(!cl.stat) {
        cl.stat <- TRUE
        cl.stat.weight.barriers <- TRUE
      }
    }
  } else cl.stat.nUni <- NULL
  
  
  # Specific cluster statistics
  if("cl.stat.wMI" %in% names(input.args)) { # dgarol
    cl.stat.wMI <- input.args$cl.stat.wMI
    stopifnot(is.logical(cl.stat.wMI))
    if(!cl.stat && cl.stat.wMI) cl.stat <- TRUE
  } else cl.stat.wMI <- FALSE      # can be long
  if("cl.stat.KL" %in% names(input.args)) { # dgarol
    cl.stat.KL <- input.args$cl.stat.KL
    stopifnot(is.logical(cl.stat.KL))
    if(!cl.stat && cl.stat.KL) cl.stat <- TRUE
  } else cl.stat.KL <- FALSE
  if("cl.stat.TE" %in% names(input.args)) { # dgarol
    cl.stat.TE <- input.args$cl.stat.TE
    stopifnot(is.logical(cl.stat.TE))
    if(!cl.stat && cl.stat.TE) cl.stat <- TRUE
  } else cl.stat.TE <- FALSE
  if("cl.stat.stft" %in% names(input.args)) { # dgarol
    cl.stat.stft <- input.args$cl.stat.stft
    stopifnot(is.logical(cl.stat.stft))
    if(!cl.stat && cl.stat.stft) cl.stat <- TRUE
  } else cl.stat.stft <- FALSE     # it is plotting too much
  
  # Entropy deeper exploration
  if("cl.stat.entropy" %in% names(input.args)) { # dgarol
    cl.stat.entropy <- input.args$cl.stat.entropy
    stopifnot(is.logical(cl.stat.entropy))
    if(!cl.stat && cl.stat.entropy) cl.stat <- TRUE
  } else cl.stat.entropy <- FALSE
  
  # Plotting
  if("plot.cl.stat" %in% names(input.args)) { # dgarol
    plot.cl.stat <- input.args$plot.cl.stat
    stopifnot(is.logical(plot.cl.stat))
    if(!cl.stat && plot.cl.stat) cl.stat <- TRUE
  } else plot.cl.stat <- FALSE
  
  # --------------------------- functions
  scale <- function(x) {
    scl <- (max(progind$Time)-1)/(max(x)-min(x))
    return(1+(x-min(x))*scl)
  }
  # refine <- function(vec) {
  #   ## Attach to the peak in the cut function
  #   for (jj in 1:length(vec)) {
  #     vec <- round(vec)
  #     tmp <- which.max(kin[(vec[jj]-round(sd.kin/2)):(vec[jj]+round(sd.kin/2))]) 
  #     vec[jj] <- tmp+vec[jj]-round(sd.kin/2)-1
  #   }
  #   vec <- unique(c(1,vec,cstored))
  #   if(length(which(vec<=sd.kin))>1) vec <- vec[-c(2:max(which(vec<sd.kin)))]
  #   if(length(which((cstored-vec)<=sd.kin)>1)) vec <- c(vec[-which((cstored-vec)<sd.kin)],cstored)
  #   while (any(diff(vec)<sd.kin)) {
  #     vec.tmp <- vec
  #     dst <- 0
  #     start <- NULL
  #     for(ii in 2:(length(vec)-1)){
  #       dst <- dst+vec[ii+1]-vec[ii]
  #       if(dst<sd.kin) {
  #         start <- c(start,ii)
  #         next
  #       } else if (length(start)==1) {
  #         start <- c(start,ii)
  #         ll <- start[which.min(kin[vec[start]])]
  #         vec.tmp <- vec.tmp[-match(vec[ll],vec.tmp)]
  #       } else if (length(start)>1) {
  #         start <- c(start,ii)
  #         ll <- which.max(kin[vec[start]])
  #         lst <- start[-ll]
  #         vec.tmp <- vec.tmp[-match(vec[lst], vec.tmp)]
  #       }
  #       start <- NULL
  #       dst <- 0
  #     }
  #     vec <- vec.tmp
  #   }
  #   vec <- vec[which(vec != 1)]
  #   vec <- vec[which(vec != cstored)]
  #   if(length(vec)<3) warning("Vector completely erased in refine function")
  #   return(vec)
  # }
  # true.peaks <- function(idxx, values, up=TRUE) {
  #   ## To select the effective peaks or minima, it return a selection of idxx 
  #   rl <- rle(diff(idxx))
  #   if(any(rl$values==1)) {
  #     disc <- NULL
  #     st <- c(0, cumsum(rl$lengths))[which(rl$values==1)] + 1
  #     end <- cumsum(rl$lengths)[which(rl$values==1)]
  #     for(ii in seq_along(st)) {
  #       ## max() means that in case two consecutive index have the same values, I pick the one on the right (.. don't know if it's the best choice honestly)
  #       if(up) idx.max <- max(which(values[idxx[st[ii]:(end[ii]+1)]]==max(values[idxx[st[ii]:(end[ii]+1)]])))
  #       else idx.max <- max(which(values[idxx[st[ii]:(end[ii]+1)]]==min(values[idxx[st[ii]:(end[ii]+1)]])))
  #       disc <- c(disc, idxx[st[ii]:(end[ii]+1)][-idx.max])
  #     }
  #     return(idxx[-match(disc, idxx)])
  #   } else return(idxx)
  # }
  # myHell <- function(x, y){
  #   return(sqrt(1-sum(sqrt(x*y))))
  # }
  # myKL <- function(x, y){ # It is defined afterwards with a rough symmetrization - dgarol
  #   xn <- x[-which(x==0 | y==0)]
  #   yn <- y[-which(x==0 | y==0)]
  #   return(sum(xn*log(xn/yn)))
  # }
  # ------------------------ end functions
  
  ## INPUT FILE 
  if(!silent) cat("Reading PROGIDX file...\n")
  if(is.data.frame(data)){
    if(ncol(data) != 20) stop('Column number must be 20 as it is the default format.')
    if(!local.cut) {
      progind <- data.frame(data[, c(1, 3, 4)])
      colnames(progind) <- c("PI", "Time", "Cut")
    } else {
      foo <- data.frame(data[, c(1, 3, 10, 12)])
      progind <- data.frame(PI=foo[[1]], Time=foo[[2]], Cut=(foo[[3]]+foo[[4]])/2)
      rm(foo)
    }
    if(!data.out.it) {
      data.out <- NULL # was NULL but for my cases I need it always evan with silent
      warning('Putting the silent mode, no input data will be outputted!')
    }else{
      data.out <- data
    }
  } else {
    if(!local.cut) {
      progind <- data.frame(fread(data, showProgress=FALSE)[, c(1, 3, 4)])
      colnames(progind) <- c("PI", "Time", "Cut")
    } else {
      foo <- data.frame(fread(data, showProgress=FALSE)[, c(1, 3, 10, 12)])
      progind <- data.frame(PI=foo[[1]], Time=foo[[2]], Cut=(foo[[3]]+foo[[4]])/2)
      rm(foo)
    }
    if("time.series" %in% names(input.args)){
      time.series <- input.args$time.series
      progind$Time <- as.vector(unlist(fread(time.series, showProgress=FALSE)[,1]))
      if(!silent) cat("Time series given by", time.series, "\n")
    }
    data.out <- data
  }
  if(nrow(progind) < 80) stop('It is impossible to recognize basins with less than 80 snapshots in the sapphire table.')
  cstored <- dim(progind)[1]
  
  ## Dynamic parameters
  # perc <- 0.90 ## Percentage to find optimal nbiny
  # xidx <- 0.5 ## Joining consecutive cells on a single raw
  # K <- 2 ## Parameter of the Haat wavelet
  # dpeaks.dyn <- 3 ## Parameter of the maximum search algorithm (window size)
  # nsample <- 50 ## Number of samples in the distribution of reshuffled Hellinger distances 
  # cutjoin <- 50 ## Number of null joining attempts required to quit the procedure 
  # conf.lev <- 0.005 ## Cut on pvalues of the Grubb Test
  # ## Kinetic parameters
  # wsize <- round(2*cstored/nx) + ((round(2*cstored/nx)+1) %% 2)  ## To make it odd
  # dpeaks.kin <- ceiling(wsize/2)+ceiling(wsize/2)%%2+1 
  # thr.ratio <- 0.05 ## Parameter of second data cleaning  
  # ## Matching parameters
  # lx <- round(cstored/nx)
  # sd.kin <- 1*lx
  # sd.dyn <- 1*lx
  
  
  #######################################################################
  #### final calculations for scores - dgarol
  #######################################################################
  if(dbg_basins_recognition) browser()
  # functions
  # ----------------------------------------------------
  myTE <- function(x, y, emb, sym = T, sd = 0.001){
    nx <- x + stats::rnorm(n = .lt(x), mean = 0, sd = sd)
    ny <- y + stats::rnorm(n = .lt(x), mean = 0, sd = sd)
    xy <- TransferEntropy::computeTE(nx, ny, embedding = emb, k = 2, safetyCheck = T)
    yx <- TransferEntropy::computeTE(ny, nx, embedding = emb, k = 2, safetyCheck = T)
    if(sym) return((unlist(xy) + unlist(yx))/2)
    else return(unlist(xy))
  }
  myShEn <- function(x){
    x2 <- replace(x = x, list = which(x == 0), values = 1)
    return(-sum(x2*log(x2)))
  }
  
  myKL <- function(x, y, sym = T){
    xn <- x[-which(x==0 | y==0)]
    yn <- y[-which(x==0 | y==0)]
    if(sym) return((sum(xn*log(xn/yn)) + sum(yn*log(yn/xn)))/2)
    else return(sum(xn*log(xn/yn)))
  }
  
  # creates the densities and the counts over bins (only ny is relevant) - doing hist2d
  dens_histCounts <- function(progind, breaks, nx, ny = nx){
    hist.internal <- hist2d(matrix(c(progind[,1], progind[,2]), ncol=2, nrow=.lt(progind$PI)), nbins=c(nx, ny), show=FALSE)
    cnts <- matrix(0, nrow=ny, ncol=(.lt(breaks) + 1))
    for (i in 1:ncol(cnts)) {
      if (i==1) {
        ncls <- 1
        ncle <- which.min(abs(hist.internal$x - breaks[i]))
      } else if (i == (.lt(breaks) + 1) ) {
        ncls <- which.min(abs(hist.internal$x - breaks[i-1])) 
        ncle <- nx
      } else {
        ncls <- which.min(abs(hist.internal$x - breaks[i-1]))
        ncle <- which.min(abs(hist.internal$x - breaks[i]))
      }
      ## ColSums doesn't work if ncls==ncle
      for (j in 1:ny) cnts[j, i] <- sum(hist.internal$counts[c(ncls:ncle), j])
    }
    dens <- apply(cnts, 2, function(x) x/sum(x)) # this can break if you change the number of bins
    # dens <- apply(cnts, 2, function(x) { # probably useless if we put the histogram in the function!
    #   if(sum(x) != 0) return(x/sum(x))
    #   else return(rep(0, .lt(x)))
    #   })
    if(nrow(dens) != ny) dens <- t(dens) # safety check
    
    return(list('density' = dens, 'counts' = cnts))
  }
  
  # once defined a split (e.g. 50 parts) it uses a certain number of slides (e.g. 10) to calculate the MI
  sliding_MI <- function(progind, n_breaks = 50, n_slides = 10, nx, ny = nx){
    lpi <- .lt(progind$PI)
    brks <- floor(seq(1, lpi, length.out = n_breaks))
    span_MI <- sapply(X = floor(seq(1, brks[2]-1, length.out = n_slides)), FUN = function(x){
      br <- floor(seq(x, lpi, by = brks[2]))
      dhc <- dens_histCounts(progind, breaks = br, nx = nx, ny = ny)
      mi.out <- list(.normalize(sapply(seq(.lt(br)), function(j) infotheo::mutinformation(dhc$counts[, j], dhc$counts[, j+1]))))
      return(mi.out)
    })
    return(c(t(do.call(cbind, span_MI))))
  }
  # ----------------------------------------------------
  
  # color exploration
  # require(RColorBrewer)
  # RColorBrewer::display.brewer.all()
  # ncolor = 7; pie(x = rep(1, ncolor), labels = 1:ncolor, col = RColorBrewer::brewer.pal(n = ncolor, name = 'Set1'))
  # ncolor = 7; pie(x = rep(1, ncolor), labels = 1:ncolor, col = RColorBrewer::brewer.pal(n = ncolor, name = 'Dark2'))
  
  # using the division to calculate the densities and the counts
  lpi <- .lt(progind$PI)
  stopifnot(cl.stat.nBreaks >= 0, cl.stat.nBreaks <= floor(lpi/2))
  
  # calculating the weights for the barriers 
  if(cl.stat.weight.barriers) {
    if(is.null(cl.stat.nUni)) n_unif <- 40
    else n_unif <- cl.stat.nUni
    stopifnot(n_unif > 0, n_unif < lpi/ 2)
    
    # SBR based mutual information
    dhc_sbr <- dens_histCounts(progind, breaks = breaks, nx = nx, ny = ny) # breaks are the barriers points on x
    sbr_cnts <- dhc_sbr$counts    
    MI_sbr <- sapply(seq(.lt(breaks)), function(j) infotheo::mutinformation(sbr_cnts[, j], sbr_cnts[, j+1]))
    MI_sbr <- .normalize(MI_sbr)
    
    # Uniform division based mutual information
    # n_unif <- round(lpi/50) 
    # n_unif <- 10
    # n_unif <- c(5,10,15,20,25,30,40,50,60)
    # n_unif <- seq(5, 150, 5)
    if(!silent) cat('Calculating the barrier statistics using ') 
    expand_MI_uni <- array(0, dim = lpi)
    for(nun in n_unif){
      if(!silent) cat(nun, ' ')
      brks_uni <- floor(seq(1, lpi, length.out = nun))[-c(1,nun)]
      dhc_uni <- dens_histCounts(progind, breaks = brks_uni, nx = nx, ny = ny) # breaks are the barriers points on x
      uni_cnts <- dhc_uni$counts    
      MI_uni <- sapply(seq(nun-2), function(j) infotheo::mutinformation(uni_cnts[, j], uni_cnts[, j+1]))
      MI_uni <- .normalize(MI_uni)
      expand_MI_uni <- expand_MI_uni + stats::approx(seq(nun), c(0, 1 - MI_uni, 0), n = lpi)$y
    }
    if(!silent) cat('divisions for the MI ratio. \n')
    expand_MI_uni <- expand_MI_uni/.lt(n_unif)
    # expand_MI_uni <- .normalize(expand_MI_uni)
    # plot(x = seq(1,80000, length.out = .lt(MI_uni)), y = 1-MI_uni, type = 'l'); lines(expand_MI_uni, col = 'red')
    # MI_sbr <- .normalize(MI_sbr, xmax = max(c(MI_sbr, MI_uni)), xmin = min(c(MI_sbr, MI_uni)))
    # MI_uni <- .normalize(MI_uni, xmax = max(c(MI_sbr, MI_uni)), xmin = min(c(MI_sbr, MI_uni)))
    
    # hard wired options for the basin weights
    MI_comb <- cl.stat.MI_comb
    denat <- cl.stat.denat 
    denat.MI <- cl.stat.denat.MI
    # Calculating the kinetical cut and the uniform MI denaturation if necessary
    # xr1 <- c(margin:(cstored-round(cstored*0.99)))
    kin.pl <- .normalize(-log(cutf / cstored)) # xr1 impose a selection of 99% of the snapshots to avoid the initial inf
    
    # normalizing the results to avoid parabolic effects (it depends on the grade of the poly!)
    if(!is.null(denat) && denat == 'poly_interpolation'){
      if(!is.null(denat.MI)){
        if(denat.MI != -1) kinpl <- .denaturate(yyy = kin.pl, xxx = seq(lpi), polydeg = denat.MI, plotit = F)
        else kinpl <- .denaturate(yyy = kin.pl, xxx = seq(lpi), polydeg = 7, plotit = F)
        if(denat.MI != -1) miuni <- .denaturate(yyy = expand_MI_uni, xxx = seq(lpi), polydeg = denat.MI, plotit = F)
        else miuni <- 1 - expand_MI_uni
        # MI_ratio <- (kin.pl[breaks] + miuni[breaks]) / 2
      }else{
        kinpl <- .denaturate(yyy = kin.pl, xxx = seq(lpi), polydeg = 7, plotit = F)
        miuni <- expand_MI_uni
      }
    }else if(!is.null(denat) && denat == 'process_subtraction'){
      kinpl <- .normalize(kin)
      if(!is.null(denat.MI)) {
        if(denat.MI != -1) miuni <- .denaturate(yyy = expand_MI_uni, xxx = seq(lpi), polydeg = denat.MI, plotit = F)
        else miuni <- 1 - expand_MI_uni
      }else miuni <- expand_MI_uni
    }else{
      kinpl <- kin.pl
      miuni <- expand_MI_uni
    }
    
    # select the final barrier_weight
    if(MI_comb == 'mean') MI_ratio <- (MI_sbr + miuni[breaks]) / 2 
    else if(MI_comb == 'multip') MI_ratio <- MI_sbr * miuni[breaks]
    else if(MI_comb == 'kin_MI') MI_ratio <- (kinpl[breaks] + miuni[breaks]) / 2
    else if(MI_comb == 'kin') MI_ratio <- kinpl[breaks]
    else if(MI_comb == 'MI') MI_ratio <- miuni[breaks]
    if(plot.cl.stat){
      
      # plot the found histogram Mutual Information
      plot(miuni, type = 'l', ylim = c(0,1), xlim = c(0,lpi))                             # from uniform sampling (average of different x-divisions)
      points(breaks, MI_sbr, pch ='o', col = 'green', cex = 2)                                    # from SBR division
      points(breaks, miuni[breaks], pch = 20, col = 'darkblue', cex = 2)                  # height of the point on barriers using uniform sampling
      
      # selecting the barriers and 
      which_barrier <- breaks > lpi/n_unif[.lt(n_unif)] & breaks < lpi - lpi/n_unif[.lt(n_unif)]
      points(breaks[which_barrier], MI_ratio[which_barrier], pch ='+', col = 'red', cex = 2)
      points(progind$PI, progind$Time/lpi, pch='.', cex = 2)
      lines(kinpl, col = 'darkred')
      # if(cl.stat.wMI) lines(seq(1, lpi, length.out = .lt(sl_MI)), sl_MI, col = 'darkgreen')
      
      # what is the relation with the kinetic curve?
      # plot((kin.pl + expand_MI_uni)/2, type = "l")
      # points(breaks, (kin.pl[breaks] + expand_MI_uni[breaks]) / 2, col = 'darkblue')
    }
  }else{
    MI_comb <- cl.stat.MI_comb # should be NULL as the others
    denat <- cl.stat.denat 
    denat.MI <- cl.stat.denat.MI
  }
  
  # sliding window MI
  if(cl.stat.nBreaks == 0) nbr <- .lt(breaks)
  else nbr <- cl.stat.nBreaks
  if(cl.stat.wMI) {
    if(!silent) cat('Calculating the slided MI using 10 sub-division of the uniform divisions. \n')
    sl_MI <- sliding_MI(progind = progind, n_breaks = nbr, n_slides = 10, nx = nx, ny = ny)
  }
  
  if(cl.stat.nBreaks == 0) tobrk <- breaks
  else tobrk <- floor(seq(1, lpi, length.out = cl.stat.nBreaks))
  dhc <- dens_histCounts(progind, breaks = tobrk, nx = nx, ny = ny) # breaks are the barriers points on x
  dens <- dhc$density
  cnts <- dhc$counts
  
  # Analysis of it - rolling means/max/min
  # require(RcppRoll)
  # plot(sl_MI, type = 'l')
  # lines(x = c(1, .lt(sl_MI)), y = rep(mean(sl_MI), 2))
  # lines(RcppRoll::roll_mean(sl_MI, n = 18), col = 'red')
  # lines(RcppRoll::roll_sd(sl_MI, n = 30), col = 'blue')
  # lines(RcppRoll::roll_min(sl_MI, n = 20), col = 'blue')
  # lines(RcppRoll::roll_max(sl_MI, n = 20), col = 'blue')
  # lines(abs(RcppRoll::roll_mean(sl_MI, n = 10) - RcppRoll::roll_mean(sl_MI, n = 60)), col = 'blue')
  
  # short time Fourier Transform
  if(cl.stat.stft){
    if(!silent) cat('Calculating the short time Fourier Transform using a window of 30 and an increment of 1. \n')
    short_time_FT <- e1071::stft(sl_MI, win = 30, inc = 1)
    histOfStft <- apply(short_time_FT$values, 1, sum)
    histOfStft <- .normalize(histOfStft)
    if(plot.cl.stat){
      plot(short_time_FT)
      # plot(apply(short_time_FT$values, 2, sum), type = 'l', xlab = 'Freq')
      plot(histOfStft, type = 'l')
      lines(abs(diff(abs(diff(histOfStft))))^2, type = 'l')
      lines(abs(diff(histOfStft))^2, type = 'l')
    }
    
    # peaks finding
    # peakss <- RcppRoll::roll_mean(abs(diff(histOfStft))^2, 10)
    if(!silent) cat('Finding the peaks using a standard span of 80. \n')
    peakss <- abs(diff(histOfStft))^2
    TFpeakss <- splus2R::peaks(peakss, span = 80)
    if(plot.cl.stat){
      plot(peakss, type = 'l')
      lines(TFpeakss, type = 'l')
    }
    
    # now we try to project it back to the original points (.lt -> 441)
    # .lt(TFpeakss)
    xpeaks <- seq(1, lpi, length.out = .lt(TFpeakss))[TFpeakss]
  }
  
  # good_areas_for_thresholds # deprecated
  # gaft <- 1 - RcppRoll::roll_min(sl_MI, n = 20)
  # plot(x= seq(1, lpi - (.lt(sl_MI) - .lt(gaft))*lpi/.lt(sl_MI), length.out = .lt(gaft)), y=gaft, type = 'l')
  # segments(x0=breaks_save, y0 = 0, y1 = 1)
  
  # calculating the standards and not
  ent <- .normalize(apply(dens, 2, myShEn)) # this is a value PER cluster
  distHell <- sapply(seq(.lt(tobrk)), function(j) myHell(dens[, j], dens[, j+1]))
  distMI <- sapply(seq(.lt(tobrk)), function(j) infotheo::mutinformation(cnts[, j], cnts[, j+1])); distMI <- .normalize(distMI)
  if(cl.stat.KL) { distSymKL <- sapply(seq(.lt(tobrk)), function(j) myKL(dens[, j], dens[, j+1], T)); distSymKL <- .normalize(distSymKL) }
  if(cl.stat.KL) { distKL <- sapply(seq(.lt(tobrk)), function(j) myKL(dens[, j], dens[, j+1], F)); distKL <- .normalize(distKL) }
  if(cl.stat.TE) { distSymTE <- sapply(seq(.lt(tobrk)), function(j) myTE(dens[, j], dens[, j+1], emb = 3, sym = T)); distSymTE <- .normalize(distSymTE) }
  
  # defining the centers/ini/end of the splits
  center_cl <- diff(c(1, tobrk, lpi))/2 + c(1, tobrk)
  ini_points <- c(1, tobrk); end_points <- c(1, tobrk) + diff(c(1, tobrk, lpi))
  
  # Plotting the barriers, breaks and statistics on the breaks 
  if(plot.cl.stat && F){ # I will put it back in a moment
    gg <- ggplot() + theme_classic() + xlab('Progress Index') + ylab('Barrier score') +
      geom_line(aes(tobrk, distHell, col = as.factor(1)), size = 0.5) +
      geom_segment(aes(x = ini_points, xend = end_points, y = ent, yend = ent, col = as.factor(2)), size = 4) +
      geom_line(aes(tobrk, distMI, col = as.factor(3)), size = 2) +
      geom_point(aes(progind$PI, progind$Time/lpi), size = 0.1) + 
      geom_vline(aes(xintercept=c(1, tobrk, lpi)), size = 0.1) + # MI and entropy calculations 
      geom_vline(aes(xintercept=c(1, breaks, lpi)), size = 0.5)  # Found barriers from SBR
    lbls <- c('Hellinger', 'Entropy', 'MI')
    if(cl.stat.KL){
      gg <- gg + geom_line(aes(tobrk, distKL, col = as.factor(4)), size = 0.5) +
        geom_line(aes(tobrk, distSymKL, col = as.factor(5)), size = 0.5)
      lbls <- c(lbls, 'KL', 'symKL')
    }
    if(cl.stat.TE){
      gg <- gg + geom_line(aes(tobrk, distSymTE, col = as.factor(6)), size = 0.5)
      lbls <- c(lbls, 'symTE')
    }
    if(cl.stat.wMI){
      gg <- gg + geom_line(aes(seq(1, lpi, length.out = .lt(sl_MI)), sl_MI, col = as.factor(7)), size = 1) + 
        geom_hline(aes(yintercept = mean(sl_MI)))
      lbls <- c(lbls, 'WinMI')
    }
    if(cl.stat.stft){
      gg <- gg + geom_vline(aes(xintercept=c(xpeaks, lpi), col = as.factor(8)), size = 0.2) # no need of new splits... better to score the bsr split (cohesion - MI) 
      lbls <- c(lbls, 'stFT barriers')
    }
    if(cl.stat.weight.barriers){
      gg <- gg + geom_segment(aes(x = breaks - round(lpi/30), xend = breaks + round(lpi/30), 
                                  y = MI_sbr*expand_MI_uni[breaks], yend = MI_sbr*expand_MI_uni[breaks]), col = 'red', size = 1.5)
    }
    gg <- gg + scale_color_manual(name = "Method", labels = lbls, values = RColorBrewer::brewer.pal(n = 8, name = 'Dark2')) +
      guides(color = guide_legend(override.aes = list(size=5)))
    print(gg)
  }
  
  # Entropy calculations
  if(cl.stat.entropy){
    if(!silent) cat('Calculating various entropy scores. This feature remains for visualization only. \n')
    # ent <- .normalize(apply(dens, 2, myShEn)/apply(dens, 2, function(x) mean(x)))
    ent <- .normalize(apply(dens, 2, myShEn)) # this is a value PER cluster
    # ent1 <- .normalize(apply(dens, 2, myShEn)/diff(c(1, tobrk, lpi))) # the size is not really relevant
    ent1 <- .normalize(apply(dens, 2, myShEn) * apply(cnts, 2, function(x) sum(x)/sum(x != 0)))
    ent2 <- .normalize(apply(dens, 2, myShEn)/apply(dens, 2, function(x) max(x)-min(x))) 
    ent3 <- .normalize(apply(dens, 2, myShEn)/apply(dens, 2, function(x) sd(x)))
    ent4 <- .normalize(apply(cnts, 2, function(x) sd(x)))
    ent5 <- .normalize(apply(dens, 2, function(x) mean(x[which(x!=0)])))
    
    # Check for NAs
    if(anyNA(ent) || anyNA(ent1) || anyNA(ent2) || anyNA(ent3) || anyNA(ent4) || anyNA(ent5))
      stop('We found NAs in the entropy calculations. you should reduce or re consider statistical binning.')
    
    # Plotting entropy calculations
    if(plot.cl.stat){
      gg <- ggplot() + theme_classic() + xlab('Progress Index') + ylab('Barrier score') +
        geom_segment(aes(x = ini_points, xend = end_points, y = ent, yend = ent, col = as.factor(1)), size = 3, alpha = 0.5) +
        geom_segment(aes(x = ini_points, xend = end_points, y = ent1, yend = ent1, col = as.factor(2)), size = 3, alpha = 0.5) +
        geom_segment(aes(x = ini_points, xend = end_points, y = ent2, yend = ent2, col = as.factor(3)), size = 3, alpha = 0.5) +
        geom_segment(aes(x = ini_points, xend = end_points, y = ent3, yend = ent3, col = as.factor(4)), size = 3, alpha = 0.5) +
        geom_segment(aes(x = ini_points, xend = end_points, y = ent4, yend = ent4, col = as.factor(5)), size = 3, alpha = 0.5) +
        geom_segment(aes(x = ini_points, xend = end_points, y = ent5, yend = ent5, col = as.factor(6)), size = 3, alpha = 0.5) +
        geom_point(aes(progind$PI, progind$Time/lpi), size = 0.1) +  
        geom_vline(aes(xintercept=c(1, tobrk, lpi)), size = 0.1) + # MI and entropy calculations 
        geom_vline(aes(xintercept=c(1, breaks, lpi)), size = 0.5) + # Found barriers from SBR
        scale_color_manual(name = "Entropy", 
                           labels = c("Classic", "/size_cl", "/max-min", "/max", "sd(cnts)", "mean(dens!=0)"), 
                           values = RColorBrewer::brewer.pal(n = 6, name = 'Dark2')) +
        guides(color = guide_legend(override.aes = list(size=5)))
      print(gg)    
    }
  }
  
  # Final assignment for output
  statistics <- NULL
  ele <- 1
  if(cl.stat.nBreaks == 0){ # case in which the statistic comes from SBR splits
    tab.st <- cbind(tab.st, 'Hellinger' =  c(-1, distHell), 'Entropy' = ent, 'MI' = c(-1, distMI))
    if(cl.stat.KL) tab.st <- cbind(tab.st, 'symKL' = c(-1, distSymKL))
    if(cl.stat.KL) tab.st <- cbind(tab.st, 'KL' = c(-1, distKL))
    if(cl.stat.TE) tab.st <- cbind(tab.st, 'symTE' = c(-1, distSymTE))
  }else{ # otherwise
    statistics <- list('Hellinger' =  distHell, 'Entropy' = ent, 'MI' = distMI); ele <- ele + 3
    if(cl.stat.KL) { statistics[[ele]] <- distSymKL; names(ele)[ele] <- 'symKL'; ele <- ele + 1 }
    if(cl.stat.KL) { statistics[[ele]] <- distKL; names(ele)[ele] <- 'KL'; ele <- ele + 1 }
    if(cl.stat.TE) { statistics[[ele]] <- distSymTE; names(ele)[ele] <- 'symTE'; ele <- ele + 1 }
  }
  if(cl.stat.wMI) { statistics[[ele]] <- sl_MI; names(ele)[ele] <- 'winMI'; ele <- ele + 1 }
  if(cl.stat.stft) { statistics[[ele]] <- xpeaks; names(ele)[ele] <- 'stFT'; ele <- ele + 1 }
  if(cl.stat.weight.barriers) {
    if(!silent) cat('Keeping only the barriers which are reasonably far from borders.\n')
    select_non_border <- breaks > lpi/n_unif[.lt(n_unif)] & breaks < lpi - lpi/n_unif[.lt(n_unif)]
    if(!silent) cat('Discarded the following barrier indexes (from left / also on the stats out):', seq(1, .lt(breaks))[!select_non_border], '\n')
    forout <- MI_ratio
    forout[!select_non_border] <- -1
    tab.st <- cbind(tab.st, 'barWeight' = c(-1, forout))
  }

  if(plot){
    if(new.dev) dev.new(width=15, height=10)
    save_par <- par()
    par(mgp=c(0, 0.4, 0))
    par(ps=6)
    par(mar = c(3.5, 0, 1, 3), oma = c(2, 4, 2, 2))
    margin <- round(cstored*0.99)
    cx <- 2.2
    sc <- 1.8
    xr1 <- c((cstored-margin):margin)
    xr <- c(1, cstored)
    yr <- range(progind$Time)*sc
    plot(0, 0, main="", xlim=xr, ylim=yr, xlab="", ylab="", axes=FALSE, frame.plot=FALSE, type="n")
    kin.pl <- -log(cutf / cstored)[xr1]
    xx.lab <- c(1,round(breaks),cstored)
    axis(1, at=xx.lab, tck=.01, cex.axis=1.8)
    axis(3, labels=rep("", .lt(xx.lab)), at=xx.lab, tck=.01)
    mtext("Progress Index", side=1, line=1.5, cex=cx )
    yy.lab1 <- format(c(min(kin.pl), min(kin.pl)+(max(kin.pl)-min(kin.pl))*c(1:3)/3), digits=2)
    axis(2, labels=yy.lab1, at=scale(as.numeric(yy.lab1)), las=3, tck=.01, cex.axis=cx)
    mtext(expression("ln(("*italic(tau["SA"]+tau["AS"])*")/2)"), at=max(progind$Time)/2, side=2, line=1.8, cex=cx)
    yy.lab2 <- round(c(1, c(1:5)/5*max(progind$Time)))
    axis(4, labels=rep("", .lt(yy.lab2)), at=max(progind$Time)+round(yy.lab2*(sc-1)), las=2, tck=.01, hadj=-0.6, col="red")
    mtext(yy.lab2, side=4, las=2, line=0.2, at=max(progind$Time)+yy.lab2*(sc-1), col="red", cex=1.8)
    mtext("Time", at=max(progind$Time)*(sc+1)/2, side=4, line=2.4, cex=cx, col="red")
    
    if (cstored>10000) { lst <- seq(1, cstored, by=5)
    } else lst <- c(1:cstored)
    points(progind$PI[lst], max(progind$Time[lst])+progind$Time[lst]*(sc-1), cex=0.005, col="tomato1")
    lines(xr1, scale(kin.pl), lwd=1, col="black")
    lines(xr1, scale(kin[xr1]), lwd=0.8, col="forestgreen") #the one without the parabol
    if(!is.null(cl.stat.denat) && cl.stat.denat == 'poly_interpolation') { # my simple fit (dgarol)
      if(!is.null(denat.MI)) den_kin <- .denaturate(-log(cutf / cstored), seq(cstored), polydeg = denat.MI, plotit = FALSE)
      else den_kin <- .denaturate(-log(cutf / cstored), seq(cstored), polydeg = 7, plotit = FALSE)
      lines(xr1, scale(den_kin[xr1]), lwd = 0.8, col = 'darkblue')
    }
    if(only.kin) abline(v=breaks, lwd=0.7, col="orange4")
    else if(!match) { # we want to plot all the vertical line also to see the effects on the plot
      ## Not-matched vdyn
      abline(v=vdyn, lwd=0.4, col="blue")
      ## Not-matched vkin
      abline(v=vkin[-match(brk.mtc, vkin)], lwd=0.4, col="orange4")
      ## Matched breaks
      abline(v=brk.mtc, lwd=0.7, col="black")
    } else abline(v=breaks, lwd=0.7, col="black")
    if(cl.stat.weight.barriers){
      select_non_border <- breaks > lpi/n_unif[.lt(n_unif)] & breaks < lpi - lpi/n_unif[.lt(n_unif)]
      # if(!silent) cat('Discarded the following barrier indexes (from left):', seq(1, .lt(breaks))[!select_non_border], '\n')
      # select_non_border <- -c(1, .lt(breaks))
      points(breaks[select_non_border], MI_ratio[select_non_border]*max(yr), pch ='+', col = 'red', cex = 5)
      abline(h=mean(MI_ratio)*max(yr), lwd=1.1, col= 'grey')
    }
    suppressWarnings(par(save_par))
  }
  invisible(list(tab.st=tab.st, nbins=c(nx,ny), seq.st=seq.st[order(progind$Time)], statistics = statistics, call=call, data.out = data.out))
}

