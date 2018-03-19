#' @title Identifying basins in the SAPPHIRE plot
#' @description
#'     \code{basins_recognition} uses the information provided by the SAPPHIRE plot in order to identify free energy basins, 
#'     thus leading to a discretized trajectory. In particular, it inspects both the kinetic annotation and times of occurence
#'     (called dynamic annotation hereafter) of the progress index and, at the end, it matches (or merges) the resulting state barriers 
#'     coming from each of the single analysis.  
#'     The analysis of the dynamical trace is driven by an underlying 2-D histogram whose number of bins has to be given either manually 
#'     or automatically (this option is available only for number of bins on the y-axis). These values has to be chosen appropriately because 
#'     they provide the minimum size that a state should have in order to be identified.
#'     Regarding the kinetic annotation, the identification of the states is equivalent to discerning its main peaks. This task is performed 
#'     by means of a smoothing filter, either moving average or Savitzky Golay, and a naive algorithm to find the local maxima. The window size 
#'     of the filter and of the peak search algorithm is equal to the number of snapshots divided by \code{nx}.
#'     Eventually, the user can decide on either selecting only the partitions in common to both the annotation (matching) or retaining them all (merging).
#'     
#' @param data Name of the PROGIDX_<...>.dat file or Data frame with three columns containing, in order, the (sorted) progress index, the relative time indices and the cut function. 
#' @param nx Number of bins on the x-axis of the 2-D histogram.
#' @param ny Number of bins on the y-axis of the 2-D histogram. Default to nx.
#' @param ny.aut Logical value indicating whether a suitable number of bins on the y-axis has to be identified automatically. Default value is \code{FALSE}.
#' @param local.cut Logical value indicating whether the localized cut function has to be used (see references). Default value is \code{FALSE}.
#' @param match Logical value indicating whether results from kinetic and dynamic annotation have to be either matched or merged. Default value is \code{FALSE}               
#' @param dyn.check The number of times the statistical check has to run on the partitions. 1 by default, 2 usually is enough. 
#' @param avg.opt Smoothing filter in the kinetic annotation analysis
#'       \itemize{
#'            \item "\code{movav}" for moving average filter
#'            \item "\code{SG}" for savitzkyGolay filter
#'       }
#' @param plot A logical value indicating whether to display the SAPPHIRE plot with the resulting partitions or not. Black partitions 
#' are the matched ones, blue ones derive only from the dynamic analysis and orange ones only from the kinetic analysis. The green curve is the kinetic annotation (black curve) 
#' where the parabolic shape has been subtracted, i.e. the actual curve used for the peaks identification. Default value is \code{FALSE}.
#' @param new.dev A logical value indicating whether a new window/device has to plotted or not when \code{plot=TRUE}.
#' @param out.file A logical value indicating whether to write an output file with the state sequence ordered with respect to the PI. Default value is \code{TRUE}.
#' @param silent A logical value indicating whether the function has to remain silent or not. Default value is \code{FALSE}.
#' @param ...
#'      \itemize{
#'          \item "\code{time.series}" File name. If specified, it substitutes the time series of the PROGIDX_<..> file with the one provided by the file.
#'          \item "\code{only.kin}" Logical value indicating whether it should be performed only the analysis of the kinetic annotation, thus disregarding completely 
#'          the dynamic trace.  
#'          \item "\code{pol.degree}" Degree of the Savitzky-Golay filter in case it was chosen. Default to 2.
#'          \item "\code{cl.stat}" If \code{TRUE} it will activate the cluster analysis with Mutual Information, Hellinger distance and Shannon Entropy
#'          as defaults.
#'          \item "\code{cl.stat.weight.barriers}" If \code{TRUE} (standard once cl.stat is active) it will calculate two different Mutual Information: 
#'          one based on the clusters founds and one based on uniform subsampling (done using \code{cl.stat.nBreaks} or the number of cluster founds using 
#'          SBR multiplied for 10). The second one will be inverted to find the discontinuities and the barriers will be weighted simply multiplying the two values 
#'          found on the specific progress index data point.  
#'          \item "\code{cl.stat.nBreaks}" Integer. This variable defines the number of splits for the analysis. If set to 0 it will use the breaks found in the SBR.
#'          \item "\code{cl.stat.denat}" Logical. \code{TRUE} if should the kinetic annotation and the MI based on uniform suvdivision be rescaled for parabolic 
#'          artifacts. The polynomial fit is by default 7 and 12.
#'          \item "\code{cl.stat.MI_comb}" This value is the representative of how the various calculation in barrier weighting are combined. In particular the options 
#'          are the following:
#'          \itemize{
#'            \item "\code{mean}" Mean between the MI based on the SBR and the MI based on uniform subdivisions
#'            \item "\code{multip}" Multiplication between the MI based on the SBR and the MI based on uniform subdivisions
#'            \item "\code{kin_ann}" Mean between the kinetic annotaton and the MI based on the SBR and the MI based on uniform subdivisions
#'            \item "\code{kin}" Only kinetic annotation value on the barrier
#'            \item "\code{ann}" Only the MI based on the SBR and the MI based on uniform subdivisions
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
#' @importFrom graphics abline axis legend lines mtext par points text
#' @importFrom stats as.ts coef filter lm nls sd weighted.mean
#' @importFrom utils head tail
#' @importFrom gplots hist2d
#' @importFrom outliers grubbs.test
#' @importFrom prospectr movav savitzkyGolay
#' @importFrom splus2R peaks
#' @importFrom grDevices dev.new
#' @importFrom data.table fread fwrite
#' @importFrom RColorBrewer brewer.pal
#' @importFrom TransferEntropy computeTE
#' @importFrom infotheo mutinformation
#' @importFrom e1071 stft
# @importFrom RcppRoll roll_mean
#' 
#' @export basins_recognition

basins_recognition <- function(data, nx, ny=nx, ny.aut=FALSE, local.cut=FALSE, match=FALSE, dyn.check=1, 
                               avg.opt=c("movav", "SG"), plot=FALSE, new.dev=TRUE, out.file=TRUE, silent=FALSE, ...) {
  call <- match.call()
  
  if(!is.character(data) && !is.data.frame(data)) stop("data must be a string or a data frame")
  if(is.character(data) && (!all(grepl("PROGIDX", data)) && !all(grepl("REPIX", data)))) stop("Please provide a data name starting with 'PROGIDX' or 'REPIX'" )
  if((nx %% 1) != 0) stop("nx must be an integer")
  if((ny %% 1) != 0) stop("ny must be an integer")
  if((nx < 7 || ny < 7)) stop('nx and ny must be > 7')
  if(!is.logical(local.cut)) stop("local.cut must be a logical value")
  if(!is.logical(match)) stop("match must be a logical value")
  if(!(dyn.check %in% c(1:10))) stop("dyn.check value not valid")
  if(!is.logical(plot)) stop("plot must be a logical value")
  if(!is.logical(new.dev)) stop("new.dev must be a logical value")
  if(!is.logical(out.file)) stop("out.file must be a logical value")

  # Extra arguments checks
  input.args <- list(...)
  avail.extra.arg <- c("pol.degree", "only.kin", "time.series",
                       "cl.stat", "plot.cl.stat", 'cl.stat.entropy', 'cl.stat.weight.barriers',
                       'cl.stat.stft', 'cl.stat.TE', 'cl.stat.KL', 'cl.stat.wMI',
                       'cl.stat.nUni', 'cl.stat.nBreaks', 'cl.stat.MI_comb', 'cl.stat.denat')
  
  if(!is.null(names(input.args)) && any(!(names(input.args) %in% avail.extra.arg))) 
    warning('There is a probable mispelling in one of the inserted variables. Please check the available extra input arguments.')
  if("only.kin" %in% names(input.args)) {
    only.kin <- input.args$only.kin
    if(only.kin) match <- TRUE
  } else only.kin <- FALSE
  
  # cluster statistics check - dgarol
  cl.stat.wMI <- FALSE      # can be long
  cl.stat.KL <- FALSE
  cl.stat.TE <- FALSE
  cl.stat.stft <- FALSE     # it is plotting too much
  cl.stat.weight.barriers <- FALSE
  cl.stat.MI_comb <- 'kin_ann'
  cl.stat.denat <- FALSE
  
  if("cl.stat" %in% names(input.args)) { # dgarol
    cl.stat <- input.args$cl.stat
    cl.stat.weight.barriers <- TRUE
    stopifnot(is.logical(cl.stat))
  } else cl.stat <- FALSE
  
  # Checking the number of breaks for the statistics (0 is the standard barriers)
  if("cl.stat.nBreaks" %in% names(input.args)) { # dgarol
    cl.stat.nBreaks <- input.args$cl.stat.nBreaks
    stopifnot(.isSingleInteger(cl.stat.nBreaks))
    if(!cl.stat && cl.stat.nBreaks) cl.stat <- TRUE
  } else cl.stat.nBreaks <- 0
  
  # checking the uniform sampling nsplits
  if("cl.stat.denat" %in% names(input.args)) { # dgarol
    cl.stat.denat <- input.args$cl.stat.denat
    stopifnot(is.logical(cl.stat.denat))
    if(!cl.stat && cl.stat.denat) cl.stat <- TRUE
  }
  
  if("cl.stat.MI_comb" %in% names(input.args)) { # dgarol
    cl.stat.MI_comb <- input.args$cl.stat.MI_comb
    stopifnot(!is.character(cl.stat.MI_comb))
    stopifnot(!(cl.stat.MI_comb %in% c('ann', 'kin', 'mean', 'kin_ann', 'multip')))
    if(!cl.stat && cl.stat.MI_comb) cl.stat <- TRUE
  } 
  
  if("cl.stat.nUni" %in% names(input.args)) { # dgarol
    cl.stat.nUni <- input.args$cl.stat.nUni
    stopifnot(any(sapply(a, function(x) x%%1) != 0))
    if(!cl.stat && cl.stat.nUni) cl.stat <- TRUE
  } else cl.stat.nUni <- NULL
  
  # Checking the weight of the barriers
  if("cl.stat.weight.barriers" %in% names(input.args)) { # dgarol
    cl.stat.weight.barriers <- input.args$cl.stat.weight.barriers
    stopifnot(is.logical(cl.stat.weight.barriers))
    if(!cl.stat && cl.stat.weight.barriers) cl.stat <- TRUE
  }
  # Specific cluster statistics
  if("cl.stat.wMI" %in% names(input.args)) { # dgarol
    cl.stat.wMI <- input.args$cl.stat.wMI
    stopifnot(is.logical(cl.stat.wMI))
    if(!cl.stat && cl.stat.wMI) cl.stat <- TRUE
  }
  if("cl.stat.KL" %in% names(input.args)) { # dgarol
    cl.stat.KL <- input.args$cl.stat.KL
    stopifnot(is.logical(cl.stat.KL))
    if(!cl.stat && cl.stat.KL) cl.stat <- TRUE
  }
  if("cl.stat.TE" %in% names(input.args)) { # dgarol
    cl.stat.TE <- input.args$cl.stat.TE
    stopifnot(is.logical(cl.stat.TE))
    if(!cl.stat && cl.stat.TE) cl.stat <- TRUE
  }
  if("cl.stat.stft" %in% names(input.args)) { # dgarol
    cl.stat.stft <- input.args$cl.stat.stft
    stopifnot(is.logical(cl.stat.stft))
    if(!cl.stat && cl.stat.stft) cl.stat <- TRUE
  }
  
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
  
  # other checks
  avg.opt.arg <- c("movav", "SG")
  if(!(avg.opt[1] %in% avg.opt.arg)) stop("Average option not valid")
  if(avg.opt[1]=="SG") {
    if(!("pol.degree" %in% names(input.args))) {
      if(!silent) cat("SG but pol.degree not specified, set to default value 2\n")
      pol.degree <- 2
    } else if(!(input.args$pol.degree %in% c(2:6))) {
      warning("Degree of the polynomial not valid, set to default value 2")
      pol.degree <- 2 
    } else pol.degree <- input.args$pol.degree
  }

  # --------------------------- functions
  as.real <- function(x) {
    return(as.double(x))
  }
  scale <- function(x) {
    scl <- (max(progind$Time)-1)/(max(x)-min(x))
    return(1+(x-min(x))*scl)
  }
  refine <- function(vec) {
    ## Attach to the peak in the cut function
    for (jj in 1:length(vec)) {
      vec <- round(vec)
      tmp <- which.max(kin[(vec[jj]-round(sd.kin/2)):(vec[jj]+round(sd.kin/2))]) 
      vec[jj] <- tmp+vec[jj]-round(sd.kin/2)-1
    }
    vec <- unique(c(1,vec,cstored))
    if(length(which(vec<=sd.kin))>1) vec <- vec[-c(2:max(which(vec<sd.kin)))]
    if(length(which((cstored-vec)<=sd.kin)>1)) vec <- c(vec[-which((cstored-vec)<sd.kin)],cstored)
    while (any(diff(vec)<sd.kin)) {
      vec.tmp <- vec
      dst <- 0
      start <- NULL
      for(ii in 2:(length(vec)-1)){
        dst <- dst+vec[ii+1]-vec[ii]
        if(dst<sd.kin) {
          start <- c(start,ii)
          next
        } else if (length(start)==1) {
          start <- c(start,ii)
          ll <- start[which.min(kin[vec[start]])]
          vec.tmp <- vec.tmp[-match(vec[ll],vec.tmp)]
        } else if (length(start)>1) {
          start <- c(start,ii)
          ll <- which.max(kin[vec[start]])
          lst <- start[-ll]
          vec.tmp <- vec.tmp[-match(vec[lst], vec.tmp)]
        }
        start <- NULL
        dst <- 0
      }
      vec <- vec.tmp
    }
    vec <- vec[which(vec != 1)]
    vec <- vec[which(vec != cstored)]
    if(length(vec)<3) warning("Vector completely erased in refine function")
    return(vec)
  }
  true.peaks <- function(idxx, values, up=TRUE) {
    ## To select the effective peaks or minima, it return a selection of idxx 
    rl <- rle(diff(idxx))
    if(any(rl$values==1)) {
      disc <- NULL
      st <- c(0, cumsum(rl$lengths))[which(rl$values==1)] + 1
      end <- cumsum(rl$lengths)[which(rl$values==1)]
      for(ii in seq_along(st)) {
        ## max() means that in case two consecutive index have the same values, I pick the one on the right (.. don't know if it's the best choice honestly)
        if(up) idx.max <- max(which(values[idxx[st[ii]:(end[ii]+1)]]==max(values[idxx[st[ii]:(end[ii]+1)]])))
        else idx.max <- max(which(values[idxx[st[ii]:(end[ii]+1)]]==min(values[idxx[st[ii]:(end[ii]+1)]])))
        disc <- c(disc, idxx[st[ii]:(end[ii]+1)][-idx.max])
      }
      return(idxx[-match(disc, idxx)])
    } else return(idxx)
  }
  myHell <- function(x, y){
    return(sqrt(1-sum(sqrt(x*y))))
  }
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
  }
  if(nrow(progind) < 80) stop('It is impossible to recognize basins with less than 80 snapshots in the sapphire table.')
  cstored <- dim(progind)[1]
  
  ## Dynamic parameters
  perc <- 0.90 ## Percentage to find optimal nbiny
  xidx <- 0.5 ## Joining consecutive cells on a single raw
  K <- 2 ## Parameter of the Haat wavelet
  dpeaks.dyn <- 3 ## Parameter of the maximum search algorithm (window size)
  nsample <- 50 ## Number of samples in the distribution of reshuffled Hellinger distances 
  cutjoin <- 50 ## Number of null joining attempts required to quit the procedure 
  conf.lev <- 0.005 ## Cut on pvalues of the Grubb Test
  ## Kinetic parameters
  wsize <- round(2*cstored/nx) + ((round(2*cstored/nx)+1) %% 2)  ## To make it odd
  dpeaks.kin <- ceiling(wsize/2)+ceiling(wsize/2)%%2+1 
  thr.ratio <- 0.05 ## Parameter of second data cleaning  
  ## Matching parameters
  lx <- round(cstored/nx)
  sd.kin <- 1*lx
  sd.dyn <- 1*lx
  
  
  #########&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&###############
  ############################&&&&&& DYNAMICS &&&&&&&###############################
  #########&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&###############
  
  if(!only.kin) {
    if(!silent) cat("Analysis of the dynamic trace...\n")
    
    ##################################################################################
    ## NBIN Y IDENTIFICATION
    ##################################################################################
    
    if(ny.aut) {
      idx <- 0
      nbin <- 0
      lthbin <- NULL
      seqbin <- seq(from=10, to=min(2000,cstored/10), by=10)
      for (nbin in seqbin) {
        lth <- NULL
        br <- seq(from=1, to=cstored, length.out=(nbin+1))
        for (i in 1:nbin) {
          smp <- progind[,2][round(br[i]):round(br[i+1])]
          lth[i] <- max(smp)-min(smp)
        }
        idx <- idx+1
        lthbin[idx] <- sd(lth)
      }
      a.0 <- max(lthbin)
      idx.fin <- which(lthbin==a.0)
      expmodel.0 <- lm(log(a.0-lthbin)[-idx.fin] ~ seqbin[-idx.fin])
      expmodel <- try(nls(lthbin ~ a-b*exp(seqbin*c), start=list(a=a.0, b=exp(coef(expmodel.0)[1]), c=coef(expmodel.0)[2] ) ), silent=TRUE)
      fracmodel <- try(nls(lthbin ~ a*(seqbin/(l+seqbin)), start=list(a=max(lthbin), l=(max(lthbin)-min(lthbin)/2))), silent=TRUE)
      
      if(is.list(expmodel)) {
        theor.exp <- coef(expmodel)["a"]-coef(expmodel)["b"] * exp(coef(expmodel)["c"] * seqbin)
        X2.exp.rid <- sum((theor.exp-lthbin)^2/theor.exp)
      }
      if(is.list(fracmodel)) {
        theor.frac <- coef(fracmodel)[1]*(seqbin/(coef(fracmodel)[2]+seqbin))
        X2.frac.rid <- sum((theor.frac-lthbin)^2/theor.frac)
      }
      if(is.list(fracmodel) && is.list(expmodel)) {
        if (X2.exp.rid < X2.frac.rid) {
          if(!silent) cat("Exponential model chosen\n")
          can <- (1/coef(expmodel)["c"])*log((coef(expmodel)["a"]/coef(expmodel)["b"]) *(1-perc))
        } else {
          if(!silent) cat("Hyperbola model chosen\n")
          can <- coef(fracmodel)[2]*perc/(1-perc) 
        }
        ny <- as.numeric(round(can))
      } else if (sum(c(is.list(expmodel), is.list(fracmodel)))==1) {
        if (is.list(expmodel)) {
          if(!silent) cat("Exponential model chosen, hyperbola not available\n")
          can <- (1/coef(expmodel)["c"])*log((coef(expmodel)["a"]/coef(expmodel)["b"]) *(1-perc))
        } else {
          if(!silent) cat("Hyperbola model chosen, exponential not available\n")
          can <- coef(fracmodel)[2]*perc/(1-perc) 
        }
        ny <- as.numeric(round(can))
      } else {
        if(!silent) cat("Non linear fit doesn't work... ny set to nx default value.\n")
        ny <- 500
      }
    }
    
    ##################################################################################
    ## 2-D HISTOGRAM
    ##################################################################################
    
    if(!silent) cat("Nbins on y is", ny, "\n")
    hist <- hist2d(matrix(c(progind[,1],progind[,2]), ncol=2, nrow=cstored), nbins=c(nx,ny), show=FALSE)
    
    ##################################################################################
    ## STRETCHES CREATION
    ##################################################################################
    ## Joining cells through density criterion: the larger values of xids the larger tolerance 
    joinx <- hist$counts
    joinx[,] <- 0
    for (j in 1:ny) {
      for(i in 1:(nx-1)) {
        if(hist$counts[i,j] < 1) next #mettere un minimo?
        if((hist$counts[i,j] >= (hist$counts[i+1,j]*(1-xidx))) & (hist$counts[i,j] <= hist$counts[i+1,j]))  joinx[i,j] <- 1
        if((hist$counts[i+1,j] >= (hist$counts[i,j]*(1-xidx))) & (hist$counts[i+1,j] < hist$counts[i,j]))  joinx[i,j] <- 1
        ## Joining if separated by an empty bin as well
        if (i==(nx-1)) next
        if((hist$counts[i,j] >= (hist$counts[i+2,j]*(1-xidx))) & (hist$counts[i,j] <= hist$counts[i+2,j]))  joinx[c(i,i+1),j] <- 1
        if((hist$counts[i+2,j] >= (hist$counts[i,j]*(1-xidx))) & (hist$counts[i+2,j] < hist$counts[i,j]))  joinx[c(i,i+1),j] <- 1
      }
    }
    
    ####################################################################################
    ## BIGG PARTITION
    ###############################################################################
    ## ovlap is the measure of overlaps of the ministretchs along the x-axis
    ovlap <- as.vector(rowSums(joinx))
    ovlap.der <- c(diff(ovlap),0) + c(0, diff(ovlap))
    interval <- c(mean(ovlap.der)-2*sd(ovlap.der), mean(ovlap.der)+2*sd(ovlap.der))
    idxx.up <- which(ovlap.der>interval[2])
    idxx.down <- which(ovlap.der<interval[1])
    if(.lt(idxx.up)>1) idxx.up.cl <- true.peaks(idxx.up, ovlap.der, up=TRUE)
    else idxx.up.cl <- NULL
    if(.lt(idxx.down)>1) idxx.down.cl <- true.peaks(idxx.down, ovlap.der, up=FALSE)
    else idxx.down.cl <- NULL
    if(is.null(idxx.up.cl) & is.null(idxx.down.cl)) bigg.idx <- NULL
    else bigg.idx <- sort(c(idxx.up.cl, idxx.down.cl))
    bigg.brk <- hist$x[bigg.idx]
    
    ## source("./Rfunctions/SBR_images_functions.R")
    ## stretches.crosses(progind$Time, rawset, hist, joinx, ovlap, ovlap.der, bigg.idx)
    
    
    ####################################################################################
    ## FILL RAWSET
    ###############################################################################
    rawset <- data.frame(min=rep(0,ny), max=rep(0,ny), center=rep(0,ny), wth1=rep(0,ny), wth2=rep(0,ny))
    for (j in 1:ny) {
      temp <- which(joinx[,j]==1)
      if (length(temp)==0) next      
      minr <- temp[1]
      maxr <- temp[length(temp)]+1
      if((maxr-minr)<= 3) rawset[j,c(1:3)] <- 0
      else  rawset[j,c(1:3)] <- c(hist$x[minr], hist$x[maxr], weighted.mean(hist$x[minr:maxr], hist$counts[c(minr:maxr),j]))
    }
    dum <- cbind(rawset$max-rawset$center, rawset$center-rawset$min)
    rawset$wth1 <- (apply(dum, 1, min)/apply(dum, 1, max))*(1-((rawset$max-rawset$min)/cstored))
    rawset$wth2 <- (apply(dum, 1, min)/apply(dum ,1, max)) / (rawset$max-rawset$min)
    rawset$wth1[which(is.na(rawset$wth1))] <- 0
    rawset$wth2[which(is.na(rawset$wth2))] <- 0
    if(any(!is.finite(unlist(rawset)))) stop("Error in rawset")
    
    rm("joinx")
    
    #################################################################################
    ## WEIGHTED SUMS FUNCTIONS
    ##############################################################################
    ## Forward
    sumwr <- function(rawset, first, meanopt, wthopt) {
      rawsetsort <- rawset[order(rawset$center),]
      if (first=="min") rawfirst <- rawsetsort$min
      if (first=="max") rawfirst <- rawsetsort$max
      if (first=="center") rawfirst <- rawsetsort$center
      if (meanopt=="min") rawmean <- rawsetsort$min
      if (meanopt=="max") rawmean <- rawsetsort$max
      if (meanopt=="center") rawmean <- rawsetsort$center
      if (wthopt==1) rawwth <- rawsetsort$wth1
      if (wthopt==2) rawwth <- rawsetsort$wth2
      idrs <- match(FALSE, rawsetsort$center==0)
      if (is.na(idrs)) idrs <- 1
      sm <- NULL
      for (i in seq(nx)) {
        idxlist <- which(rawfirst[idrs:ny] < hist$x.breaks[i+1]) + idrs - 1 
        if(length(idxlist)==0) sm[i] <- 0
        else sm[i] <- sum(rawwth[idxlist])
      }
      return(sm)
    }
    
    #################################################################################
    ## Backward
    backsumwr <- function(rawset, first, meanopt, wthopt) {
      rawsetsort <- rawset[order(rawset$center),]
      if (first=="min") rawfirst <- rawsetsort$min
      if (first=="max") rawfirst <- rawsetsort$max
      if (first=="center") rawfirst <- rawsetsort$center
      if (meanopt=="min") rawmean <- rawsetsort$min
      if (meanopt=="max") rawmean <- rawsetsort$max
      if (meanopt=="center") rawmean <- rawsetsort$center
      if (wthopt==1) rawwth <- rawsetsort$wth1
      if (wthopt==2) rawwth <- rawsetsort$wth2
      idrs <- match(FALSE,rawsetsort$center==0)
      if (is.na(idrs)) idrs <- 1
      sm <- NULL
      for (i in nx:1) {
        idxlist <- which(rawfirst[idrs:ny] > hist$x.breaks[i]) + idrs - 1
        if(length(idxlist)==0) sm[i] <- 0
        else sm[i] <- sum(rawwth[idxlist])
      }
      return(sm)
    }
    
    #########################################################################################
    ##FILTER SECTION
    ########################################################################################
    Haarfilt <- c(rep(1,K),1,rep(-1,K))
    
    #################################################################################
    ## Forward Filter with Haar wavelet
    MaxHaar <- function(rawset, first, meanopt, wthopt) {
      filt <- filter(as.ts(sumwr(rawset, first, meanopt, wthopt)), Haarfilt, method="convolution", sides=2)
      filt <- c(rep(0,K),filt[c((K+1):(length(filt)-K))],rep(0,K))  #Filt with 0 instead of NA
      xbr <- hist$x[which(peaks(filt,dpeaks.dyn))]
      return(c(1,xbr,cstored))
    }
    #################################################################################
    ## Backward Filter with Haar wavelet
    BackMaxHaar <- function(rawset, first, meanopt,wthopt) {
      filt <- filter(as.ts(rev(backsumwr(rawset, first, meanopt, wthopt))), Haarfilt, method="convolution", sides=2)
      filt <- c(rep(0,K),filt[c((K+1):(length(filt)-K))],rep(0,K))  #Filt with 0 instead of NA
      xbr <- hist$x[which(rev(peaks(filt, dpeaks.dyn)))]
      return(c(1,xbr,cstored))
    }
    
    #################################################################################
    ##CHOICE OF THE BREAKS
    #################################################################################
    breaks.max <- MaxHaar(rawset, "max","min",1) 
    breaks.min <- BackMaxHaar(rawset, "min","max",1)
    
    rm("rawset")
    
    ###############################################################################
    ##HARD BREAKS: Joining selected breaks.min and breaks.max with res
    ###############################################################################
    if(.lt(breaks.max)==2 && .lt(breaks.min)>2) {
      breaks.tot <- breaks.min
    } else if (.lt(breaks.max)>2 && .lt(breaks.min)==2) {
      breaks.tot <- breaks.max
    } else if (.lt(breaks.max)==2 && .lt(breaks.min)==2) {
      breaks.tot <- c(1, cstored)
    } else {
      sep <- NULL
      idx <- 0
      selbreaks.max <- breaks.max[-c(1,.lt(breaks.max))]
      selbreaks.min <- breaks.min[-c(1,.lt(breaks.min))]
      softbreaks.max <- selbreaks.max
      softbreaks.min <- selbreaks.min
      for (i in 1:.lt(selbreaks.max)) {
        selcell.max <- which(hist$x==selbreaks.max[i])
        if (i==.lt(selbreaks.max)) selcell2.max <- .lt(hist$x)
        else selcell2.max <- which(hist$x==selbreaks.max[i+1])
        for (j in seq_along(selbreaks.min)) {
          selcell.min <- which(hist$x==selbreaks.min[j])
          if (selcell.max==selcell.min) {
            idx <- idx+1
            sep[idx] <- hist$x[selcell.max]
            softbreaks.max <- softbreaks.max[-which(softbreaks.max==selbreaks.max[i])]
            softbreaks.min <- softbreaks.min[-which(softbreaks.min==selbreaks.min[j])]
            next
          }
          if ((selcell.max+1) == selcell.min) {
            if (!is.na(match(hist$x[selcell.max], sep))) {
              ## If already in sep it's just removed
              softbreaks.min <- softbreaks.min[-which(softbreaks.min==selbreaks.min[j])]
              next
            }
            idx <- idx+1
            sep[idx] <- hist$x[selcell.max] #Sep border included in the trailing partition -> =blue
            softbreaks.max <- softbreaks.max[-which(softbreaks.max==selbreaks.max[i])]
            softbreaks.min <- softbreaks.min[-which(softbreaks.min==selbreaks.min[j])]
            next
          }
          if ((selcell.max+2) == selcell.min & (selcell2.max-selcell.max) >= 3) {
            ## if are distant two cells and the subsequent one is far enough
            if (!is.na(match(hist$x[selcell.max],sep))) {
              softbreaks.min <- softbreaks.min[-which(softbreaks.min==selbreaks.min[j])]
              next
            }
            idx <- idx+1
            sep[idx] <- hist$x[selcell.min-1]
            softbreaks.max <- softbreaks.max[-which(softbreaks.max==selbreaks.max[i])]
            softbreaks.min <- softbreaks.min[-which(softbreaks.min==selbreaks.min[j])]
            next
          }
          ## if (selcell.max[i]>selcell.min+3) break    
        }
      }
      breaks.tot <- sort(c(1, cstored, sep, softbreaks.max, softbreaks.min))
    }
    
    #################################################################################
    ## UNIFYING with BIGG Partition 
    ################################################################################
    
    ## It works in pathological cases (i.e. NULL vecto)
    ## if(!silent) cat("Number of bigg is", .lt(bigg.brk), "\n")
    breaks.tot <- unique(sort(c(breaks.tot, bigg.brk)))
    ## if(!silent) cat("Number of matched is", .lt(bigg.brk), "\n")
    
    if(any(breaks.tot==hist$x[nx])) breaks.tot <- breaks.tot[-which(breaks.tot==hist$x[nx])]
    if(any(breaks.tot==hist$x[1])) breaks.tot <- breaks.tot[-which(breaks.tot==hist$x[1])]
    
    #################################################################################
    ## HISTOGRAM of each PARTITION 
    ################################################################################
    if(.lt(breaks.tot)==2) {
      brk.dyn <- breaks.tot
    } else {
      ## Breaks.Tot
      for(iii in seq(dyn.check)) {  
        brkjy <- matrix(rep(0,(.lt(breaks.tot)-1)*ny), nrow=ny, ncol=(.lt(breaks.tot)-1))
        for (i in 1:(.lt(breaks.tot)-1)) {
          if (i==1) {
            ncls <- 1
            ncle <- which(hist$x==breaks.tot[i+1])
          } else if (i == (.lt(breaks.tot)-1) ) {
            ncls <- which(hist$x==breaks.tot[i])+1
            ncle <- nx
          } else {
            ncls <- which(hist$x==breaks.tot[i])+1
            ncle <- which(hist$x==breaks.tot[i+1])
          }
          ## ColSums doesn't work if ncls==ncle
          for (j in 1:ny) brkjy[j,i] <- sum(hist$counts[c(ncls:ncle),j])
        }
        dens <- apply(brkjy, 2, function(x) x/sum(x))
        if(nrow(dens)!=ny) dens <- t(dens)
        
        
        ########################################################################
        ##JOINING PARTITIONS METHODS
        ########################################################################
        ## Computation of Distances Hell between consecutive partitions
        
        distHell <- sapply(seq(.lt(breaks.tot)-2), function(j) myHell(dens[,j], dens[,j+1]))
        discbreaks <- NULL
        ncounts <- NULL
        flagbreak <- 0
        ll <- 0
        lstHell <- order(distHell)
        sampleHell <- rep(-1, nsample)
        for (i in lstHell) { 
          l1 <- which(hist$x==breaks.tot[i])+1 
          if (i==1) l1 <- 1
          l2 <- which(hist$x==breaks.tot[i+1])
          l3 <- which(hist$x==breaks.tot[i+2])
          if (i==.lt(lstHell)) l3 <- nx
          pr <- (l2-l1+1)/(l3-l1+1)
          if(!silent) {
            cat.str <- paste("Would-be basins", paste(round(breaks.tot[i:(i+2)]), collapse=" "))
            cat(cat.str)
            nlett <- nchar(cat.str)
          }
          ncounts <- colSums(hist$counts[c(l1:l3),])
          for (idx in seq(nsample)) {
            unif1 <- rep(0,ny)
            unif2 <- rep(0,ny)
            for (j in 1:ny) {
              if(ncounts[j]==0) next
              dum <- sample(c(0:1), ncounts[j], replace=TRUE, prob=c(pr,1-pr))
              unif1[j] <- ncounts[j]-sum(dum)
              unif2[j] <- sum(dum)
            }
            if(sum(unif1)==0 || sum(unif2)==0) sampleHell[idx] <- 1
            else sampleHell[idx] <- myHell(unif1/sum(unif1), unif2/sum(unif2))
          }
          grubbsHell <- grubbs.test(c(sampleHell, distHell[i]), type=10)
          if (distHell[i] < max(sampleHell) | grubbsHell$p.value>conf.lev) {
            if(!silent) cat(rep(" ",15+4+3*nchar(as.character(cstored))-nlett), "--> Joining partitions\n", sep="")
            flagbreak <- 0
            ll <- ll+1
            discbreaks[ll] <- breaks.tot[i+1]
          } else {
            if(!silent) cat("\n")
            flagbreak <- flagbreak+1
          }
          if (flagbreak==cutjoin) break
        }
        
        if(!is.null(discbreaks)) {
          brk.dyn <- sort(breaks.tot[-match(discbreaks, breaks.tot)])
        } else brk.dyn <- sort(breaks.tot)
        
        if(!silent) cat("Discarded", .lt(discbreaks), ":: Final number", .lt(brk.dyn), "\n")
        ## source("./Rfunctions/SBR_images_functions.R")
        ## if(iii==1) probdist.comparison(progind$Time, rawset, hist, joinx, breaks.tot, distHell, brk.dyn, new=TRUE)
        ## else probdist.comparison(progind$Time, rawset, hist, joinx, breaks.tot, distHell, brk.dyn, new=FALSE)
        breaks.tot <- brk.dyn
        if(.lt(brk.dyn) == 2) break
      }
    }
    if(!silent) cat("End of the dynamic analysis\n")
  } 
  
  ## From here onwards the unique resu.lt of this analysis is just brk.dyn
  
  
  #########&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&###############
  ############################&&&&&& KINETIC &&&&&&&###############################
  #########&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&###############
  
  if(!silent) cat("Analysis of the kinetic annotation...\n")
  parabol <- 2*progind$PI*(cstored-progind$PI)/cstored
  parabol <- replace(parabol, which(parabol==0), min(parabol[-which(parabol==0)]))
  parabol.log <- -log(parabol/cstored)
  cutf <- replace(progind$Cut, which(progind$Cut==0), min(progind$Cut[which(progind$Cut>0)]))
  kin <- -log(cutf / cstored)-parabol.log
  
  ### SMOOTHING Section
  if(avg.opt[1]=="SG") {
    if(!silent) cat("Savitzky Golay smoothing filter with degree", pol.degree,"\n")
    kin.mv <- savitzkyGolay(kin, p=pol.degree, w=wsize, m=0)
  } else {
    if(!silent) cat("Moving average smoothing filter\n")
    kin.mv <- movav(kin, w=wsize)
  }
  kin.mv <- c(kin[c(1:((wsize-1)/2))], kin.mv, tail(kin,(wsize-1)/2))
  
  ### MAXIMA search Section
  max.mv <- which(peaks(kin.mv, dpeaks.kin, strict=FALSE))
  
  if(.lt(max.mv)!=0) {
    ### First Cleaning on mv:: check separation between consecutive max
    max.mv.tmp <- max.mv
    rif <- max.mv[.lt(max.mv)]
    if(.lt(max.mv)>1){
      for (i in (.lt(max.mv)-1):1) {
        if (rif-max.mv[i]<dpeaks.kin/2) {
          max.mv.tmp <- max.mv.tmp[-match(max.mv[i],max.mv.tmp)] ##Remove the smallest
        } else rif <- max.mv[i]
      }
    }
    max.mv <- max.mv.tmp
    
    ### Second Cleaning Attempt on mv
    adj.mv <- NULL
    adj.mv[1] <- which.min(kin.mv[1:max.mv[1]])
    for (i in seq_along(max.mv) ) {
      if(.lt(max.mv)==1) next
      if (i==.lt(max.mv)) adj.mv[i+1] <- which.min(kin.mv[round(max.mv[i]):cstored])+round(max.mv[i])-1
      else adj.mv[i+1] <- which.min(kin.mv[round(max.mv[i]):round(max.mv[i+1])])+round(max.mv[i])-1 
      amax <- mean( c(kin.mv[max.mv[i]]-kin.mv[adj.mv[i]], kin.mv[max.mv[i]]-kin.mv[adj.mv[i+1]]) ) ## Lts of adjacent (closest) vertical bars
      if (i==1) amin <- kin.mv[max.mv[i+1]]-kin.mv[adj.mv[i+1]]
      else if (i==.lt(max.mv)) amin <- kin.mv[max.mv[i-1]]-kin.mv[adj.mv[i]]
      else amin <- mean(c(kin.mv[max.mv[i-1]]-kin.mv[adj.mv[i]], kin.mv[max.mv[i+1]]-kin.mv[adj.mv[i+1]]) ) 
      ## print(paste("Evaluating", i, max.mv[i], "with ratio", 100*amax/amin, "%"))
      if ( amax/amin < thr.ratio) {
        ## print(paste("Excluding", i, max.mv[i]))
        max.mv.tmp <- max.mv.tmp[-which(max.mv.tmp==max.mv[i])]
      }
    }
    brk.kin <- round(max.mv.tmp)
  } else {
    brk.kin <- NULL
  }
  
  #########&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&###############
  ############################&&&&&& MATCHING/MERGING &&&&&&&#######################
  #########&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&###############
  if(!only.kin) {
    
    if(!silent) {
      if(match) cat("Matching the results...\n")
      else  cat("Merging the results...\n")
    } 
    breaks <- NULL
    
    vkin <- brk.kin
    if(.lt(vkin)!=0) vkin <- refine(vkin)
    vdyn <- brk.dyn
    vdyn <- vdyn[which(vdyn != 1)]
    vdyn <- vdyn[which(vdyn != cstored)]
    if(.lt(vdyn)!=0) {
      vdyn <- refine(vdyn)
    }
    if(!silent) cat("Number of dynamic and kinetic breaks are respectively",.lt(vdyn), .lt(vkin), "\n")
    ##################################################################################
    ## Matching
    if(.lt(vdyn) !=0 && .lt(vkin!=0) ) {
      ll <- 0
      for (i in seq_along(vkin) ) { 
        dist <- NULL
        ## if(!silent) cat("********************************************\n")
        ## print(paste("Analyzing kin break n", i, " ::: ", vkin[i]))
        set <- c(tail(vdyn[which(vdyn<vkin[i])],1), head(vdyn[which(vdyn>=vkin[i])],1) )
        if (.lt(set)==0) break
        for (j in 1:.lt(set) ) {
          ## Naive criteria (binary decision)
          dist[j] <- abs(set[j]-vkin[i])
          ## if(!silent) cat("Comparing with ", set[j], "Distance is", dist[j],"\n")
          dyn.cand <- set[which.min(abs(set-vkin[i]))]
        }
        if (min(dist)< (sd.kin+sd.dyn)) {
          ## print(paste("Adding", vkin[i]))
          ll <- ll+1
          breaks[ll] <- vkin[i]
          vdyn <- vdyn[-match(dyn.cand,vdyn)]
        }
      }
      if(!silent) {
        cat("Number of matched partitions is",.lt(breaks), "\n")
        if(.lt(vdyn)==0) cat("All the dynamic breaks are matched\n")
        if(.lt(vkin)==.lt(breaks)) cat("All the kinetic breaks are matched\n")
      }
    } else {
      if(!match) {
        breaks <- c(vkin, vdyn)
        match <- TRUE ## To skip next session
        if(.lt(vkin)!=0) only.kin <- TRUE
      } else breaks <- NULL
    }
    #######################################################################
    ## Merging option: adding the remaining partitions 
    if(!match) {
      if(.lt(vdyn)!=0 && .lt(vkin)!=0) {
        ## Looking for residual vdyn too close to any of the vkin (rare)
        dist <- NULL
        ll <- 0
        for(i in 1:.lt(vdyn)) {
          for(j in 1:.lt(vkin)) {
            ll <- ll+1
            dist[ll] <- abs(vdyn[i]-vkin[j])
          }
        }
        ## Identifying and removing them from vdyn vector
        dist.mtx <- matrix(dist<sd.kin, nrow=.lt(vdyn), ncol=.lt(vkin), byrow=TRUE)
        if(any(dist.mtx)) {
          near <- unique(which(dist.mtx==TRUE, arr.ind=TRUE)[,1])
          vdyn <- vdyn[-near]
        }
        ## Merging residual vdyn with breaks and vkin
        brk.mtc <- breaks
        breaks <- sort(unique(c(breaks,vdyn,vkin)))
      } else {
        brk.mtc <- breaks
        breaks <- sort(unique(c(breaks,vkin)))
      }
    }
  } else if(.lt(brk.kin)!=0) breaks <- refine(brk.kin) ## If only.kin==TRUE
  else breaks <- NULL
  
  
  ######################################################################
  ## OUTPUT section
  #####################################################################
  if(is.null(breaks)) {
    if(!silent) cat("NO barriers have been found")
    seq.st <- rep(1, nrow(progind))
    if(out.file) {
      output.match <- data.frame(PI=progind$PI, Time=progind$Time, State=seq.st)
      if(is.character(data)) {
        fileout <- gsub("PROGIDX", "BASINS", strsplit(data,"/",fixed=T)[[1]][.lt(strsplit(data,"/",fixed=T)[[1]])])
        if(!silent) cat("Writing", fileout, "...\n")
        fwrite(output.match, file=fileout, sep='\t', row.names=FALSE, col.names=FALSE)
      } else {
        if(!silent) cat(paste0("Writing BASINS_", as.character(progind$Time[1]), ".dat\n"))
        fwrite(output.match, file=paste("./BASINS_", as.character(progind$Time[1]), ".dat", sep=""), sep='\t', row.names=FALSE, col.names=FALSE)
      }
    }
  } else {
    if(!silent) cat("Number of states is", .lt(breaks)+1, "\n")
    if(!silent) cat(breaks, "\n")
    vec <- sort(breaks)
    seq.st <- NULL
    for (i in 1:(.lt(vec)+1)) {
      if (i==1) ib <- 0
      else ib <- vec[i-1]
      if (i==.lt(vec)+1) fb <- cstored
      else fb <- vec[i]
      ## if(!silent) cat(i,ib,fb,fb-ib,"\n")
      seq.st <- c(seq.st,rep(i,fb-ib))
    }
    if(out.file) {
      output.match <- data.frame(PI=progind$PI, Time=progind$Time, State=seq.st)
      if(is.character(data)) {
        fileout <- gsub("PROGIDX", "BASINS", strsplit(data,"/",fixed=T)[[1]][.lt(strsplit(data,"/",fixed=T)[[1]])])
        if(!silent) cat("Writing", fileout, "...\n")
        fwrite(output.match, file=fileout, sep='\t', row.names=FALSE, col.names=FALSE)
      } else {
        if(!silent) cat(paste0("Writing BASINS_", as.character(progind$Time[1]), ".dat\n"))
        fwrite(output.match, file=paste("./BASINS_", as.character(progind$Time[1]), ".dat", sep=""), sep='\t', row.names=FALSE, col.names=FALSE)
      }
    }
  }    
  
  if(is.null(breaks)){
    tab.st <- data.frame(n_cl=1, start=1, end=cstored, .lt=cstored, type=1)
  } else {
    tab.st <- data.frame(n_cl=c(1:(.lt(breaks)+1)), start=c(1,vec+1), end=c(vec, cstored), .lt=diff(c(0,vec, cstored)), type=c(rep(NaN,.lt(breaks)),1))
    if(only.kin) tab.st$type <- c(rep(3, .lt(breaks)), 1)
        else if(match) tab.st$type <- 1
        else {
            if(.lt(brk.mtc)>0) tab.st$type[match(brk.mtc,tab.st$end)] <- 1
            if(.lt(vdyn)>0) tab.st$type[match(vdyn, tab.st$end)] <- 2
            tab.st$type[match(vkin[-match(brk.mtc, vkin)], tab.st$end)] <- 3
        }
    }

  
  #######################################################################
  #### final calculations for scores - dgarol
  #######################################################################
    browser()
    if(cl.stat && !is.null(breaks)){
      
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
        if(!silent) cat('Calculating the barrier statistics using', n_unif, 'division for the MI ratio. \n')
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
        expand_MI_uni <- array(0, dim = lpi)
        for(nun in n_unif){
          brks_uni <- floor(seq(1, lpi, length.out = nun))
          dhc_uni <- dens_histCounts(progind, breaks = brks_uni, nx = nx, ny = ny) # breaks are the barriers points on x
          uni_cnts <- dhc_uni$counts    
          MI_uni <- sapply(seq(nun), function(j) infotheo::mutinformation(uni_cnts[, j], uni_cnts[, j+1]))
          MI_uni <- .normalize(MI_uni)
          expand_MI_uni <- expand_MI_uni + stats::approx(seq(nun), 1 - MI_uni, n = lpi)$y
        }
        expand_MI_uni <- expand_MI_uni/.lt(n_unif)
        # MI_sbr <- .normalize(MI_sbr, xmax = max(c(MI_sbr, MI_uni)), xmin = min(c(MI_sbr, MI_uni)))
        # MI_uni <- .normalize(MI_uni, xmax = max(c(MI_sbr, MI_uni)), xmin = min(c(MI_sbr, MI_uni)))
        
        # hard wired options for the basin weights
        MI_comb <- cl.stat.MI_comb
        denat <- cl.stat.denat 
        
        # Calculating the kinetical cut and the uniform MI denaturation if necessary
        # xr1 <- c(margin:(cstored-round(cstored*0.99)))
        kin.pl <- .normalize(-log(cutf / cstored)) # xr1 impose a selection of 99% of the snapshots to avoid the initial inf
        
        # normalizing the results to avoid parabolic effects (it depends on the grade of the poly!)
        if(denat){
          kinpl <- .denaturate(yyy = kin.pl, xxx = seq(lpi), polydeg = 7, plotit = F)
          miuni <- .denaturate(yyy = expand_MI_uni, xxx = seq(lpi), polydeg = 12, plotit = F)
          # MI_ratio <- (kin.pl[breaks] + miuni[breaks]) / 2
        }else{
          kinpl <- kin.pl
          miuni <- expand_MI_uni
        }
        
        # select the final barrier_weight
        if(MI_comb == 'mean') MI_ratio <- (MI_sbr + miuni[breaks]) / 2 
        else if(MI_comb == 'multip') MI_ratio <- MI_sbr * miuni[breaks]
        else if(MI_comb == 'kin_ann') MI_ratio <- (kinpl[breaks] + miuni[breaks]) / 2
        else if(MI_comb == 'kin') MI_ratio <- kinpl[breaks]
        else if(MI_comb == 'ann') MI_ratio <- miuni[breaks]
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
          
          # what is the relation with the kinetic curve?
          # plot((kin.pl + expand_MI_uni)/2, type = "l")
          # points(breaks, (kin.pl[breaks] + expand_MI_uni[breaks]) / 2, col = 'darkblue')
        }
      }
        
      if(cl.stat.nBreaks == 0) tobrk <- breaks
      else tobrk <- floor(seq(1, lpi, length.out = cl.stat.nBreaks))
      dhc <- dens_histCounts(progind, breaks = tobrk, nx = nx, ny = ny) # breaks are the barriers points on x
      dens <- dhc$density
      cnts <- dhc$counts
      
      # sliding window MI
      if(cl.stat.nBreaks == 0) nbr <- .lt(breaks)
      else nbr <- cl.stat.nBreaks
      if(cl.stat.wMI) {
        if(!silent) cat('Calculating the slided MI using 10 sub-division of the uniform divisions. \n')
        sl_MI <- sliding_MI(progind = progind, n_breaks = nbr, n_slides = 10, nx = nx, ny = ny)
      }
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
      if(plot.cl.stat){
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
      if(cl.stat.weight.barriers) tab.st <- cbind(tab.st, 'barWeight' = c(-1, MI_sbr*expand_MI_uni[breaks]))
    } else statistics <- NULL
  if(is.null(breaks)) statistic <- FALSE
  if(cl.stat && is.null(breaks)) {
    if(!silent) message('The statistics were impossible to calculate as no barrier was succesfully found. The statistic variable has been set to FALSE for this reason.')
    warning('The statistics were impossible to calculate as no barrier was succesfully found. The statistic variable has been set to FALSE for this reason.')
  }
  
  ####################################################################
  ## PLOT Section - at the end if you do cl.stat
  ####################################################################
  
  if(plot && !is.null(breaks)){
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
    if(cl.stat.denat) { # my simple fit (dgarol)
      den_kin <- .denaturate(-log(cutf / cstored), seq(cstored), polydeg = 7, plotit = FALSE)
      lines(xr1, scale(den_kin[xr1]), lwd = 0.8, col = 'darkblue')
    }
    if(only.kin) abline(v=breaks, lwd=0.7, col="orange4")
    else if(!match) {
      ## Not-matched vdyn
      abline(v=vdyn, lwd=0.4, col="blue")
      ## Not-matched vkin
      abline(v=vkin[-match(brk.mtc, vkin)], lwd=0.4, col="orange4")
      ## Matched breaks
      abline(v=brk.mtc, lwd=0.7, col="black")
    } else abline(v=breaks, lwd=0.7, col="black")
    if(cl.stat.weight.barriers){
      # select_non_border <- -c(1, .lt(breaks))
      select_non_border <- seq(.lt(breaks))
      points(breaks[select_non_border], MI_ratio[select_non_border]*max(yr), pch ='+', col = 'red', cex = 5)
    }
    par(save_par)
  }
    invisible(list(tab.st=tab.st, nbins=c(nx,ny), seq.st=seq.st[order(progind$Time)], statistics = statistics, call=call))
}

