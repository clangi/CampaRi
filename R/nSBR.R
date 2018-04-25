#' @title New identifying basins in the SAPPHIRE plot
#' @description
#'     Extimation of the MI curve on the temporal annotation (dots). Using the inverted MIC in particular we can find the clearest separations between basins
#'     (barriers) by looking at the highest values. To calculate this line we randomly split the PI (x-axis) in randomly chosen uniformly sized pieces and we 
#'     calculate the y-axis histogram to compare using MIC. We invert it and we interpolate all the results.
#'     
#' @param data Name of the PROGIDX_<...>.dat file or Data frame with three columns containing, in order, the (sorted) progress index, the relative time 
#' indices and the cut function. 
#' @param ny Number of bins on the y-axis of the histogram. Finer cuts can detect smaller differencies, but they could stretch to the meaningless when data points
#' are sparse.
#' @param local.cut Logical value indicating whether the localized cut function has to be used (see references). Default value is \code{FALSE}.
#' @param comb_met Method to weight the barriers. One of the following:
#'      \itemize{
#'            \item "\code{MIC}" Use this to weight the barriers with the Maximal Information Coefficient.
#'            \item "\code{MIC_kin}" Combine the MIC score and the kinetic annotation for weighting the barriers.
#'            \item "\code{kin}" Use only kinetic annotation for weighting the barriers.
#'            \item "\code{diff_kin}" Better not to use.
#'            \item "\code{diff}" Better not to use.
#'          }
#' @param unif.splits Integer. This variable defines the number of splits for the uniform sampling. Using only one number, e.g. 40 divisions, 
#' could not be optimal for clogged data. Generally speaking, this value is useful for the MI - ratio weighting of the barriers. It is also possible to insert a 
#' vector of integer and the resulting value is an average of these expanded results. Standard is \code{c(5,10,15,20,25,30,40,50,60)} generally.
#' @param n.cluster Integer. If this value is inserted only the first \code{n.cluster - 1} barriers are shown with dotted lines on the plot. Note that no
#' preselection of the barriers is done in the final output and this must be done again in a second instance. The output contains all the possible barriers (ordered).
#' @param plot A logical value indicating whether to display the SAPPHIRE plot with the resulting partitions or not. The dark red curve indicates the MIC interpolation
#' while black dots indicate the peaks found. Here it is possible to see the effect of the spanning variable. Grey dots in the background are the temporal annotation as
#' it has been ordere by the progress index.
#' @param pk_span Integer. The spanning window that should be the minimum size of the basin. It defines the peak finding step.
#' @param return_plot Logical. This can return the plot object in the return list.
#' @param silent A logical value indicating whether the function has to remain silent or not. Default value is \code{FALSE}.
#' @param ...
#'      \itemize{
#'          \item "\code{time.series}" File name. If specified, it substitutes the time series of the PROGIDX_<..> file with the one provided by the file.
#'          \item "\code{data.out.it}" Should I output the data.frame in input?
#'      }
#'      
#' @return A list containing
#' returning_list <- list('nbins' = c(ny, ny), 'barriers' = rank.pk, 'call' = call, 'data.out' = data.out, 'plot' = ggp)
#'       \itemize{
#'         \item "\code{nbins}" c(ny, ny). This is not influential.
#'         Type=1 means that it is a matched partition, type=2 is only dynamic, type=3 is only kinetic. The type of the last state is always equal to 1.
#'         \item "\code{barriers}" Found barriers (no cut is made with the number of clusters which is only for plotting).
#'         \item "\code{call}" The matched call. 
#'         \item "\code{data.out}" The sapphire file name if it is provided, otherwise the data obj. It is used by \code{\link{score_sapphire}}.
#'         \item "\code{plot}" The plot if you wanted it returned.
#'       }   
#' 
#' @details For details regarding the SAPPHIRE plot, please refer to the relative publications \url{http://www.nature.com/articles/srep06264}. 
#' Main documentation of the original campari software \url{http://campari.sourceforge.net/documentation.html}.
#' 
#' @importFrom data.table fread fwrite
#' @importFrom RColorBrewer brewer.pal
#' @importFrom minerva mine
#' @importFrom gplots hist2d
#' @importFrom infotheo mutinformation
#' @import ggplot2
# @importFrom TransferEntropy computeTE
# @importFrom RcppRoll roll_mean
#' 
#' @export nSBR

nSBR <- function(data, ny, local.cut=FALSE, n.cluster=NULL, comb_met=c('MIC'), 
                 unif.splits = NULL, pk_span = NULL, plot=FALSE, silent=FALSE, return_plot = FALSE, ...) {
  call <- match.call()
  
  if(!is.character(data) && !is.data.frame(data)) stop("data must be a string or a data frame")
  if(is.character(data) && (!all(grepl("PROGIDX", data)) && !all(grepl("REPIX", data)))) stop("Please provide a data name starting with 'PROGIDX' or 'REPIX'" )
  if(!is.character(comb_met[1])) stop('comb_met must be a string')
  if(!(comb_met[1] %in% c('MIC', 'MIC_kin', 'kin', 'diff', 'diff_kin'))) stop("comb_met must be in c('MIC', 'MIC_kin', 'kin', 'diff', 'diff_kin')")
  if(!is.null(pk_span) && !.isSingleInteger(pk_span)) stop('pk_span must be a single integer')
  if(!is.null(n.cluster) && !.isSingleInteger(n.cluster)) stop('pk_span must be a single integer')
  if((ny %% 1) != 0) stop("ny must be an integer")
  if(ny < 7) stop('ny must be > 7')
  if(!is.logical(local.cut)) stop("local.cut must be a logical value")
  if(!is.logical(return_plot)) stop("return_plot must be a logical value")
  if(!is.logical(plot)) stop("plot must be a logical value")
  if(!is.logical(silent)) stop("plot must be a logical value")
  if(!is.null(unif.splits)) stopifnot(all(sapply(unif.splits, function(x) x%%1) == 0))
  # Extra arguments checks
  input.args <- list(...)
  avail.extra.arg <- c("time.series",
                       # 'cs.entropy', # DEPRECATED
                       # 'cs.stft', 'cs.TE', 'cs.KL', 'cs.wMI',
                       # 'cs.nUni', 'cs.nBreaks', 'cs.MI_comb', 'cs.denat', 'cs.denat.MI', 
                       'dbg_nSBR', 'data.out.it')
  
  if(!is.null(names(input.args)) && any(!(names(input.args) %in% avail.extra.arg))){
    warning('There is a probable mispelling in one of the inserted variables. Please check the available extra input arguments.')
    if(!silent) cat('!!!!!!!!! We found the following variables without a father (not between our extra input arguments) !!!!!!!!!!\n')
    if(!silent) cat(names(input.args)[!(names(input.args) %in% avail.extra.arg)], '\n')
  }
  
  # dbg_nSBR for stats - dgarol
  if("dbg_nSBR" %in% names(input.args)) { # dgarol
    dbg_nSBR <- input.args$dbg_nSBR
    stopifnot(is.logical(dbg_nSBR))
  } else dbg_nSBR <- FALSE
  
  
  if("data.out.it" %in% names(input.args)) { # dgarol
    data.out.it <- input.args$data.out.it
    stopifnot(is.logical(data.out.it))
  } else data.out.it <- FALSE
  
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
  
  # Making the kinetic notation great again
  parabol <- 2*progind$PI*(cstored-progind$PI)/cstored
  parabol <- replace(parabol, which(parabol==0), min(parabol[-which(parabol==0)]))
  parabol.log <- -log(parabol/cstored)
  cutf <- replace(progind$Cut, which(progind$Cut==0), min(progind$Cut[which(progind$Cut>0)]))
  kin <- -log(cutf / cstored)-parabol.log
  kin[kin < 0] <- 0
  kin <- .normalize(kin)
  
  if(dbg_nSBR) browser()

  # using the division to calculate the densities and the counts
  lpi <- .lt(progind$PI)
  
  # calculating the weights for the barriers 
  if(is.null(unif.splits)) n_unif <- c(5,10,15,20,25,30,40,50,60)
  else n_unif <- unif.splits
  stopifnot(n_unif > 0, n_unif < lpi/ 2)
  
  
  # Uniform division based mutual information
  # n_unif <- round(lpi/50) 
  # n_unif <- 10
  # n_unif <- seq(50, 350, 15)
  if(!silent) cat('Calculating the barrier statistics using ') 
  e_MI_uni <- array(0, dim = lpi)
  e_MIC_uni <- array(0, dim = lpi)
  e_MAS_uni <- array(0, dim = lpi)
  e_MEV_uni <- array(0, dim = lpi)
  e_MCN_uni <- array(0, dim = lpi)
  e_MICR2_uni <- array(0, dim = lpi)
  for(nun in n_unif){
    if(!silent) cat(nun, ' ')
    brks_uni <- floor(seq(1, lpi, length.out = nun))[-c(1,nun)]
    dhc_uni <- .dens_histCounts(progind, breaks = brks_uni, nx = ny, ny = ny) # breaks are the barriers points on x
    uni_cnts <- dhc_uni$counts    
    MI_uni <- sapply(seq(nun-2), function(j) infotheo::mutinformation(uni_cnts[, j], uni_cnts[, j+1]))
    minerva.out <- sapply(seq(nun-2), function(j) minerva::mine(x = uni_cnts[, j], y = uni_cnts[, j+1]))
    MI_uni <- .normalize(MI_uni)
    MIC_uni <- .normalize(unlist(minerva.out[1,]))
    MAS_uni <- .normalize(unlist(minerva.out[2,]))
    MEV_uni <- .normalize(unlist(minerva.out[3,]))
    MCN_uni <- .normalize(unlist(minerva.out[4,]))
    MICR2_uni <- .normalize(unlist(minerva.out[5,]))
    e_MI_uni <- e_MI_uni + stats::approx(seq(nun), c(0, 1 - MI_uni, 0), n = lpi)$y
    e_MIC_uni <- e_MIC_uni + stats::approx(seq(nun), c(0, 1 - MIC_uni, 0), n = lpi)$y
    e_MAS_uni <- e_MAS_uni + stats::approx(seq(nun), c(0, 1 - MAS_uni, 0), n = lpi)$y
    e_MEV_uni <- e_MEV_uni + stats::approx(seq(nun), c(0, 1 - MEV_uni, 0), n = lpi)$y
    e_MCN_uni <- e_MCN_uni + stats::approx(seq(nun), c(0, 1 - MCN_uni, 0), n = lpi)$y
    e_MICR2_uni <- e_MICR2_uni + stats::approx(seq(nun), c(0, MICR2_uni, 0), n = lpi)$y
  }
  if(!silent) cat('divisions for the MI ratio. \n')
  e_MI_uni <- e_MI_uni/.lt(n_unif)
  e_MIC_uni <- e_MIC_uni/.lt(n_unif)
  e_MAS_uni <- e_MAS_uni/.lt(n_unif)
  e_MEV_uni <- e_MEV_uni/.lt(n_unif)
  e_MCN_uni <- e_MCN_uni/.lt(n_unif)
  e_MICR2_uni <- e_MICR2_uni/.lt(n_unif)
  
  
  convMIC_cov <- .normalize(convolve(x = 1:lpi, y = e_MIC_uni))
  convMIC <- .cut0(e_MIC_uni - convMIC_cov)
  df.plot <- cbind('PI' = 1:lpi,
                   'MI' = e_MI_uni, 
                   'MIC' = e_MIC_uni, 
                   'MAS' = e_MAS_uni, 
                   'MEV' = e_MEV_uni, 
                   'MCN' = e_MCN_uni, 
                   'MICR2' = e_MICR2_uni,
                   'PI_points_x' = progind$PI, 
                   'PI_points_y' = progind$Time/lpi,
                   'kin' = kin,
                   'convDiff' = convMIC,
                   'conv' = convMIC_cov,
                   'diff' = .normalize(c(0,diff(e_MIC_uni))),
                   'diffcut' = .normalize(.cut0(c(0,diff(e_MIC_uni)))),
                   'ddiffcut' = .normalize(.cut0(-c(0,diff(.normalize(.cut0(c(0,diff(e_MIC_uni)))))))),
                   'ddcutkin' = .normalize(.cut0(-c(0,diff(.normalize(.cut0(c(0,diff(e_MIC_uni + kin)))))))) # wtf
                   )
  
  # setting up the stats
  df.plot <- as.data.frame(df.plot)
  if(comb_met[1] == 'MIC') wtp <- df.plot$MIC
  else if(comb_met[1] == 'kin') wtp <- kin
  else if(comb_met[1] == 'MIC_kin') wtp <- (df.plot$MIC + kin )/ 2 # wtp what to plot (the ending peak finding subject)
  else if(comb_met[1] == 'diff_kin') wtp <- (df.plot$ddiffcut + kin )/ 2
  else if(comb_met[1] == 'diff') wtp <- df.plot$ddiffcut
  else stop('comb_met[1] not in permitted ones. This here should not happen.')
  
  if(is.null(pk_span)) span <- round(lpi/50) # number of snapshots in which only one peak
  else span <- pk_span
  stopifnot(span > 3, span < lpi/2)
  
  # finding peaks and ordering 
  pky <- array(0, lpi)
  where.pk <- peaks(wtp, span = span)
  pky <- as.numeric(where.pk)[where.pk]
  pkx <- 1:lpi
  pkx <- pkx[where.pk]
  rank.pk <- order(wtp[where.pk], decreasing = T)
  rank.pk <- pkx[rank.pk]
  df.pp <- cbind(pkx, pky)
  
  if(is.null(n.cluster)) nclu <- 15
  else nclu <- n.cluster
  nba <- nclu - 1 # cluster and barrier number 
  rnk_ts <- rank.pk[1:(min(max(15, nclu), .lt(rank.pk)))] # rank to show (number of peaks to show)
  
  dhc <- .dens_histCounts(progind, breaks = rank.pk, nx = ny, ny = ny) # breaks are the barriers points on x
  dens <- dhc$density
  ent <- .normalize(apply(dens, 2, .myShEn))
  
  
  # if(nrow(df.plot) > 50000) df.plot <- df.plot[unique(round(seq(1, lpi, length.out = 50000))),] # it is not clear for the barriers
  ggp <- ggplot(data = df.plot, mapping = aes(x = PI)) + theme_classic() + ylab('UniDiv score') +
         geom_point(mapping = aes(x = PI_points_x, y = PI_points_y), size = 0.8, col = 'grey') +
         # geom_line(mapping = aes(y = MI), col = 'blue') +
         geom_line(mapping = aes(y = MIC), color = 'darkred') +
         # geom_line(mapping = aes(y = MAS)) +
         # geom_line(mapping = aes(y = MEV)) +
         # geom_line(mapping = aes(y = MCN)) +
         # geom_line(mapping = aes(y = convDiff), col = 'purple') +
         # geom_line(mapping = aes(y = MICR2), col = 'green') +
         # geom_line(mapping = aes(y = ddiffcut, col = 'ddiffcut')) +
         # geom_line(mapping = aes(y = ddcutkin), col = 'yellow') +
         # geom_line(mapping = aes(y = kin, col = 'kin'))  + 
         geom_point(data = df.pp, mapping = aes(y = pky, x = pkx)) +
         geom_vline(xintercept = rank.pk[1:nba], col = 'black', size = 1, linetype="dotted") +
         geom_point(data = data.frame(a = rnk_ts, b = wtp[rnk_ts]), 
                   mapping = aes(x = a, y = b),
                   shape = 3, col = 'darkblue', size = 5, stroke = 2.5) + xlab('Progress Index') + ylab(paste0('I', comb_met))
  #+ #deeppink3
    # scale_colour_manual(name="Statistic", values=c("MIC"="darkred","ddiffcut"="lightgreen","kin"="black"))
  if(plot) print(ggp)
  
  # microbenchmarking MIC MI
  if(F){
    mbm <- microbenchmark::microbenchmark(
      minerva.out = sapply(seq(nun-2), function(j) minerva::mine(x = uni_cnts[, j], y = uni_cnts[, j+1])),
      MI_uni = sapply(seq(nun-2), function(j) infotheo::mutinformation(uni_cnts[, j], uni_cnts[, j+1])), 
      times = 100
    )
    ggplot2::autoplot(mbm)  
  }

  
  # FINAL OUTPUT
  returning_list <- list('nbins' = c(ny, ny), 'barriers' = rank.pk, 'call' = call, 'data.out' = data.out)
  if(return_plot) returning_list <- list('nbins' = c(ny, ny), 'barriers' = rank.pk, 'call' = call, 'data.out' = data.out, 'plot' = ggp)
  invisible(returning_list)
}


#########################################################
#          Functions ext
#########################################################
# myTE <- function(x, y, emb, sym = T, sd = 0.001){
#   nx <- x + stats::rnorm(n = .lt(x), mean = 0, sd = sd)
#   ny <- y + stats::rnorm(n = .lt(x), mean = 0, sd = sd)
#   xy <- TransferEntropy::computeTE(nx, ny, embedding = emb, k = 2, safetyCheck = T)
#   yx <- TransferEntropy::computeTE(ny, nx, embedding = emb, k = 2, safetyCheck = T)
#   if(sym) return((unlist(xy) + unlist(yx))/2)
#   else return(unlist(xy))
# }
.myShEn <- function(x){
  x2 <- replace(x = x, list = which(x == 0), values = 1)
  return(-sum(x2*log(x2)))
}

# myKL <- function(x, y, sym = T){
#   xn <- x[-which(x==0 | y==0)]
#   yn <- y[-which(x==0 | y==0)]
#   if(sym) return((sum(xn*log(xn/yn)) + sum(yn*log(yn/xn)))/2)
#   else return(sum(xn*log(xn/yn)))
# }

# creates the densities and the counts over bins (only ny is relevant) - doing hist2d
.dens_histCounts <- function(progind, breaks, nx, ny = nx){
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
# sliding_MI <- function(progind, n_breaks = 50, n_slides = 10, nx, ny = nx){
#   lpi <- .lt(progind$PI)
#   brks <- floor(seq(1, lpi, length.out = n_breaks))
#   span_MI <- sapply(X = floor(seq(1, brks[2]-1, length.out = n_slides)), FUN = function(x){
#     br <- floor(seq(x, lpi, by = brks[2]))
#     dhc <- .dens_histCounts(progind, breaks = br, nx = nx, ny = ny)
#     mi.out <- list(.normalize(sapply(seq(.lt(br)), function(j) infotheo::mutinformation(dhc$counts[, j], dhc$counts[, j+1]))))
#     return(mi.out)
#   })
#   return(c(t(do.call(cbind, span_MI))))
# }

# --------------------------- functions - fra
# .scale <- function(x) {
#   scl <- (max(progind$Time)-1)/(max(x)-min(x))
#   return(1+(x-min(x))*scl)
# }
# ------------------------ end functions

# ------------------------------------------------------------------------------------------------------------------------  
# DENATURATION - useless
# 
# # hard wired options for the basin weights
# MI_comb <- cs.MI_comb
# denat <- cs.denat 
# denat.MI <- cs.denat.MI
# # Calculating the kinetical cut and the uniform MI denaturation if necessary
# # xr1 <- c(margin:(cstored-round(cstored*0.99)))
# kin.pl <- .normalize(-log(cutf / cstored)) # xr1 impose a selection of 99% of the snapshots to avoid the initial inf
# 
# # normalizing the results to avoid parabolic effects (it depends on the grade of the poly!)
# if(!is.null(denat) && denat == 'poly_interpolation'){
#   if(!is.null(denat.MI)){
#     if(denat.MI != -1) kinpl <- .denaturate(yyy = kin.pl, xxx = seq(lpi), polydeg = denat.MI, plotit = F)
#     else kinpl <- .denaturate(yyy = kin.pl, xxx = seq(lpi), polydeg = 7, plotit = F)
#     if(denat.MI != -1) miuni <- .denaturate(yyy = e_MI_uni, xxx = seq(lpi), polydeg = denat.MI, plotit = F)
#     else miuni <- 1 - e_MI_uni
#     # MI_ratio <- (kin.pl[breaks] + miuni[breaks]) / 2
#   }else{
#     kinpl <- .denaturate(yyy = kin.pl, xxx = seq(lpi), polydeg = 7, plotit = F)
#     miuni <- e_MI_uni
#   }
# }else if(!is.null(denat) && denat == 'process_subtraction'){
#   kinpl <- .normalize(kin)
#   if(!is.null(denat.MI)) {
#     if(denat.MI != -1) miuni <- .denaturate(yyy = e_MI_uni, xxx = seq(lpi), polydeg = denat.MI, plotit = F)
#     else miuni <- 1 - e_MI_uni
#   }else miuni <- e_MI_uni
# }else{
#   kinpl <- kin.pl
#   miuni <- e_MI_uni
# }
# 

# ------------------------------------------------------------------------------------------------------------------------  
# MIXING THE ANNOTATION AND THE MI
# select the final barrier_weight
# if(MI_comb == 'mean') MI_ratio <- (MI_sbr + miuni[breaks]) / 2 
# else if(MI_comb == 'multip') MI_ratio <- MI_sbr * miuni[breaks]
# else if(MI_comb == 'kin_MI') MI_ratio <- (kinpl[breaks] + miuni[breaks]) / 2
# else if(MI_comb == 'kin') MI_ratio <- kinpl[breaks]
# else if(MI_comb == 'MI') MI_ratio <- miuni[breaks]
# if(plot){
#   
#   # plot the found histogram Mutual Information
#   plot(miuni, type = 'l', ylim = c(0,1), xlim = c(0,lpi))                             # from uniform sampling (average of different x-divisions)
#   points(breaks, MI_sbr, pch ='o', col = 'green', cex = 2)                                    # from SBR division
#   points(breaks, miuni[breaks], pch = 20, col = 'darkblue', cex = 2)                  # height of the point on barriers using uniform sampling
#   
#   # selecting the barriers and 
#   which_barrier <- breaks > lpi/n_unif[.lt(n_unif)] & breaks < lpi - lpi/n_unif[.lt(n_unif)]
#   points(breaks[which_barrier], MI_ratio[which_barrier], pch ='+', col = 'red', cex = 2)
#   points(progind$PI, progind$Time/lpi, pch='.', cex = 2)
#   lines(kinpl, col = 'darkred')
#   # if(cs.wMI) lines(seq(1, lpi, length.out = .lt(sl_MI)), sl_MI, col = 'darkgreen')
#   
#   # what is the relation with the kinetic curve?
#   # plot((kin.pl + e_MI_uni)/2, type = "l")
#   # points(breaks, (kin.pl[breaks] + e_MI_uni[breaks]) / 2, col = 'darkblue')
# }
# 