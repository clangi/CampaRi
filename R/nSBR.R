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
#'          \item "\code{cs.nBreaks}" Integer. This variable defines the number of splits for the analysis. If set to 0 it will use the breaks found in the SBR.
#'          \item "\code{cs.denat}" This value can be set to \code{"process_subtraction"} or \code{"poly_interpolation"} and it is defining the removal of 
#'          parabolic artifacts in the kinetic trace. The polynomial fit is by default 7 and 12 in degree for the kinetic annotation and uniform MI curves. The process 
#'          option istead is referring to the simulated process which is behind the kin ann. This option can be set to \code{TRUE} for using the process way.
#'          \item "\code{cs.denat.MI}" Integer. If NULL its default is 7 for poly_interpolation of the kinetic ann but nothing is done on the MI score. 
#'          \item "\code{cs.MI_comb}" This value is the representative of how the various calculation in barrier weighting are combined. In particular the options 
#'          are the following:
#'          \itemize{
#'            \item "\code{mean}" Mean between the MI based on the SBR and the MI based on uniform subdivisions
#'            \item "\code{multip}" Multiplication between the MI based on the SBR and the MI based on uniform subdivisions
#'            \item "\code{kin_MI}" Mean between the kinetic annotaton and the MI based on the SBR and the MI based on uniform subdivisions
#'            \item "\code{kin}" Only kinetic annotation value on the barrier
#'            \item "\code{MI}" Only the MI based on the SBR and the MI based on uniform subdivisions
#'          }
#'          \item "\code{cs.nUni}" Integer. This variable defines the number of splits for the uniform sampling. If not set the algorithm will use 40 divisions which 
#'          could not be optimal for clogged data. Generally speaking, this value is useful for the MI - ratio weighting of the barriers. It is also possible to insert a 
#'          vector of integer and the resulting value is an average of these expanded results. We advice to use \code{c(5,10,15,20,25,30,40,50)} generally.
#'          \item "\code{cs.entropy}" Logical. Calculate Shannon entropy.
#'          \item "\code{cs.stft}" Logical. Calculate short time fourier transform.
#'          \item "\code{cs.TE}" Logical. Calculate symmetric Transfer Entropy.
#'          \item "\code{cs.KL}" Logical. Calculate symmetric and non-symmetric Kullback-Leibler divergence.
#'          \item "\code{cs.wMI}" Logical. Use the number of breaks (\code{cs.nBreaks}) to find the Mutual Information for 10 slided divisions.
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

nSBR <- function(data, ny, local.cut=FALSE, comb_met= c('MIC', 'MIC_kin', 'kin'), pk_span = NULL, plot=FALSE, silent=FALSE,...) {
  call <- match.call()
  
  if(!is.character(data) && !is.data.frame(data)) stop("data must be a string or a data frame")
  if(is.character(data) && (!all(grepl("PROGIDX", data)) && !all(grepl("REPIX", data)))) stop("Please provide a data name starting with 'PROGIDX' or 'REPIX'" )
  if((ny %% 1) != 0) stop("ny must be an integer")
  if(ny < 7) stop('ny must be > 7')
  if(!is.logical(local.cut)) stop("local.cut must be a logical value")
  if(!is.logical(plot)) stop("plot must be a logical value")
  
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
  stopifnot(cs.nBreaks >= 0, cs.nBreaks <= floor(lpi/2))
  
  # calculating the weights for the barriers 
  if(is.null(cs.nUni)) n_unif <- 40
  else n_unif <- cs.nUni
  stopifnot(n_unif > 0, n_unif < lpi/ 2)
  
  
  # Uniform division based mutual information
  # n_unif <- round(lpi/50) 
  # n_unif <- 10
  n_unif <- c(5,10,15,20,25,30,40,50,60)
  n_unif <- seq(50, 350, 15)
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
    dhc_uni <- dens_histCounts(progind, breaks = brks_uni, nx = ny, ny = ny) # breaks are the barriers points on x
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
  # wtp <- (df.plot$ddiffcut + kin )/ 2 # wtp what to plot (the ending peak finding subject)
  # wtp <- df.plot$ddiffcut
  wtp <- kin
  nclu <- 4; nba <- nclu - 1 # cluster and barrier number 
  span <- 5000 # number of snapshots in which only one peak
  
  # finding peaks and ordering 
  pky <- array(0, lpi)
  where.pk <- peaks(wtp, span = span)
  pky <- as.numeric(where.pk)[where.pk]
  pkx <- 1:lpi
  pkx <- pkx[where.pk]
  rank.pk <- order(wtp[where.pk], decreasing = T)
  rank.pk <- pkx[rank.pk]
  df.pp <- cbind(pkx, pky)
  rnk_ts <- rank.pk[1:(min(15, .lt(rank.pk)))] # rank to show (number of peaks to show)
  
  # if(nrow(df.plot) > 50000) df.plot <- df.plot[unique(round(seq(1, lpi, length.out = 50000))),] # it is not clear for the barriers
  ggp <- ggplot(data = df.plot, mapping = aes(x = PI)) + theme_classic() + ylab('UniDiv score') +
         geom_point(mapping = aes(x = PI_points_x, y = PI_points_y), size = 0.8, col = 'grey') +
         # geom_line(mapping = aes(y = MI), col = 'blue') +
         geom_line(mapping = aes(y = MIC, col = 'MIC')) +
         # geom_line(mapping = aes(y = MAS)) +
         # geom_line(mapping = aes(y = MEV)) +
         # geom_line(mapping = aes(y = MCN)) +
         # geom_line(mapping = aes(y = convDiff), col = 'purple') +
         # geom_line(mapping = aes(y = MICR2), col = 'green') +
         geom_line(mapping = aes(y = ddiffcut, col = 'ddiffcut')) +
         # geom_line(mapping = aes(y = ddcutkin), col = 'yellow') +
         geom_line(mapping = aes(y = kin, col = 'kin'))  + 
         geom_point(data = df.pp, mapping = aes(y = pky, x = pkx)) +
         geom_vline(xintercept = rank.pk[1:nba], col = 'black', size = 1, linetype="dotted") +
         geom_point(data = data.frame(a = rnk_ts, b = wtp[rnk_ts]), 
                   mapping = aes(x = a, y = b),
                   shape = 3, col = 'deeppink3', size = 5, stroke = 4) +
    scale_colour_manual(name="Statistic", values=c("MIC"="darkred","ddiffcut"="lightgreen","kin"="black"))
  print(ggp)
  # microbenchmarking MIC MI
  if(F){
    mbm <- microbenchmark::microbenchmark(
      minerva.out = sapply(seq(nun-2), function(j) minerva::mine(x = uni_cnts[, j], y = uni_cnts[, j+1])),
      MI_uni = sapply(seq(nun-2), function(j) infotheo::mutinformation(uni_cnts[, j], uni_cnts[, j+1])), 
      times = 100
    )
    ggplot2::autoplot(mbm)  
  }
  
  
  
  if(F){
    
    # SBR based mutual information - it needs breaks
    dhc_sbr <- dens_histCounts(progind, breaks = breaks, nx = ny, ny = ny) # breaks are the barriers points on x
    sbr_cnts <- dhc_sbr$counts    
    MI_sbr <- sapply(seq(.lt(breaks)), function(j) infotheo::mutinformation(sbr_cnts[, j], sbr_cnts[, j+1]))
    MI_sbr <- .normalize(MI_sbr)
    
    # sliding window MI
    if(cs.nBreaks == 0) nbr <- .lt(breaks)
    else nbr <- cs.nBreaks
    if(cs.wMI) {
      if(!silent) cat('Calculating the slided MI using 10 sub-division of the uniform divisions. \n')
      sl_MI <- sliding_MI(progind = progind, n_breaks = nbr, n_slides = 10, nx = ny, ny = ny)
    }
    
    if(cs.nBreaks == 0) tobrk <- breaks
    else tobrk <- floor(seq(1, lpi, length.out = cs.nBreaks))
    dhc <- dens_histCounts(progind, breaks = tobrk, nx = ny, ny = ny) # breaks are the barriers points on x
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
    if(cs.stft){
      if(!silent) cat('Calculating the short time Fourier Transform using a window of 30 and an increment of 1. \n')
      short_time_FT <- e1071::stft(sl_MI, win = 30, inc = 1)
      histOfStft <- apply(short_time_FT$values, 1, sum)
      histOfStft <- .normalize(histOfStft)
      if(plot){
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
      if(plot){
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
    if(cs.KL) { distSymKL <- sapply(seq(.lt(tobrk)), function(j) myKL(dens[, j], dens[, j+1], T)); distSymKL <- .normalize(distSymKL) }
    if(cs.KL) { distKL <- sapply(seq(.lt(tobrk)), function(j) myKL(dens[, j], dens[, j+1], F)); distKL <- .normalize(distKL) }
    if(cs.TE) { distSymTE <- sapply(seq(.lt(tobrk)), function(j) myTE(dens[, j], dens[, j+1], emb = 3, sym = T)); distSymTE <- .normalize(distSymTE) }
    
    # defining the centers/ini/end of the splits
    center_cl <- diff(c(1, tobrk, lpi))/2 + c(1, tobrk)
    ini_points <- c(1, tobrk); end_points <- c(1, tobrk) + diff(c(1, tobrk, lpi))
    
    # Plotting the barriers, breaks and statistics on the breaks 
    if(plot && F){ # I will put it back in a moment
      gg <- ggplot() + theme_classic() + xlab('Progress Index') + ylab('Barrier score') +
        geom_line(aes(tobrk, distHell, col = as.factor(1)), size = 0.5) +
        geom_segment(aes(x = ini_points, xend = end_points, y = ent, yend = ent, col = as.factor(2)), size = 4) +
        geom_line(aes(tobrk, distMI, col = as.factor(3)), size = 2) +
        geom_point(aes(progind$PI, progind$Time/lpi), size = 0.1) + 
        geom_vline(aes(xintercept=c(1, tobrk, lpi)), size = 0.1) + # MI and entropy calculations 
        geom_vline(aes(xintercept=c(1, breaks, lpi)), size = 0.5)  # Found barriers from SBR
      lbls <- c('Hellinger', 'Entropy', 'MI')
      if(cs.KL){
        gg <- gg + geom_line(aes(tobrk, distKL, col = as.factor(4)), size = 0.5) +
          geom_line(aes(tobrk, distSymKL, col = as.factor(5)), size = 0.5)
        lbls <- c(lbls, 'KL', 'symKL')
      }
      if(cs.TE){
        gg <- gg + geom_line(aes(tobrk, distSymTE, col = as.factor(6)), size = 0.5)
        lbls <- c(lbls, 'symTE')
      }
      if(cs.wMI){
        gg <- gg + geom_line(aes(seq(1, lpi, length.out = .lt(sl_MI)), sl_MI, col = as.factor(7)), size = 1) + 
          geom_hline(aes(yintercept = mean(sl_MI)))
        lbls <- c(lbls, 'WinMI')
      }
      if(cs.stft){
        gg <- gg + geom_vline(aes(xintercept=c(xpeaks, lpi), col = as.factor(8)), size = 0.2) # no need of new splits... better to score the bsr split (cohesion - MI) 
        lbls <- c(lbls, 'stFT barriers')
      }
      gg <- gg + geom_segment(aes(x = breaks - round(lpi/30), xend = breaks + round(lpi/30), 
                                  y = MI_sbr*e_MI_uni[breaks], yend = MI_sbr*e_MI_uni[breaks]), col = 'red', size = 1.5)
      gg <- gg + scale_color_manual(name = "Method", labels = lbls, values = RColorBrewer::brewer.pal(n = 8, name = 'Dark2')) +
        guides(color = guide_legend(override.aes = list(size=5)))
      print(gg)
    }
    
    # Entropy calculations
    if(cs.entropy){
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
      if(plot){
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
  }
  
  
  # Final assignment for output
  if(F){
    statistics <- NULL
    ele <- 1
    if(cs.nBreaks == 0){ # case in which the statistic comes from SBR splits
      tab.st <- cbind(tab.st, 'Hellinger' =  c(-1, distHell), 'Entropy' = ent, 'MI' = c(-1, distMI))
      if(cs.KL) tab.st <- cbind(tab.st, 'symKL' = c(-1, distSymKL))
      if(cs.KL) tab.st <- cbind(tab.st, 'KL' = c(-1, distKL))
      if(cs.TE) tab.st <- cbind(tab.st, 'symTE' = c(-1, distSymTE))
    }else{ # otherwise
      statistics <- list('Hellinger' =  distHell, 'Entropy' = ent, 'MI' = distMI); ele <- ele + 3
      if(cs.KL) { statistics[[ele]] <- distSymKL; names(ele)[ele] <- 'symKL'; ele <- ele + 1 }
      if(cs.KL) { statistics[[ele]] <- distKL; names(ele)[ele] <- 'KL'; ele <- ele + 1 }
      if(cs.TE) { statistics[[ele]] <- distSymTE; names(ele)[ele] <- 'symTE'; ele <- ele + 1 }
    }
    if(cs.wMI) { statistics[[ele]] <- sl_MI; names(ele)[ele] <- 'winMI'; ele <- ele + 1 }
    if(cs.stft) { statistics[[ele]] <- xpeaks; names(ele)[ele] <- 'stFT'; ele <- ele + 1 }
    if(!silent) cat('Keeping only the barriers which are reasonably far from borders.\n')
    select_non_border <- breaks > lpi/n_unif[.lt(n_unif)] & breaks < lpi - lpi/n_unif[.lt(n_unif)]
    if(!silent) cat('Discarded the following barrier indexes (from left / also on the stats out):', seq(1, .lt(breaks))[!select_non_border], '\n')
    forout <- MI_ratio
    forout[!select_non_border] <- -1
    tab.st <- cbind(tab.st, 'barWeight' = c(-1, forout))
  }

  invisible(list('nbins' = c(ny, ny), 
                 'barriers' = rank.pk, 
                 # 'seq.st' = seq.st[order(progind$Time)], 
                 'call' = call, 
                 'data.out' = data.out))
}




#########################################################
#          Functions ext
#########################################################
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

# --------------------------- functions
.scale <- function(x) {
  scl <- (max(progind$Time)-1)/(max(x)-min(x))
  return(1+(x-min(x))*scl)
}
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