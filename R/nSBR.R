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
#' @param random_picks Insert the number of times the barriers should be randomly recalculated for having a random baseline.
#' @param ann If random_picks is inserted you need to score the the random barriers against an annotation!!
#' @param ...
#'      \itemize{
#'          \item "\code{time.series}" File name. If specified, it substitutes the time series of the PROGIDX_<..> file with the one provided by the file.
#'      }
#'      
#' @return A list containing
#' \itemize{
#'         \item "\code{nbins}" c(ny, ny). This is not influential.
#'         \item "\code{barriers}" Found barriers (no cut is made with the number of clusters which is only for plotting).
#'         \item "\code{plot}" The plot if you wanted it returned.
#'         \item "\code{rnd.picks}" This is an object containing the results of the random.picks, i.e. a list contatining:
#'         \itemize{
#'                 \item "\code{ggp.bar}" ggplot containing the barriers found and the scores found with the random procedure.
#'                 \item "\code{ggp.dist.rnd.bar}"  ggplot containing the estimate (in z-scores) of the distance from random distribution (showed in light blue).
#'                 \item "\code{rnd.scores.nas}" Number of NAs in the randomly picked scores.
#'                 \item "\code{rnd.scores}" data.frame with all the scores found.
#'                 \item "\code{max.rnd.nmi}" score_sapphire obj on the best random pick.
#'                 \item "\code{max.rnd.bar}" data.frame with the best parriers found randomly and their score (xm = barrier points, ym = score).
#'                 \item "\code{ref.score}" score_sapphire obj with the standard best score optimization procedure.
#'                 \item "\code{z.scores}" Found zscore between optimal barrier and random picks distribution.
#'                 \item "\code{std.n.bar}" Total number of barriers found without imposing the number of clusters in optimization procedure.
#'                 \item "\code{pval_zsc}" The p-value calculated using the 
#'                 \item "\code{pval_95per}" Total number of barriers found without imposing the number of clusters in optimization procedure.
#'         }
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

nSBR <- function(data, ny, local.cut = FALSE, n.cluster = NULL, comb_met = c('MIC'),                        # fundamental vars
                 unif.splits = NULL, pk_span = NULL,                                                        # algorithm details
                 plot = FALSE, silent = FALSE, return_plot = FALSE,                                         # plots and prints
                 random_picks = NULL, ann = NULL,                                                           # randomization of the barriers and comparison
                 ...) { 
  
  # Standard input checks  
  allowed_mets <- c('MI', 'MIC', 'MAS', 'MEV', 'MCN', 'MICR2', 'MIC_kin', 'kin', 'diff', 'convDiff', 'conv')
  if(!is.character(data) && !is.data.frame(data)) stop("data must be a string or a data frame")
  if(is.character(data) && (!all(grepl("PROGIDX", data)) && !all(grepl("REPIX", data)))) stop("Please provide a data name starting with 'PROGIDX' or 'REPIX'" )
  if(!is.character(comb_met[1])) stop('comb_met must be a string')
  if(!(comb_met[1] %in% allowed_mets)) stop(paste0("comb_met must be in ", paste(allowed_mets, collapse = " ")))
  if(!is.null(pk_span) && !.isSingleInteger(pk_span)) stop('pk_span must be a single integer')
  if(!is.null(n.cluster) && !.isSingleInteger(n.cluster)) stop('n.cluster must be a single integer')
  if(!is.null(random_picks) && !.isSingleInteger(random_picks)) stop('random_picks must be a single integer') 
  if(!is.null(random_picks) && is.null(ann)) stop('To test the random pick of variables we need to know the real annotation. Please insert ann variable.')
  
  if((ny %% 1) != 0) stop("ny must be an integer")
  if(ny < 7) stop('ny must be > 7')
  if(!is.logical(local.cut)) stop("local.cut must be a logical value")
  if(!is.logical(return_plot)) stop("return_plot must be a logical value")
  if(!is.logical(plot)) stop("plot must be a logical value")
  if(!is.logical(silent)) stop("plot must be a logical value")
  if(!is.null(unif.splits)) stopifnot(all(sapply(unif.splits, function(x) x%%1) == 0))
  
  # Extra arguments checks
  input.args <- list(...)
  avail.extra.arg <- c("time.series", 'dbg_nSBR')
  
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
  
  lpi <- .lt(progind$PI)
  
  # Making the kinetic notation great again
  parabol <- 2*progind$PI*(lpi-progind$PI)/lpi
  parabol <- replace(parabol, which(parabol==0), min(parabol[-which(parabol==0)]))
  parabol.log <- -log(parabol/lpi)
  cutf <- replace(progind$Cut, which(progind$Cut==0), min(progind$Cut[which(progind$Cut>0)]))
  kin <- -log(cutf / lpi)-parabol.log
  kin[kin < 0] <- 0
  kin <- .normalize(kin)
  
  # Calculating the weights for the barriers 
  if(is.null(unif.splits)) n_unif <- c(5,10,15,20,25,30,40,50,60)
  else n_unif <- unif.splits
  stopifnot(n_unif > 0, n_unif < lpi/ 2)

  # Calculating the curves
  if(!silent) cat('Calculating the barrier statistics using ') 
  if('MI' %in% comb_met) e_MI_uni <- array(0, dim = lpi)
  if(any(comb_met %in% allowed_mets)){
    e_MIC_uni <- array(0, dim = lpi)
    e_MAS_uni <- array(0, dim = lpi)
    e_MEV_uni <- array(0, dim = lpi)
    e_MCN_uni <- array(0, dim = lpi)
    e_MICR2_uni <- array(0, dim = lpi)
  }
  
  # Main loop on the divisions
  for(nun in n_unif){
    if(!silent) cat(nun, ' ')
    
    # Calculating the histogram
    brks_uni <- floor(seq(1, lpi, length.out = nun))[-c(1,nun)]
    dhc_uni <- .dens_histCounts(progind, breaks = brks_uni, nx = ny, ny = ny) # breaks are the barriers points on x
    uni_cnts <- dhc_uni$counts
    
    if('MI' %in% comb_met) MI_uni <- sapply(seq(nun-2), function(j) infotheo::mutinformation(uni_cnts[, j], uni_cnts[, j+1]))
    if(any(comb_met %in% allowed_mets)) minerva.out <- sapply(seq(nun-2), function(j) minerva::mine(x = uni_cnts[, j], y = uni_cnts[, j+1]))
    if('MI' %in% comb_met) {
      MI_uni <- .normalize(MI_uni)
      e_MI_uni <- e_MI_uni + stats::approx(seq(nun), c(0, 1 - MI_uni, 0), n = lpi)$y
    }
    if(any(comb_met %in% allowed_mets)){
      MIC_uni <- .normalize(unlist(minerva.out[1,]))
      MAS_uni <- .normalize(unlist(minerva.out[2,]))
      MEV_uni <- .normalize(unlist(minerva.out[3,]))
      MCN_uni <- .normalize(unlist(minerva.out[4,]))
      MICR2_uni <- .normalize(unlist(minerva.out[5,]))
      e_MIC_uni <- e_MIC_uni + stats::approx(seq(nun), c(0, 1 - MIC_uni, 0), n = lpi)$y
      e_MAS_uni <- e_MAS_uni + stats::approx(seq(nun), c(0, 1 - MAS_uni, 0), n = lpi)$y
      e_MEV_uni <- e_MEV_uni + stats::approx(seq(nun), c(0, 1 - MEV_uni, 0), n = lpi)$y
      e_MCN_uni <- e_MCN_uni + stats::approx(seq(nun), c(0, 1 - MCN_uni, 0), n = lpi)$y
      e_MICR2_uni <- e_MICR2_uni + stats::approx(seq(nun), c(0, MICR2_uni, 0), n = lpi)$y
    }
  }
  if(!silent) cat('divisions for the MI ratio. \n')
  df.main <- data.frame('PI' = 1:lpi, 
                        'PI_points_x' = progind$PI, 
                        'PI_points_y' = progind$Time/lpi,
                        'kin' = kin)
  if('MI' %in% comb_met) {
    e_MI_uni <- e_MI_uni/.lt(n_unif)
    df.main <- cbind(df.main, 'MI' = e_MI_uni)
  }
  if(any(comb_met %in% allowed_mets)){
    e_MIC_uni <- e_MIC_uni/.lt(n_unif)
    e_MAS_uni <- e_MAS_uni/.lt(n_unif)
    e_MEV_uni <- e_MEV_uni/.lt(n_unif)
    e_MCN_uni <- e_MCN_uni/.lt(n_unif)
    e_MICR2_uni <- e_MICR2_uni/.lt(n_unif)
    convMIC_cov <- .normalize(stats::convolve(x = 1:lpi, y = e_MIC_uni))
    convMIC <- .cut0(e_MIC_uni - convMIC_cov)
    df.main <- cbind(df.main,
                     'MIC' = e_MIC_uni, 
                     'MIC_kin' = (e_MIC_uni + kin) / 2,
                     'MAS' = e_MAS_uni, 
                     'MEV' = e_MEV_uni, 
                     'MCN' = e_MCN_uni, 
                     'MICR2' = e_MICR2_uni,
                     'convDiff' = convMIC,
                     'conv' = convMIC_cov,
                     'diff' = .normalize(c(0,diff(e_MIC_uni))))
                     # 'diffcut' = .normalize(.cut0(c(0,diff(e_MIC_uni)))),
                     # 'ddiffcut' = .normalize(.cut0(-c(0,diff(.normalize(.cut0(c(0,diff(e_MIC_uni)))))))),
                     # 'ddcutkin' = .normalize(.cut0(-c(0,diff(.normalize(.cut0(c(0,diff(e_MIC_uni + kin))))))))) # wtf
  }
  
  # setting up the stats
  df.main <- as.data.frame(df.main)
  if(comb_met[1] == 'MIC') wtp <- df.main$MIC
  else if(comb_met[1] == 'MI') wtp <- df.main$MI
  else if(comb_met[1] == 'MAS') wtp <- df.main$MAS
  else if(comb_met[1] == 'MEV') wtp <- df.main$MEV
  else if(comb_met[1] == 'MCN') wtp <- df.main$MCN
  else if(comb_met[1] == 'MICR2') wtp <- df.main$MICR2
  else if(comb_met[1] == 'convDiff') wtp <- df.main$convDiff
  else if(comb_met[1] == 'conv') wtp <- df.main$conv
  else if(comb_met[1] == 'MIC_kin') wtp <- df.main$MIC_kin # wtp what to plot (the ending peak finding subject)
  # else if(comb_met[1] == 'diff_kin') wtp <- (df.main$ddiffcut + kin )/ 2
  else if(comb_met[1] == 'diff') wtp <- df.main$ddiffcut
  else if(comb_met[1] == 'kin') wtp <- kin
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
  df.pp <- data.frame(pkx, pky)
  
  # Setting the number of clusters to select first tot barriers
  if(is.null(n.cluster)) nclu <- 15
  else nclu <- n.cluster
  nbar.selected <- (min(max(15, nclu), .lt(rank.pk))) # this is only showing the specific number of crosses
  if(nbar.selected < 2) stop('Found less than 2 clusters (i.e. less than 1 barrier). Try to decrease the peak span.')
  nba <- nclu - 1 # cluster and barrier number 
  rnk_ts <- rank.pk[1:nbar.selected] # rank to show (number of peaks to show)
  
  if(F){ # entropy consideration for heuristic future
    dhc <- .dens_histCounts(progind, breaks = rank.pk, nx = ny, ny = ny) # breaks are the barriers points on x
    dens <- dhc$density
    ent <- .normalize(apply(dens, 2, .myShEn))
  }
  
  # Color handling for curves
  col_vec <- c('darkred', RColorBrewer::brewer.pal(max(.lt(comb_met), 3), name = 'Set1')[-1]) # first color is a red-violet for which I prefer darkred
    
  # if(nrow(df.main) > 50000) df.main <- df.main[unique(round(seq(1, lpi, length.out = 50000))),] # it is not clear for the barriers and not nec (it is fast)
  tpl <- ggplot(data = df.main, mapping = aes_string(x = 'PI')) + theme_classic() + ylab('UniDiv score') +
         geom_point(mapping = aes_string(x = 'PI_points_x', y = 'PI_points_y'), size = 0.8, col = 'grey') 
  
  # Adding the specific curves
  for(met.i in 1:.lt(comb_met)) tpl <- tpl + geom_line(mapping = aes_string(y = comb_met[met.i]), col = col_vec[met.i])
  
  # if('MI' %in% comb_met) tpl <- tpl + geom_line(mapping = aes(y = MI), col = 'darkred') 
  # if('MIC' %in% comb_met) tpl <- tpl + geom_line(mapping = aes(y = MIC), col = 'darkred') 
  # if('MAS' %in% comb_met) tpl <- tpl + geom_line(mapping = aes(y = MAS), col = 'darkred') 
  # if('MEV' %in% comb_met) tpl <- tpl + geom_line(mapping = aes(y = MEV), col = 'darkred') 
  # if('MCN' %in% comb_met) tpl <- tpl + geom_line(mapping = aes(y = MCN), col = 'darkred') 
  # if('MICR2' %in% comb_met) tpl <- tpl + geom_line(mapping = aes(y = MICR2), col = 'darkred') 
  # if('convDiff' %in% comb_met) tpl <- tpl + geom_line(mapping = aes(y = convDiff), col = 'darkred') 
  # if('conv' %in% comb_met) tpl <- tpl + geom_line(mapping = aes(y = conv), col = 'darkred') 
  # if('MIC_kin' %in% comb_met) tpl <- tpl + geom_line(mapping = aes(y = MIC_kin), col = 'darkred') 
  # if('diff' %in% comb_met) tpl <- tpl + geom_line(mapping = aes(y = diff), col = 'darkred') 
  # if('kin' %in% comb_met) tpl <- tpl + geom_line(mapping = aes(y = kin), col = 'darkred') 
  
  if(dbg_nSBR) browser()
  # Final add of barriers and similaria
  tpl <- tpl + geom_point(data = df.pp, mapping = aes_string(y = 'pky', x = 'pkx')) +
               geom_vline(xintercept = rank.pk[1:nba], col = 'black', size = 1, linetype="dotted") +
               geom_point(data = data.frame(a = rnk_ts, b = wtp[rnk_ts]), mapping = aes_string(x = 'a', y = 'b'),
                          shape = 3, col = 'darkblue', size = 5, stroke = 2.5) + 
               xlab('Progress Index') + ylab(paste0('I', comb_met[1])) + 
               theme(axis.title.y.left = element_text(colour = col_vec[1]))
  # tpl + scale_colour_manual(name="Statistic", values=c("MIC"= 'darkred',"MI"="lightgreen")) # does not work with this setting
  # tpl <- tpl + scale_y_continuous(sec.axis = sec_axis(~ . *3000 , name = "Temporal progress"), limits = c(0, 1))
  # tpl + theme(axis.title.y.right = element_text(color = 'darkgray'),
  #             axis.line.y.right = element_blank(),
  #             axis.ticks.y.right = element_line(color = 'darkgray'),
  #             axis.text.y.right = element_text(color = 'darkgray'))
  
  
  
  # randomization procedure
  if(!is.null(random_picks)){
    if(!silent) cat('Random barrier picks started. Calculation made with the following number of uniformly random picks:', random_picks, '\n')
    
    # Checks on random_picks
    if(random_picks < 10 ) stop('Use more than 10 random_picks to have some interesting result. Please.')

    # main random pick function
    rnd.obj <- .rnd_bar_estimation(pi.tab = data, ann = ann, ncl = nclu, span = span, ny = ny, 
                                   unifsplits = n_unif, random_trials = random_picks, comb_met = comb_met,
                                   silent = silent)
  }
  
  # microbenchmarking MIC MI
  # if(F){
  #   mbm <- microbenchmark::microbenchmark(
  #     minerva.out = sapply(seq(nun-2), function(j) minerva::mine(x = uni_cnts[, j], y = uni_cnts[, j+1])),
  #     MI_uni = sapply(seq(nun-2), function(j) infotheo::mutinformation(uni_cnts[, j], uni_cnts[, j+1])), 
  #     times = 100
  #   )
  #   ggplot2::autoplot(mbm)  
  # }

  # FINAL OUTPUT
  if(plot && is.null(random_picks)) print(tpl)
  if(plot && !is.null(random_picks)) print(rnd.obj$ggp.bar)
  returning_list <- list('nbins' = c(ny, ny), 'barriers' = rank.pk)
  if(return_plot) returning_list <- list('nbins' = c(ny, ny), 'barriers' = rank.pk, 'plot' = tpl)
  else if(return_plot && !is.null(random_picks)) returning_list <- list('nbins' = c(ny, ny), 'barriers' = rank.pk, 'plot' = tpl, 'rnd.picks' = rnd.obj)
  else if(!is.null(random_picks)) returning_list <- list('nbins' = c(ny, ny), 'barriers' = rank.pk, 'rnd.picks' = rnd.obj)
  invisible(returning_list)
}


#########################################################
#          Functions ext
#########################################################
# major function to estimate random barrier variability
.rnd_bar_estimation <- function(pi.tab, ann, ncl = 4, span = 1000, ny = 50, 
                                unifsplits = seq(5, 100, 8), random_trials = 500, comb_met = c('MIC'),
                                silent = F){
  # small init
  nba <- ncl - 1
  
  # for the reference score
  ref_nSBR <- CampaRi::nSBR(data = pi.tab, n.cluster = ncl, 
                            comb_met = comb_met,
                            unif.splits = unifsplits,  
                            pk_span = span, ny = ny, plot = F, 
                            silent = T, return_plot = T)
  babar <- ref_nSBR$barriers
  if(!silent) cat('Number of barriers found with optimization:', length(babar), '\n')
  ref_sc <- CampaRi::score_sapphire(the_sap = pi.tab, ann = ann, manual_barriers = babar[1:nba], silent = T)
  if(!silent) cat('Found the following score with optimization:', ref_sc$score.out, '\n')
  if(is.na(ref_sc$score.out)) stop('Found NA in the score.')
  if(is.na(ref_sc$score.out)) ref_sc$score.out <- 0 # no more necessary
  
  # calculating the repetition 
  repe <- lapply(X = 1:random_trials, FUN = function(x){ 
    if(random_trials > 99 && x%%round(random_trials/100) == 0) cat(round(x*100.0/random_trials,2),'%\r')
    ba <- .rand_pick_min(ar = seq(2, length(ann) - 2), min.dist = span, n.picks = nba)
    # ba <- round(runif(nba, 2, length(ann) - 2))
    l1 <- CampaRi::score_sapphire(the_sap = pi.tab, ann = ann, manual_barriers = ba, silent = T)
    return(list('sc' = l1$score.out, 'ba' = ba))
  })
  
  # split the results accordingly
  sco <- sapply(repe, function(x) return(x$sc))
  sco.nas <- sum(is.na(sco))
  # sco[is.na(sco)] <- 0 # no more needed
  best_ba <- sapply(repe, function(x) return(x$ba))
  
  # calculating rough zscore
  zcs <- (ref_sc$score.out - mean(sco)) / sd(sco)
  pval_zsc <- 2*stats::pnorm(-abs(zcs))
  pval_95per <- stats::quantile(c(sco),  probs = 0.95)

  # use the classic plot to see the points (barrier points) which exceded the combination found.
  df_pl <- data.frame('barr' = c(best_ba), 'scr' = rep(sco, each=nba))
  # find the two maxima (red points)
  dhe_point <- data.frame('xm' = df_pl[df_pl[, 'scr'] == max(df_pl[, 'scr']), 1], 
                          'ym' = df_pl[df_pl[, 'scr'] == max(df_pl[, 'scr']), 2])
  
  # add layers to the plot
  gtemp1 <- ref_nSBR$plot + geom_point(data = df_pl, aes_string(x = 'barr', y = 'scr'), col = 'green4', size = 0.65) +
    geom_point(data = dhe_point, aes_string(x = 'xm', y = 'ym'), col = 'red', size = 2.5) + 
    geom_hline(aes_string(yintercept = ref_sc$score.out), size = 1, col = 'darkblue') + 
    scale_y_continuous(sec.axis = sec_axis(~ . , name = "NMI"), limits = c(0, 1)) +
    theme(axis.title.y.right = element_text(color = 'green4', margin = margin(l = 5)),
          axis.line.y.right = element_blank(),
          axis.ticks.y.right = element_blank(),
          axis.text.y.right = element_text(color = 'green4', margin = margin(l = -16)))
  #+ ylab('') #+ 
  # annotate('text', x = -10000, y = 0.25, label = TeX('\text{IMIC / } {\\color{DarkGreen} \text{NMI}}'), angle = 90) # fail
  max_rnd_nmi <- CampaRi::score_sapphire(the_sap = pi.tab, ann = ann, manual_barriers = dhe_point$xm, silent = T)
  if(!silent) cat('Found the following score using random procedure:', max_rnd_nmi$score.out, '\n')
  
  # showing how many are better!
  df_pl_better <- df_pl[df_pl[, 'scr'] > ref_sc$score.out,]
  if(nrow(df_pl_better)>0){
    if(!silent) cat('Found', nrow(df_pl_better), 'better barriers.\n')
  }else{
    if(!silent) cat('Nothing better!\n')
  }
  
  # plotting the simple final result
  gtemp2 <- ggplot() + geom_density(aes_string(x = 'sco'), fill = 'lightblue') + theme_classic() + 
    geom_vline(aes_string(xintercept = ref_sc$score.out), col = 'darkred') + #ggtitle(paste0('dataset ', n_to_anal, ' randomization of the barriers 500 times')) +
    xlab('Score') + ylab('Density')
  
  # extracting exact positions for the text annotations
  ggbuild <- ggplot_build(gtemp2)
  yeight_text <- ggbuild$layout$panel_scales_y[[1]]$range$range
  xwid_text <- ggbuild$layout$panel_scales_x[[1]]$range$range
  gtemp2 <- gtemp2 + annotate('text', x = ref_sc$score.out + 0.03*diff(xwid_text), y = diff(yeight_text)/2, angle = -90,
                              label = paste0('Z-score: ', round(zcs, digits = 2)))
  
  
  # final return
  invisible(list('ggp.bar' = gtemp1, 'ggp.dist.rnd.bar' = gtemp2, 'rnd.scores.nas' = sco.nas, 
                 'rnd.scores' = df_pl, 'max.rnd.nmi' = max_rnd_nmi, 'max.rnd.bar' = dhe_point,
                 'ref.score' = ref_sc, 'z.scores' = zcs, 'pval_zsc' = pval_zsc, 'pval_95per' = pval_95per,
                 'std.n.bar' = length(babar)))
}

# function for random picking values from an array with a certain minimum distance between the values
.rand_pick_min <- function(ar, min.dist, n.picks){
  stopifnot(is.numeric(min.dist), 
            is.numeric(n.picks), n.picks%%1 == 0)
  if(length(ar)/n.picks < min.dist) 
    stop('The number of picks exceeds the maximum number of divisions that the array allows which is: ', 
         floor(length(ar)/min.dist))
  picked <- array(NA, n.picks)
  copy <- ar
  for (i in 1:n.picks) {
    # if(length(copy) < 1) print(i)
    # if(length(copy) < 1) print(length(copy))
    # if(length(copy) < 1) print(copy)
    stopifnot(length(copy) > 0)  
    picked[i] <- sample(copy, 1)
    copy <- copy[ abs(copy - picked[i]) >= min.dist ]
  }
  return(picked)
}

# simple uniform splitting function (it returns a list)
.split_by <- function(ar, by) return(split(ar, ceiling(seq_along(ar)/by)))

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
    for (j in 1:ny) cnts[j, i] <- sum(hist.internal$counts[c(ncls:ncle), j]) # here is where the y is summed up
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
# # xr1 <- c(margin:(lpi-round(lpi*0.99)))
# kin.pl <- .normalize(-log(cutf / lpi)) # xr1 impose a selection of 99% of the snapshots to avoid the initial inf
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