#' @title Identifying basins in the SAPPHIRE plot
#' @description
#'     \code{basins_recognition} uses the information provided by the SAPPHIRE plot in order to identify free energy basins, 
#'     thus leading to a discretized trajectory. In particular, it inspects both the kinetic annotation and times of occurence
#'     (called dynamic annotation hereafter) of the progress index and, at the end, it matches (or merges) the resulting state barriers coming from each of the single analysis.  
#'     The analysis of the dynamical trace is driven by an underlying 2-D histogram whose number of bins has to be given either manually or automatically (this option is available only for number of bins on the y-axis). These values has to be chosen appropriately because they provide the minimum size that a state should have in order to be identified.
#'     Regarding the kinetic annotation, the identification of the states is equivalent to discerning its main peaks. This task is performed by means of a smoothing filter, either moving average or Savitzky Golay, and a naive algorithm to find the local maxima. The window size of the filter and of the peak search algorithm is equal to the number of snapshots divided by \code{nx}.
#'     Eventually, the user can decide on either selecting only the partitions in common to both the annotation (matching) or retaining them all (merging).
#'     
#' @param data Name of the PROGIDX_<...>.dat file or Data frame with three columns containing, in order, the (sorted) progress index, the relative time indices and the cut function. 
#' @param nx Number of bins on the x-axis of the 2-D histogram.
#' @param ny Number of bins on the y-axis of the 2-D histogram. Default to nx.
#' @param ny.aut Logical value indicating whether a suitable number of bins on the y-axis has to be identified automatically. Default value is \code{FALSE}.
#' @param local.cut Logical value indicating whether the localized cut function has to be used (see references). Default value is \code{FALSE}.
#' @param match Logical value indicating whether results from kinetic and dynamic annotation have to be either matched or merged. Default value is \code{FALSE}               
#' @param avg.opt Smoothing filter in the kinetic annotation analysis
#'       \itemize{
#'            \item "\code{movav}" for moving average filter
#'            \item "\code{SG}" for savitzkyGolay filter
#'       }
#' @param plot A logical value indicating whether to display the SAPPHIRE plot with the resulting partitions or not. Black partitions are the matched ones, blue ones derive only from the dynamic analysis and orange ones only from the kinetic analysis. Default value is \code{FALSE}.
#' @param out.file A logical value indicating whether to write an output file with the state sequence ordered with respect to the PI. Default value is \code{TRUE}.
#' @param silent A logical value indicating whether the function has to remain silent or not. Default value is \code{FALSE}
#' @param ...
#'      \itemize{
#'        \item "\code{time.series}" File name. If specified, it substitutes the time series of the PROGIDX_<..> file with the one provided by the file.
#'        \item "\code{only.kin}" Logical value indicating whether it should be performed only the analysis of the kinetic annotation, thus disregarding completely the dynamic trace.  
#'        \item "\code{pol.degree}" Degree of the Savitzky-Golay filter in case it was chosen. Default to 2.
#'      }
#'      
#' @return A list containing
#'       \itemize{
#'         \item "\code{tab.st}" Data frame containing the boundaries of each state, their lengths and the type of the right boundary. Type=1 means that it is a matched partition, type=2 is only dynamic, type=3 is only kinetic. The type of the last state is always equal to 1.
#'         \item "\code{nbins}" 2-D vector containing number of bins on x-axis and y-axis.
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
#' @importFrom gplots hist2d
#' @importFrom distrEx HellingerDist
#' @importFrom outliers grubbs.test
#' @importFrom prospectr movav savitzkyGolay
#' @importFrom splus2R peaks
#' @importFrom grDevices dev.new
#' @importFrom sfsmisc is.whole 
#' @importFrom distr DiscreteDistribution
#' @importFrom data.table fread fwrite
#' @export basins_recognition

basins_recognition <- function(data, nx, ny=nx, ny.aut=FALSE, local.cut=FALSE, match=FALSE, avg.opt=c("movav", "SG"), plot=FALSE, out.file=TRUE, silent=FALSE, ...) {
    call <- match.call()
    
    if(!is.character(data) && !is.data.frame(data)) stop("data must be a string or a data frame")
    if(is.character(data) && !all(grepl("PROGIDX", data))) stop("Please provide a data name starting with 'PROGIDX'" )
    if(!is.whole(nx)) stop("nx must be an integer")
    if(!is.whole(ny)) stop("ny must be an integer")
    if(!is.logical(local.cut)) stop("local.cut must be a logical value")
    if(!is.logical(match)) stop("match must be a logical value")
    if(!is.logical(plot)) stop("plot must be a logical value")
    if(!is.logical(out.file)) stop("out.file must be a logical value")
    avg.opt.arg <- c("movav", "SG")
    if(!(avg.opt[1] %in% avg.opt.arg)) stop("Average option not valid")

    input.args <- list(...)
    avail.extra.arg <- c("pol.degree","only.kin","time.series")
    if(!is.null(names(input.args)) && any(!(names(input.args) %in% avail.extra.arg))) 
        warning('There is a probable mispelling in one of the inserted variables. Please check the available extra input arguments.')
    if("only.kin" %in% names(input.args)) {
        only.kin <- input.args$only.kin
        if(only.kin) match <- TRUE
    } else only.kin <- FALSE
    if(avg.opt[1]=="SG") {
        if(!("pol.degree" %in% names(input.args))) {
            if(!silent) cat("SG but pol.degree not specified, set to default value 2\n")
            pol.degree <- 2
        } else if(!(input.args$pol.degree %in% c(2:6))) {
            warning("Degree of the polynomial not valid, set to default value 2")
            pol.degree <- 2 
        } else pol.degree <- input.args$pol.degree
    }

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

    ## INPUT FILE 
    if(!silent) cat("Reading PROGIDX file...\n")
    if(is.data.frame(data)){
        progind <- data
        colnames(progind) <- c("PI", "Time", "Cut")
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
                }
                else {
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
        rawsetsort <- rawset[order(rawset$center),]

#################################################################################
        ## WEIGHTED SUMS FUNCTIONS
##############################################################################
        ## Forward
        sumwr <- function(first,meanopt,wthopt) {
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
            for (i in 1:nx) {
                idxlist <- which(rawfirst[idrs:ny] < hist$x.breaks[i+1]) + idrs - 1 
                if(length(idxlist)==0) sm[i] <- 0
                else sm[i] <- sum(rawwth[idxlist])
            }
            return(sm)
        }

#################################################################################
        ## Backward
        backsumwr <- function(first,meanopt,wthopt) {
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
        MaxHaar <- function(first,meanopt,wthopt) {
            filt <- filter(as.ts(sumwr(first,meanopt,wthopt)), Haarfilt, method="convolution", sides=2)
            filt <- c(rep(0,K),filt[c((K+1):(length(filt)-K))],rep(0,K))  #Filt with 0 instead of NA
            xbr <- hist$x[which(peaks(filt,dpeaks.dyn))]
            return(c(1,xbr,cstored))
        }
#################################################################################
        ## Backward Filter with Haar wavelet
        BackMaxHaar <- function(first,meanopt,wthopt) {
            filt <- filter(as.ts(rev(backsumwr(first,meanopt,wthopt))), Haarfilt, method="convolution", sides=2)
            filt <- c(rep(0,K),filt[c((K+1):(length(filt)-K))],rep(0,K))  #Filt with 0 instead of NA
            xbr <- hist$x[which(rev(peaks(filt, dpeaks.dyn)))]
            return(c(1,xbr,cstored))
        }

#################################################################################
        ##CHOICE OF THE BREAKS
#################################################################################
        breaks.max <- MaxHaar("max","min",1) 
        breaks.min <- BackMaxHaar("min","max",1)
       
###############################################################################
        ##HARD BREAKS: Joining selected breaks.min and breaks.max with res
###############################################################################
        if(lt(breaks.max)==2 && lt(breaks.min)>2) {
            breaks.tot <- breaks.min
        } else if (lt(breaks.max)>2 && lt(breaks.min)==2) {
            breaks.tot <- breaks.max
        } else if (lt(breaks.max)==2 && lt(breaks.min)==2) {
            breaks.tot <- NULL
        } else {
            sep <- NULL
            idx <- 0
            selbreaks.max <- breaks.max[-c(1,length(breaks.max))]
            selbreaks.min <- breaks.min[-c(1,length(breaks.min))]
            softbreaks.max <- selbreaks.max
            softbreaks.min <- selbreaks.min
            for (i in 1:length(selbreaks.max)) {
                selcell.max <- which(hist$x==selbreaks.max[i])
                if (i==length(selbreaks.max)) selcell2.max <- length(hist$x)
                else selcell2.max <- which(hist$x==selbreaks.max[i+1])
                for (j in 1:length(selbreaks.min)) {
                    selcell.min <- which(hist$x==selbreaks.min[j])
                    if (selcell.max==selcell.min) {
                        idx <- idx+1
                        sep[idx] <- hist$x[selcell.max]
                        softbreaks.max <- softbreaks.max[-which(softbreaks.max==selbreaks.max[i])]
                        softbreaks.min <- softbreaks.min[-which(softbreaks.min==selbreaks.min[j])]
                        next
                    }
                    if ((selcell.max+1) == selcell.min) {
                        if (!is.na(match(hist$x[selcell.max],sep))) {
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
            breaks.tot <- sort(c(1,cstored,sep,softbreaks.max,softbreaks.min))
        }
#################################################################################
        ## HISTOGRAM of each PARTITION and Max density 
################################################################################
        if(is.null(breaks.tot)) {
            brk.dyn <- NULL
        } else {
            ## Breaks.Tot
            brkjy.tot <- matrix(rep(0,(length(breaks.tot)-1)*ny), nrow=ny, ncol=(length(breaks.tot)-1))
            for (i in 1:(length(breaks.tot)-1)) {
                if (i==1) {
                    ncls <- 1
                    ncle <- which(hist$x==breaks.tot[i+1])
                }
                else if (i == (length(breaks.tot)-1) ) {
                    ncls <- which(hist$x==breaks.tot[i])+1
                    ncle <- nx
                }
                else {
                    ncls <- which(hist$x==breaks.tot[i])+1
                    ncle <- which(hist$x==breaks.tot[i+1])
                }
                for (j in 1:ny) {
                    brkjy.tot[j,i] <- sum(hist$counts[c(ncls:ncle),j])
                }
            }
            dens.tot <- brkjy.tot
            for (i in 1:(length(breaks.tot)-1)) {
                ## dens.tot[,i] <-  brkjy.tot[,i]/(breaks.tot[i+1]-breaks.tot[i])
                dens.tot[,i] <-  brkjy.tot[,i]/(sum(brkjy.tot[,i]))
            }

########################################################################
            ##JOINING PARTITIONS METHODS
########################################################################
            ## Computation of Distances Hell and Kolm between consecutive partitions
            distHell.tot <- NULL
            for(j in 1:(length(breaks.tot)-2)) {
                prova1 <- DiscreteDistribution(supp = c(1:ny) , prob=dens.tot[,j])
                prova2 <- DiscreteDistribution(supp = c(1:ny) , prob=dens.tot[,j+1])
                distHell.tot[j] <- HellingerDist(prova1,prova2)
            }

########################################################################################
            ##MAIN JOINING procedure: comparison with training uniform samples
########################################################################################

            lstHell.tot <- sort.int(distHell.tot, index.return=TRUE)$ix
            discbreaks.tot <- NULL
            ncounts <- NULL
            flagbreak <- 0
            ll <- 0
            for (i in lstHell.tot) { 
                sampleHell.tot <- NULL
                l1 <- which(hist$x==breaks.tot[i])+1 
                if (i==1) l1 <- 1
                l2 <- which(hist$x==breaks.tot[i+1])
                l3 <- which(hist$x==breaks.tot[i+2])
                if (i==length(lstHell.tot)) l3 <- nx
                pr <- (l2-l1+1)/(l3-l1+1)
                if(!silent) {
                    cat.str <- paste("Would-be basins", paste(round(breaks.tot[i:(i+2)]), collapse=" "))
                    cat(cat.str)
                    nlett <- nchar(cat.str)
                }
                for (j in 1:ny) ncounts[j] <- sum(hist$counts[c(l1:l3),j])
                for (idx in 1:nsample) {
                    unif1 <- rep(0,ny)
                    unif2 <- rep(0,ny)
                    for (j in 1:ny) {
                        if(ncounts[j]==0) next
                        dum <- sample(c(0:1), ncounts[j], replace=TRUE, prob=c(pr,1-pr))
                        unif1[j] <- ncounts[j]-sum(dum)
                        unif2[j] <- sum(dum)
                    }
                    dens1 <- unif1/sum(unif1)
                    dens2 <- unif2/sum(unif2)
                    part1 <- DiscreteDistribution(supp = c(1:ny) , prob=dens1)
                    part2 <- DiscreteDistribution(supp = c(1:ny) , prob=dens2)
                    sampleHell.tot[idx] <- HellingerDist(part1,part2)
                }
                grubbsHell.tot <- grubbs.test(c(sampleHell.tot,distHell.tot[i]), type=10)
                foo <- c(sampleHell.tot,distHell.tot[i])
                if (distHell.tot[i] < max(sampleHell.tot) | grubbsHell.tot$p.value>conf.lev) {
                    if(!silent) cat(rep(" ",15+4+3*nchar(as.character(cstored))-nlett), "--> Joining partitions\n", sep="")
                    flagbreak <- 0
                    ll <- ll+1
                    discbreaks.tot[ll] <- breaks.tot[i+1]
                }
                else {
                    if(!silent) cat("\n")
                    flagbreak <- flagbreak+1
                }
                if (flagbreak==cutjoin) break
            }

            if(!is.null(discbreaks.tot)) {
                brk.dyn <- sort(breaks.tot[-match(discbreaks.tot, breaks.tot)])
            } else brk.dyn <- sort(breaks.tot)
        }
        if(!silent) cat("End of the dynamic analysis\n")
    } 


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

    if(lt(max.mv)!=0) {
### First Cleaning on mv:: check separation between consecutive max
        max.mv.tmp <- max.mv
        rif <- max.mv[length(max.mv)]
        if(length(max.mv)>1){
            for (i in (length(max.mv)-1):1) {
                if (rif-max.mv[i]<dpeaks.kin/2) {
                    max.mv.tmp <- max.mv.tmp[-match(max.mv[i],max.mv.tmp)] ##Remove the smallest
                }
                else rif <- max.mv[i]
            }
        }
        max.mv <- max.mv.tmp

### Second Cleaning Attempt on mv
        adj.mv <- NULL
        adj.mv[1] <- which.min(kin.mv[1:max.mv[1]])
        for (i in 1:length(max.mv) ) {
            if(length(max.mv)==1) next
            if (i==length(max.mv)) adj.mv[i+1] <- which.min(kin.mv[round(max.mv[i]):cstored])+round(max.mv[i])-1
            else adj.mv[i+1] <- which.min(kin.mv[round(max.mv[i]):round(max.mv[i+1])])+round(max.mv[i])-1 
            amax <- mean( c(kin.mv[max.mv[i]]-kin.mv[adj.mv[i]], kin.mv[max.mv[i]]-kin.mv[adj.mv[i+1]]) ) ## Lengths of adjacent (closest) vertical bars
            if (i==1) amin <- kin.mv[max.mv[i+1]]-kin.mv[adj.mv[i+1]]
            else if (i==length(max.mv)) amin <- kin.mv[max.mv[i-1]]-kin.mv[adj.mv[i]]
            else amin <- mean(c(kin.mv[max.mv[i-1]]-kin.mv[adj.mv[i]], kin.mv[max.mv[i+1]]-kin.mv[adj.mv[i+1]]) ) 
            ## print(paste("Evaluating", i, max.mv[i], "with ratio", 100*amax/amin, "%"))
            if ( amax/amin < thr.ratio) {
                ## print(paste("Excluding", i, max.mv[i]))
                max.mv.tmp <- max.mv.tmp[-which(max.mv.tmp==max.mv[i])]
            }
        }
        brk.kin <- round(max.mv.tmp)
    }
    else {
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
        if(lt(vkin)!=0) vkin <- refine(vkin)
        vdyn <- brk.dyn
        vdyn <- vdyn[which(vdyn != 1)]
        vdyn <- vdyn[which(vdyn != cstored)]
        if(lt(vdyn)!=0) {
            vdyn <- refine(vdyn)
        }
        if(!silent) cat("Number of dynamic and kinetic breaks are respectively",length(vdyn), length(vkin), "\n")
##################################################################################
        ## Matching
        if(lt(vdyn) !=0 && lt(vkin!=0) ) {
            ll <- 0
            for (i in 1:length(vkin) ) { 
                dist <- NULL
                ## if(!silent) cat("********************************************\n")
                ## print(paste("Analyzing kin break n", i, " ::: ", vkin[i]))
                set <- c(tail(vdyn[which(vdyn<vkin[i])],1), head(vdyn[which(vdyn>=vkin[i])],1) )
                if (length(set)==0) break
                for (j in 1:length(set) ) {
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
                cat("Number of matched partitions is",lt(breaks), "\n")
                if(lt(vdyn)==0) cat("All the dynamic breaks are matched\n")
                if(lt(vkin)==lt(breaks)) cat("All the kinetic breaks are matched\n")
            }
        } else {
            if(!match) {
                breaks <- c(vkin, vdyn)
                match <- TRUE ## To skip next session
                if(lt(vkin)!=0) only.kin <- TRUE
            } else breaks <- NULL
        }
#######################################################################
        ## Merging option: adding the remaining partitions 
        if(!match) {
            if(lt(vdyn)!=0 && lt(vkin)!=0) {
                ## Looking for residual vdyn too close to any of the vkin (rare)
                dist <- NULL
                ll <- 0
                for(i in 1:length(vdyn)) {
                    for(j in 1:length(vkin)) {
                        ll <- ll+1
                        dist[ll] <- abs(vdyn[i]-vkin[j])
                    }
                }
                ## Identifying and removing them from vdyn vector
                dist.mtx <- matrix(dist<sd.kin, nrow=length(vdyn), ncol=length(vkin), byrow=TRUE)
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
    }
    else if(lt(brk.kin)!=0) breaks <- refine(brk.kin) ## If only.kin==TRUE
    else breaks <- NULL


######################################################################
    ## OUTPUT section
#####################################################################
    if(is.null(breaks)) {
        if(!silent) cat("NO barriers have been found")
        if(out.file) {
            output.match <- data.frame(PI=progind$PI, Time=progind$Time, State=rep(1, nrow(progind)))
            if(is.character(data)) {
                if(!silent) cat("Writing", gsub("PROGIDX", "BASINS", strsplit(data,"/",fixed=T)[[1]][length(strsplit(data,"/",fixed=T)[[1]])]), "...\n")
                fwrite(output.match, file=gsub("PROGIDX", "BASINS", data), sep='\t', row.names=FALSE, col.names=FALSE)
            } else {
                if(!silent) cat("Writing PROGIDX_", as.character(progind$Time[1]))
                fwrite(output.match, file=paste("./BASINS_", as.character(progind$Time[1]), ".dat", sep=""), sep='\t', row.names=FALSE, col.names=FALSE)
            }
        }
    } else {
        if(!silent) cat("Number of states is", length(breaks)+1, "\n")
        if(!silent) cat(breaks, "\n")

        vec <- sort(breaks)
        state <- NULL
        for (i in 1:(length(vec)+1)) {
            if (i==1) ib <- 0
            else ib <- vec[i-1]
            if (i==length(vec)+1) fb <- cstored
            else fb <- vec[i]
            ## if(!silent) cat(i,ib,fb,fb-ib,"\n")
            state <- c(state,rep(i,fb-ib))
        }
        if(out.file) {
            output.match <- data.frame(PI=progind$PI, Time=progind$Time, State=state)
            if(is.character(data)) {
                if(!silent) cat("Writing", gsub("PROGIDX", "BASINS", strsplit(data,"/",fixed=T)[[1]][length(strsplit(data,"/",fixed=T)[[1]])]), "...\n")
                fwrite(output.match, file=gsub("PROGIDX", "BASINS", data), sep='\t', row.names=FALSE, col.names=FALSE)
            } else {
                if(!silent) cat("Writing PROGIDX_", as.character(progind$Time[1]))
                fwrite(output.match, file=paste("./BASINS_", as.character(progind$Time[1]), ".dat", sep=""), sep='\t', row.names=FALSE, col.names=FALSE)
            }
        }
    }    

####################################################################
    ## PLOT Section
####################################################################

    if(plot == TRUE){
        dev.new(width=15, height=10)
        par(mgp=c(0, 0.4, 0))
        par(ps=6)
        par(mar = c(3.5, 0, 1, 3), oma = c(2, 4, 2, 2))
        margin <- round(cstored*0.99)
        cx <- 2.2
        sc <- 1.8
        xr1 <- c(margin:(cstored-margin))
        xr <- c(1,cstored)
        yr <- range(progind$Time)*sc
        plot(0, 0, main="", xlim=xr, ylim=yr, xlab="", ylab="", axes=FALSE, frame.plot=FALSE, type="n")
        kin.pl <- -log(cutf / cstored)[xr1]
        xx.lab <- c(1,round(breaks),cstored)
        axis(1, at=xx.lab, tck=.01, cex.axis=1.8)
        axis(3, labels=rep("", length(xx.lab)), at=xx.lab, tck=.01)
        mtext("Progress Index", side=1, line=1.5, cex=cx )
        yy.lab1 <- format(c(min(kin.pl), min(kin.pl)+(max(kin.pl)-min(kin.pl))*c(1:3)/3), digits=2)
        axis(2, labels=yy.lab1, at=scale(as.numeric(yy.lab1)), las=3, tck=.01, cex.axis=cx)
        mtext(expression("ln(("*italic(tau["SA"]+tau["AS"])*")/2)"), at=max(progind$Time)/2, side=2, line=1.8, cex=cx)
        yy.lab2 <- round(c(1, c(1:5)/5*max(progind$Time)))
        axis(4, labels=rep("", length(yy.lab2)), at=max(progind$Time)+round(yy.lab2*(sc-1)), las=2, tck=.01, hadj=-0.6, col="red")
        mtext(yy.lab2, side=4, las=2, line=0.2, at=max(progind$Time)+yy.lab2*(sc-1), col="red", cex=1.8)
        mtext("Time", at=max(progind$Time)*(sc+1)/2, side=4, line=2.4, cex=cx, col="red")

        if (cstored>10000) { lst <- seq(1, cstored, by=5)
        } else lst <- c(1:cstored)
        points(progind$PI[lst], max(progind$Time[lst])+progind$Time[lst]*(sc-1), cex=0.005, col="tomato1")
        lines(xr1, scale(kin.pl), lwd=1, col="black")
        if(only.kin) abline(v=breaks, lwd=0.7, col="orange4")
        else if(!match) {
            ## Not-matched vdyn
            abline(v=vdyn, lwd=0.4, col="blue")
            ## Not-matched vkin
            abline(v=vkin[-match(brk.mtc, vkin)], lwd=0.4, col="orange4")
            ## Matched breaks
            abline(v=brk.mtc, lwd=0.7, col="black")
        } else abline(v=breaks, lwd=0.7, col="black")
    }
    
    if(is.null(breaks)){
        tab.st <- data.frame(n=1, start=1, end=cstored, length=cstored, type=1)
    } else {
        tab.st <- data.frame(n=c(1:(length(breaks)+1)), start=c(1,vec+1), end=c(vec, cstored), length=diff(c(0,vec, cstored)), type=c(rep(NaN,length(breaks)),1))
        if(only.kin) tab.st$type <- c(rep(3, length(breaks)), 1)
        else if(match) tab.st$type <- 1
        else {
            if(lt(brk.mtc)>0) tab.st$type[match(brk.mtc,tab.st$end)] <- 1
            if(lt(vdyn)>0) tab.st$type[match(vdyn, tab.st$end)] <- 2
            tab.st$type[match(vkin[-match(brk.mtc, vkin)], tab.st$end)] <- 3
        }
    }
    invisible(list(tab.st=tab.st, nbins=c(nx,ny), call=call))
    
}

