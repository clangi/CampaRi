#' @title Identifying basins in the SAPPHIRE plot
#' @description
#'     \code{basins_recognition} uses the information provided by the SAPPHIRE plot in order to identify free energy basins, 
#'     thus leading to a discretized trajectory. In particular, it inspects both the kinetic annotation and times of occurence
#'     (called dynamic annotation hereafter) of the progress index and, at the end, it matches (or merges) the resulting state barriers coming from each of the single analysis.  
#'     The analysis of the dynamical trace is driven by an underlying 2-D histogram whose number of bins has to be given either manually or automatically (this option is available only for number of bins on the y-axis). These values has to be chosen appropriately because they provide the minimum size that a state should have in order to be identified.
#'     Regarding the kinetic annotation, the identification of the states is equivalent to discerning its main peaks. This task is performed by means of a smoothing filter, either moving average or Savitzky Golay, and a naive algorithm to find the local maxima.
#'     Finally, the user can decide on either selecting only the partitions in common to both the annotation (matching) or retaining them all (merging).
#'     
#' @param file Name of the PROGIDX_<...>.dat file 
#' @param nx Number of bins on the x-axis of the 2-D histogram
#' @param ny.aut Logical value indicating whether a suitable number of bins on the y-axis has to be identified automatically. Default value is \code{TRUE}
#' @param match Logical value indicating whether results from kinetic and dynamic annotation have to be either matched or merged.                
#' @param avg.opt Smoothing filter in the kinetic annotation analysis
#'       \itemize{
#'            \item "\code{movav}" for moving average filter
#'            \item "\code{SG}" for savitzkyGolay filter
#'       }
#' @param plot A logical value indicating whether to plot results or not. Default value if \code{FALSE} 
#' @param ...
#'      \itemize{
#'        \item "\code{ny}" Number of bins on the y-axis of the 2-D histogram if \code{ny.aut=FALSE}. Default value is nx.
#'        \item "\code{pol.degree}" degree of the Savitzky-Golay filter
#'      }
#'      
#' @return None
#' 
#' @examples
#' 
#' adjl <- mst_from_trj(trj = matrix(rnorm(1000), nrow = 100, ncol = 10))
#' ret<-gen_progindex(adjl = adjl)
#' gen_annotation(ret_data = ret, local_cut_width = 10)
#' \dontrun{
#' basins_recognition("PROGIDX_000000000001.dat")
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
#' @importFrom is.whole sfsmisc
#' @importFrom distr DiscreteDistribution
#' @export basins_recognition

library(distrEx) #HellingerDist
library(gplots) #hist2d
library(outliers) #grubbs.test
library(RcppArmadillo) #loading required package
library(prospectr) #movav, savitzkyGolay
library(splus2R) #peaks 

file <- "/home/fcocina/Desktop/Campari/GitHub/feature_selection/BPTIdata/PROGIDXSAMPLE_000000020521.dat"
ptm <- proc.time()

#basins_recognition <- function(file, nx=500, ny.aut=TRUE, match=TRUE, avg.opt="movav", plot=FALSE, ...) {

    ## if(!is.character(file)) stop("file must be a string")
    ## if(!is.whole(nx)) stop("nx must be an integer")
    ## if(!is.logical(match)) stop("match must be a logical value")
    ## if(!is.logical(plot)) stop("plot must be a logical value")
    ## avg.opt.arg <- c("movav", "SG")
    ## if(!(avg.opt %in% avg.opt.arg)) stop("Average option not valid")

    ## input.args <- list(...)
    ## avail.extra.arg <- c("pol.degree","ny")
    ## if(!is.null(names(input.args)) && any(!(names(input.args) %in% avail.extra.arg))) 
    ##     warning('There is a probable mispelling in one of the inserted variables. Please check the available extra input arguments.')
    ## if(!ny.aut) {
    ##     if(!("ny" %in% names(input.args))) {
    ##         warning("ny not specified, set to default value nx")
    ##         ny <- nx
    ##     } else ny <- input.args$ny
    ## }
    ## if(avg.opt=="SG") {
    ##     if(!("pol.degree" %in% names(input.args))) {
    ##         print("SG but pol not specified")
    ##         pol.degree <- 2
    ##     } else if(!(input.args$pol.degree %in% c(2:6))) {
    ##         warning("Degree of the polynomial not valid, set to default value 2")
    ##         pol.degree <- 2 
    ##     } else pol.degree <- input.args$pol.degree
    ## }

    ## TO BE REMOVED
    nx <- 500
    plot <- TRUE
    match <- FALSE
    avg.opt <- "movav"
    ny.aut <- TRUE

    ## print(nx)
    ## print(plot)
    ## print(avg.opt)
    ## if(any("pol.degree" %in% ls())) print(pol.degree)

    as.real <- function(x) {
        return(as.double(x))
    }
    scale <- function(x) {
        mn <- min(x)
        scl <- (cstored-1)/(max(x)-min(x))
        return(1+(x-mn)*scl)
    }
    refine <- function(vec) {
        for (jj in 1:length(vec)) {
            vec <- round(vec)
            tmp <- which.max(kin[(vec[jj]-round(sd.kin/2)):(vec[jj]+round(sd.kin/2))]) 
            vec[jj] <- tmp+vec[jj]-round(sd.kin/2)-1
        }
        vec <- unique(c(1,vec,cstored))
        vec.tmp <- vec
        idx <- which(diff(vec)<sd.kin)
        for (ii in idx) {
            if (ii == 1) {
                vec.tmp <- vec.tmp[-match(vec[2],vec.tmp)]
            } else if (ii == (length(vec)-1)) {
                vec.tmp <- vec.tmp[-match(vec[length(vec)-1],vec.tmp)]
            } else {
                kk <- which.min(c(kin[vec[ii]], kin[vec[ii+1]]))-1
                vec.tmp <- vec.tmp[-match(vec[ii+kk],vec.tmp)]
            }
        }
        vec.tmp <- vec.tmp[which(vec.tmp != 1)]
        vec.tmp <- vec.tmp[which(vec.tmp != cstored)]
        return(vec.tmp)
    }

    ## INPUT FILE 
    progind <- read.table(file)[, c(1, 3, 4)]
    colnames(progind) <- c("PI", "Time", "Cut")
    cstored <- dim(progind)[1]

    ## Dynamic parameters
    perc <- 0.90 ## Percentage to find optimal nbiny
    xidx <- 0.5 ## Joining consecutive cells on a single raw
    K <- 2 ## Parameter of the Haat wavelet
    dpeaks.dyn <- 3 ## Parameter of the maximum search algorithm (window size)
    nsample <- 50 ## Number of samples in the distribution of reshuffled Hellinger distances 
    cutjoin <- 10 ## Number of null joining attempts required to quit the procedure 
    conf.lev <- 0.005 ## Cut on pvalues of the Grubb Test
    ## Kinetic parameters
    wsize <- round(2*cstored/nx) + ((round(2*cstored/nx)+1) %% 2)  ## To make it odd
    dpeaks.kin <- ceiling(wsize/2) 
    thr.ratio <- 0.05 ## Parameter of second data cleaning  
    ## Matching parameters
    lx <- round(cstored/nx)
    sd.kin <- 1*lx
    sd.dyn <- 1*lx
    
#########&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&###############
############################&&&&&& DYNAMICS &&&&&&&###############################
#########&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&###############

    datap <- matrix(c(progind[,1],progind[,2]), ncol=2, nrow=cstored)

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
        expmodel <- nls(lthbin ~ a-b*exp(seqbin*c), start=list(a=a.0, b=exp(coef(expmodel.0)[1]), c=coef(expmodel.0)[2] ) )
        fracmodel <- nls(lthbin ~ a*(seqbin/(l+seqbin)), start=list(a=max(lthbin), l=(max(lthbin)-min(lthbin)/2)))

        ## Chi Square Test
        gdl.exp <- length(seqbin)-3
        gdl.frac <- length(seqbin)-2
        theor.exp <- coef(expmodel)["a"]-coef(expmodel)["b"] * exp(coef(expmodel)["c"] * seqbin)
        theor.frac <- coef(fracmodel)[1]*(seqbin/(coef(fracmodel)[2]+seqbin))
        X2.exp.rid <- sum((theor.exp-lthbin)^2/theor.exp)/gdl.exp
        X2.frac.rid <- sum((theor.frac-lthbin)^2/theor.frac)/gdl.frac
        print(paste("X2rid exp is", X2.exp.rid, "X2rid frac is", X2.frac.rid))
        can <- coef(fracmodel)[2]*perc/(1-perc)
        if (X2.exp.rid < X2.frac.rid) {
            print("Exponential model chosen")
            can <- (1/coef(expmodel)["c"])*log((coef(expmodel)["a"]/coef(expmodel)["b"]) *(1-perc))
        } else {
            print("Hyperbola model chosen")
            can <- coef(fracmodel)[2]*perc/(1-perc) 
        }
        ny <- as.numeric(round(can))
    }
    
    ## dev.new()
    ## plot(seqbin, lthbin, xlab="nbin y", ylab="StdDev Stretches",frame.plot=FALSE, type="p",  lwd=1, pch=16, lty=1, cex=0.7)
    ## curve(coef(expmodel)["a"]-coef(expmodel)["b"] * exp(coef(expmodel)["c"] * x), lwd=2, col="darkgreen", add=TRUE)
    ## curve(coef(fracmodel)[1]*(x/(coef(fracmodel)[2]+x)), lwd=2, col="blue", add=TRUE)
    ## axis(1, at=round(can), lwd=2, col="red")
    ## #text(can, max(lthbin)*1.2, paste(sep="",as.character(perc*100),"%"))
    ## abline(v=can, lwd=2, col="red")

##################################################################################
    ## 2-D HISTOGRAM
##################################################################################

    print(paste("Nbins on y is", ny))
    hist <- hist2d(datap, nbins=c(nx,ny), show=FALSE)

##################################################################################
    ## STRETCHES CREATION
##################################################################################
    # Joining cells through density criterion: the larger values of xids the larger tolerance 
    joinx <- hist$counts
    joinx[,] <- 0
    for (j in 1:ny) {
        for(i in 1:(nx-1)) {
            if(hist$counts[i,j] < 1) next #mettere un minimo?
            if((hist$counts[i,j] >= (hist$counts[i+1,j]*(1-xidx))) & (hist$counts[i,j] <= hist$counts[i+1,j]))  joinx[i,j] <- 1
            if((hist$counts[i+1,j] >= (hist$counts[i,j]*(1-xidx))) & (hist$counts[i+1,j] < hist$counts[i,j]))  joinx[i,j] <- 1
                                        #Joining if separated by an empty bin as well
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
        else  rawset[j,c(1:3)] <- c(hist$x[minr], hist$x[maxr], sum(hist$counts[c(minr:maxr),j]*hist$x[minr:maxr])/sum(hist$counts[c(minr:maxr),j]))
    }
    dum <- data.frame(x1=rawset$max-rawset$center, x2=rawset$center-rawset$min)
    rawset$wth1 <- (pmin(dum$x1,dum$x2)/pmax(dum$x1,dum$x2))*(1-((rawset$max-rawset$min)/cstored))
    rawset$wth2 <- (pmin(dum$x1,dum$x2)/pmax(dum$x1,dum$x2)) / (rawset$max-rawset$min)
    rawset$wth1[which(is.na(rawset$wth1))] <- 0
    rawset$wth2[which(is.na(rawset$wth2))] <- 0
    rawsetsort <- rawset[order(rawset$center),]

#################################################################################
    ## WEIGHTED SUMS FUNCTIONS
##############################################################################
    #Forward
    sumwr <- function(first,meanopt,wthopt) {
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
        for (i in 1:nx) {
            idxlist <- which(rawfirst[idrs:ny] < hist$x.breaks[i+1]) + idrs - 1 
            if(length(idxlist)==0) sm[i] <- 0
            else sm[i] <- sum(rawwth[idxlist]*rawmean[idxlist])/sum(rawwth[idxlist])
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
            else sm[i] <- sum(rawwth[idxlist]*(cstored-rawmean[idxlist]))/sum(rawwth[idxlist])
        }
        return(sm)
    }

#########################################################################################
###FILTER SECTION
########################################################################################
    matchfilt <- c(rep(1,K),1,rep(-1,K))

#################################################################################
### Forward Filter with Haar wavelet
    MaxHaar <- function(first,meanopt,wthopt) {
        filt <- filter(as.ts(sumwr(first,meanopt,wthopt)), matchfilt, method="convolution", sides=2)
        filt <- c(rep(0,K),filt[c((K+1):(length(filt)-K))],rep(0,K))  #Filt with 0 instead of NA
        xbr <- hist$x[which(peaks(filt,dpeaks.dyn))]
        return(c(1,xbr,cstored))
    }
#################################################################################
### Backward Filter with Haar wavelet
    BackMaxHaar <- function(first,meanopt,wthopt) {
        filt <- filter(as.ts(rev(backsumwr(first,meanopt,wthopt))), matchfilt, method="convolution", sides=2)
        filt <- c(rep(0,K),filt[c((K+1):(length(filt)-K))],rep(0,K))  #Filt with 0 instead of NA
        xbr <- hist$x[which(rev(peaks(filt, dpeaks.dyn)))]
        return(c(1,xbr,cstored))
    }

#################################################################################
###CHOICE OF THE BREAKS
#################################################################################
    breaks.max <- MaxHaar("max","min",1) 
    breaks.min <- BackMaxHaar("min","max",1)

###############################################################################
    ##HARD BREAKS: Joining selected breaks.min and breaks.max with res
###############################################################################
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

#################################################################################
#### HISTOGRAM of each PARTITION and Max density 
################################################################################
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
        #dens.tot[,i] <-  brkjy.tot[,i]/(breaks.tot[i+1]-breaks.tot[i])
        dens.tot[,i] <-  brkjy.tot[,i]/(sum(brkjy.tot[,i]))
    }

########################################################################
    ##JOINING PARTITIONS METHODS
########################################################################
    #Computation of Distances Hell and Kolm between consecutive partitions
    distHell.tot <- NULL
    for(j in 1:(length(breaks.tot)-2)) {
        prova1 <- DiscreteDistribution(supp = c(1:ny) , prob=dens.tot[,j])
        prova2 <- DiscreteDistribution(supp = c(1:ny) , prob=dens.tot[,j+1])
        distHell.tot[j] <- HellingerDist(prova1,prova2)
    }

########################################################################################
###MAIN JOINING procedure: comparison with training uniform samples
########################################################################################
    lstHell.tot <- sort.int(distHell.tot, index.return=TRUE)$ix
    discbreaks.tot <- NULL
    flagbreak <- 0
    ll <- 0
    for (i in lstHell.tot) { 
        sampleHell.tot <- NULL
        l1 <- which(hist$x==breaks.tot[i])+1 
        if (i==1) l1 <- 1
        l2 <- which(hist$x==breaks.tot[i+1])
        l3 <- which(hist$x==breaks.tot[i+2])
        if (i==length(lstHell.tot)) l3 <- nx
        print(paste("Dyn Break Number ", i, ", Borders", round(breaks.tot[i]),round(breaks.tot[i+1]),round(breaks.tot[i+2]))) 
        for (idx in 1:nsample) {
            unif1 <- matrix(rep(0,ny*(l2-l1+1)), nrow=ny)
            unif2 <- matrix(rep(0,ny*(l3-l2)), nrow=ny)
            for (j in 1:ny) {
                ncounts <- sum(hist$counts[c(l1:l3),j])
                if(ncounts==0) next
                dum <- sample(1:(l3-l1+1), ncounts, replace=TRUE)
                dumhist <- hist(dum, breaks=seq(from=0.5,to=(l3-l1+1)+0.5,by=1), plot=FALSE)
                unif1[j,] <- dumhist$counts[c(1:(l2-l1+1))]
                unif2[j,] <- dumhist$counts[c((l2-l1+2):(l3-l1+1))]
            }
            dens1 <- rowSums(unif1)/sum(unif1)
            dens2 <- rowSums(unif2)/sum(unif2)
            part1 <- DiscreteDistribution(supp = c(1:ny) , prob=dens1)
            part2 <- DiscreteDistribution(supp = c(1:ny) , prob=dens2)
            sampleHell.tot[idx] <- HellingerDist(part1,part2)
        }
        grubbsHell.tot <- grubbs.test(c(sampleHell.tot,distHell.tot[i]), type=10)
                                        #print(paste(grubbsHell.tot$p.value))
        if (distHell.tot[i] < max(sampleHell.tot) | grubbsHell.tot$p.value>conf.lev) {
            print(paste("Joining partitions"))
            flagbreak <- 0
            ll <- ll+1
            discbreaks.tot[ll] <- breaks.tot[i+1]
        }
        else flagbreak <- flagbreak+1
        if (flagbreak==cutjoin) break
                                        #print(" ")
    }

    brk.dyn <- sort(breaks.tot[-match(discbreaks.tot, breaks.tot)])

#########&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&###############
############################&&&&&& KINETIC &&&&&&&###############################
#########&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&###############

    parabol <- 2*progind$PI*(cstored-progind$PI)/cstored
    parabol.log <- -log(parabol/cstored)
    red <- parabol.log
    kin <- -log(progind$Cut / cstored)-red
    kin <- replace(kin, which(abs(kin)==Inf), rep(tail(kin,2)[1], length(which(abs(kin)==Inf))) )

### SMOOTHING Section
    if(avg.opt=="SG") {
        kin.mv <- savitzkyGolay(kin, p=pol.degree, w=wsize, m=0)
    } else kin.mv <- movav(kin, w=wsize)
    kin.mv <- c(kin[c(1:((wsize-1)/2))], kin.mv, tail(kin,(wsize-1)/2))

### MAXIMA search Section
    max.mv <- which(peaks(kin.mv, dpeaks.kin, strict=FALSE))

### First Cleaning on mv:: check separation between consecutive max
    max.mv.tmp <- max.mv
    rif <- max.mv[length(max.mv)]
    for (i in (length(max.mv)-1):1) {
        if (rif-max.mv[i]<dpeaks.kin/2) {
            max.mv.tmp <- max.mv.tmp[-match(max.mv[i],max.mv.tmp)] ##Remove the smallest
        }
        else rif <- max.mv[i]
    }
    max.mv <- max.mv.tmp

### Second Cleaning Attempt on mv
    adj.mv <- NULL
    adj.mv[1] <- which.min(kin.mv[1:max.mv[1]])
    for (i in 1:length(max.mv) ) {
        if (i==length(max.mv)) adj.mv[i+1] <- which.min(kin.mv[round(max.mv[i]):cstored])+round(max.mv[i])-1
        else adj.mv[i+1] <- which.min(kin.mv[round(max.mv[i]):round(max.mv[i+1])])+round(max.mv[i])-1 
        amax <- mean( c(kin.mv[max.mv[i]]-kin.mv[adj.mv[i]], kin.mv[max.mv[i]]-kin.mv[adj.mv[i+1]]) )
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

#########&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&###############
############################&&&&&& MATCHING &&&&&&&###############################
#########&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&###############


################################################################################
    ## MATCHING/MERGING SECTION
################################################################################
    breaks <- NULL
    ll <- 0

    vkin <- brk.kin
    vkin <- refine(vkin)

    vdyn <- brk.dyn
    vdyn <- vdyn[which(vdyn != 1)]
    vdyn <- vdyn[which(vdyn != cstored)]
    vdyn <- refine(vdyn)

##################################################################################
    ## Matching
    for (i in 1:length(vkin) ) { 
        dist <- NULL
        ## cat("********************************************\n")
        ## print(paste("Analyzing kin break n", i, " ::: ", vkin[i]))
        set <- c(tail(vdyn[which(vdyn<vkin[i])],1), head(vdyn[which(vdyn>=vkin[i])],1) )
        if (length(set)==0) break
        for (j in 1:length(set) ) {
            ## Naive criteria (binary decision)
            dist[j] <- abs(set[j]-vkin[i])
            ## cat("Comparing with ", set[j], "Distance is", dist[j],"\n")
            dyn.cand <- set[which.min(abs(set-vkin[i]))]
        }
        if (min(dist)< (sd.kin+sd.dyn)) {
            ## print(paste("Adding", vkin[i]))
            ll <- ll+1
            breaks[ll] <- vkin[i]
            vdyn <- vdyn[-match(dyn.cand,vdyn)]
        }
    }

#######################################################################
    ## Merging option 
    if(!match) {
        dist <- NULL
        ll <- 0
                                        #Looking for residual vdyn too close to any of the vkin (it should happen rarely)
        for(i in 1:length(vdyn)) {
            for(j in 1:length(vkin)) {
                ll <- ll+1
                dist[ll] <- abs(vdyn[i]-vkin[j])
            }
        }
                                        #Identifying and removing them from vdyn vector
        dist.mtx <- matrix(dist<sd.kin, nrow=length(vdyn), ncol=length(vkin), byrow=TRUE)
        if(any(dist.mtx)) {
            near <- unique(which(dist.mtx==TRUE, arr.ind=TRUE)[,1])
                                        #print(paste("removing",near,vdyn[near]))
            vdyn <- vdyn[-near]
        }
                                        #Merging residual vdyn with breaks and vkin
        brk.mtc <- breaks
        breaks <- sort(unique(c(breaks,vdyn,vkin)))
    }

######################################################################
    ## OUTPUT section
#####################################################################
    print("Breaks are ")
    print(breaks)

    vec <- sort(breaks)
    state <- NULL
    for (i in 1:(length(vec)+1)) {
        if (i==1) ib <- 0
        else ib <- vec[i-1]
        if (i==length(vec)+1) fb <- cstored
        else fb <- vec[i]
        ## cat(i,ib,fb,fb-ib,"\n")
        state <- c(state,rep(i,fb-ib))
    }
    output.match <- data.frame(PI=progind$PI, Time=progind$Time, State=state)
    write.table(output.match, file=gsub("PROGIDX", "BASINS", file), row.names=FALSE, quote=FALSE)

####################################################################
    ## PLOT Section
####################################################################

    if(plot == TRUE){
        dev.new(width=15, height=10)
        par(mgp=c(0, 0.4, 0))
        par(ps=6)
        par(mar = c(3.5, 0, 1, 3), oma = c(2, 4, 2, 2))
        margin <- 100
        cx <- 2.2
        xr1 <- c(margin:(cstored-margin))
        xr <- c(1,cstored)
        yr <- range(progind$Time)

        plot(0, 0, main="", xlim=xr, ylim=yr, xlab="", ylab="", axes=FALSE, frame.plot=FALSE, type="n")

        kin.pl <- -log(progind$Cut[xr1] / cstored)
        xx.lab <- c(1,round(breaks),cstored)
        axis(1, at=xx.lab, tck=.01, cex.axis=1.8)
        axis(3, labels=rep("", length(xx.lab)), at=xx.lab, tck=.01)
        mtext("Progress Index", side=1, line=2.4, cex=cx )

        yy.lab1 <- format(c(min(kin.pl), min(kin.pl)+(max(kin.pl)-min(kin.pl))*c(1:3)/3), digits=2)
        axis(2, labels=yy.lab1, at=scale(as.numeric(yy.lab1)), las=3, tck=.01, cex.axis=cx)
        mtext(expression("ln(("*italic(tau["SA"]+tau["AS"])*")/2)"), side=2, line=1.8, cex=cx)

        yy.lab2 <- c(1, c(1:5)/5*cstored)
        axis(4, labels=rep("", length(yy.lab2)), at=yy.lab2, las=2, tck=.01, hadj=-0.6, col="red")
        mtext(yy.lab2, side=4, las=2, line=0.2, at=yy.lab2, col="red", cex=1.8)
        mtext("Time", side=4, line=2.1, cex=cx, col="red")

        points(progind$PI, progind$Time, cex=0.05, col="red")
        lines(xr1, scale(kin.pl), lwd=1, col="black")
        if(!match) {
            ## Not-matched vdyn
            abline(v=vdyn, lwd=0.4, col="blue")
            ## Not-matched vkin
            abline(v=vkin[-match(brk.mtc, vkin)], lwd=0.4, col="orange4")
            ## Matched breaks
            abline(v=brk.mtc, lwd=0.7, col="black")
        } else abline(v=breaks, lwd=0.7, col="black")
    }

print(proc.time()-ptm)

#}


