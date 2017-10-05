#' @title Building a Markov State Model
#' @description 
#'    \code{msm} uses the discretized trajectory (state sequence) to infer the transition matrix of the process. Further values are extracted, like the eigenvalues, eigenvectors, non-Markovian flux and metastability. If requested, a Chapman-Kolmogorv test is carried out to check the Markovianity of model.
#' @param seq It can be provided (1) A file name containing only one row with the state sequence, (2) a file type BASINS_<...>, i.e. the output file of the \code{basins_recognition} function, (3) the state sequence itself, (4) the count matrix.
#' @param lag The lag time in terms of number of snapshots. Default to 1.
#' @param tm.opt How to infer the transition matrix. Using \code{"symm"} we derive the transition matrix simply symmetrizing the count matrix (and normalizing it). \code{"mle"} uses the maximum probability estimator for a reversible transition matrix (see references).
#' @param eig.plot Logical value indicating whether to plot the eigenvalues spectrum or not. Default to FALSE.
#' @param CK.test Logical value indicatig whether to perform the Chapman-Kolmogorov test or not. The number of points computed is 10. We recall that, if a count matrix is provided as \code{seq}, the test cannot be performed. Default to FALSE.
#' @param CK.plot Logical value indicating whether to plot the results of the CK test spectrum or not. Default to TRUE.
#' @param setA In case \code{CK.test} is true, this vector indicates the list of states that assemble the reference (macro)state A used in the CK test. Default to the largest state. 
#' @param silent A logical value indicating whether the function has to remain silent or not. Default value is \code{FALSE}.
#' @param dev.new.eig Logical value indicating whether, in case an eigenvalues plot is requested, a new device has to open or not. Default to TRUE.
#' @param dev.new.CK Same as above but regarding the CK test plot. Default to TRUE.
#'
#' @return A list containing
#'    \itemize{
#'       \item \code{TM} The transition matrix
#'       \item \code{cnt} The count matrix 
#'       \item \code{values} The eigenvalues of the transition matrix
#'       \item \code{left} The left eigenvectors of the transition matrix, normalized as in [1].
#'       \item \code{right} The right eigenvectors of the transition matrix, normalized as in [1].
#'       \item \code{Q} The metastability, i.e. the trace of the transition matrix.
#'       \item \code{flux} Fluxes relative to each timescale. see ref. [2].
#'       \item \code{mtp} Mean transition time between different states (see reference [2])
#'       \item \code{pMSM} A vector of probabilities to be at set A, starting from A itself, aftera certain number of steps, calculated by means of the transition matrix (see reference [1]).
#'       \item \code{pMD} A vector of probabilities to be at set A, starting from A itself, aftera certain number of steps, calculated directly from the trajectory (see reference [1]).
#'       \item \code{epsMD} The statistical errors associated to the \code{pMD} (see reference [1]).
#'       \item \code{chi} The reduced chi-square value of the CK test.
#'       \item \code{n.st} The number of states. 
#'       \item \code{seq.st} The state sequence.
#'       \item \code{tab.st} A data.frame containing the state number and its relative weigth.
#'       \item \code{cstored} Number of snapshots of the trajectory. 
#'       \item \code{call} The matched call of the function.
#'}
#'
#' @examples 
#' \dontrun{
#'     msm("BASINS_1.dat", 100, tm.opt="mle", CK.test=TRUE, setA=c(1,2))
#' }
#'
#' @details The main reference for the Markov state model creation is [1][Prinz, J.H., et al. "Markov models of molecular kinetics: Generation and validation." \emph{The Journal of chemical physics} 134.17 (2011): 174105]. [2] [Noe' F., et. al. "Hierarchical analysis of conformational dynamics in biomolecules: Transition metwork of metastable states". \emph{The Journal of Chemical Physics} 126, (2007): 155102.]
#'
#' @importFrom graphics abline axis legend lines mtext par points text
#' @importFrom stats as.ts coef filter lm nls sd weighted.mean
#' @importFrom utils head tail
#' @importFrom data.table fread fwrite
#' @importFrom Hmisc errbar
#' @importFrom expm %^%
#' @export msm


msm <-  function(seq, lag=1, tm.opt=c("symm", "mle"), eig.plot=FALSE, CK.test=FALSE, CK.plot=TRUE, setA=NULL, silent=FALSE, dev.new.eig=TRUE, dev.new.CK=TRUE)  {
    call <- match.call()
    cnt.mode <- FALSE
    
    lt <- function(x) {
        return(length(x))
    }
    symm <- function(x) {
        return((x+t(x))/2)
    }
    tscale <- function(x) {
        return(-lag/log(x))
    }

    if(!is.numeric(lag)) stop("Please provide an integer value of lag")
    if(lt(tm.opt)==2 || !(tm.opt %in% c("symm", "mle"))) {
        warning("tm.opt value absent or ill-defined, set to default value \"symm\"", immediate.=T)
        tm.opt <- "symm"
    }

    if(is.character(seq)) {
        data.bas <- data.frame(fread(seq))
        cstored <- nrow(data.bas)
        if(ncol(data.bas)!=1) {
            ##answer <- readline(paste("File", seq, "has", ncol(data.bas), "columns, \nis it a BASINS type file? (y|n)\t"))
            if(grepl("BASIN", seq)) {
                colnames(data.bas) <- c("PI", "Time", "State")
                seq.st <- data.bas[order(data.bas$Time),]$State
            } else {
                answer <- readline(paste("File", seq, "has", ncol(data.bas), "columns, \n which one is the state sequence? ")) 
                seq.st <- data.bas[, as.numeric(answer)]
            }
        } else {
            seq.st <- data.bas[,]
        }
    } else if(!is.matrix(seq)) {
        seq.st <- unlist(seq)
    } else if(nrow(seq)==ncol(seq)) {
        cnt.mode <- TRUE
    } else stop("Input \"seq\" not recognized")

    if(!cnt.mode) {
        n.st <- lt(unique(seq.st))
        cstored <- lt(seq.st)
        tmp <- rle(sort(seq.st))
        tab.st <- data.frame(n=tmp$values, length=tmp$lengths)
        if(CK.test) {
            if(!all((setA) %in% tab.st$n)) stop("State in setA not existing")
            if(is.null(setA)) {
                setA <- tab.st$n[which.max(tab.st$length)[1]]
                if(!silent) cat(paste("setA is not defined.\n  Default to the largest state", setA, "\n"))
            }
        }
        if(!all(tab.st$n==c(1:n.st))) stop("State are not numbered from 1 to n.st")
    } else {
        cnt <- seq
        cnt.v <- rowSums(cnt)
        n.st <- nrow(cnt)
        if(CK.test) {
            warning("CK test cannot be performed if you provide only the count matrix")
            CK.test <- FALSE
        }
    }

    if(!silent) cat("Number of states", n.st, "\n")    
    
    
    
########################################################################
    if(!cnt.mode) {
        ## Count Matrices each one for a different lag time
        ## Overlapping window count method
        cnt <- matrix(0, nrow=n.st, ncol=n.st)
        for (i in 1:(cstored-lag)) {
            row <- seq.st[i]
            col <- seq.st[i+lag]
            cnt[row,col] <- cnt[row,col]+1
        }
        cnt.v <- rowSums(cnt)
        if(!all(is.finite(cnt))) stop("ERROR in Window Count Matrices")

        ## Cleaning of the cnt; very rare, only if a state appear solely in the last lag snaps
        st.del <- which(cnt.v==0)
        if(lt(st.del)!=0) {
            warning(paste("State", st.del, "empty in the count matrix"), immediate.=T)
            cnt <- cnt[-st.del,-st.del]
            n.st <- nrow(cnt)
        }
        cnt.v <- rowSums(cnt)
    }
    
#######################################################################
    ## Naive Transition Matrix (Simple)
    if(tm.opt=="symm") {
        tm <- t(apply(symm(cnt), 1, function(x) x/sum(x) ))
        if(!all(is.finite(tm))) stop("ERROR in simple TM")
    } else if (tm.opt=="mle") {
        xcnt <- 2*symm(cnt)
        xcnt.v <- rowSums(xcnt)

        ## (2) Main Algorithm
        tm <- matrix(NaN, nrow=n.st, ncol=n.st)
        for(fl in 1:20) {
            ## (2.1)
            for (i in 1:n.st) {   
                if ((cnt.v[i]-cnt[i,i]) > 0) {
                    xcnt[i,i] <- cnt[i,i]*(xcnt.v[i]-xcnt[i,i])/(cnt.v[i]-cnt[i,i])
                }
                else xcnt[i,i] <- 0
            }
            xcnt.v <- rowSums(xcnt)
            ## (2.2)
            for (i in 1:(n.st-1)) {  
                for (j in (i+1):n.st) {
                    a <- cnt.v[i]-cnt[i,j]+cnt.v[j]-cnt[j,i]
                    b <- cnt.v[i]*(xcnt.v[j]-xcnt[i,j]) + cnt.v[j]*(xcnt.v[i]-xcnt[i,j])-(cnt[i,j]+cnt[j,i])*(xcnt.v[i]+xcnt.v[j]-2*xcnt[i,j])
                    c <- -(cnt[i,j]+cnt[j,i])*(xcnt.v[i]-xcnt[i,j])*(xcnt.v[j]-xcnt[i,j])
                    if (a==0) {  ## Internal update OK
                        xcnt[i,j] <- 0
                        xcnt[j,i] <- 0
                    } else {     ## Internal update OK
                        xcnt[i,j] <- (-b+sqrt(b^2-(4*a*c)))/(2*a)
                        xcnt[j,i] <- xcnt[i,j]
                    }
                }
            }
            ## Second update
            xcnt.v <- rowSums(xcnt)
            ## (2.3) Updating Transition matrix
            for (i in 1:n.st) {
                if (xcnt.v[i] != 0) tm[i,] <- xcnt[i,]/xcnt.v[i]
                else tm[i,] <- rep(0,n.st)
            }
        } 
        if(!all(is.finite(tm))) stop("ERROR in MLE TM")
        if(any(is.na(tm)) | any(tm<0)) stop("ERROR in MLE TM")
    }

##############################################################################

    ## Eigenvalues and stationary distibution
    ##     tmp <- eigen(t(tm)) ## Left eigenvectors
    ##     if (!(max(Re(tmp$values))<1.01 & max(Re(tmp$values))>0.99)) stop("ERROR in MLE Eigenvalues    eig.val <- Re(sort(L$values, decreasing=TRUE))
    ## ")
    ##     eigs <- Re(sort(tmp$values,decreasing=TRUE))
    ##     eig.vec <- Re(tmp$vectors[,order(Re(tmp$values), decreasing=TRUE)])
    ##     eig.vec.right <- Re(eigen(tm)$vectors[,order(Re(tmp$values), decreasing=TRUE)])
    ##     stat <- Re(tmp$vectors[,which.max(Re(tmp$values))])/Re(sum(tmp$vectors[,which.max(Re(tmp$values))]))
    ##     if (sum(stat)<0) stat <- -stat

    ## R <- eigen(tm) ## Right Eigenvectors
    L <- eigen(t(tm)) ## L eigenvectors
    if (!(max(Re(L$values))<1.01 & max(Re(L$values))>0.99)) stop("ERROR in MLE Eigenvalues")
    eig <- Re(sort(L$values, decreasing=TRUE))
    ## stat <- Re(L$vectors[,which.max(Re(L$values))])/Re(sum(L$vectors[,which.max(Re(L$values))]))
    ## if (sum(stat)<0) stat <- -stat
    L <- Re(L$vectors[,order(Re(L$values), decreasing=TRUE)])
    ## R <- Re(R$vectors[,order(Re(R$values), decreasing=TRUE)])

    L <- apply(L, 2, function(x) {
        if(sum(x)<0) return(-x)
        else return(x)
    })
    stat <- L[,1]/sum(L[,1])
    if (sum(stat)<0) stat <- -stat
    ## R <- apply(R, 2, function(x) {
    ##     if(sum(x)<0) return(-x)
    ##     else return(x)
    ## })
    ## for(kk in 1:nrow(R)) if(sum(R[,kk]*L[,kk]) < 0) R[,kk] <- -R[,kk]


    ## Left eigenvectors properly normalized
    A.L <- apply(L, 2, function(x) sqrt(sum(x^2/stat)))
    L.tmp <- apply(L, 1, function(x) x/A.L)
    if(identical(dim(L.tmp), dim(L))) L <- t(L.tmp)
    else L <- t(L.tmp)
    rm("L.tmp")

    ## Right eigenvectors properly normalized
    ## A.R <- apply(R, 2, function(x) sqrt(sum(x^2*stat)))
    ## R.tmp <- apply(R, 1, function(x) x/A.R)
    ## if(identical(dim(R.tmp), dim(R))) R <- t(R.tmp)
    ## else R <- t(R.tmp)
    ## rm("R.tmp")

    R <- apply(L, 2, function(x) x/stat)
    if(!identical(dim(R),dim(L))) R <- t(R)

    if(abs(mean(sapply(c(1:nrow(L)), function(i) sum(L[,i]*R[,i])))-1) > 0.0001) stop("Error in left and right eigenvalues")

#########################################################################
    ## Calculating fluxes ...
    flux <- colSums(L^2)

##############################################################################
    ## Calculating the mean transition time
    mtp <- lag/sum(as.vector(stat %*% (tm-(diag(nrow(tm))*diag(tm)))))
    
##############################################################################
    if(CK.test) {
        ## State tab
        wA <- replace(stat, c(1:n.st)[-setA], 0)/sum(stat[setA])

        ## Calculating p_MSM(A,A)  
        pMSM <- NULL
        np.CK <- 10
        for (i in 1:np.CK) {
            tmp <- wA %*% (tm %^% i)
            pMSM[i] <- sum(tmp[setA])
        }

        ## Additional cnts
        cnt.CK <- list(cnt)
        if(!silent) cat("Computing the count matrices for the CK test...\n")
        for (ll in 2:np.CK ) {  
            set <- c(1:(cstored-ll*lag)) 
            cnt.CK[[ll]] <- matrix(0, nrow=n.st, ncol=n.st)
            for (i in 1:lt(set)) {
                row <- seq.st[set[i]]
                col <- seq.st[set[i]+lag*ll]
                cnt.CK[[ll]][row,col] <- cnt.CK[[ll]][row,col]+1
            }
        }
        cnt.v <- lapply(cnt.CK, rowSums)
        if(!all(sapply(cnt.CK, function(x) all(is.finite(x))))) stop("ERROR in Window Count Matrices")
        ## Cleaning
        st.del <- unique(sort(unlist(lapply(cnt.v, function(x) which(x==0)))))
        if(lt(st.del)!=0) {
            warning(paste("State", st.del, "empty in the count matrix during CK test"), immediate.=T)
            for(ll in c(1:np.CK)) cnt.CK[[ll]] <- cnt.CK[[ll]][-st.del,-st.del]
            n.st <- unique(sapply(cnt.CK, nrow))
            cnt.v <- lapply(cnt.CK, rowSums)
        }
        
        ## Calculating p_MD(A,A), and eps_MD(A,A)
        pMD <- NULL    
        for (ll in 1:np.CK) {
            if (lt(setA)>1) pi <- rowSums(cnt.CK[[ll]][,setA])/rowSums(cnt.CK[[ll]])
            else pi <- cnt.CK[[ll]][,setA]/rowSums(cnt.CK[[ll]])
            pMD[ll] <- sum(wA[setA]*pi[setA])  ## Which wA should I use?
        }
        epsMD <- NULL
        for (ll in 1:np.CK) epsMD[ll] <- sqrt(ll*lag*(pMD[ll]-pMD[ll]^2)/sum(cnt.v[[ll]][setA]) )
        
        ## Calculating chi2 deviations
        chi <- sum(((pMD-pMSM)/epsMD)^2)/np.CK 
    }

    #############################################################################
    ## Plot Eigenvalues
    if(eig.plot) {
        if(dev.new.eig) dev.new(width=7, height=7)
        ng <- c(1:min(20, n.st))
        xr <- range(ng)
        yr <- c(0,max(eig))
        plot(ng, eig[ng], main="Eigenvalues", xlim=xr, ylim=yr, xlab="", ylab="Eigenvalues", axes=TRUE, col="blue", pch=19)
        for(i in ng) lines(rep(i,2), c(0, eig[i]))
        mtext(text=paste("First", lt(ng), "eigenvalues out of", lt(eig)), side=1, line=2)
    }

    ## #############################################################################
    ## ## Plot Chapman-Kolmogorov Test
    if(CK.test && CK.plot) {
        if(dev.new.CK) dev.new(width=9,height=7)
        par(mar=c(5.1, 4.1, 2.1, 2.1))
        xr <- range(lag*c(1:np.CK))
        ##yr <- range(c(pMD+epsMD, pMD-epsMD, pMSM))
        ## yr <- c(0, 1)
        yr <- c(max(0, min(c(pMD-epsMD, pMSM))),1)
        main <- NULL
        main2 <- NULL
        for (i in 1:lt(setA)) main <- paste(main,as.character(setA[i]), sep=" ")

        ## plot(0, 0, main=paste("Chapman-Kolmogorv Test :: setA={", main,"}, lag=", lag,  sep=""), xlim=xr, ylim=yr, xlab="Lag Time", ylab="Probability ", axes=TRUE, frame.plot=TRUE, type="n")
        plot(0, 0, main="", xlim=xr, ylim=yr, xlab="Lag Time", ylab="Probability ", axes=TRUE, frame.plot=TRUE, type="n")
        points(lag*c(1:np.CK), pMSM, type="o", cex=1.5, pch=16, col="red")
        errbarHmisc <- get(x="errbar", pos="package:Hmisc")
        errbarHmisc(lag*c(1:np.CK), pMD, pMD+epsMD, pMD-epsMD, pch=4, lwd=0.8, cex=1.5, col='black', add=TRUE)
        legend((2/3)*(np.CK*lag), yr[2]*0.97, legend=c(paste("MSM lag=", lag), "MD"), col=c("red", "black"), bty="n", pch=c(16,4), cex=1.2)
        lab <- NULL
        text((2/3)*(np.CK*lag), yr[2]*0.99, labels=paste("Chi2 =", chi), font=2)
    }

    if(!CK.test) {
        chi <- NaN
        pMSM <- NaN
        pMD <- NaN
        epsMD <- NaN
    }
    if(cnt.mode) {
        seq.st <- NaN
        tab.st <- NaN
        cstored <- NaN
    }
    
    if(!silent) cat("Done.\n")
    invisible(list(TM=tm, cnt=cnt, values=eig, left=L, right=R, Q=sum(diag(tm)), flux=flux, mtp=mtp, pMSM=pMSM, pMD=pMD, epsMD=epsMD, chi=chi, n.st=n.st, seq.st=seq.st, tab.st=tab.st, cstored=cstored, call=call))

}
