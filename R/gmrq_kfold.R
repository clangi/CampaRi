#' @title GMRQ score from McGibbon RT et al. (2015)
#' @description
#'      Function to calculate the GMRQ score. Input are the multidim trajectory,vector of GMRQ values, an ATOMICS lag
#'      
#' @param traj trajectory matrix.
#' @param gmrq number of eigenvalues.
#' @param lag number of lags in the Markov state model (MSM)
#' @param kfolds cross-validation number of folds
#' @param clust.method clustering method. Available is kmeans sbr and nsbr
#' @param tm.opt definition for matrix creation (C and S)
#' @param clust.param the cluster parameters to optimize - it is a matrix with parameters on the columns and 
#' different set of them in the rows (it loops on the rows)
#' @param preproc vars for the campari run
#' @param chunks if \code{TRUE} it splits consecutive chunks of data. It is needed for time-series data.
#' @param dir.out directory for files
#' @param plot if \code{TRUE} plots.
#' @param return_plot if \code{TRUE} return the plot.
#' @param silent if \code{TRUE} print less strings.
#' @param ... mainly debugging options
#' @return A list containing tot
#' @details For details regarding the SAPPHIRE plot, please refer to the relative publications \url{http://www.nature.com/articles/srep06264}. 
#' Main documentation of the original campari software \url{http://campari.sourceforge.net/documentation.html}.
#' 
#' @importFrom data.table fread fwrite
# @importFrom dbscan dbscan
#' @importFrom ClusterR MiniBatchKmeans
#' @importFrom ClusterR predict_KMeans
#' @importFrom MASS ginv
#' @importFrom fields rdist
#' @import ggplot2
#' 
#' @export gmrq.kfold
gmrq.kfold <- function(traj, gmrq=c(2:6), lag=1, kfolds=5, clust.method=c("kmeans", "sbr", "nsbr", "tbc", "dbscan"), tm.opt = NULL,
                       clust.param=NULL, preproc = NULL, chunks=FALSE, dir.out=FALSE, plot = FALSE, return_plot = FALSE, silent=FALSE, ...) {
  ###################################################################################
  ## Function to calculate the GMRQ score. Input are the multidim trajectory,     ###
  ## vector of GMRQ values, an ATOMICS lag,                                       ###
  ###################################################################################
  
  avail_tm.opt <- c("symm", "mle")
  avail_clust.method <- c("kmeans", "sbr", "nsbr", "tbc", "dbscan")
  stopifnot(is.logical(chunks), is.logical(plot), is.logical(return_plot), is.logical(dir.out), is.logical(silent))
  stopifnot(.isSingleInteger(kfolds), .isSingleInteger(lag), all(sapply(gmrq, .isSingleInteger)))
  stopifnot(any(gmrq > 0), kfolds > 2)
  if(is.null(tm.opt)) tm.opt <- 'mle'
  if(!is.null(tm.opt) && !tm.opt %in% avail_tm.opt) stop('tm.opt must be between c("symm", "mle")')
  if(!is.null(clust.method) && !clust.method %in% avail_clust.method) stop('clust.method must be between  c("kmeans", "sbr", "nsbr")')
  if(is.null(clust.param)) stop("Provide the clusters parameters")
  if(is.character(dir.out)) if(!grepl(clust.method, dir.out)) stop("Out file names and directory doesn't coincide")
  tot <- rep(list(data.frame(matrix(NaN, ncol=1+3*kfolds, nrow=nrow(clust.param)))), .lt(gmrq))
  names(tot) <- paste0("gmrq = ", gmrq)
  for(ii in seq_along(tot)) colnames(tot[[ii]]) <- c(rep("nstates", kfolds), rep("train_score", kfolds), rep("test_score", kfolds))
  
  # Extra arguments checks
  input.args <- list(...)
  avail.extra.arg <- c('dbg_gmrq.kfold')
  
  if(!is.null(names(input.args)) && any(!(names(input.args) %in% avail.extra.arg))){
    warning('There is a probable mispelling in one of the inserted variables. Please check the available extra input arguments.')
    if(!silent) cat('!!!!!!!!! We found the following variables without a father (not between our extra input arguments) !!!!!!!!!!\n')
    if(!silent) cat(names(input.args)[!(names(input.args) %in% avail.extra.arg)], '\n')
  }
  
  if("dbg_gmrq.kfold" %in% names(input.args)) dbg_gmrq.kfold <- input.args$dbg_gmrq.kfold else dbg_gmrq.kfold <- FALSE
  
  # if(clust.method=="tbc") {
  #   birchmulti <- 14
  #   clust.param <- data.frame(preproc=rep(clust.param[[1]], birchmulti))
  # }
  
  # grid search on the parameters
  for(ll in seq(nrow(clust.param))) {
    if(!silent) cat("\n-----------------------------------------------\n")
    if(!silent) cat("Analysis started with the following parameters and using", clust.method[1], "method: \n")
    if(!silent) {
      if(ncol(clust.param) == 1) cat(colnames(clust.param), ': ', clust.param[ll,], '\n')
      else print(clust.param[ll,])
    }
    # chunking
    if(chunks) folds <- cut(seq(nrow(traj)), breaks=kfolds, labels=FALSE)
    else folds <- cut(sample(nrow(traj)), breaks=kfolds, labels=FALSE)
    
    # kfolding cross-validation of methods
    for(i in seq(kfolds)) {
      if(!silent) cat("\nAnalysing fold", i, "out of", kfolds,'folds.\n')
      # sets init
      sets <- list(train=NaN, test=NaN)
      if(clust.method[1]=="kmeans") {
        #########################################################################
        if(!identical(names(clust.param), c("k"))) stop("Provide the right parameters for the chosen method")
        if(!silent) cat("Kmeans clustering using" , clust.param$k[ll], 'clusters.')
        km <- ClusterR::MiniBatchKmeans(traj[which(folds!=i),], clust.param$k[ll], num_init=2, initializer="kmeans++")
        sets$train <- ClusterR::predict_KMeans(traj[which(folds!=i),], km$centroids)
        sets$test <- ClusterR::predict_KMeans(traj[which(folds==i),], km$centroids)
      } else if(clust.method[1]=="sbr")  {
        #########################################################################
        if(ll==1) {
          if(!silent) cat('Creating SAPPHIRE table from', kfolds - 1, '/', kfolds, 'parts of data-set using', i, 'fold...')
          .pi_profile(traj[which(folds!=i),], preproc = preproc, silent = silent)
          file.rename(from = 'PROGIDX_000000000001.dat', to = paste0('PROGIDX_000000000001_fold',i,'.dat'))
          if(!silent) cat('done. Saved to file ', paste0('PROGIDX_000000000001_fold',i,'.dat'))
        }
        # pi.file <- grep(list.files(), pattern = "PROGIDX_", value = T)[1] # taking only the first one!!
        pi.file <- paste0('PROGIDX_000000000001_fold',i,'.dat') # taking only the first one!!
        stopifnot(!is.na(pi.file), file.exists(pi.file))
        tmp <- CampaRi::basins_recognition(pi.file, nx=clust.param$nx[ll], ny=clust.param$nx[ll],
                                                  match=FALSE, dyn.check=2, avg.opt="SG", out.file=F, silent=T)
        sets$train <- tmp$seq.st
        centers <- tmp$tab.st$centers
        sets$test <- .predict.basin(traj[which(folds==i),], traj[which(folds!=i),], centers, sets$train, silent)
      } else if(clust.method[1]=="tbc")  {
        #########################################################################
        stop('tbc not available')
        # if(!identical(names(clust.param), c("preproc"))) stop("Provide the right parameters for the chosen method")
        # if(ll==1) {
        #   cat("TBC clustering on", as.character(clust.param[[1]])[ll], "\n")
        #   ## browser()
        #   tbc_clustering(traj[which(folds!=i),], preproc=as.character(clust.param[[1]])[ll], nfold=i)
        # }
        # ## browser()
        # ncstruct <- as.numeric(system(paste0("awk \'{print NF}\' ", paste0("./BPTIdata/gmrq/tbc/training_set/STRUCT_CLUSTERING_fold", i, ".dat"), " | sort -nu"), intern =T))
        # if(ncstruct != birchmulti) warning("supposed birchmulti different from the actual number of columns", immediate.=T)
        # sets$train <- fread(paste0("./BPTIdata/gmrq/tbc/training_set/STRUCT_CLUSTERING_fold", i, ".dat"), select=ll, data.table=F)[[1]]
        # logfile <- paste0("./BPTIdata/gmrq/tbc/training_set/training_", as.character(clust.param[[1]])[ll], "_fold", i, ".log")
        # centers <- extract.centers(logfile, ll) ## Return a list with cluster centres
        # sets$test <- predict.tbc(traj[which(folds==i),], traj[which(folds!=i),], centers)
      } else if(clust.method[1]=="dbscan")  {
        #########################################################################
        stop('dbscan not available')
        
        # if(ll==1 && i==1) cat("DBSCAN clustering\n")
        # if(!identical(names(clust.param), c("eps", "minpts"))) stop("Provide the right parameters for the chosen method")
        # browser()
        # db <- dbscan::dbscan(traj[which(folds!=i),], eps=clust.param$eps[ll], minPts=clust.param$minpts[ll], borderPoints=TRUE) 
        # sets$train <- db$cluster + 1
        # sets$test <- as.numeric(predict(db, traj[which(folds!=i),], newdata=traj[which(folds==i),])) + 1
        # if(.lt(unique(sets$test)) == 1) warning("most likely test set of noise..")
      } else if(clust.method[1]=="nsbr")  {
        #########################################################################
        if(ll==1) {
          if(!silent) cat('Creating SAPPHIRE table from', kfolds - 1, '/', kfolds, 'parts of data-set using', i, 'fold...')
          .pi_profile(traj[which(folds!=i),], preproc = preproc, silent = silent)
          file.rename(from = 'PROGIDX_000000000001.dat', to = paste0('PROGIDX_000000000001_fold',i,'.dat'))
          if(!silent) cat('done. Saved to file ', paste0('PROGIDX_000000000001_fold',i,'.dat'), '\n')
        }
        # pi.file <- grep(list.files(), pattern = "PROGIDX_", value = T)[1] # taking only the first one!!
        pi.file <- paste0('PROGIDX_000000000001_fold',i,'.dat') # taking only the first one!!
        stopifnot(!is.na(pi.file), file.exists(pi.file))
        if(is.list(clust.param$unif.splits[ll])) uspli <- clust.param$unif.splits[[ll]] # if I use the param as I(rep(list(c(1,2))), 5)
        else uspli <- clust.param$unif.splits[ll]                                       # case for single number usually
        tmp <- CampaRi::nSBR(pi.file, ny=clust.param$ny[ll], unif.splits=uspli, n.cluster = clust.param$n.cluster[ll],
                             pk_span=clust.param$pk_span[ll], comb_met=clust.param$comb_met[ll], 
                             silent=T, force_correct_ncl = F, return_ordered_predicted = T)
        
        sets$train <- tmp$predicted_div
        centers <- tmp$centers
        sets$test <- .predict.basin(traj[which(folds==i),], traj[which(folds!=i),], centers, sets$train, silent)
        # if(.lt(unique(sets$test)) < 2) browser()
      } else stop("Method not valid for clustering")
      
      # test clustering output
      if(any(is.na(sets))) stop("Clustering assignement of training and test set failed")
      ns <- max(sets$train) # -1 ?
      if(!silent) cat("Found number of states:", ns, '\n')
      # if(check) train.msm <- .create.trans(sets$train, ns, n.eigen=gmrq, lag=lag, tm.opt="mle")[c(1,2,3,4)]
      # else {
      rr <- try(train.msm <- .create.trans(sets$train, ns, n.eigen=gmrq, lag=lag, tm.opt=tm.opt, silent = silent)[c(3,4)])
      if(class(rr)=="try-error") browser()
      # }
      test.msm <- try(.create.trans(sets$test, ns, n.eigen=gmrq, lag=lag, tm.opt=tm.opt, silent = silent))
      if(class(test.msm)=="try-error") browser()
      s <- test.msm$stat
      for(nn in seq_along(gmrq)) {
        if(gmrq[nn] > ns) {
          tot[[nn]][ll, seq(i, by=kfolds, length.out=3)] <- rep(-1, 3)  
          # stop('Number of states is less than number of gmrq!!!')
        } else {
          train.score <- sum(train.msm$eig[seq(gmrq[nn])])
          a <- train.msm$right.vec[,seq(gmrq[nn])]
          test.score <- .tr(t(a) %*% diag(s) %*% test.msm$tm %*% a %*% MASS::ginv(t(a) %*% diag(s) %*% a))
          # if(check) {
          #   cck <- .tr(t(a) %*% diag(train.msm$stat) %*% train.msm$tm %*% a %*% MASS::ginv(t(a) %*% diag(train.msm$stat) %*% a))
          #   if(abs(cck-train.score)> 10^(-10)) {
          #     stop(paste("Score problem::: Nstates", ns, "Fold and gmrq ", i, gmrq[nn], ", Check is  ", cck-train.score))
          #   }
          # }
          tot[[nn]][ll, seq(i, by=kfolds, length.out=3)] <- c(ns, train.score, test.score)
        }
        # if(!silent) cat('Final scores found for gmrq', nn, 'in train and test sets are:', train.score, test.score, '\n')
        if(!silent) cat('GMRQ:', gmrq[nn], 'Scores (train, test):', train.score, test.score, '\n')
      }
    }
  }
  ## browser()
  ## In case write on file the results
  if(is.character(dir.out)) {
    for(nn in seq_along(gmrq)) {
      fout <- paste0(dir.out, "gmrq", gmrq[nn], '_lag', lag, '.dat')
      data.table::fwrite(tot[[nn]], file=fout, sep="\t", col.names=TRUE)
      if(!silent) cat("Written", fout, "\n")
    }
  }
  # if(dbg_gmrq.kfold) browser()
  ## Output 
  if(plot | return_plot) {
    gg <- .plot_gmrq.kfold(tot, gmrq, kfolds, skip_the_neg = TRUE, silent = silent)
  }
  if(plot) print(gg) 
  if(return_plot) invisible(list('plot' = gg, 'tot' = tot))
  else invisible(list('tot' = tot))
}

# -------------------------------------------------------------------------------------------- internal functions
.plot_gmrq.kfold <- function(tot, gmrq, kfolds, skip_the_neg = TRUE, silent = FALSE){
  dt_fin <- data.frame()
  for(i in 1:.lt(gmrq)) {
    if(any(tot[[i]]==-1, na.rm = T)) {
      if(!silent) warning('Skipping the gmrq with negative values (i.e. there were less clusters than gmrq).')
      if(skip_the_neg) next
    }
    nst <- tot[[i]][,1:kfolds]
    which_sel <- nst[,1] != -1
    nst <- nst[which_sel, ]
    train <- tot[[i]][which_sel,(kfolds+1):(2*kfolds)]
    test <- tot[[i]][which_sel,(2*kfolds+1):(3*kfolds)]
    xl <- range(nst)
    yl <- range(cbind(train, test))
    dt <- cbind(rowMeans(nst), rowMeans(train), apply(train, 1, sd), rowMeans(test), apply(test, 1, sd))
    # if(order.nst)
    dt <- cbind(as.data.frame(dt),  rep(gmrq[i], nrow(nst)))
    dt2 <- data.frame()
    for(i in sort(unique(dt[,1]))){
      dt2 <- rbind(dt2, apply(dt[dt[,1] == i,], 2, mean))
    }
    colnames(dt2) <- c('n_states', 'train_m', 'train_sd', 'test_m', 'test_sd', 'gmrq')
    dt_fin <- rbind(dt_fin, dt2)
  }
  dt_fin$gmrq <- as.factor(dt_fin$gmrq)
  pd <- "identity" # move them to the left and right # so the error bars dont overlap
  gg <- ggplot(dt_fin, aes(x=n_states)) + theme_classic() + 
    geom_errorbar(aes(ymin=train_m - train_sd, ymax=train_m + train_sd), width=.1, position=pd, col = 'blue', alpha = 0.3) +
    geom_line(aes(y = train_m, col = gmrq), position=pd) +
    geom_point(aes(y = train_m), position=pd, size=1, col = 'blue') +
    geom_errorbar(aes(ymin=test_m - test_sd, ymax=test_m + test_sd), width=.1, position=pd, col = 'red', alpha = 0.3) +
    geom_line(aes(y = test_m, col = gmrq), position=pd) +
    geom_point(aes(y = test_m), position=pd, size=1, col = 'red') + 
    xlab('State') + ylab('GMRQ') + scale_color_manual(label = gmrq, values = RColorBrewer::brewer.pal(length(gmrq) + 1, 'Greens')[2:(length(gmrq)+1)])
  return(gg)
}

.pi_profile <- function(traj, preproc = NULL, silent = FALSE){
  
  stopifnot(!is.null(preproc))
  
  if(!is.null(preproc$FMCSC_CPROGINDSTART)) {
    if(!silent) warning('preproc$FMCSC_CPROGINDSTART coerced to 1 for file management.')
  }
  preproc[['FMCSC_CPROGINDSTART']] <- 1
  do.call(run_campari, c('trj'=list(traj), 
                         'print_status'=FALSE,
                         'multi_threading'=TRUE,
                         'run_in_background'=FALSE,
                         'return_log'=TRUE,
                         'silent' = TRUE,
                         preproc))
}


.create.trans <- function(seq.st, n.st, n.eigen, lag=1, tm.opt="symm", silent = FALSE) {
  ## Count Matrices each one for a different lag time
  cnt <- matrix(0, nrow=n.st, ncol=n.st)
  ## browser()
  for (i in 1:(.lt(seq.st)-lag)) {
    row <- seq.st[i]
    col <- seq.st[i+lag]
    cnt[row,col] <- cnt[row,col]+1
  }
  cnt.v <- rowSums(cnt)
  cnt.h <- colSums(cnt)
  if(any(cnt.h == 0) || any(cnt.v == 0)){
    if(!silent) warning('In the test or train set there is one cluster which is never visited. For this reason symm mode is coerced')
    tm.opt <- "symm"
  }
  
  if(!all(is.finite(cnt))) stop("ERROR in Window Count Matrices")
  ## Naive Transition Matrix (Simple)
  if(tm.opt=="symm") {
    tm <- t(apply(.symm(cnt), 1, function(x) {
      xx <- sum(x)
      if(xx == 0) xx <- 1
      return(x/sum(xx))
      }))
    ## ADDITION
    if(any(is.na(t))) tm[which(is.na(tm), arr.ind=TRUE)] <- 
        if(!all(is.finite(tm))) stop("ERROR in simple TM")
  } else if (tm.opt=="mle") {
    xcnt <- 2*.symm(cnt)
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
    if(any(is.na(tm)) | any(tm<0)) stop("ERROR in MLE TM")
    ## if(!all(is.finite(tm))) stop("ERROR in MLE TM")
  }
  right <- eigen(tm)     ## Right Eigenvectors
  left <- eigen(t(tm))   ## Left eigenvectors
  if (!(max(Re(left$values))<1.01 & max(Re(left$values))>0.99)) stop("ERROR in MLE Eigenvalues")
  eig.val <- Re(sort(left$values, decreasing=TRUE))
  stat <- Re(left$vectors[,which.max(Re(left$values))])/Re(sum(left$vectors[,which.max(Re(left$values))]))
  if (sum(stat)<0) stat <- -stat
  R <- Re(right$vectors[,order(Re(right$values), decreasing=TRUE)])
  A.R <- apply(R, 2, function(x) sqrt(sum(x^2*stat)))
  R.tmp <- apply(R, 1, function(x) x/A.R)
  if(identical(dim(R.tmp), dim(R))) R <- t(R.tmp)
  else R <- t(R.tmp)
  return(list(tm=tm, stat=stat, right.vec=R, eig=eig.val[c(1:max(n.eigen))], bella=R.tmp))
}

# .extract.centers <- function(log, nst) {
#   ## print.cluster.summary(log[1])
#   ## Removing the output
#   cents <- print.cluster.details(log, nst)[-c(1,2)]
#   cents <- trimws(rev(rev(cents)[-c(1:7)]))
#   cents <- sapply(strsplit(trimws(cents), "  ", fixed=TRUE), function(x) {
#     x <- x[which(x != "")]
#     as.integer(x)[3]
#   })
#   return(cents)
# }

.predict.basin <- function(test, train, cents, train.seq, silent, force_more1_state = FALSE) {
  distmat <- fields::rdist(test, train[cents,]) # distance between test and centers
  cl <- apply(distmat, 1, which.min)            # which is the minimum between the two senters
  test.seq <- train.seq[cents[cl]]              # return the cluster value for that cluster
  if(force_more1_state){
    uts <- unique(test.seq)
    if(length(uts) < 2) {
      if(!silent) warning('Insufficient variety in the test set. Every snap was in a single basin. reduce the number of kfold or avoid chunking.')
      test.seq[1] <- train.seq[cents[cents != cents[uts]]][1]
    }
  }
  return(test.seq)
}