#' @title writing the pdb file
#' @description
#'      pwriting the pdb file(X)
#'
#' @param fl file location
#' @param n_fold Number of links contractions (folds)
#' Default: \code{0}
#'
#'
#' @return tree: degree list, connectivity matrix and weights
#'
#' @export 
write.pdb.d<-function(x,base_name,round=FALSE,digit=4){
  if(!is.numeric(x)) stop("Numeric input must be supplied.")
  
  if(file.exists(paste0(base_name,".pdb"))) system(paste0("rm ",paste0(base_name,".pdb")))
  if(file.exists(paste0(base_name,".in"))) system(paste0("rm ",base_name,".in"))
  if(is.null(nrow(x))&&is.atomic(x)) {
    nr <- 1
    nc <- length(x)
    x <- matrix(x,nrow = nr, ncol = nc)
  } else {
    nr <- nrow(x)
    nc <- ncol(x)
  }
  if(nc%%3 != 0) stop("The number of row in input must be divisible for 3 for atom-like representation.")
  if(nchar(nr) > 4) stop("Unfeasible dimensions of the input matrix. It would get collapsed on other character space in the pdb.")
  space <- c()
  for(n in seq(nc/3)){
    cat("CL-", sep = "\n", file = paste0(base_name,".in"), append = TRUE)
  }
  cat("END", sep = "\n", file = paste0(base_name,".in"), append = TRUE)
  for (i in seq(nr)){
    cat(paste0("MODEL        ",i), file = paste0(base_name,".pdb"), sep = "\n", append=TRUE)
    for(j in seq(to = (nc),by = 3)){
#       if(!(((j-1)/3)%%1==0)) next
#       if(((j+2)/3)<10) space<-" " 
#       else space<-""
      space[1] <- paste0(rep(" ", 6 - nchar((j+2)/3)), collapse = '')
      space[2] <- paste0(rep(" ", 1 - (if(x[i,j]<0) 1 else 0)), collapse = '')
      space[3] <- paste0(rep(" ", 1 - (if(x[i,j+1]<0) 1 else 0)), collapse = '')
      space[4] <- paste0(rep(" ", 1 - (if(x[i,j+2]<0) 1 else 0)), collapse = '')
      if(round) cat(paste0("ATOM ",space[1], (j+2)/3,"  ","CL  CL-",space[1], (j+2)/3, "      ", space[2],
                           format(round(x[i,j],digits = digit), nsmall = digit),"  ", space[3],
                           format(round(x[i,j+1],digits = digit), nsmall = digit),"  ", space[4],
                           format(round(x[i,j+2],digits = digit), nsmall = digit)), file=paste0(base_name,".pdb"), sep='\n', append=TRUE)
      else cat(paste0("ATOM ",space[1], (j+2)/3,"  ","CL  CL-",space[1],(j+2)/3, "      ", space[2],
                      format(x[i,j],nsmall=3),"  ", space[3],
                      format(x[i,j+1],nsmall=3),"  ", space[4],
                      format(x[i,j+2],nsmall=3)),file=paste0(base_name,".pdb"),sep='\n',append=TRUE)
    }
    cat("ENDMDL",file = paste0(base_name,".pdb"),sep = "\n",append=TRUE)
  }
}