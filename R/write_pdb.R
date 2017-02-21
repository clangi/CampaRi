#' @title Writing a pdb file from trajectory matrix in R
#' @description
#'      \code{write.pdb.d} writes to file a pdb file with the trajectory matrix in input. 
#'      The number of variables (number of columns) must be divisible by 3 to mantain the atomic 3D structure.
#'
#' @param trj Time series in a matrix shape (also data.frame numeric). The number of variables (nrow) must be divisible by 3 for atom-like pdb writing.
#' @param base_name File base name.
#' @param round \code{TRUE} will truncate the output following the \code{digit} variable.
#' @param digit Number of digits that will be kept from truncation (\code{round}).
#' @details For details, please refer to the main documentation of the original campari software \url{http://campari.sourceforge.net/documentation.html}.
#' @seealso 
#' \code{\link{campari}}
#' @examples 
#' write.pdb.d(matrix(rnorm(900), nrow = 100, ncol = 9), base_name = "rnorm_trj")
#' @export 
write.pdb.d<-function(trj, base_name, round=FALSE, digit=4, dim_check = TRUE){
  if(!is.matrix(trj)){
    if(!is.data.frame(trj)) stop('trj input must be a matrix or a data.frame')
    trj <- as.matrix(trj)
    if(!is.numeric(trj)) stop("Numeric input must be supplied.")
  }
  if(file.exists(paste0(base_name, ".pdb"))) system(paste0("rm ", paste0(base_name, ".pdb")))
  if(file.exists(paste0(base_name, ".in"))) system(paste0("rm ", base_name, ".in"))
  nr <- nrow(trj)
  nc <- ncol(trj)
  # checking if the input has not been inverted in dimensions
  if(dim_check){
    if(nc > nr) stop("The number of columns (variable dimensions) are more than the number of snapshots (rows). 
                     To continue set the dim_check flag off (FALSE).")
  }
  if(nc%%3 != 0) stop("The number of row in input must be divisible for 3 for atom-like representation.")
  if(nchar(nc) > 4) stop("Unfeasible dimensions of the input matrix. It would get collapsed on other character space in the pdb.")
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
      space[2] <- paste0(rep(" ", 1 - (if(trj[i,j]<0) 1 else 0)), collapse = '')
      space[3] <- paste0(rep(" ", 1 - (if(trj[i,j+1]<0) 1 else 0)), collapse = '')
      space[4] <- paste0(rep(" ", 1 - (if(trj[i,j+2]<0) 1 else 0)), collapse = '')
      if(round) cat(paste0("ATOM ",space[1], (j+2)/3,"  ","CL  CL-",space[1], (j+2)/3, "      ", space[2],
                           format(round(trj[i,j],digits = digit), nsmall = digit),"  ", space[3],
                           format(round(trj[i,j+1],digits = digit), nsmall = digit),"  ", space[4],
                           format(round(trj[i,j+2],digits = digit), nsmall = digit)), file=paste0(base_name,".pdb"), sep='\n', append=TRUE)
      else cat(paste0("ATOM ",space[1], (j+2)/3,"  ","CL  CL-",space[1],(j+2)/3, "      ", space[2],
                      format(trj[i,j],nsmall=3),"  ", space[3],
                      format(trj[i,j+1],nsmall=3),"  ", space[4],
                      format(trj[i,j+2],nsmall=3)),file=paste0(base_name,".pdb"),sep='\n',append=TRUE)
    }
    cat("ENDMDL",file = paste0(base_name,".pdb"),sep = "\n",append=TRUE)
  }
}