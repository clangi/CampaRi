#' @title Writing a pdb file from trajectory matrix in R
#' @description
#'      \code{write.pdb.d} writes to file a pdb file with the trajectory matrix in input.
#'      The number of variables (number of columns) must be divisible by 3 to mantain the atomic 3D structure.
#'
#' @param trj Time series in a matrix shape (also data.frame numeric). The number of variables (nrow) must be divisible by 3 for atom-like pdb writing.
#' @param method default method will use bio3d support to write pdb file (other methods are "automatic") but it is buggy
#' @param file_name File name of the output pdb.
#' @param round \code{TRUE} will truncate the output following the \code{digit} variable.
#' @param digit Number of digits that will be kept from truncation (\code{round}).
#' @param dim_check \code{TRUE} to check if number of columns is higher than number of rows and stop consequently.
#' @param ... Extra variables for bio3d function write.pdb.
#' 
#' @details For details, please refer to the main documentation of the original campari software \url{http://campari.sourceforge.net/documentation.html}.
#' @seealso
#' \code{\link{campari}}
#' 
#' @importFrom bio3d write.pdb
#' @export write.pdb.d
#' 
#' 

write.pdb.d<-function(trj, file_name, method = "bio3d", round=FALSE, digit=4, dim_check = TRUE, ...){
  # input checks
  if(! method %in% c('bio3d', 'automatic_safe', 'automatic_unsafe'))
    stop('not supported methods')
  
  if(!is.character(file_name))
    stop('file_name must be a character.')
  
  if(!is.matrix(trj)){
    if(!is.data.frame(trj)) stop('trj input must be a matrix or a data.frame')
    trj <- as.matrix(trj)
    if(!is.numeric(trj)) stop("Numeric input must be supplied.")
  }
  
  # main functions
  if(method == 'bio3d'){
    message('Selected bio3d handling of PBD. Please check the documentation of bio3d::write.pdb and add all the necessary ')
    bio3d::write.pdb(file = file_name, xyz = trj, ...) # todo
  }else if(method == 'automatic_safe'){
    pdbfile <- file(file_name, open="w")
    
    sst <- sprintf("MODEL %8i\nATOM      1  CL  CL-     1    %8.3f%8.3f%8.3f\nATOM      2  CL  CL-     2    %8.3f%8.3f%8.3f\nATOM      3  CL  CL-     3    %8.3f%8.3f%8.3f\nATOM      4  CL  CL-     4    %8.3f%8.3f%8.3f\nATOM      5  CL  CL-     5    %8.3f%8.3f%8.3f\nENDMDL",
                   1:dd[1], as.double(wd2[, 1]), as.double(wd2[, 2]), as.double(wd2[, 3]), as.double(wd2[, 4]), as.double(wd2[, 5]), 
                   as.double(wd2[, 6]), as.double(wd2[, 7]), as.double(wd2[, 8]), as.double(wd2[, 9]), as.double(wd2[, 10]),
                   as.double(wd2[, 11]), as.double(wd2[, 12]), as.double(wd2[, 13]), as.double(wd2[, 14]), as.double(wd2[, 15]))
    
    write(sst, file=pdbfile)
    close(pdbfile)
    
  }else if(method == 'automatic_unsafe'){
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
    if(any(trj)< 0) stop("Negative values in trj are not supported yet.")
    
      
    space <- c()
    for(n in seq(nc/3)){
      cat("CL-", sep = "\n", file = paste0(base_name,".in"), append = TRUE)
    }
    cat("END", sep = "\n", file = paste0(base_name,".in"), append = TRUE)
    for (i in seq(nr)){
      cat(paste0("MODEL        ",i), file = paste0(base_name,".pdb"), sep = "\n", append=TRUE)
      if(i == (nr/100)) message("1% done...")
      if(((i*100)/nr)%%10 == 0) message(((i*100)/nr),"% done...")
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
}
