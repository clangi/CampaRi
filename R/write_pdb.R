#' @title Writing a pdb file from trajectory matrix in R
#' @description
#'      \code{write.pdb.d} writes to file a pdb file with the trajectory matrix in input.
#'      The number of variables (number of columns) must be divisible by 3 to mantain the atomic 3D structure.
#'
#' @param trj Time series in a matrix shape (also data.frame numeric). The number of variables (nrow) must be divisible by 3 for atom-like pdb writing.
#' @param method default method will use automatic_safe support to write pdb file (other methods are 'bio3d', 'automatic_unsafe'). With bio3d extra options will be passed to write.pdb function.
#' The first one uses an old system with printf. The other one is fragile for large dimensions and negative values.
#' @param file_name File name of the output pdb (e.g. BPTI.pdb).
#' @param seq_name Sequence name of the output .in file (principally used by CAMPARI).
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
#' @export write_pdb
#' 
#' 

write_pdb<-function(trj, file_name, seq_name=NULL, method = "automatic_safe", round=FALSE, digit=4, dim_check = TRUE, ...){
  # input checks
  if(! method %in% c('bio3d', 'automatic_safe', 'automatic_unsafe'))
    stop('not supported methods')
  
  if(!is.character(file_name))
    stop('file_name must be a character.')
  if(!is.null(seq_name) && !is.character(seq_name))
    stop('seq_name must be a character.')
  
  if(is.null(seq_name))
    seq_name <- "seq.in"
    
  if(!is.matrix(trj)){
    if(!is.data.frame(trj)) stop('trj input must be a matrix or a data.frame')
    trj <- as.matrix(trj)
    if(!is.numeric(trj)) stop("Numeric input must be supplied.")
  }
  if(ncol(trj)%%3 != 0)
    stop('The trj matrix must be a multiple of 3 to write a pdb file correctly formatted (even if it is not molecular data).')
  
  # main functions
  if(method == 'bio3d'){
    message('Selected bio3d handling of PBD. Please check the documentation of bio3d::write.pdb and add all the necessary variables to this function (or use it directly')
    for(i in 1:nrow(trj))
      bio3d::write.pdb(file = file_name, xyz = trj[i,], append = TRUE, ...) # todo
  }else if(method == 'automatic_safe'){
    message('Selected sprintf support! it will be printed to file in append using a general scheme with CL CL- types. For more specifications please use the bio3d option.')
    pdbfile <- file(file_name, open="w")
    for(u in 1:nrow(trj)){
      sst <- sprintf("MODEL %8i\n", as.integer(u))
      sst2 <- NULL
      for(u2 in 1:(ncol(trj)/3))
        sst2 <- paste0(sst2, sprintf("ATOM  %5i  CL  CL- %5i    %8.3f%8.3f%8.3f\n", as.integer(u2), as.integer(u2), trj[u, u2], trj[u, u2 + 1], trj[u, u2 + 2]))
      write(paste0(sst, sst2, "ENDMDL"), file=pdbfile, append = TRUE)
    }
    close(pdbfile)
    
  }else if(method == 'automatic_unsafe'){
    if(file.exists(file_name)) system(paste0("rm ", file_name))
    if(file.exists(seq_name)) system(paste0("rm ", seq_name))
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
      cat("CL-", sep = "\n", file = seq_name, append = TRUE)
    }
    cat("END", sep = "\n", file = seq_name, append = TRUE)
    for (i in seq(nr)){
      cat(paste0("MODEL        ",i), file = file_name, sep = "\n", append=TRUE)
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
                             format(round(trj[i,j+2],digits = digit), nsmall = digit)), file=file_name, sep='\n', append=TRUE)
        else cat(paste0("ATOM ",space[1], (j+2)/3,"  ","CL  CL-",space[1],(j+2)/3, "      ", space[2],
                        format(trj[i,j],nsmall=3),"  ", space[3],
                        format(trj[i,j+1],nsmall=3),"  ", space[4],
                        format(trj[i,j+2],nsmall=3)), file=file_name, sep='\n',append=TRUE)
      }
      cat("ENDMDL",file = file_name,sep = "\n",append=TRUE)
    }
  }
}
