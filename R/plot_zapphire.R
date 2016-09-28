#' @title samsahara
#' @description
#'      plotting from file(X)
#'
#' @param fl file location
#' @param n_fold Number of links contractions (folds)
#' Default: \code{0}
#'
#' @examples
#' none
#'
#' @return tree: degree list, connectivity matrix and weights
#'
#' @export 
#' @import ggplot2
zap_ggplot<-function(sap_file, write=F, folderPlot = "plots/", 
                       timeline=T, basin_call=F, 
                       ann_trace = NULL, ann_trace_ret = F,
                       title = "no title"){
  # check on title
  if(!is.character(title)) stop("title var must be a string")
  # check on output folder
  if(file.exists(folderPlot)&write) print(paste0(folderPlot," already exixts. Posting plots there."))
  else if(write){
    dir.create(folderPlot)
    cat(paste0(folderPlot," created in order to contain my plots."))
  }
  # loading data
  pin <- read.table(sap_file, skip=0)
  dp <- dim(pin)
  ann_tr <- array("NA",dim = dp[1])
  
  if(is.null(ann_trace)){
    message("Annotation trace not selected. It will be considered bepartite along the timeline.")
    cat("Half random mode selected for the trace annotation.")
    ann_tr[pin[,3]>=dp[1]/2 & ann_tr == "NA"]<-"olivedrab3"
    ann_tr[pin[,3] < dp[1]/2 & ann_tr == "NA"] <- "royalblue3"
  }else if(is.numeric(ann_trace)&&length(ann_trace)==1){ 
    message("Only 5 shades of grey are possible for the 'number' option of ann_trace.
            If you inserted more than 6 it will be truncated. Please consider manual color insertion.")
    if(ann_trace>5) ann_trace = 5
    ann_tr[pin[,3]<dp[1]/ann_trace] <- "gray1"
    for(i in 1:(ann_trace-1)) ann_tr[pin[,3]<dp[1]*(i+1)/ann_trace 
                                 & pin[,3]>=dp[1]*(i)/ann_trace 
                                 & ann_tr == "NA"] <- paste0("gray",floor(100/ann_trace)*i)
  }else if(is.character(ann_trace)&&length(ann_trace)==dp[1]){
      ann_tr <- array(ann_trace)
  }else{
      stop("check the input of ann_trace or read the documentation. It is neither a number nor a color array")
    }
  
  # Set range of x and y values for the plot:
  Nsnap<-dp[1]
  xx = seq(from=1, by=1, to=Nsnap)
  ymin = 0
  ymax = -log(pin[,4]/Nsnap)
  ymax = ymax[!is.infinite(ymax)&!is.na(ymax)]
  ymax = max(ymax)
  yr = range(c(ymin, ymax))
  xr = range(c(1, Nsnap))
  
  # initial creation of the plot
  gg <- ggplot(data = pin, mapping = aes(x = xx, y = -log((pin[,4]/Nsnap)))) +
    xlab("Progress Index") + ylab("Annotation") + ggtitle(title)
  
  # plotting the trace and the timeline
  gg <- gg + geom_segment(aes(xx, y = rep(ymax*3/4,length(xx)),
                              xend = xx, yend = rep(ymax*3/4+ymax/8, length(xx))), 
                          col = ann_tr)
  if(timeline)gg <- gg + geom_point(aes(xx,y=(pin[,3]*1.0*ymax*1/5)/dp[1]-1/10),col=ann_tr,size=0.01)
  gg <- gg + geom_line(color="black") +
    geom_point(mapping = aes(x=xx,y=2.5 - (1./3.)*log((pin[,10] + pin[,12]) / Nsnap)), color="red", size=0.01)
  if(basin_call) gg <- gg + 
    geom_text(data = data.frame(), aes(Nsnap/4, ymax-1*ymax/14, label = "Basin 1")) +
    geom_text(data = data.frame(), aes(Nsnap*3/4, ymax-1*ymax/14, label = "Basin 2"))
  # p + annotate("rect", xmin = 3, xmax = 4.2, ymin = 12, ymax = 21,
  #              alpha = .2)
  plot(gg)
  if(ann_trace_ret)return(ann_tr)
}
