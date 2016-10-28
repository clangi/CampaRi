#' @title plotting
#' @description
#'      plotting from file(X)
#'
#' @param fl file location
#' @param n_fold Number of links contractions (folds)
#' Default: \code{0}
#'
#'
#' @return tree: degree list, connectivity matrix and weights
#'
#' @export 
#' @import ggplot2
zap_ggplot<-function(sap_file=NULL, sap_table=NULL,write=F, folderPlot = "plots/", 
                     timeline=T, subtitle = NULL, basin_call=F, local_cut=T,
                     ann_trace = F, ann_trace_ret = F, background_height = NULL, 
                     ann_names_L = NULL,ann_names_R = NULL,
                     title = "no title"){
  # check on title
  if(!is.character(title)) stop("title var must be a string")
  # check on output folder
  if(file.exists(folderPlot)&write) print(paste0(folderPlot," already exixts. Posting plots there."))
  else if(write){
    dir.create(folderPlot)
    cat(paste0(folderPlot," created in order to contain my plots."))
  }
  # loading data - sapphire table
  if(is.null(sap_table)&&!is.null(sap_file)) pin <- read.table(sap_file)
  else if(is.null(sap_file)&&!is.null(sap_table)) pin <- sap_table
  else stop("Sapphire table needed in input. Check the documentation")
  dp <- dim(pin)
  ann_tr <- array("NA",dim = dp[1])
  if(!is.logical(ann_trace)&&is.numeric(ann_trace)&&length(ann_trace)!=1)
    nrow_an_tr <- nrow(ann_trace)
  else if(length(ann_trace)!=dp[1])
    nrow_an_tr <- 1
  else
    nrow_an_tr <- NULL
  
  #checking the trace input
  if(is.null(nrow_an_tr)&&length(ann_trace)!=dp[1]&&length(ann_trace)!=1&&!is.null(ann_trace)&&!is.logical(ann_trace)) 
    stop("The annotation trace must be eighter a number vector (or matrix) with length = input trj eighter
         a single value (1-10,T/F). ") 

  if(!is.null(ann_trace)&&!is.logical(ann_trace)&&(any(!sapply(ann_trace,is.numeric)) ||
     max(ann_trace)>10)) stop("For manual insertion of the trace use numbers 1-10 for each value (also more than one row)")
  if(!is.null(nrow_an_tr)) max_an_tr <- max(ann_trace)
  
  
#Main ann_trace constructor
  if(is.null(ann_trace)||
     (is.logical(ann_trace)&&ann_trace)||
     (is.numeric(ann_trace)&&length(ann_trace)==1&&ann_trace==2)){
    message("Annotation trace not selected. It will be considered bepartite along the timeline.")
    cat("Half random mode selected for the trace annotation. First half will be light grey")
    ann_tr[pin[,3]>=dp[1]/2 & ann_tr == "NA"]<-"gray75"
    ann_tr[pin[,3] < dp[1]/2 & ann_tr == "NA"] <- "gray30"
    nrow_an_tr <- 1
  }else if(is.numeric(ann_trace)&&length(ann_trace)==1){ 
    message("Only 10 shades of grey are possible for the 'number' option of ann_trace.
            If you inserted more than 10 it will be truncated. Please consider manual color insertion.")
    if(ann_trace>10) ann_trace = 10
    ann_tr[pin[,3]<dp[1]/ann_trace] <- "gray1"
    for(i in 1:(ann_trace-1)) ann_tr[pin[,3]<dp[1]*(i+1)/ann_trace 
                                 & pin[,3]>=dp[1]*(i)/ann_trace 
                                 & ann_tr == "NA"] <- paste0("gray",floor(100/ann_trace)*i)
    nrow_an_tr <- 1
  }else if(nrow_an_tr==1){
    ann_tr <- sapply(ann_trace,FUN = function(x){
      paste0("gray",floor(100/max_an_tr)*x)
    })
  }else if(nrow_an_tr>1){
    ann_tr <- array("NA", dim = dim(ann_trace))
    for(i in 1:nrow_an_tr){
      ann_tr[i,] <- sapply(ann_trace[i,],FUN = function(x){
        paste0("gray",floor(100/max_an_tr)*x)
      })
    }
  }else if(!ann_trace){
    warning("ann_trace = F silenced the annotation trace.")
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
  
  # initial creation of the plot
  gg <- ggplot(data = pin, mapping = aes(x = xx, y = -log((pin[,4]/Nsnap)))) +
    xlab("Progress Index") + ylab("Annotation") 
    # theme_bw() +
    # theme(panel.grid.minor = element_line(colour="gray80"))
  
  #Trace height from the top. This is the 0-16 parts out of ymax
  if(!is.null(background_height)&&is.numeric(background_height)&&length(background_height)==1){
    if(background_height>14||background_height<0){
      warning("Inserted background height too small or too big.")
      background_height <- 12
    }
    tr_init <- background_height
    if(background_height>8)
      main_col<-"darkblue"
    else
      main_col <- "dodgerblue"  
  }else if(!is.null(background_height)&&is.character(background_height)&&background_height=="full"){
    if(timeline)
      tr_init <- 4
    else
      tr_init <- 0
    main_col <- "dodgerblue"
  }else{
    tr_init <- 12
    main_col<-"darkblue"
  }
  # plotting the trace and the timeline
  if(!is.logical(ann_trace)&&nrow_an_tr==1){
    gg <- gg + geom_segment(aes(xx, y = rep(ymax*3/4,length(xx)),
                                xend = xx, yend = rep(ymax*3/4+ymax/8, length(xx))),
                            col = ann_tr)
  } else if(!is.logical(ann_trace)){
    for(i in 0:(nrow_an_tr-1))
      gg <- gg + geom_segment(x=xx, y = rep(ymax*((tr_init+((i*(16-tr_init))/nrow_an_tr))/16),length(xx)),
                                  xend = xx, yend = rep(ymax*((tr_init+(((i+1)*(16-tr_init))/nrow_an_tr))/16), length(xx)), 
                              col = ann_tr[(i+1),])
  }
  # annotation names LEFT
  if(!is.null(ann_names_L)&&is.character(ann_names_L)&&length(ann_names_L)==nrow_an_tr){
    for(i in 0:(nrow_an_tr-1))
      gg <- gg + annotate("text",x = -(length(xx)/24)*nchar(ann_names_L[i+1])/2, 
                          y = ymax*((tr_init + (((i+0.5)*(16-tr_init))/nrow_an_tr))/16), label = ann_names_L[i+1])
  }else if(!is.null(ann_names_L)){
    stop('The annotation names have not been inserted correctly')
  }
  # annotation names RIGHT
  if(!is.null(ann_names_R)&&is.character(ann_names_R)&&length(ann_names_R)==nrow_an_tr){
    for(i in 0:(nrow_an_tr-1))
      gg <- gg + annotate("text",x = length(xx)+(length(xx)/24)*nchar(ann_names_R[i+1])/2, 
                          y = ymax*((tr_init + (((i+0.5)*(16-tr_init))/nrow_an_tr))/16), label = ann_names_R[i+1])
  }else if(!is.null(ann_names_R)){
    stop('The annotation names have not been inserted correctly')
  }
  #basic annotation
  gg <- gg + geom_line(color=main_col,size=0.2)
  
  #timeline at the bottom
  if(timeline&&!is.logical(ann_trace)&&nrow_an_tr==1) {
    gg <- gg + geom_point(aes(x=xx,y=(pin[,3]*1.0*ymax*1/5)/dp[1]-1/10),col=ann_tr,size=0.01)
    }else if(timeline&&!is.logical(ann_trace)){
      gg <- gg + geom_point(aes(x=xx,y=(pin[,3]*1.0*ymax*1/5)/dp[1]-1/10),col=rep("black",length(xx)),size=0.01)
    }
  
  #local cut
  if(local_cut) gg <- gg + geom_point(mapping = aes(x=xx,y=2.5 - (1./3.)*log((pin[,10] + pin[,12]) / Nsnap)), 
                                      color="red3", size=0.1)
  
  #basin call
  if(basin_call) gg <- gg + 
    geom_text(data = data.frame(), aes(Nsnap/4, ymax-1*ymax/14, label = "Basin 1")) +
    geom_text(data = data.frame(), aes(Nsnap*3/4, ymax-1*ymax/14, label = "Basin 2"))
  # p + annotate("rect", xmin = 3, xmax = 4.2, ymin = 12, ymax = 21,
  #              alpha = .2)
  if(!is.null(subtitle)&&is.character(subtitle))
    gg <- gg + ggtitle(bquote(atop(.(title), atop(italic(.(subtitle)), "")))) 
  else
    gg <- gg + ggtitle(title)
  plot(gg)
  if(ann_trace_ret) return(ann_tr)
}