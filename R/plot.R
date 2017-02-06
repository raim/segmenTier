### plotting segmentation results and segments


### PLOTTING RESULTS - from genomeBrowser 20161102

## from lib TeachingDemos; used in plotFeatures
shadowtext <- function(x, y=NULL, labels, col='white', bg='black',
                       theta= seq(pi/4, 2*pi, length.out=8), r=0.1, ... ) {
  
  xy <- xy.coords(x,y)
  xo <- r*strwidth('A')
  yo <- r*strheight('A')
  
  for (i in theta) {
    text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ... )
  }
  text(xy$x, xy$y, labels, col=col, ... )
}

## plot features as arrows; kept private since official version
## of this function is maintained in package genomeBrowser
segment.plotFeatures <- function(data, coors, types, strand,
                         typord=FALSE, cuttypes=FALSE,
                         names=FALSE, legend=FALSE, axis1=FALSE, ylab=NA,
                         args,line=1,ycx=1,tcx=1,tpos=NULL,
                         columns=c(name="name",chr="chr",strand="strand",
                           start="start",end="end",type="type",color="color")) {

  ## parse arguments, which will also override other 
  if ( !missing(args) ) {    
    for ( i in 1:length(args) )       
      assign(names(args)[i],args[[i]])
  }

  ## get all available types
  if ( missing(types) ) 
    types <- unique(data[,columns["type"]])

  ## use only features from the requested strand
  if ( !missing(strand) )
    data <- data[data[,columns["strand"]]==strand,]

  ## coordinates
  chr <- as.character(coors[1])
  xlim <- sort(coors[2:3])

  ## cut to xlim
  coor <- cbind(as.numeric(data[,columns["start"]]),
                as.numeric(data[,columns["end"]]))
    
  ## find and filter features overlapping with requested coors
  lovl <- apply(coor,1,max) >= xlim[1] & apply(coor,1,max) <= xlim[2]
  rovl <- apply(coor,1,min) >= xlim[1] & apply(coor,1,min) <= xlim[2]
  insd <- lovl & rovl
  span <- apply(coor,1,min) <= xlim[1] & apply(coor,1,max) >= xlim[2]
  ovl <- ( lovl | rovl ) | span
  ## index features to determine where to plot names
  if ( names )
    data <- cbind(data, lovl=lovl, rovl=rovl, insd=insd, span=span)

    ## filter overlapping features
    ## to avoid factor confusion!
    if ( "chr"%in%names(columns) ) {
        chrs <- as.character(data[,columns["chr"]])
        ovl <- chrs == chr & ovl
    }
    feat <- data[ovl,, drop=FALSE]

    ## get new coordinates
    coor <- cbind(as.numeric(feat[,columns["start"]]),
                  as.numeric(feat[,columns["end"]]))

    ## order start/end by strand (arrow direction!)
    if ( "strand"%in%names(columns) ) {
        str <- as.character(feat[,columns["strand"]])
        ## TODO draw different arrow for strand=.
        coor[str==".",] <-
            t(apply(coor[str==".",,drop=FALSE],1,
                    function(x) sort(x,decreasing=FALSE)))
        coor[str=="+"|str=="+1",] <-
            t(apply(coor[str=="+"|str=="+1",,drop=FALSE],1,
                    function(x) sort(x,decreasing=FALSE)))
        coor[str=="-"|str=="-1",] <-
            t(apply(coor[str=="-"|str=="-1",,drop=FALSE],1,
                    function(x) sort(x,decreasing=TRUE)))
    }

  ## TODO: redirect to plotFeatureBlocks if too many features

  ## only take types present in the current set
  ## NOTE: this was an if(missing(types))elseif before 20160212: 
  ## thus setting types overruled cuttypes: was this necessary sometimes?
  if ( cuttypes ) 
    types <- types[types%in%feat[,columns["type"]]]

  ## ylimits without typord
  typey <- NULL
  ylim <- c(-.1,1.1)
  ## generate y-order and ylim, unless y is set by types
  ## TODO: generate a matrix of feature overlaps to select y-value and ylim
  if ( typord & (nrow(coor)>0 & length(types)>0) ) {
    offset <- 1/length(types)
    typey <- 1-offset*0:(length(types)-1)
    names(typey) <- types
    ylim <- c(min(typey)-offset, max(typey)+offset)
  }
  ## plot
  #old.par <- par(no.readonly = TRUE)
  #on.exit(par(old.par))
  plot(0,xlim=xlim,ylim=ylim,col=NA,axes=FALSE,ylab=ylab,xlab=NA)
  offset <- - 1/nrow(feat)
  y <- .5
  for ( i in order(rowMeans(coor)) ) {
      
    typ <- as.character(feat[i,columns["type"]])
    if ( !typ %in% types ) next
    yt <- ifelse(typord, typey[typ], y) 
    col <- 1
    if ( "color"%in%names(columns) ) # take feature color if present
        if ( columns["color"] %in% colnames(feat) )
            col <- as.character(feat[i,columns["color"]])
    
    
    x <- coor[i,]
    if ( x[1]!=x[2] )
      arrows(x0=x[1], y0=yt, col=col,y1=yt, x1=x[2], lwd=2, length=.075)
    points(x=x[1], y=yt, col=col, pch=3, lwd=2)

      ## plot names
      if ( names ) {
          if ( feat[i,"insd"] ) x <- mean(x)
          else if ( feat[i,"span"] ) x <- mean(xlim)
          else if ( feat[i,"lovl"] ) x <- max(x)
          else if ( feat[i,"rovl"] ) x <- min(x)
          shadowtext(x=x,y=yt,labels=feat[i,columns["name"]],
                     col=col,bg="white",r=.15,cex=tcx,pos=tpos)
      }
    if ( typord ) next
    y <- y + offset
    if ( y > 0.9 & offset> 0 ) { offset <- -offset; y <- 1 }
    if ( y < 0.1 & offset< 0 ) { offset <- -offset; y <- 0 + offset/2 }
  }
  if ( typord & (nrow(coor)>0 & length(types)>0) ) 
    axis(2,at=typey,labels=names(typey),las=2,line=line,tcl=.3, mgp=c(2,.5,0),
         cex.axis=ycx)
  if ( axis1 ) axis(1,line=-1)

  ## silently return y-position of types (list of plotted features)
  features <- typey #rownames(feat[feat[,columns["type"]] %in% types,])
}

## plot 3D genomeData as a heatmap; kept private since official version
## of this function is maintained in package genomeBrowser
segment.plotHeat <- function(data, coors, orig.idx, breaks, colors, orig.col, ylab="", axes=FALSE, axis2=FALSE, colnorm=FALSE, chrS, coor, args) {
  
  ## parse arguments, which will also override other 
  if ( !missing(args) ) {    
    for ( i in 1:length(args) )       
      assign(names(args)[i],args[[i]])
  }

  chr <- coors[1]
  xlim <- sort(coors[2:3])

  ## split off coordinate columns, if not separately supplied
  ## TODO: use as alternative to chrS
  ## TODO: clean this up a bit (chrS vs. coor vs. data)
  firstcol <- 1
  if ( sum(c("chr","coor")%in%colnames(data))==2 )
    firstcol <- 3

  ## cut to xlim - TODO: consolidate with coor in plotFeatures and plotBlocks
  if ( missing(chrS) ) { ## search indices - SLOW
    idx <- which(as.numeric(data[,"chr"]) == chr &
                 (as.numeric(data[,"coor"]) >= xlim[1] &
                  as.numeric(data[,"coor"]) <= xlim[2]))
  } else { ## take direct index - ONLY FOR EXPANDED and ORDERED DATA
    start <- xlim[1]+chrS[chr]
    end <- xlim[2]+chrS[chr]
    if ( start < 1 ) start <- 1
    if ( end > nrow(data) ) end <- nrow(data) ## TODO: cut by chromosome!
    idx <- start:end
  }
  
  dat <- as.matrix(data[idx,,drop=FALSE])
  if ( !missing(orig.idx) ) orig.idx <- orig.idx[idx] 
  if ( missing(chrS) )
    x <- dat[,"coor"]
  else
    x <- xlim[1]:xlim[2] # start:end - chrS[chr]
  ##  x <- coor[idx,"coor"]
  dat <- as.matrix(dat[,firstcol:ncol(data),drop=FALSE])
  y <- 1:ncol(dat)

  if ( length(x)==0 | !any(!is.na(dat)) ) {  ## empty?
    plot(1,xlim=xlim,col=NA,axes=FALSE,xlab="chromosome position",ylab=ylab)
  }else{

    if ( colnorm ) { # between 0 and 1
      dat <- t(apply(dat,1,function(x) {
        if ( sum(is.na(x))==length(x) ) return(x)
        mx <- max(x,na.rm=TRUE)
        mn <- min(x,na.rm=TRUE)
        if ( mx==mn ) rep(0.5,length(x)) else ((x-mn)/(mx-mn))}))
    }
    
    if ( missing(breaks) ) {
      if ( !missing(colors) ) nbreaks <- length(colors)+1
      else nbreaks <- 101
      breaks <- seq(min(dat,na.rm=TRUE),max(dat,na.rm=TRUE),length.out=nbreaks)
    }
    nbreaks <- length(breaks)
    if ( missing(colors) )
      colors <- gray(seq(1,0,length.out=nbreaks-1))
    
    
    p1 <- image(x=x,y=y,z=dat, axes=FALSE, ylab=ylab, xlim=xlim,
                xlab="chromosome position",breaks=breaks,col=colors)
  }
  ## if present/requested, plot locations of original data
  ## (or abitrary discrete info on positions)
  if ( !missing(orig.idx) ) {
    if ( missing(orig.col) ) orig.col <- 2
    # cat(paste("plotting original data locations\n"))
    axis(1,at=x[which(orig.idx)],labels=NA,tcl=.2,col=NA,col.ticks=orig.col,lwd.ticks=1)
  }
  if ( axes ) {
    axis(1)
    axis(2,las=2)
  } else if ( axis2 )
    axis(2,las=2)
}

### PLOTTING SEGMENTATION INTERNAL STRUCTURES 




## plot a matrix as seen in R; kept private since official version
## of this function is maintained in package segmenTools
image_matrix <- function(dat, text, text.col, axis=1:2, axis1.col, axis2.col, ...) {

    ## reverse columns and transpose
    if ( nrow(dat)>1 )
        imgdat <- t(apply(dat, 2, rev))
    else
        imgdat <- t(dat)
    image(x=1:ncol(dat), y=1:nrow(dat), z=imgdat, axes=FALSE, ...)

    ## add text
    if ( !missing(text) ) {
        if ( missing(text.col) )
            text.col <- rep(1, length(c(text)))
        text(x=rep(1:ncol(dat),nrow(dat)), y=rep(nrow(dat):1,each=ncol(dat)),
             paste(t(text)),col=t(text.col))
    }
    
    ## add axes
    ## TODO : handle axes=FALSE
    if ( !missing(axis) ) {
        if ( 1 %in% axis ) 
            if ( !missing(axis1.col) ) # colored ticks
                for ( i in 1:ncol(dat) )
                    axis(1, at=i,colnames(dat)[i],
                         col.axis=axis1.col[i], col=axis1.col[i],
                         las=2, cex.axis=1.5, lwd=2)
            else
                axis(1, at=1:ncol(dat), labels=colnames(dat), las=2)
        if ( 2 %in% axis )
            if ( !missing(axis2.col) ) # colored ticks
                for ( i in 1:nrow(dat) )
                        axis(2, at=nrow(dat)-i+1, rownames(dat)[i],
                             col.axis=axis2.col[i], col=axis2.col[i],
                             las=2, cex.axis=1.5, lwd=2)
            else
                axis(2, at=nrow(dat):1, rownames(dat),las=2)        
    }
} 


#' plot the processed time-series object returned from
#' \code{\link{processTimeseries}}.
#' @param x the time-series object returned by
#' \code{\link{processTimeseries}}
#' @param plot a string vector indicating the values to be plotted;
#' `total': plot of the total signal, summed over
#' the time-points, and indicating the applied threshold \code{low.thresh};
#' note that the total levels may have been transformed (e.g. by \code{log_1}
#' or \code{ash}) depending on the arguments \code{trafo} and \code{dc.trafo}
#' in \code{\link{processTimeseries}}; `timeseries': plot the complete
#' time-series as a heatmap, where time is plotted bottom-up on the y-axis
#' @param ... currently unused additional arguments to plot
#'@export
plot.timeseries <- function(x, plot=c("total","timeseries"), ...) {
    
    tset <- x
    
    ## get time-series data
    ts <- tset$ts # incl. all trafos and zeros set to NA
    ts[tset$zero.vals,] <- NA
    
    tot <- tset$tot # total of the time-series

    ## mock "chromosome" coordinates
    N <- nrow(ts)
    coors <- c(chr=1,start=1,end=N) 

    ## timeseries heatmap colors
    colors0 <- rev(grDevices::gray.colors(100)) 
    colors0[1] <- "#FFFFFF" ## replace minimal by white

    ## settings
    low.thresh <- tset$settings$low.thresh
    logged <- tset$settings$trafo != "raw" # plot with log if not already done
    
    if ( "total" %in% plot ) {
        plot(1:N,tot,log=ifelse(logged,"","y"),
             type="l",lwd=2,axes=FALSE,ylab=NA,xlab=NA)
        graphics::polygon(x=c(1,1,N,N),
                          y=c(min(tot,na.rm=TRUE),rep(low.thresh,2),
                              min(tot,na.rm=TRUE)),col="#00000055",border=NA)
        graphics::abline(h=low.thresh,col="#000000BB")
        lines(1:N,tot)
        axis(2);
        axis(1)
        graphics::mtext("total signal", 2, 2)
    }
    if ( "timeseries" %in% plot ) {
        segment.plotHeat(ts,coors=coors,chrS=0,colors=colors0, colnorm=TRUE)
        axis(2,at=1:ncol(ts))
        axis(1)
        graphics::mtext("time points", 2, 2)
    }
}

#' plot the clustering object returned by \code{\link{clusterTimeseries}}
#' @param x a set of clusterings as returned by
#' \code{\link{clusterTimeseries}}
#' @param k a numeric or string vector indicating the clusterings to be plotted;
#' specifically the column numbers or names in the matrix of clusterings
#' in \code{cset$clusters}; if missing all columns will be plotted
#' and the calling code must take care of properly assigning \code{par(mfcol)}
#' or \code{layout} for the plot
#' @param xaxis optinally x-values to use as x-axis (e.g. to reflect absolute
#' chromosomal coordinates)
#' @param ... currently unused additional arguments to plot
#'@export
plot.clustering <- function(x, k, xaxis, ...) {

    cset <- x
    
    ## cluster sorting via Ccc (cluster-cluster correlation)
    if ( !"sorting" %in% names(cset) )
      cset <- sortClusters(cset, verb=1)
    ## cluster colors
    if ( !"colors" %in% names(cset) )
      cset <- colorClusters(cset)

    ## plotting all: layout or mfcol must be set from outside;
    ## or a specific k chosen to only plot the k'th clustering
    if ( missing(k) ) 
        k <- 1:ncol(cset$clusters)
    if ( length(k)>1 )
        par(mfcol=c(length(k),1))
    for ( i in k ) {
        seq <- as.character(cset$clusters[,i])
        cls.srt <- cset$sorting[[i]]
        if ( missing(xaxis) )
          xaxis <- 1:length(seq)
        y <- 1:length(cls.srt)
        names(y) <- cls.srt
        cols <- cset$colors[[i]]
        ## plot original clustering
        plot(xaxis,y[seq],axes=FALSE,xlab="",ylab=NA,
             col=cols[seq],cex=1,pch=16)
        axis(2, at=y, labels=names(y), las=2)
        graphics::mtext("cluster", 2, 2)
    }
}

#' plot the final segmentation objects returned by
#' \code{\link{segmentClusters}} and \code{\link{segmentCluster.batch}}
#' @param x a set of segmentations as returned by
#' \code{\link{segmentClusters}} and \code{\link{segmentCluster.batch}}
#' @param types a string vector indicating segment types to plot (a subset of
#' \code{sset$ids}; defaults to all in \code{sset$ids})
#' @param xaxis optional x-values to use as x-axis (e.g. to reflect absolute
#' chromosomal coordinates)
#' @param plot string list indicating which data should be plotted;
#' `segments': plot segments as arrows; `S1' plot the scoring vectors
#' \code{s(i,j,c} for all \code{c}; `S' plot the derivative of
#' matrix \code{S(i,c)} for all \code{c}
#' @param ... currently unused additional arguments to plot
#'@export
plot.segments <- function(x, types, xaxis, plot=c("segments", "S", "S1"), ...) {

    sset <- x
    
    if ( missing(types) ) 
        types <- rownames(sset$settings)

    ## one plot for all segments 
    if ( "segments" %in% plot ) {

        ## mock "chromosome" coordinates
        N <- sset$N
        coors <- c(chr=1,start=1,end=N) 

        segs <- sset$segments
        columns <- c(name="ID", type="type", start="start", end="end",
                     color="color")
        ## filter allsegs by segments for the current clustering
        ypos <- segment.plotFeatures(segs, types=types,
                                     coors=coors, typord=TRUE,cuttypes=TRUE,
                                     ylab="", names=FALSE,columns=columns,
                                     tcx=.5)
        axis(1)
        ## plot fuse tag
        fuse <- segs[segs[,"fuse"],]
        points(fuse[,"start"], ypos[fuse[,"type"]], col="black",
               pch=4, lwd=1, cex=1.5)
    }

    ## colors for S1 heatmap
    colors0 <- rev(grDevices::gray.colors(100)) 
    colors0[1] <- "#FFFFFF" ## replace minimal by white

    if ( any(c("S","S1") %in% plot) ) {

        ## get matrices, cluster sorting and colors
        SK <- sset$SK[types]
        sk.srt <- sset$sorting[types]
        sk.col <- sset$colors[types]

        ## one plot for each segmentation!
        for ( j in 1:length(SK) ) {

            ## cluster sorting
            srt <- as.character(sk.srt[[j]])
            ## cluster coloring; add black for nuissance
            sgcols <- c("#000000", sk.col[[j]][srt])
            srt <- c("0",srt)

            ## plot S1 as heatmap
            if ( "S1" %in% plot ) {
                S1 <- t(SK[[j]]$S1[,rev(srt)]) # sort by cluster sorting!
                S1 <- S1/apply(S1,1,mean)
                ## TODO: add from segmenTools
                image_matrix(S1,axis=2, col=colors0, ylab="cluster")
            }
            if ( !"S" %in% plot ) next
            
            ## add alpha
            sgcols <- paste(sgcols,"EE",sep="") 
            names(sgcols) <- srt
            ## only show clusters that actually produced a segment
            tp <- sset$segments[,"type"]%in%names(SK)[j]
            sgcols[!names(sgcols)%in%as.character(sset$segments[tp,"CL"])] <- NA

            ## get matrix and sort according to cluster sorting
            ## TODO: plot by segment; highlight winning segment!!
            S <- SK[[j]]$S[,srt]
            dS <- apply(S,2,function(x) c(0,diff(x)))

            ## x-axis: x can be passed to use real coordinates
            if ( missing(xaxis) )
              xaxis <- 1:nrow(S)
            xlim <- range(xaxis)             

            xrng <- stats::quantile(xaxis,c(.05,.95))
            xidx <- which(xaxis>xrng[1]&xaxis<xrng[2]) #x%in%xrng[1]:xrng[2]
            ylim <- stats::quantile(ash(dS[xidx,]),c(0,1))
            plot(1,ylim=ylim,xlim=xlim,ylab=expression(ash(Delta~S["i,C"])))
            lines(xaxis,ash(dS[,1]),lwd=7,col="#00000015") # NUI: BACKGROUND 
            lines(xaxis,ash(dS[,1]),lwd=1,lty=3,col="#00000099") 
            graphics::matplot(xaxis, ash(dS), type="l",
                              lty=1, lwd=1, add=TRUE, col=sgcols)
            graphics::mtext(names(SK)[j], side=2 , line=4, las=2)
        }
    }
}

#' plot all objects from the segmentation pipeline, i.e. the processed
#' time-series, the clustering, the internal scoring matrices and
#' the final segments
#' @param tset the time-series object returned by
#' \code{\link{processTimeseries}}
#' @param cset a set of clusterings as returned by
#' \code{\link{clusterTimeseries}}
#' @param sset a set of segmentations as returned by
#' @param plot.matrix include the internal scoring matrices in the plot
#' @param mai margins of invidual plots, see \code{par}
#' \code{\link{segmentClusters}} and \code{\link{segmentCluster.batch}}
#'@export
plotSegmentation <- function(tset, cset, sset, plot.matrix=FALSE,
                             mai=c(.01,1.5,.01,.01)) {

    nsg <- length(sset$ids)# total number of segmentations
    nk <- length(cset$ids) # number of clusterings
    spk <- nsg/nk # segmentations per clustering
    ## number of plots
    ## each clustering can have multiple segmentations; plot each
    ## 2 for time-series; and for each clustering 2 (clustering and segments),
    ## plus S1 & S 
    nplots <- 2 + nk * (2 + ifelse(plot.matrix, 2*spk, 0))
    par(mfcol=c(nplots,1), xaxs="i", mai=mai)

    ## TIME-SERIES PLOT UTILITY: plot both the total signal (optionally used
    ## for threshold) and a heatmap of the time-series
    plot(tset, plot=c("total","timeseries"))
    
    ## CLUSTERING PLOT UTILITY: 
    ## NOTE that clusterings are sorted (by their similarity matrix `Ccc`)
    ## and colored along a color-wheel
    for ( k in 1:ncol(cset$clusters) ) {
        ## plot clustering
        plot(cset, k)
        ## SEGMENTATION PLOT UTILITY
        ## plot all segments, S1(i,c), S(i,c) for this clustering
        kid <- cset$ids[k]
        types <- rownames(sset$settings)[sset$settings[,"K"] %in% kid]
        plot <- "segments"
        if ( plot.matrix ) plot <- c("segments","S","S1")
        plot(sset, plot=plot, types=types) 
    }
}
