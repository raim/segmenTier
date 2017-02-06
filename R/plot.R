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

#' plot features as arrows
#'@export
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

#' plot 3D genomeData as a heatmap
## NOTE: fast indexing via chrS (cumulative chrosome lengths)
## is possible only for fully expanded 'data')
#'@export
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



#' plot segmentation data based on data returned by
#' high-level wrappers, incl. S(c,i), K(c,i) and segments
#' @details This is mostly used for testing, where a small data set
#' is segmented by a high-level wrapper for multiple scoring functions
#' and parameters.
#' @param scrR a list containing multiple segmentations for (1) different
#' scoring functions, (2) different parameters of the forward-step (dyn.prog.)
#' and (3) different parameters of the back-tracing step.
#' @param seq the sequence of clusters to be segments
#' @param ts an optional time-series from for which the clustering was
#' calculated
#' @param tot an optional total signal from the timeseries in \code{ts}
#' @param out.file optional out.file name (w/o file extension) to which a
#' png will be plotted
#' @param use.log plot the total data (\code{tot}) with logged y-axis
#' @param add.plots in the plot layout, leave these rows for additional
#' external plots
#' @param verb level of verbosity; 0: no output, 1: progress messages
#' @return Returns the file.name if out.file was specified.
#' @export
plotSegments <- function(cset, ts, scrR,tot, out.file, use.log=FALSE,
                         add.plots=0, verb=2) {
  
    if ( verb>0 )
        cat(paste("plotting",
                  ifelse(missing(out.file),"results",out.file),"\t",date(),"\n"))

    ## if nuissance cluster is present, increase
    ## coloring by 1
    seq <- cset$clusters[,1]
    if ( 0 %in% seq ) {
        cols <- sort(unique(seq))
        names(cols) <- cols
        cols <- cols+1
    }
    
    ## how many rows?
    ## time-series and clustering
    rws <- as.numeric(!missing(tot)) + as.numeric(!missing(ts)) + 1 + add.plots
    scores <- names(scrR)
    scores <- scores[scores!="SM"]
    for ( scor in scores ) {
        multis <- names(scrR[[scor]])
        multis <- multis[multis!="SM"]
        rws <- rws + 1 + 2*(length(multis)) # one for each S(c,i) 
    }
    x <- 1:length(seq)
    xlim <- c(min(x)-1,max(x)+1)
    ## plot results
    if ( !missing(out.file) ) {
        file.name <- paste(out.file, ".png",sep="")
        png(file.name, width=6,height=rws*.8,res=200,units="in")
    }
    orig.par <- par(c("mfcol","mai","mgp"))
    par(mfcol=c(rws,1),mai=c(.01,.5,.01,.01),mgp=c(1.7,.5,0),xaxs="i")
    
    ## plot original data
    if( !missing(tot) ) {
        plot(x,tot,log=ifelse(use.log,"","y"),type="l",lwd=2,axes=FALSE,xlim=xlim)
        axis(1);axis(2);
    }
    if( !missing(ts) ) {
        #require("colorRamps")
        colors0 <- rainbow(100) #matlab.like(100)  ## read count
        colors0[1] <- "#FFFFFF" ## replace minimal by white
        ts.nrm <- t(apply(ts,1,function(x) {
            if ( sum(is.na(x))==length(x) ) return(x)
            mx <- max(x,na.rm=TRUE)
            mn <- min(x,na.rm=TRUE)
            if ( mx==mn ) rep(NA,length(x)) else ((x-mn)/(mx-mn))}))
        image(ts.nrm,col=colors0,axes=FALSE,xlim=xlim)
    }
    
    ## plot original clustering
    plot(x,seq,axes=FALSE,xlab="i",ylab="cluster",
         col=cols[as.character(seq)],cex=1,pch=16,xlim=xlim)
    axis(2)
    points(which(seq==0),seq[seq==0],pch=16,cex=1)
    
    ## plot dyn.prog results
    for ( scor in scores ) {
        multS <- scrR[[scor]]
        multi <- names(multS)
        multi <- multi[multi!="SM"]
        ## Scoring Matrix S(c,i) - should be the same over multi
        S <- multS[[multi[1]]]$SK$S # total score S(c,i)
        matplot(S,axes=FALSE,xlab="i",ylab="score S(i,c)",type="l",lty=1,lwd=1,
                col=cols[1:ncol(S)],xlim=xlim)
        axis(1)#,at=x,labels=seq,cex.axis=.6);axis(2)
        ##axis(3,at=x,cex.axis=.6,tcl=.05,mgp=c(0,-.75,0))
        legend("left",legend=scor,bty="n")
        
        for ( mult in multi ) {
            ## back-tracing matrix K(c,i) - should be the same over multib
            K <- multS[[mult]]$SK$K
            matplot(K,axes=FALSE,xlab="i",ylab=paste(multS[[mult]]$SK$mink,"k"),
                    type="l",lty=1,lwd=1, col=cols[1:ncol(K)],xlim=xlim)
            legend("left",legend=paste("scoring:",mult),bty="n")
            axis(1)#,at=x,labels=seq,cex.axis=.6);axis(2)
            
            ## plot segments
            multib <- names(multS[[mult]])
            multib <- multib[!multib%in% c("SK")]
            
            yl <- length(multS[[mult]])-1
            par(mgp=c(2,.5,0))
            plot(NA,ylim=c(0,yl+1),axes=FALSE,
                 ylab="back-tracing",xlab=NA,xlim=xlim)
            
            for ( k in 1:length(multib) ) {
                multb <- multib[k]
                segments <- multS[[mult]][[multb]]
                if ( nrow(segments)>0 )
                    arrows(x0=segments[,2],y0=0+k,x1=segments[,3],
                           col=cols[as.character(segments[,1])],
                           length=0.075,code=3,angle=90,lwd=2)
            }    
            axis(1)#,at=x,labels=seq,cex.axis=.6)
            axis(2,at=1:length(multib),labels=multib,las=2)
        }
    }
    par(orig.par)
    if ( !missing(out.file) )  {
        dev.off()
        if ( verb>0 )
            cat(paste("plotted\t", file.name, date(), "\n"))
        return(file.name)
    }
    if ( verb>0 ) cat(paste("done\t", date(), "\n"))
}

### USING CLASSES
#' plot the processed time-series object returned from
#' \code{\link{processTimeseries}}.
#' @param tset the time-series object returned by
#' \code{\link{processTimeseries}
#' @param plot a string vector indicating the values to be plotted;
#' that is: `total' for a plot of the total signal, summed over
#' the time-points, and indicating the applied threshold \code{low.thresh};
#' note that the total levels may have been transformed (e.g. by \code{log_1}
#' or \code{ash}) depending on the arguments \code{trafo} and \code{dc.trafo}
#' in \code{\link{processTimeseries}
#'@export
plot.tset <- function(tset, plot=c("total","timeseries")) {

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

    ##totlog <- tset$settings$trafo!="raw" | tset$settings$dc.trafo!="raw"
    low.thresh <- tset$settings$low.thresh
    
    if ( "total" %in% plot ) {
        plot(1:N,tot,log="",type="l",lwd=2,axes=FALSE,ylab=NA,xlab=NA)
        polygon(x=c(1,1,N,N),y=c(min(tot,na.rm=TRUE),rep(low.thresh,2),
                               min(tot,na.rm=TRUE)),col="#00000055",border=NA)
        abline(h=low.thresh,col="#000000BB")
        lines(1:N,tot)
        axis(2);
        axis(1)
        mtext("total signal", 2, 2)
    }
    if ( "timeseries" %in% plot ) {
        segment.plotHeat(ts,coors=coors,chrS=0,colors=colors0, colnorm=TRUE)
        axis(2,at=1:ncol(ts))
        axis(1)
        mtext("time points", 2, 2)
    }
}

#' plot the clustering of each time-series; clusters will be sorted
#' and then colored by their similarity (see \code{\link{sortClusters}}
#' and \code\link{colorClusters}}
#' @param cset a set of clusterings as returned by
#' \code{\link{clusterTimeseries}}
#' @param k a numeric vector indicating the clusterings to be plotted;
#' specifically the column number in the matrix of clusterings
#' in \code{cset$clusters}; if missing all columns will be plotted
#' and the calling code must take care of properly assigning \code{par(mfcol)}
#' or \code{layout} for the plot
#' @param x optinally x-values can be passed here to use real (chromosomal)
#' coordinates for the x-axis
plot.cset <- function(cset, k, x) {
    
  
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
        if ( missing(x) )
          x <- 1:length(seq)
        y <- 1:length(cls.srt)
        names(y) <- cls.srt
        cols <- cset$colors[[i]]
        ## plot original clustering
        plot(x,y[seq],axes=FALSE,xlab="",ylab="cluster",
             col=cols[seq],cex=1,pch=16)
        axis(2, at=y, labels=names(y), las=2)
    }
}

plot.sset <- function(sset, types, x, plot=c("segments", "S", "S1")) {

    if ( missing(types) ) {
        types <- rownames(sset$settings)
    }
      
    if ( "segments" %in% plot ) {
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
    
    colors0 <- rev(grDevices::gray.colors(100)) 
    colors0[1] <- "#FFFFFF" ## replace minimal by white

    if ( any(c("S","S1") %in% plot) ) {
        SK <- sset$SK[types]
        sk.srt <- sset$sorting[types]
        sk.col <- sset$colors[types]
        #if ( missing(k) ) 
            k <- 1:length(SK)
        #if ( length(k)>1 )
        #    par(mfcol=c(length(k),1))
        for ( j in k ) {

            ## sorting
            srt <- as.character(sk.srt[[j]])
            ## coloring; add black for nuissance
            sgcols <- c("#000000", sk.col[[j]][srt])
            srt <- c("0",srt)

            ## plot S1 as heatmap
            if ( "S1" %in% plot ) {

                S1 <- t(SK[[j]]$S1[,rev(srt)])
                S1 <- S1/apply(S1,1,mean)
                image_matrix(S1,axis=2, col=colors0, ylab="cluster")
            }
            
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
            if ( missing(x) )
              x <- 1:nrow(S)
            xlim <- range(x) 
            

            ##cat(paste(paste(range(ash(dS)),collapse="-"),"\n"))
            xrng <- stats::quantile(x,c(.05,.95))
            xidx <- which(x>xrng[1]&x<xrng[2]) #x%in%xrng[1]:xrng[2]
            ylim <- stats::quantile(ash(dS[xidx,]),c(0,1))
            plot(1,ylim=ylim,xlim=xlim,ylab=expression(ash(Delta~S["i,C"])))
            lines(x,ash(dS[,1]),lwd=7,col="#00000015") # NUI: BACKGROUND GRAY
            lines(x,ash(dS[,1]),lwd=1,lty=3,col="#00000099") # NUI: BACKGROUND GRAY
            matplot(x, ash(dS), type="l", lty=1, lwd=1, add=TRUE, col=sgcols)
            graphics::mtext(names(SK)[j], side=2 , line=4, las=2)

            
        }
    }
}
