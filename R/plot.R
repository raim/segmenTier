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

#' converts the vector form of the triangular scoring matrix
#' into a matrix
#' @param SV vector form of the triangular scoring matrix as returned
#' by \code{\link{calculateScoringMatrix}}
#'@export
fillScoringMatrix <- function(SV) {
    ## L <- (N+1)*N/2
    L <- length(SV)
    N <- -.5 + sqrt(.25 + 2*L)
    SM <- matrix(NA,nrow=N, ncol= N)
    getIdx <- function(i,k) ((i + 1) * i / 2 + k)
    for ( i in 0:(N-1) ) {
        idx <- getIdx(i,i)
        SM[i+1,i+1] <- SV[idx+1] 
        for ( k in (i-1):0 ) {
            if ( k<0 ) next
            idx <- getIdx(i,k)
            SM[k+1,i+1] <- SV[idx+1] 
        }
    }
    SM
}

#' plot the scoring function matrices as a heatmap
#' @param SML a list with scoring function matrices for each cluster,
#' as provided by \code{\link{calculateScoringMatrix}}
#' @param seq the original cluster sequence, optional for axis labeling
#' @param score name of the used scoring function
#' @param out.file if supplied the scoring matrices will be plotted to
#' individual png files named <out.file>_<number>.png
#' @param verb level of verbosity; 0: no output, 1: progress messages
#' @export
plotScoring <- function(SML, seq, score, out.file, verb=2) {
    files <- NULL
    for ( c in 1:length(SML) ) {
        if ( !missing(out.file) ) {
            file.name <- paste(out.file,"_",c,".png",sep="")
            files <- c(files, file.name)
            png(file.name,width=5,height=5,res=200,units="in")
        }
        ## expand scoring function vector to triangular matrix
        SM <- fillScoringMatrix(SML[[c]])
        image(x=1:nrow(SM),y=1:nrow(SM),z=SM,axes=FALSE,
              main=paste("scoring function:",
                         ifelse(missing(score),"",score)),
              ylab=NA,xlab=NA)
        if ( !missing(seq) ) {
            axis(2,at=1:nrow(SM),labels=seq,cex.axis=.5,las=2)
            axis(3,at=1:nrow(SM),labels=seq,cex.axis=.5,las=2)
        }
        legend("bottomright",paste("cluster", c),cex=2,bty="n")
        
        if ( !missing(out.file) )
            dev.off()
        else {
            cat(paste("plotted cluster ", c,
                      ". please enter to proceed to the next plot"))
            scan()
        }
    }
    if ( !missing(out.file) ) {
        if ( verb>0 )
            cat(paste("plotted", paste(files,collapse=" ; "), date(), "\n"))
        res <- return(files)
    }
}

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
plotSegments <- function(scrR, seq, ts, tot, out.file, use.log=FALSE,
                         add.plots=0, verb=2) {
  
    if ( verb>0 )
        cat(paste("plotting",
                  ifelse(missing(out.file),"results",out.file),"\t",date(),"\n"))

    ## if nuissance cluster is present, increase
    ## coloring by 1
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

