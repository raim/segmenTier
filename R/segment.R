#' segmenTier : cluster-based segmentation
#' from a sequential clustering
#'@author Rainer Machne, Peter F. Stadler
#'@docType package
#'@name segmenTier
#'@section Dependencies: The package strictly depends on \code{package::Rcpp},
#' and \code{package:parallel} allows to signficantly speed up scoring
#' function matrix calculations. All other dependencies are
#' usually present in a basuc installation (\code{graphics}, \code{grDevices})).
#'@importFrom Rcpp evalCpp
#'@importFrom parallel mclapply
#'@importFrom graphics image axis par plot points lines legend arrows matplot
#'@importFrom  grDevices png dev.off rainbow
#'@useDynLib segmenTier
NULL # this just ends the global package documentation



### DYNAMIC PROGRAMMING BASED SEGMENTATION OF A CLUSTERING
### implemented by Rainer Machne, hopefully
### as conceived by Peter F. Stadler

### NOTE: these functions are tested in $GENBRO/src/segment_test.R !

## PROBLEM
## find optimal segments k:i in a sequence of clusters
##      k      i     N
## 1333533335134424413   
## SCORING (dyn.prog.) - dynamically fill the matrix S(i,c):
## for i=1..N   // position
##   for c in 1:4 // clusters 1..5
##     // find the segmentation points k at the maximal sum of
##     // S(k-1,c')    , which is the end of the best prev. local segment c'!=c
##     //   PLUS
##     // score(k,i,c) , which is a measure of cluster c enrichment in k..i
##     S[i,c] <- max_c' max_k S(k-1,c') + score(k,i,c)
##     // and store the used k for back-tracing
##     K[i,c] <- k  // the k which delivered S[i,c]
## BACKTRACING: 
## ... and find segments by back-tracing the cluster c and position
## k which had delivered the maximal score S(c,i) at each i
## find max score in S(i,c) at i=N
## i = N
## while(i != 0)
##   for all k<i und color c'!= c:
##      if S(i,c) == S(k-1,c') + score[k,c,i] )
##      break;
##   output:  interval k bis i mit Farbe c
##   i <- k - 1 // nextmax: search next non-decreasing S(i,c)
##   c <- c'




### FUNCTIONS

### UTILS
warn <- function(w, warnings,verb=FALSE) {
  if (verb) cat(w)
  c(warnings,w)
}

### HIGH-LEVEL WRAPPERS

## TODO: high-level wrapper that takes a time-series as input
## and clusters the time-series  calling segmentClusters
## additionally reports data medians/centers for all segments
## and cluster centers; can be used as input for clusterSegments
## to finally cluster all segments genome-wide
segmentData <- function() {}

## TODO: high-level wrapper that segments genome data into primary domains,
## and sub-divides these into coherent segments based on clustering
## of the data and a dynamic programming algo; ....
clusterSegments <- function() {}

#' segmenTier's main wrapper interface, calculates segments from a
#' clustering sequence.
#' @param seq a clustering sequence. The only strict requirement is that
#' nuissance clusters (which will not be segmented) have to be numeric or
#' character "0" (zero).
#' @param csim cluster-cluster or position-cluster similarity
#' matrix, for scoring functions ccor and icor, respectively
#' @param score the scoring function to be used: "ccor", "icor" or "cls"
#' @param M minimal sequence length; Note, that this is not a strict
#' cut-off but defined as a penalty that must be "overcome" by good score.
#' @param Mn minimal sequence length for nuissance cluster, Mn<M will allow
#' shorter distances between segments; only used in scoring functions
#' "ccor" and "icor" 
#' @param a an additional penalty only used for pure cluster-based
#' scoring w/o cluster similarity measures in scoring function "cls"
#' @param nui the similarity score to be used for nuissance clusters in the
#' cluster similarity matrices
#' @param nextmax go backwards while score is increasing before openening a
#' new segment, default is TRUE
#' @param multi handling of multiple k with max. score in forward phase,
#' either "min" (default) or "max"
#' @param multib handling of multiple k with max. score in back-trace phase,
#' either "min" (default), "max" or "skip"
#' @param ncpu number of available cores (CPUs), passed to
#' \code{parallel::mclapply} by \code{\link{calculateScoringMatrix}}
#' @param verb level of verbosity, 0: no output, 1: progress messages
#' @param save.mat store the scoring function matrix SM or the back-tracing
#' matrix K by adding "SM" and "SK" to the string vector save.mat; useful
#' in testing stage or for debugging or illustration of the algorithm;
#' see \code{\link{plotScoring}}
#' @details This is the main R wrapper function for the segmentation algorithm.
#' It takes a sequence of clusterings and returns segments of
#' consistent clusters. It runs the dynamic programing algorithm for
#' a selected scoring function and an according cluster similarity matrix,
#' followed by the  back-tracing step to find segments.
#' Some more details of the algorithm can be tuned, but these usually
#' have little effect on real-life data sets.
#' @return Returns a list containing the main result ("segments"), "warnings"
#' from the dynamic programing and back-tracing phases, and optionally (see
#' option \code{save.mat} the score function matrices \code{SM}, the
#' total score matrix \code{S(i,c)} and the backtracing matrix \code{K(i,c)}.
#' The main result structure "segments" is a 3-column matrix, where column 1
#' is the cluster assignment and colums 2 and 3 are start and end position
#' of the segments.
#' @export
segmentClusters <- function(seq, csim, score="ccor", M=175, Mn=20, a=2, nui=1,
                            nextmax=TRUE, multi="min",multib="min", 
                            ncpu=1, verb=1, save.mat="") {
    
    ## 1: set up sequence and data
    N <- length(seq)
    seqr <- seq

    ## 1a: map to internal 1:K clustering:
    map <- sort(unique(seqr)) # clusters
    map <- map[map!=0] ##  nuissance cluster
    names(map) <- map
    map[] <- 1:length(map)
  
    seqr <- map[as.character(seq)]
    seqr[seq==0] <- 0      # original nuissance clusters
    seqr[is.na(seqr)] <- 0 # replace NA by nuissance
        
    ## 1b: add nuissance clusteer if present:
    ## columns and rows are be added to the similarity matrices, using
    ## nui and -nui as "correlations";
    ## cluster index increased by +1, and the nuissance cluster will be "1"
    ## throughout further processing!

    nui.present <- FALSE
    if ( 0 %in% seqr & score!="cls" ) {
        
        nui.present <- TRUE
        ## increase clustering by +1
        ## nuissance cluster will be cluster 1 
        seqr <- seqr + 1
        ## 1c: add data for nuissance cluster to matrices cr and P
        if ( score=="icor" ) {
            ## cor(i,c) - similarity of position i to cluster medians
            ## reduce passed matrix to actually present clusters
            csim <- csim[,as.numeric(names(map))]
            ## add nuissance cluster
            csim <- cbind(rep(-nui,N),csim)
            csim[seqr==1,] <- -nui
            csim[seqr==1,1] <- nui
        }
        if ( score %in% c("ccor","xcor") ) {
            ## cor(c,c) - similarity of cluster medians
            ## reduce passed matrix to actually present clusters
            csim <- csim[as.numeric(names(map)),as.numeric(names(map)),
                         drop=FALSE] # allow single cluster? or abort here?
            ## add nuissance cluster
            csim <- rbind(rep(-nui,nrow(csim)+1),
                          cbind(rep(-nui,nrow(csim)), csim))
            csim[1,1] <- nui
        }
        ## experimental: testing scaling of correlations
        if ( score == "xcor" ) {
            csim <- csim^3 # testing scaling of correlations
            csim[csim<0] <- csim[csim<0] -1 # punish neg.cor
        }
        map <- c('0'=1, map  + 1)
    }

    ## for the pure cluster segmentation pass par. a to ccSMcls
    if ( score=="cls" ) csim <- a 
    
    ## get clusters
    C <- sort(unique(seqr))
    C <- C[C!=0] # rm nuissance - should only be there for "cls"

    ## 2: generate scoring function matrices (TODO: lists, to save mem)
    if ( verb>0 ) 
      cat(paste("scoring function", score, "\t", date(), "\n"))
    SM <- calculateScoringMatrix(seqr, C=C, score=score, M=M, Mn=Mn,
                                 csim=csim, ncpu=ncpu)

    ## 3: calculate total scoring S(i,c) and backtracing K(i,c)
    if ( verb>0 )
        cat(paste("\ttotal score with", multi, "\t", date(), "\n"))
    SK <- calculateTotalScore(seq=seqr, C=C, SM=SM, multi=multi)
    
    ## 4: back-tracing to generate segments
    if ( verb>0 )
        cat(paste("\tbacktracing with", multib, "\t", date(), "\n"))
    seg <- backtrace(S=SK$S, K=SK$K, multib=multib, nextmax=nextmax, verb=verb)

    ## map back segments to original
    remap <- as.numeric(names(map))
    seg$segments[,1] <- remap[seg$segments[,1]]

    ## rm nuissance segments
    seg$segments <- seg$segments[seg$segments[,1]!=0,,drop=FALSE]
    
    ## add matrices if requested!
    ## ... can be used for plotting or re-analysis
    if ( "SM" %in% save.mat ) seg$SM <- SM
    if ( "SK" %in% save.mat ) seg$SK <- SK
    
    return(seg)
    
}

## PLOTTING RESULTS

#' plot the scoring function matrices as a heatmap
#' @param SM a list with scoring function matrices for each cluster,
#' see \code{\link{calculateScoringMatrix}}
#' @param seq the original cluster sequence, optional for axis labeling
#' @param out.file if supplied the scoring matrices will be plotted to
#' individual png files named <out.file>_<number>.png
#' @param verb level of verbosity; 0: no output, 1: progress messages

#' @export
plotScoring <- function(SM, seq, out.file, verb=2) {
    files <- NULL
    for ( c in 1:length(SM) ) {
        if ( !missing(out.file) ) {
            file.name <- paste(out.file,"_",c,".png",sep="")
            files <- c(files, file.name)
            png(file.name,
                width=5,height=5,res=200,units="in")
        }
        image(x=1:nrow(SM[[c]]),y=1:nrow(SM[[c]]),z=SM[[c]],axes=FALSE,main=c)
        if ( !missing(seq) ) {
            axis(2,at=1:nrow(SM[[c]]),labels=seq,cex.axis=.5,las=2)
            axis(3,at=1:nrow(SM[[c]]),labels=seq,cex.axis=.5,las=2)
        }
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
 
    ## plot results
    if ( !missing(out.file) ) {
        file.name <- paste(out.file, ".png",sep="")
        png(file.name, width=6,height=rws*.8,res=200,units="in")
    }
    par(mfcol=c(rws,1),mai=c(.01,.5,.01,.01),mgp=c(1.7,.5,0),xaxs="i")
    
    ## plot original data
    if( !missing(tot) ) {
        plot(x,tot,log=ifelse(use.log,"","y"),type="l",lwd=2,axes=FALSE)
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
        image(ts.nrm,col=colors0,axes=FALSE)
    }
    
    ## plot original clustering
    plot(x,seq,axes=FALSE,xlab="i",ylab="cluster",
         col=cols[as.character(seq)],cex=1,pch=16)
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
                col=cols[1:ncol(S)])
        axis(1)#,at=x,labels=seq,cex.axis=.6);axis(2)
        ##axis(3,at=x,cex.axis=.6,tcl=.05,mgp=c(0,-.75,0))
        legend("left",legend=scor,bty="n")
        
        for ( mult in multi ) {
            ## back-tracing matrix K(c,i) - should be the same over multib
            K <- multS[[mult]]$SK$K
            matplot(K,axes=FALSE,xlab="i",ylab=paste(multS[[mult]]$SK$mink,"k"),
                    type="l",lty=1,lwd=1, col=cols[1:ncol(K)])
            legend("left",legend=paste("scoring:",mult),bty="n")
            axis(1)#,at=x,labels=seq,cex.axis=.6);axis(2)
            
            ## plot segments
            multib <- names(multS[[mult]])
            multib <- multib[!multib%in% c("SK")]
            
            yl <- length(multS[[mult]])-1
            par(mgp=c(2,.5,0))
            plot(NA,xlim=c(1,length(seq)),ylim=c(0,yl+1),axes=FALSE,
                 ylab="back-tracing",xlab=NA)
            
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
    if ( !missing(out.file) )  {
        dev.off()
        if ( verb>0 )
            cat(paste("plotted\t", file.name, date(), "\n"))
        return(file.name)
    }
    if ( verb>0 ) cat(paste("done\t", date(), "\n"))
}


### SEGMENTATION BY DYNAMIC PROGRAMMING - MAIN FUNCTIONS
## NOTE that most of the work is done in segment.cpp
## using the Rcpp interface to C++

## SCORING MATRIX GENERATION
## TODO: instead generate a list of length seq, each
## containing a vector or kmin:i
## TODO: choose a reasonable kmin, and fill up rest only
## if requested
## NOTE: in parallel mode preschedule FALSE will take longer
## to collect data from forked processes, but this avoids
## the error "long vectors not supported yet: fork.c:378"

#' generate scoring function matrices for each cluster
#' @param seq a sequence of cluster assignments
#' @param C optional cluster sorting
#' @param score the scoring function to be used: ccor, icor or cls
#' @param M minimal sequence length; Note, that this is not a strict
#' cut-off but defined as a penalty that must be "overcome" by good score.
#' @param Mn minimal sequence length for nuissance cluster, Mn<M will allow
#' shorter distances between segments; only used in scoring functions
#' "ccor" and "icor" 
#' @param csim cluster-cluster or position-cluster similarity
#' matrix, for scoring functions ccor and icor, respectively
#' @param ncpu number of available cores (CPUs), passed to
#' \code{parallel::mclapply}
#' @param preschedule \code{parallel::mclapply} option that currently
#' requires to be set to FALSE to avoid an error in data collection from
#' parallel processes
#' @return Returns the scoring function matrices \code{SM} for all clusters
#' in the sequence \code{seq}.
#' @export
calculateScoringMatrix <- function(seq, C, score="ccor", M, Mn, csim, 
                                   ncpu=1, preschedule=FALSE) {
  if ( missing(C) ) 
    C <- sort(unique(seq)) # get clusters

  ## get requested scoring functions - from segment.cpp !
  ## NOTE: csim is an integer for scoring function "cls"; and
  ## a matrix NxC for "icor"; and CxC for "cor"
  score <- paste("ccSM",score,sep="")
  getMat <- get(score,mode="function") 

  if ( ncpu<2 ) { # single CPU
      return(lapply(C, function(c) getMat(seq, c, M, Mn, csim)))
  } else { # multiple CPUs: use package parallel!
      
      ## TODO: make parLapply work
      options(warn=0)
      x <- parallel::mclapply(C, function(c) getMat(seq, c, M, Mn, csim),
                              mc.cores=ncpu,mc.preschedule=preschedule)
      options(warn=0)
      return(x)
  }
}


#' back-tracing : collect clustered segments from the scoring function matrix
#' @param S matrix S, containing the local scores
#' @param K matrix K, containing the position k used for score maximization
#' @param multib if multiple k produce the maximal score, take either the
#' shortest k ("max") or the longest k ("min"); if multib is set to "skip"
#' the next unique k will be searched
#' @param nextmax proceed backwards while score is increasing before
#' openening a new segment
#' @param verb print messages
#' @export
backtrace <- function(S, K, multib, nextmax=FALSE, verb=TRUE) {

    segments <- NULL
    i <- nrow(S)
    warnings <- NULL
    
    while ( i>0 ) {
        
        ##  WHICH cluster had max S at i?
        c <- which(S[i,]==max(S[i,]))
        
        ## SEARCH END
        ## search next max(S[i,c]) over i--
        if ( nextmax ) {
            while( i>0 & sum(S[i,c] <= S[i-1,c])==length(c) ) 
              i <- i - 1
            ##  which cluster had max S at i?
            c <- which(S[i,]==max(S[i,]))
        }
        
        ## WHICH k was used?
        k <- K[i,c]
        
        ## handle multiple max scores?
        if ( length(c)>1 ) {
            w <- paste(i,"back-trace warning:", length(c),
                       "c:", paste(c,collapse=";"),
                       "with max_k:", paste(k,collapse=";"),
                       "\ttaking", multib, "\n")
            warnings <- warn(w,warnings,verb=verb)
            if ( multib=="skip" ) {
                i <- i-1
                next 
            }
            ## max: shortest k, min: longest k
            km <- get(multib,mode="function")(k)
            c <- c[which(k==km)] 
            k <- km
        }
        if ( length(c)>1 ) {
            w <- paste(i,"back-trace warning STILL:", length(c), "c:",
                       paste(c,collapse=";"),"\n") 
            warnings <- warn(w,warnings,verb=verb)
            c <- c[1] #paste(c,collapse=";")
        }
        
        ## ignore k==i segments
        ## NOTE this should only happen for obsolete k<=i
        if ( i %in% k ) {
            w<- paste(i,"back-trace warning:",sum(k%in%i),"0 length segments\n")
            warnings <- warn(w,warnings,verb=verb)
            c <- c[k!=i]
            k <- k[k!=i]
        }
        if ( length(k)==0 ) {
            w <- paste(i,"back-trace warning: no k left\n")
            warnings <- warn(w,warnings,verb=verb)
            i <- i-1
            next
        }
        
        segments <- rbind(segments, c(c,k,i))
        i <- k -1
    }
    segments <- segments[order(as.numeric(segments[,2])),,drop=FALSE]
    list(segments=segments,warnings=warnings)
}
    
    
