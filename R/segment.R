#' segmenTier : cluster-based segmentation
#' from a sequential clustering
#'@author Rainer Machne \email{raim@tbi.univie.ac.at}, Peter F. Stadler \email{studla@bioinf.uni-leipzig.de}
#'@docType package
#'@name segmenTier
#'@section Dependencies: The package strictly depends only on
#' \code{\link[Rcpp:Rcpp]{Rcpp}}.
#' Package \code{\link[parallel:parallel]{parallel}} allows to
#' speed up scoring function matrix calculations when more then 1 cores
#' (CPUs) are available. All other dependencies are usually present in a
#' basic installation (\code{stats}, \code{graphics}, \code{grDevices})).
#' @references
#'   @bibliography segmenTier.bib
#'@importFrom Rcpp evalCpp
#'@importFrom stats qt sd
#'@importFrom parallel mclapply
#'@importFrom graphics image axis par plot matplot points lines legend arrows strheight strwidth text
#'@importFrom grDevices png dev.off rainbow gray xy.coords
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

### MESSAGE UTILS

## nicer timestamp
time <- function() format(Sys.time(), "%Y%m%d %H:%M:%S")
## messages
msg <- function(x) cat(x, file=stdout()) # until piping is implemented
## stored "warnings" (actually detailed messages)
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
#' matrix, for scoring functions "ccor" and "icor", respectively
#' @param csim.scale exponent to scale similarity matrices, must be odd
#' to maintain negative correlations!
#' @param cset alternatively to arguments \code{seq} and \code{csim}, a
#' set of clusterings as returned by \code{\link{clusterTimeseries}} can
#' be provided; this requires the additional argument \code{k} to select
#' the kth clustering from the set
#' @param k the kth clustering of argument \code{cset} will be used
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
#' \code{\link[parallel:mclapply]{parallel::mclapply}} by
#' \code{\link{calculateScoringMatrix}}
#' @param verb level of verbosity, 0: no output, 1: progress messages
#' @param save.matrix store the total score matrix \code{S(i,c)} and the
#' backtracing matrix \code{K(i,c)}; useful in testing stage or for
#' debugging or illustration of the algorithm; see \code{\link{plotScoring}}
#' @details This is the main R wrapper function for the segmentation algorithm.
#' It takes a sequence of clusterings and returns segments of
#' consistent clusters. It runs the dynamic programing algorithm for
#' a selected scoring function and an according cluster similarity matrix,
#' followed by the  back-tracing step to find segments.
#' Some more details of the algorithm can be tuned, but these usually
#' have little effect on real-life data sets.
#' @return Returns a list containing the main result ("segments"), "warnings"
#' from the dynamic programing and back-tracing phases, and optionally (see
#' option \code{save.matrix}) (\code{results$SK}) the total score matrix
#' \code{S(i,c)} and the backtracing matrix \code{K(i,c)}.
#' The main result structure "segments" is a 3-column matrix, where column 1
#' is the cluster assignment and colums 2 and 3 are start and end position
#' of the segments.
#' @export
segmentClusters <- function(seq, csim, csim.scale=1,
                            cset, k=1,
                            score="ccor",
                            M=175, Mn=20, a=-2, nui=1,
                            nextmax=TRUE, multi="max",multib="max", 
                            ncpu=1, verb=1, save.matrix=FALSE) {

    stime <- as.numeric(Sys.time()) 
    
    ## cluster set from clusterTimeseries
    if ( !missing(cset) ) {
        seq <- cset$clusters[,k]
        if ( score=="ccor" ) csim <- cset$Ccc[[k]]
        if ( score=="icor" ) csim <- cset$Pci[[k]]
    }
    
    ## 1: set up sequence and data
    N <- length(seq)
    seqr <- seq

    ## 1a: map to internal 1:K clustering:
    ## TODO: allow character clusters!?
    map <- sort(unique(seqr)) # clusters
    map <- map[map!=0] ##  nuissance cluster
    names(map) <- map
    map[] <- 1:length(map)
  
    seqr <- map[as.character(seq)]
    seqr[seq==0] <- 0      # original nuissance clusters
    seqr[is.na(seqr)] <- 0 # replace NA by nuissance
        
    ## set-up similarity matrix for ccls
    ## internally 'ccor' is used, and we set up the
    ## cluster-cluster similarity function (matrix) here
    if ( score=="ccls" ) {
        L <- length(unique(seqr))
        csim <- matrix(a, nrow=L, ncol=L) # Delta(C,D!=C) = a
        diag(csim) <- 1 # Delta(C,C) = 1
        nui <- -a
    }
    
    ## 1b: add nuissance cluster if present:
    ## columns and rows are be added to the similarity matrices, using
    ## nui and -nui as "correlations";
    ## cluster index increased by +1, and the nuissance cluster will be "1"
    ## throughout further processing!

    nui.present <- FALSE
    if ( 0 %in% seqr ) {
        
        nui.present <- TRUE
        ## increase clustering by +1
        ## nuissance cluster will be cluster 1 
        seqr <- seqr + 1 ## TODO: get rid of this, avoid correction in .cpp
        
        if ( score=="icor" ) {
            ## cor(i,c) - similarity of position i to cluster medians
            ## reduce passed matrix to actually present clusters
            csim <- csim[,as.numeric(names(map))]
            ## add nuissance cluster
            csim <- cbind(rep(-nui,N),csim)
            csim[seqr==1,] <- -nui
            csim[seqr==1,1] <- nui
        }
        if ( score %in% c("ccor") ) {
            ## cor(c,c) - similarity of cluster medians
            ## reduce passed matrix to actually present clusters
            csim <- csim[as.numeric(names(map)),as.numeric(names(map)),
                         drop=FALSE] # allow single cluster? or abort here?
            ## add nuissance cluster
            csim <- rbind(rep(-nui,nrow(csim)+1),
                          cbind(rep(-nui,nrow(csim)), csim))
            csim[1,1] <- nui
        }
        map <- c('0'=1, map  + 1)
    }
    
    ## get clusters
    C <- sort(unique(seqr))
    ##C <- C[C!=0] # rm nuissance - should only be there for "cls"

    ## scale similarity matrix!
    sgn <- sign(csim) # store sign
    csim <- csim^csim.scale # scale matrix
    ## if exponent is even (checking within machine tolerance)
    ## the sign is re-added
    if ( csim.scale %% 2 < .Machine$double.eps^0.5 )
      csim <- sgn*csim # restore sign
           #warning("csim.scale should be odd: ", csim.scale)
 
    #SM <- calculateScoringMatrix(seqr, C=C, score=score, M=M, Mn=Mn,
    #                             csim=csim, ncpu=ncpu)

    ## 2: calculate total scoring S(i,c) and backtracing K(i,c)
    if ( verb>0 ) {
        cat(paste("Scoring matrix\t", time(), "\n",sep=""))
        cat(paste("parameters\t",paste("function:", score,
                                        "; scale:", csim.scale,
                                        "; max/min:", multi,sep=""),
                  "\n",sep=""))
    }
    ## TODO: handle Mn in scoring functions
    ## add official nuissance cluster

    #SK<- calculateTotalScore_test(seq=seqr,C=C,SM=SM,
    #                              csim=csim,M=M,Mn=Mn,multi=multi)
    SK<- calculateScore(seq=seqr, C=C, score=score, csim=csim,
                        M=M, Mn=Mn, multi=multi)

    ## 4: back-tracing to generate segments
    if ( verb>0 ) {
        cat(paste("Backtracing\t", time(), "\n",sep=""))
        cat(paste("parameters\t", paste("multib:",multib,sep=""), "\n",sep=""))
    }
    seg <- backtrace(S=SK$S, K=SK$K, multib=multib, nextmax=nextmax, verb=verb)

    ## remap: map back to original cluster names
    remap <- as.numeric(names(map))
    ## re-name to original clusters if stored
    #if ( "SM" %in% save.mat )
    #    names(SM) <- remap[as.numeric(names(SM))]
    ## map back segments to original
    seg$segments[,1] <- remap[seg$segments[,1]]

    ## rm nuissance segments
    seg$segments <- seg$segments[seg$segments[,1]!=0,,drop=FALSE]
    
    ## add matrices if requested!
    ## ... can be used for plotting or re-analysis
    #if ( "SM" %in% save.mat ) seg$SM <- SM
    if ( save.matrix ) seg$SK <- SK

    seg$csim <- csim
    class(seg) <- "segments"
    
    if ( verb>0 ) 
        cat(paste("Done at  \t", time(), "\n",sep=""))

    if ( verb>0 ) {
        etime <- as.numeric(Sys.time())
        cat(paste("elapsed\t", round((etime-stime)), " sec\n",sep=""))
    }
    return(seg)
    
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



## TODO: move this to .cpp as well, then the whole algo is available in C++
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
        
        ## FIND SEGMENT END
        ## Note, that this determines the segment's cluster assigment,
        ## unless multiple clusters deliver the maximal score
        
        ##  WHICH cluster(s) had max S at i?
        c <- which(S[i,]==max(S[i,]))
        
        ## search next max(S[i,c]) over i--
        if ( nextmax ) {
            while( i>0 & sum(S[i,c] <= S[i-1,c])==length(c) ) 
              i <- i - 1
            ##  which cluster had max S at i?
            c <- which(S[i,]==max(S[i,]))
        }

        ## FIND SEGMENT START
        ## WHICH k was used?

        ## get the k that was used for the maximum score S(i,c)
        k <- K[i,c]

        ## ... and handle multiple max scores from several clusters
        if ( length(c)>1 ) {
            w <- paste(i,"back-trace warning:", length(c),
                       "clusters with maximal score; c:", paste(c,collapse=";"),
                       "with max_k:", paste(k,collapse=";"), ":",
                       ifelse(multib=="skip", paste("skipping breakpoint",i),
                              ifelse(multib=="max",
                                     "taking shorter segment",
                                     "taking longer segment")), "\n")
                                     
            warnings <- warn(w,warnings,verb=verb)
            if ( multib=="skip" ) {
                i <- i-1
                next 
            }
            ## multib: max is shortest k, min is longest k
            ## i.e. max delivers the shorter segment, min the longer segment 
            km <- get(multib,mode="function")(k)
            c <- c[which(k==km)] 
            k <- km
        }
        ## if we still can't find one cluster 
        if ( length(c)>1 ) {
            w <- paste(i,"back-trace warning STILL:", length(c),
                       "clusters with maximal score; c:",
                       paste(c,collapse=";"),"\n") 
            warnings <- warn(w,warnings,verb=verb)
            c <- c[1] #paste(c,collapse=";")
        }
        
        ## ignore k==i segments
        ## NOTE this should only happen for obsolete k<=i instead of k<i
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
    if ( !is.null(segments) )
        colnames(segments) <- c("CL","start","end")
    list(segments=segments,warnings=warnings)
}
    
    
### DATA SET DOC

#' transcriptome time series data from a region encompassing
#' four genes and regulatory upstream non-coding RNA
#'
#' @name tsd
#' @docType data
NULL
