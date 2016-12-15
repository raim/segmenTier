### TIME-SERIES CLUSTERING PARAMETERS

## get Discrete Fourier Transformation
get.fft <- function(x) {
    n <- floor(ncol(x)/2) +1 ## Nyquist-freq
    fft <- t(stats::mvfft(t(x)))[,1:n]
    colnames(fft) <- c("DC",as.character(1:(n-1)))
    fft
}

## asinh trafo: alternative to log
ash <- function(x) log(x+sqrt(x^2+1))
## log trafo handling zeros by adding 1
log_1 <- function(x) log(x+1)

## moving average
ma <- function(x,n=5){stats::filter(x,rep(1/n,n), sides=2)}

#' process a time-series apt for the \code{\link{segmenTier}}
#' clustering wrapper \code{\link{clusterTimeseries}}
#' @param ts the timeseries as a matrix, where columns are the timepoints
#' and rows individual measurements (e.g., genomic positions for transcriptome
#' data)
#' @param smooth use stats' package \code{link[stats:smooth]{smooth}} to
#' smooth timeseries before processing
#' @param trafo prior data transformation, pass any function name, e.g.,
#' "log", or the package functions "asinh" (\code{ln(x + sqrt(x^2+1))}) or
#' "log_1" for (\code{ln(ts+1)}) 
#' @param low.thresh use this threshold to cut-off data, which will be
#' added to the absent/nuissance cluster later
#' @param use.fft use the Discrete Fourier Transform of the data
#' @param dft.range a vector of integers, giving the components of the
#' Discrete Fourier Transform to be used where 1 is the first component (DC)
#' corresponding to the mean value, and 2:n are the higher components
#' correspondong to 2:n full cycles in the data
#' @param use.snr use a scaled amplitude, where each component of the
#' Discrete Fourier Transform is divided by the mean of all other components,
#' which is similar to a signal-to-noise ratio (SNR)
#' @param dc.trafo data transformation for the first (DC) component of
#' the DFT, pass any function name, e.g., "log", or the package functions
#' "asinh" (\code{ln(x + sqrt(x^2+1))}) or "log_1" for (\code{ln(ts+1)}) 
#' @details This function exemplifies the processing of an oscillatory
#' transcriptome time-series data as used in the establishment of this
#' algorithm and the demo \code{segment_test}. As suggested by Machne & Murray
#' (PLoS ONE 2012) and Lehmann et al. (BMC Bioinformatics 2014) a Discrete
#' Fourier Transform of time-series data allows to cluster time-series by
#' their change pattern.
#' @references 
#'   @cite Machne2012 Lehmann2013
#'@export
processTimeseries <- function(ts, trafo="raw",
                              use.fft=TRUE, dc.trafo="raw", dft.range=2:7,
                              use.snr=TRUE, low.thresh=1, 
                              smooth=FALSE, keep.zeros=FALSE) {

    ## processing ID - this will be inherited to clusters
    ## and from there to segment ID and type
    processing <- trafo # ifelse(trafo=="identity","raw",trafo)
    if ( use.fft )
      processing <- paste(processing,
                          paste("dft",paste(range(dft.range),collapse="-"),
                                sep=""),
                          paste("dc",dc.trafo,sep=""),
                          ifelse(use.snr,"snr",""),
                          sep="_")
    
    tsd <-ts 
    tsd[is.na(tsd)] <- 0 # set NA to zero (will become nuissance cluster)
    zs <- apply(tsd,1,sum)==0 # remember all zeros

    ## smooth time-series
    ## TESTED, doesn't help to avoid fragmentation!
    if ( smooth ) {
        tsm <- apply(tsd,2,ma,150,FALSE)
        tsd[!zs,] <- tsm[!zs,]
    }
    
    ## transform raw data?
    ## NOTE that DFT and SNR below (use.fft) are an alternative
    ## data normalization procedure
    ## default: identity
    if ( trafo!="raw" )
        tsd <- get(trafo, mode="function")(tsd)
    
    ## get DFT
    if ( use.fft ) {
        fft <- get.fft(tsd)
        if ( !keep.zeros ) 
            fft[zs,] <- NA
        
        ## sequence length
        N <- nrow(fft)

        ## amplitude-scaling (~SNR), see Machne&Murray 2012
        if ( use.snr ) {
          amp <- abs(fft)
          snr <- fft
          for ( a in 2:ncol(fft) )
            snr[,a] <- fft[,a]/apply(amp[,-c(1,a)],1,mean)
          fft <- snr
        }
        if ( 1 %in% dft.range & dc.trafo!="raw" )
            fft[,1] <- get(dc.trafo,mode="function")(fft[,1]) #
        
        ## get low expression filter!
        tot <- Re(fft[,1]) # NOTE: DC component = rowSums(tsd)
        low <- tot < low.thresh

        ## filter selected components
        dat <- fft[,dft.range]
         
        ## get Real and Imaginary pars
        dat <- cbind(Re(dat),Im(dat))

    }else {

        dat <- tsd
        if ( !keep.zeros ) 
            dat[zs,] <- NA
        
        ## sequence length
        N <- nrow(dat)

        ## get low expression filter
        tot <- rowSums(dat,na.rm=TRUE)
        low <- rep(FALSE, nrow(dat))
        low <- tot < low.thresh
    }
    ## for plots
    #tsd[zs,] <- NA

    ## store which are NA and set to 0
    na.rows <- rowSums(is.na(dat))==ncol(dat)
    ##dat[is.na(dat)] <- 0 ## shouldn't happen?

    ## remove data rows: NA or low
    rm.vals <- na.rows | low

    ## time-series data set for clustering in clusterTimeseries
    tset <- list(dat=dat, ts=tsd, tot=tot,
                 zero.vals=zs, rm.vals=rm.vals, low.vals=low,
                 id=processing)
    class(tset) <- "timeseries"
    
    ## silent return
    tmp <- tset
}

#' simple wrapper for \code{\link[stats:kmeans]{kmeans}} clustering
#' of a time-series preprocessed by \code{\link{processTimeseries}}.
#' @param tset a timeseries processed by \code{\link{processTimeseries}}
#' @param K selected cluster numbers, the argument \code{centers}
#' of \code{\link[stats:kmeans]{kmeans}} 
#' @param iter.max the maximum number of iterations allowed in
#' \code{\link[stats:kmeans]{kmeans}}, see there
#' @param nstart initialization \code{\link[stats:kmeans]{kmeans}}:
#' "how many random sets should be chosen?", see there
#'@export
clusterTimeseries <- function(tset, K=16, iter.max=100000, nstart=100) {


    ## TODO:
    ## call recursively if multiple tsets are available
    ## and prepend tset names 

    ## get time series data
    id <- tset$id
    dat <- tset$dat
    rm.vals <- tset$rm.vals
    N <- nrow(dat)

    ## enought distinct values?
    ## TODO: issue segment based on low-filter
    ## OOR: postprocessing - extend segments into low levels?
    warn <- NULL
    if ( sum(!rm.vals)<10 ) 
        warn <- "not enough data"
    else if ( sum(!duplicated(dat[!rm.vals,]))<2 ) 
      warn <- "not enough data diversity"
    if ( !is.null(warn) ) {
        warning(warn)
        return(NULL)
    }
    
    ## CLUSTERING
    ## stored data
    clusters <- matrix(NA, nrow=nrow(dat), ncol=length(K))
    centers <- Pci <- Ccc <- rep(list(NA), length(K))
    
    usedk <- K
    for ( k in 1:length(K) ) {

        ## get cluster number K
        Kused <- min(c(K[k],sum(!duplicated(dat[!rm.vals,]))))
        cat(paste("clustering, N=",N,", K=",Kused, "\n"))
        
        ## cluster
        km <- stats::kmeans(dat[!rm.vals,], Kused, iter.max=iter.max,
                            nstart=nstart, algorithm="Hartigan-Wong")
        ## use alternative algo if this error occured
        if (km$ifault==4) {
            km <- stats::kmeans(dat[!rm.vals,], Kused,
                                iter.max=iter.max,nstart=nstart,
                                algorithm="MacQueen")
            warn <- "quick-transfer error in kmeans algorithm Hartigan-Wong, taking MacQueen"
            warning(warn)
        }
        
        ## prepare cluster sequence
        seq <- rep(0, N) ## init. to nuissance cluster 0
        seq[!rm.vals] <- km$cluster
        
        ## store which K was used, the clustering and cluster centers
        usedk[k] <- Kused
        clusters[,k] <- seq
        centers[[k]] <- km$centers
        
        ## C(c,c) - cluster X cluster cross-correlation matrix
        cr <- stats::cor(t(km$centers))

        Ccc[[k]] <- cr
        
        ## P(c,i) - position X cluster correlation
        ## NOTE: we could also calculate P(i,c) for "low" values
        ## (see processTimeseries)
        ##       but probably weaken the splits between adjacent ?
        P <- matrix(NA,nrow=N,ncol=Kused)
        P[!rm.vals,] <- clusterCor_c(dat[!rm.vals,], km$centers)

        Pci[[k]] <- P
    }
    ## count duplicate K
    if ( any(duplicated(K)) ) {
        sel <- paste(K,".1",sep="")
        cnt <- 2
        while( sum(duplicated(sel)) ) {
            sel[duplicated(sel)] <- sub("\\..*",paste(".",cnt,sep=""),
                                        sel[duplicated(sel)])
            cnt <- cnt+1
        }
        K <- sub("\\.1$","",sel)
    }
    ## name all results by K, will be used!
    colnames(clusters) <- names(centers) <-
        names(Pci) <- names(Ccc) <- paste(id,"_K",K,sep="")

    ## clustering data set for use in segmentCluster.batch 
    cset <- list(clusters=clusters, centers=centers, Pci=Pci, Ccc=Ccc,
                 K=K, usedk=usedk, warn=warn)
    class(cset) <- "clustering"

    ## silent return
    tmp <- cset
}

#' high-level wrapper for multiple runs of segmentation by
#' \code{\link{segmentClusters}} for multiple clusterings and
#' multiple segmentation parameters
#' @param cset a clustering set as returned by \code{\link{clusterTimeseries}}
#' @param csim.scale exponent to scale similarity matrices, must be odd
#' to maintain negative correlations!
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
#' @param fuse.threshold if adjacent segments are associated with clusters
#' the centers of which have a Pearson correlation \code{>fuse.threshold}
#' the field "fuse" will be set to 1 for the second segments (top-to-bottom
#' as reported)
#' @param ncpu number of available cores (CPUs), passed to
#' \code{\link[parallel:mclapply]{parallel::mclapply}} by
#' \code{\link{calculateScoringMatrix}}
#' @param short.name if TRUE (default) parameters that are not varied
#' will not be part of the segment type and ID
#' @param id if set, the default segment IDs are replaced by this
#' @param verb level of verbosity, 0: no output, 1: progress messages
#' @param save.mat store the scoring function matrix SM or the back-tracing
#' matrix K by adding "SM" and "SK" to the string vector save.mat; useful
#' in testing stage or for debugging or illustration of the algorithm;
#' see \code{\link{plotScoring}}
#' @details This is a high-level wrapper for \code{\link{segmentClusters}}
#' which allows segmentation over multiple clusterings as provided by the
#' function \code{\link{clusterTimeseries}} and over multiple segmentation
#' paramers. Specifially, parameters \code{csim.scale}, \code{score},
#' \code{M} and \code{Mn} can all be vectors.
#'@export
segmentCluster.batch <- function(cset, csim.scale=1, score="ccor",
                                 M=175, Mn=20, a=-2, nui=1,
                                 fuse.threshold=0.2,
                                 nextmax=TRUE, multi="max", multib="max",
                                 short.name=TRUE, id,
                                 ncpu=1, verb=1, save.mat="") {

        

    ## generate parameter combinations as matrix/list
    ## TODO: allow explict combinations via a list
    nk <- length(cset$K)
    nscore <- length(score)
    nscale <- length(csim.scale)
    nm <- length(M)
    nmn <- length(Mn)
    ## TODO
    ##nn <- length(nui)
    ##na <- length(a)

    ## parameter matrix
    params <- as.data.frame(matrix(NA,nrow=nk*nscore*nscale*nm,ncol=5))
    colnames(params) <- c("K","S","E","M","Mn") # clustering, scoring, exponent, M, Mn
    params[,1] <- rep(colnames(cset$clusters), each=nscore*nscale*nm*nmn)
    params[,2] <- rep(rep(score, nk), each=nscale*nm*nmn)
    params[,3] <- rep(rep(csim.scale, nk*nscore), each=nm*nmn)
    params[,4] <- rep(rep(M, nk*nscore*nscale), each=nmn)
    params[,5] <- rep(Mn, each= nk*nscore*nscale*nm)

    ## segment type name construction
    ## TODO: do this smarter? 
    typenm <- colnames(params)
    ## rm those with length==1 to keep short names
    if ( short.name ) {
        if ( nm==1 ) typenm <- typenm[-which(typenm=="M")]
        if ( nmn==1 ) typenm <- typenm[-which(typenm=="Mn")]
        if ( nscore==1 ) typenm <- typenm[-which(typenm=="S")]
        if ( nscale==1 ) typenm <- typenm[-which(typenm=="E")]
    }
    
    if ( verb>0 )
        cat(paste("CALCULATING",nrow(params),"SEGMENTATIONS\n"))

    allsegs <- NULL
    for ( i in 1:nrow(params) ) {

        sgtype <- paste(paste(typenm,params[i,typenm],sep=":"),collapse="_")
        K <- as.character(params[i,"K"])
        seq <- cset$clusters[,K]
        scr <- params[i,"S"]
        scale <- params[i,"E"]
        m <- params[i,"M"]
        mn <- params[i,"Mn"]

        if ( scr=="ccor" ) csim <- cset$Ccc[[K]]
        if ( scr=="icor" ) csim <- cset$Pci[[K]]
        if ( scr=="ccls" ) csim <- NULL

        if ( verb>0 )
            cat(paste("Calculating segment type",sgtype,";",
                      i,"of",nrow(params),"\n"))
        
        seg <-segmentClusters(seq=seq,csim=csim,csim.scale=scale,
                              score=scr,M=m,Mn=mn,nui=nui,a=a,
                              multi=multi,multib=multib,nextmax=nextmax,
                              save.mat="",verb=verb)

        ## tag adjacent segments from correlating clusters
        close <- fuseSegments(seg$segments, Ccc=cset$Ccc[[K]],
                              fuse.threshold=fuse.threshold)
        if ( sum(close)>0 & verb>0 )
          cat(paste("\t",sum(close), "segments could be fused\n"))

        ## collect results
        if ( nrow(seg$segments) > 0 ) {

            if ( missing(id) )
                id <- sgtype
            sgids <- paste(id, 1:nrow(seg$segments),sep="_")
            segs <- data.frame(ID=sgids,
                               type=rep(sgtype,length(sgids)),
                               seg$segments,
                               fuse=close)
            allsegs <- rbind(allsegs,segs)
            
            if ( verb>0 )
              cat(paste("\t", nrow(seg$segments), "added.\n"))
        } else if ( verb>0 )
          cat(paste("\tno segments.\n"))
    }
    if ( is.null(allsegs) & verb>0 )
      cat(paste("\tNO SEGMENTS FOUND, returning NULL\n"))

    allsegs
}

## tags adjacent segments if they are from correlating (>\code{fuse.thresh})
## clusters
fuseSegments <- function(segs, Ccc, fuse.threshold=.2) {

    if ( nrow(segs)==0 ) return(NULL)
    if ( nrow(segs)==1 ) return(FALSE)
    
    fuse <- rep(NA,nrow(segs))
    for ( j in 2:nrow(segs) ) 
        fuse[j] <- Ccc[segs[j,1],segs[j-1,1]]

    ## FUSE directly adjacent if clusters correlate?
    adj <- segs[2:nrow(segs),2] - segs[2:nrow(segs)-1,3] ==1
    close <- c(FALSE,adj) & fuse > fuse.threshold

    close
}
    

