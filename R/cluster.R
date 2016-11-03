### TIME-SERIES CLUSTERING PARAMETERS

## get Discrete Fourier Transformation
get.fft <- function(x) {
    n <- floor(ncol(x)/2) +1 ## Nyquist-freq
    fft <- t(stats::mvfft(t(x)))[,1:n]
    colnames(fft) <- c("DC",as.character(1:(n-1)))
    fft
}

## moving average
ma <- function(x,n=5){stats::filter(x,rep(1/n,n), sides=2)}

#' process a time-series apt for the \code{\link{segmenTier}}
#' clustering wrapper \code{\link{clusterTimeseries}}
#' @param ts the timeseries as a matrix, where columns are the timepoints
#' and rows individual measurements (e.g., genomic positions for transcriptome
#' data)
#' @param smooth use stats' package \code{link[stats:smooth]{smooth}} to
#' smooth timeseries before processing
#' @param trafo prior data transformation, either empty ("") or "log"
#' for (\code{ln(ts+1)}) or "ash" for "asinh x = log(x + sqrt(x^2+1))"
#' transformation which has less effects on extreme values
#' @param use.fft use the Discrete Fourier Transform of the data
#' @param dft.range a vector of integers, giving the components of the
#' Discrete Fourier Transform to be used where 1 is the first component (DC)
#' corresponding to the mean value, and 2:n are the higher components
#' correspondong to 2:n full cycles in the data
#' @param use.snr use a scaled amplitude, where each component of the
#' Discrete Fourier Transform is divided by the mean of all other components,
#' which is similar to a signal-to-noise ratio (SNR)
#' @param low.thresh use this threshold to cut-off data, which will be
#' added to the absent/nuissance cluster later
#' @details This function exemplifies the processing of an oscillatory
#' transcriptome time-series data as used in the establishment of this
#' algorithm and the demo \code{segment_test}. As suggested by Machne & Murray
#' (PLoS ONE 2012) and Lehmann et al. (BMC Bioinformatics 2014) a Discrete
#' Fourier Transform of time-series data allows to cluster time-series by
#' their change pattern.
#' @references 
#'   @cite Machne2012 Lehmann2013
#'@export
processTimeseries <- function(ts,
                              smooth=FALSE, trafo="",
                              use.fft=TRUE, dft.range=2:7,
                              use.snr=TRUE, low.thresh=1) {
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
    if ( trafo == "log" ) # ln
        tsd <- log(tsd+1)
    if ( trafo == "ash" ) # asinh x = log(x + sqrt(x^2+1))
        tsd <- log(tsd+sqrt(tsd^2+1))
    
    ## get DFT
    if ( use.fft ) {
        fft <- get.fft(tsd)
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
        
        ## get low expression filter!
        tot <- Re(fft[,1]) # NOTE: DC component = rowSums(tsd)
        low <- tot < low.thresh

        ## filter selected components
        dat <- fft[,dft.range]
         
        ## get Real and Imaginary pars
        dat <- cbind(Re(dat),Im(dat))

    }else {

        dat <- tsd
        dat[zs,] <- NA
        
        ## sequence length
        N <- nrow(dat)

        ## get low expression filter
        tot <- rowSums(dat,na.rm=TRUE)
        low <- rep(FALSE, nrow(dat))
        low <- tot < low.thresh
    }
    ## for plots
    tsd[zs,] <- NA

    ## store which are NA and set to 0
    na.rows <- rowSums(is.na(dat))==ncol(dat)
    ##dat[is.na(dat)] <- 0 ## shouldn't happen?

    ## remove data rows: NA or low
    rm.vals <- na.rows | low

    ## enought distinct values?
    ## TODO: issue segment based on low-filter
    ## OOR: postprocessing - extend segments into low levels?
    if ( sum(!rm.vals)<10 ) {
        cat(paste("not enough data\n"))
        next
    } 
    if ( sum(!duplicated(dat[!rm.vals,]))<2 ) {
        cat(paste("not enough diversity\n"))
        next
    }

    list(dat=dat, ts=tsd, tot=tot, zero.vals=zs, rm.vals=rm.vals, low.vals=low)
}

#' simple wrapper for \code{\link[stats:kmeans]{kmeans}} clustering
#' of a time-series preprocessed by \code{\link{processTimeseries}}.
#' @param tset a timeseries processed by \code{\link{processTimeseries}}
#' @param selected selected cluster numbers, the argument \code{centers}
#' of \code{\link[stats:kmeans]{kmeans}} 
#' @param iter.max the maximum number of iterations allowed in
#' \code{\link[stats:kmeans]{kmeans}}, see there
#' @param nstart initialization \code{\link[stats:kmeans]{kmeans}}:
#' "how many random sets should be chosen?", see there
#'@export
clusterTimeseries <- function(tset, selected=16, iter.max=100000, nstart=100) {

    dat <- tset$dat
    rm.vals <- tset$rm.vals
    N <- nrow(dat)

    ## CLUSTERING
    ## stored data
    clusters <- matrix(NA, nrow=nrow(dat), ncol=length(selected))
    centers <- Pci <- Ccc <- rep(list(NA), length(selected))
    
    usedk <- selected
    for ( k in 1:length(selected) ) {

        ## get cluster number K
        K <- min(c(selected[k],sum(!duplicated(dat[!rm.vals,]))))
        cat(paste("clustering, N=",N,", K=",K, "\n"))
        
        ## cluster
        km <- stats::kmeans(dat[!rm.vals,],K,iter.max=iter.max,nstart=nstart)
        ## use alternative algo if this error occured
        if (km$ifault==4) {
            km <- stats::kmeans(dat[!rm.vals,],K,iter.max=iter.max,nstart=nstart,
                                algorithm="MacQueen")
            cat(paste("quick-transfer error in kmeans, taking MacQueen\n"))
        }
        
        ## prepare cluster sequence
        seq <- rep(0, N) ## init. to nuissance cluster 0
        seq[!rm.vals] <- km$cluster
        
        ## store which K was used, the clustering and cluster centers
        usedk[k] <- K
        clusters[,k] <- seq
        centers[[k]] <- km$centers
        
        ## C(c,c) - cluster X cluster cross-correlation matrix
        cr <- stats::cor(t(km$centers))

        Ccc[[k]] <- cr
        
        ## P(c,i) - position X cluster correlation
        ## NOTE: we could also calculate P(i,c) for "low" values
        ## (see processTimeseries)
        ##       but probably weaken the splits between adjacent ?
        P <- matrix(NA,nrow=N,ncol=K)
        P[!rm.vals,] <- clusterCor_c(dat[!rm.vals,], km$centers)

        Pci[[k]] <- P
    }
    ## count duplicate K
    if ( any(duplicated(selected)) ) {
        sel <- paste(selected,".1",sep="")
        cnt <- 2
        while( sum(duplicated(sel)) ) {
            sel[duplicated(sel)] <- sub("\\..*",paste(".",cnt,sep=""),
                                        sel[duplicated(sel)])
            cnt <- cnt+1
        }
        selected <- sub("\\.1$","",sel)
    }
    ## name all results by K, will be used!
    colnames(clusters) <- names(centers) <-
        names(Pci) <- names(Ccc) <- paste("K",selected,sep="")

    list(clusters=clusters, centers=centers, Pci=Pci, Ccc=Ccc,
         selected=selected, usedk=usedk)
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
#' @param ncpu number of available cores (CPUs), passed to
#' \code{\link[parallel:mclapply]{parallel::mclapply}} by
#' \code{\link{calculateScoringMatrix}}
#' @param verb level of verbosity, 0: no output, 1: progress messages
#' @param save.mat store the scoring function matrix SM or the back-tracing
#' matrix K by adding "SM" and "SK" to the string vector save.mat; useful
#' in testing stage or for debugging or illustration of the algorithm;
#' see \code{\link{plotScoring}}
#' @param fuse.threshold if adjacent segments are associated with clusters
#' the centers of which have a Pearson correlation \code{>fuse.threshold}
#' the field "fuse" will be set to 1 for the second segments (top-to-bottom
#' as reported)
#' @details This is a high-level wrapper for \code{\link{segmentClusters}}
#' which allows segmentation over multiple clusterings as provided by the
#' function \code{\link{clusterTimeseries}} and over multiple segmentation
#' paramers. Specifially, parameters \code{csim.scale}, \code{score},
#' \code{M} and \code{Mn} can all be vectors.
#'@export
segmentCluster.batch <- function(cset, csim.scale=1, score="ccor",
                                 M=175, Mn=20, a=2, nui=1,
                                 fuse.threshold=0.2,
                                 nextmax=TRUE, multi="max", multib="max", 
                                 ncpu=1, verb=1, save.mat="") {

        
    allsegs <- NULL

    ## generate parameter combinations as matrix/list
    ##plst <- list(NA) ## TODO do this via  list ?
    nk <- length(cset$selected)
    nscore <- length(score)
    nscale <- length(csim.scale)
    nm <- length(M)
    nmn <- length(Mn)
    ## TODO
    ##nn <- length(nui)
    ##na <- length(a)

    ## parameter matrix
    params <- as.data.frame(matrix(NA,nrow=nk*nscore*nscale*nm,ncol=5))
    colnames(params) <- c("K","score","scale","M","Mn")
    params[,1] <- rep(colnames(cset$clusters), each=nscore*nscale*nm*nmn)
    params[,2] <- rep(rep(score, nk), each=nscale*nm*nmn)
    params[,3] <- rep(rep(csim.scale, nk*nscore), each=nm*nmn)
    params[,4] <- rep(rep(M, nk*nscore*nscale), each=nmn)
    params[,5] <- rep(Mn, each= nk*nscore*nscale*nm)

    ## segment type name construction
    ## TODO: do this smarter? 
    typenm <- colnames(params)
    ## rm those with length==1 to keep short names
    if ( nm==1 ) typenm <- typenm[-which(typenm=="M")]
    if ( nmn==1 ) typenm <- typenm[-which(typenm=="Mn")]
    if ( nscore==1 ) typenm <- typenm[-which(typenm=="score")]
    if ( nscale==1 ) typenm <- typenm[-which(typenm=="scale")]
   
    if ( verb>0 )
        cat(paste("CALCULATING",nrow(params),"SEGMENTATIONS\n"))

    allsegs <- NULL
    for ( i in 1:nrow(params) ) {

        sgtype <- paste(params[i,typenm],collapse="_")
        K <- as.character(params[i,"K"])
        seq <- cset$clusters[,K]
        scr <- params[i,"score"]
        scale <- params[i,"scale"]
        m <- params[i,"M"]
        mn <- params[i,"Mn"]

        if ( scr=="ccor" ) csim <- cset$Ccc[[K]]
        if ( scr=="icor" ) csim <- cset$Pci[[K]]
        if ( scr=="xcor" ) csim <- cset$Ccc[[K]]
        if ( scr=="cls" ) csim <- a

        if ( verb>0 )
            cat(paste("Calculating segment type",sgtype,";",
                      i,"of",nrow(params),"\n"))
        
        seg <-segmentClusters(seq=seq,csim=csim,csim.scale=scale,
                              score=scr,M=m,Mn=mn,nui=nui,
                              multi=multi,multib=multib,nextmax=nextmax,
                              save.mat="",verb=verb)

        ## tag adjacent segments from correlating clusters
        if ( nrow(seg$segments)>1 ) {
            close <- fuseSegments(seg$segments, Ccc=cset$Ccc[[K]],
                                  fuse.threshold=fuse.threshold)
            if ( sum(close)>0 & verb>0 )
                cat(paste("\t",sum(close), "segments could be fused\n"))
        }
        if ( nrow(seg$segments) > 0 ) {
            if ( nrow(seg$segments)==1 ) close <- FALSE 
            
            sgids <- paste(sgtype,1:nrow(seg$segments),sep="_")
            segs <- data.frame(ID=sgids,
                               type=rep(sgtype,length(sgids)),
                               seg$segments,
                               fuse=close)
            allsegs <- rbind(allsegs,segs)
        }
    }
    allsegs
}

## tags adjacent segments if they are from correlating (>\code{fuse.thresh})
## clusters
fuseSegments <- function(segs, Ccc, fuse.threshold=.2) {

    fuse <- rep(NA,nrow(segs))
    for ( j in 2:nrow(segs) ) 
        fuse[j] <- Ccc[segs[j,1],segs[j-1,1]]

    ## FUSE directly adjacent if clusters correlate?
    adj <- segs[2:nrow(segs),2] - segs[2:nrow(segs)-1,3] ==1
    close <- c(FALSE,adj) & fuse > fuse.threshold

    close
}
    

