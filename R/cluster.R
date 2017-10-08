

### DATA TRANSFORMATION UTILS
## get Discrete Fourier Transformation
get.fft <- function(x) {
    n <- floor(ncol(x)/2) +1 ## Nyquist-freq
    fft <- t(stats::mvfft(t(x)))[,1:n,drop=FALSE]
    if ( n==1 )
        colnames(fft) <- "DC"
    else
        colnames(fft) <- c("DC",as.character(1:(n-1)))
    fft
}
## fourier permutation
do.perm <- function(x, fft=NULL, perm, verb=0) {

    N <- ncol(x)
    if ( is.null(fft) ) fft <- get.fft(x)
    xam <- abs(fft)/N
    pvl <- matrix(0,nrow=nrow(fft), ncol=ncol(fft))
    dimnames(pvl) <- dimnames(fft)
    ## TODO: use apply and parallel!
    for ( i in 1:perm ) {
        if ( verb>0 & i%%round(perm/10)==0 )
          cat(paste(round(i/perm,2)*100,"%, "))
        ## randomize columns and get fourier
        rft <- get.fft(x[,sample(1:ncol(x))])
        ram <- abs(rft)/N
        pvl <- pvl + as.numeric(ram >= xam)
    }
    if ( verb ) cat("\n")
    Re(pvl/perm)
}

#' asinh data transformation
#'
#' The asinh trafo, an alternative to log transformation that has less
#' (compressing) effects on the extreme values (low and high values),
#' and naturally handles negative numbers and 0.
#' @param x a numeric vector
#' @export
ash <- function(x) log(x+sqrt(x^2+1))

#' log trafo handling zeros by adding 1
#' @param x a numeric vector
#' @export
log_1 <- function(x) log(x+1)

## moving average
ma <- function(x, n=5, circular=FALSE) {
    stats::filter(x,rep(1/n,n), sides=2, circular=circular)
}

# calculate 95% confidence intervals for the given
# data vector using a t-distribution
ci95 <- function(data,na.rm=FALSE) {
    if ( na.rm ) data <- data[!is.na(data)]
    n <- length(data)
    if ( n<2 ) return(NA)
    error <- qt(0.975, df=n-1) * sd(data)/sqrt(n)
    return(error)
}

## cluster/segment colors; function derived from scale_colour_hue in ggplot2
## TODO: describe better
color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}


#' Process a time-series for \code{\link{segmenTier}}.
#' 
#' Performs the requested data transformations, including a Discrete
#' Fourier Transform (DFT) of the time-series as direct input for 
#' the clustering wrapper \code{\link{clusterTimeseries}}.
#' @param ts the timeseries as a matrix, where columns are the timepoints
#' and rows individual measurements (e.g., genomic positions for transcriptome
#' data)
#' @param trafo prior data transformation, pass any function name, e.g.,
#' "log", or the package functions "ash" (\code{asinh = ln(x + sqrt(x^2+1))})
#' or "log_1" for (\code{ln(ts+1)}) 
#' @param low.thresh use this threshold to cut-off data, which will be
#' added to the absent/nuissance cluster later
#' @param perm number of permutations of the data set, to obtain
#" p-values for the oscillation
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
#' "ash" (\code{asinh= ln(x + sqrt(x^2+1))}) or "log_1" for (\code{ln(ts+1)}).
#' @param smooth.space integer, if set a moving average is calculated for
#' each time-point between adjacent data points using stats
#' package's \code{link[stats:smooth]{smooth}} with span \code{smooth.space}
#' @param smooth.time integer, if set the time-series will be smoothed
#' using stats package's \code{link[stats:filter]{filter}} to calculate a
#' moving average with span \code{smooth.time} and
#' \code{link[stats:smoothEnds]{smoothEnds}} to extrapolate smoothed first
#' and last time-points (again using span \code{smooth.time})
#' @param verb level of verbosity, 0: no output, 1: progress messages
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
                              use.fft=TRUE, dc.trafo="raw", dft.range,
                              perm=0, use.snr=TRUE, low.thresh=-Inf, 
                              smooth.space=1, smooth.time=1, verb=0) {

    if ( typeof(ts)=="list" )
        ts <- as.matrix(ts) # smoothEnds causes problems for data.frames!?
     
    
    tsd <-ts 
    tsd[is.na(tsd)] <- 0 # set NA to zero (will become nuissance cluster)
    zs <- apply(tsd,1,sum)==0 # remember all zeros

    ## smooth time-points between adjacent positions
    ## currently not used, doesn't help to avoid fragmentation!
    if ( smooth.space>1 ) {
        tsm <- apply(tsd[!zs,], 2 ,ma, smooth.space,FALSE)
        tsd[!zs,] <- tsm
    }
    ## smooth time-series
    ## currently used only in clustering final segment time series
    ## NOTE/TODO: smooth.time must be ODD for smoothEnds
    if ( smooth.time>1 ) {
        tsm <- t(apply(tsd[!zs,], 1, ma, n=smooth.time, circular=FALSE))
        ends <- stats::smoothEnds(tsd[!zs,], k=smooth.time)
        tsd[!zs,] <- tsm
        tsd[!zs,c(1,ncol(tsd))] <- ends[,c(1,ncol(tsd))]
    }

        
    ## transform raw data?
    ## NOTE that DFT and SNR below (use.fft) are an alternative
    ## data normalization procedure
    ## default: identity
    if ( trafo!="raw" )
        tsd <- get(trafo, mode="function")(tsd) # ash, log_1, etc
    
    ## get DFT
    fft <- pvl <- NULL
    if ( use.fft ) {

        ## get DFT
        tmp <- get.fft(tsd[!zs,])
        fft <- matrix(NA, ncol=ncol(tmp), nrow=nrow(tsd))
        colnames(fft) <- colnames(tmp)
        fft[!zs,] <- tmp

        ## do DFT on permuted time-series to obtain p-values
        ## TODO: include amplitude and DC scaling in permutation?
        if ( perm>0 ) { 
            tmp <- do.perm(tsd[!zs,],fft=fft[!zs,], perm, verb=verb)
            pvl <- matrix(NA, ncol=ncol(tmp), nrow=nrow(tsd))
            colnames(pvl) <- colnames(tmp)
            pvl[!zs,] <- tmp
        }
        
        ## amplitude-scaling (~SNR), see Machne&Murray 2012
        if ( use.snr ) {
            amp <- abs(fft)
          snr <- fft
          for ( a in 2:ncol(fft) )
            snr[,a] <- fft[,a]/apply(amp[,-c(1,a)],1,mean)
          fft <- snr
        }
        
        ## PREPARE DATA FOR CLUSTERING
        ## get low expression filter!
        tot <- Re(fft[,1]) # NOTE: DC component = rowSums(tsd)
        low <- tot < low.thresh

        ## DC scaling
        if ( dc.trafo!="raw" )
            fft[,1] <- get(dc.trafo,mode="function")(fft[,1]) # ash, log_1, etc

        ## filter selected components
        if ( missing(dft.range) )
            dft.range <- 1:ncol(fft)
        dat <- fft[,dft.range]
         
        ## get Real and Imaginary pars
        re <- Re(dat)
        colnames(re) <- paste("Re_",colnames(re),sep="")
        im <- Im(dat)
        colnames(im) <- paste("Im_",colnames(im),sep="")
        #if ( 1 %in% dft.range ) # rm 0 DC component from Im
        #    im <- im[,-1]
        ## filter 0 imaginary components: DC and Nyquist!
        im <- im[,apply(im,2,function(x) any(x!=0,na.rm=TRUE))]
        dat <- cbind(re,im)

    }else {

        dat <- tsd
        dat[zs,] <- NA # set zero-vals to NA

        ## get low expression filter
        tot <- rowSums(dat,na.rm=TRUE)
        low <- rep(FALSE, nrow(dat))
        low <- tot < low.thresh
    }

    ## store which are NA and set to 0
    na.rows <- rowSums(is.na(dat))==ncol(dat)
    ##dat[is.na(dat)] <- 0 ## shouldn't happen?

    ## remove data rows: NA or low
    rm.vals <- na.rows | low

    settings <- list(trafo=trafo, 
                     use.fft=use.fft,
                     dc.trafo=dc.trafo,
                     dft.range=dft.range,
                     perm=perm,
                     use.snr=use.snr,
                     low.thresh=low.thresh, 
                     smooth.space=smooth.space,
                     smooth.time=smooth.time)

    ## generate processing ID - this will be inherited to clusters
    ## and from there to segment ID and type
    processing <- paste("T:",trafo,sep="")
    if ( use.fft )
      processing <- paste(processing,"_",
                          paste("D:dft",paste(range(dft.range),collapse="-"),
                                sep=""),".",
                          paste("dc",dc.trafo,sep=""),".",
                          ifelse(use.snr,"snr","raw"),
                          sep="")

    ## time-series data set for clustering in clusterTimeseries
    tset <- list(dat=dat, ts=tsd, dft=fft, pvalues=pvl, tot=tot,
                 zero.vals=zs, rm.vals=rm.vals, low.vals=low,
                 settings=settings, id=processing)
    class(tset) <- "timeseries"
    
    ## silent return
    tmp <- tset
}

## TODO: adapt to be used in segmentation as well, is fcls@mu
## equal/similar to kmeans' 'centers'? Are Ccc and Pci calculated correctly?

#' Cluster a processed time-series with \pkg{flowClust} & \pkg{flowMerge}.
#' 
#' A wrapper for \pkg{flowClust}, clustering
#' a time-series object \code{tset} provided by \code{\link{processTimeseries}},
#' where specifically the DFT of a time-series and requested data
#' transformation were calculated. The clustering is performed
#' on the \code{tset$dat} matrix. This currently only used for
#' clustering of final segment time-series; it could in principle also
#' be used for segmentation, but that has not been tested.
#' @param tset processed time-series as provided by
#' \code{\link{processTimeseries}}
#' @param ncpu number of cores available for parallel mode of
#' \pkg{flowClust}
#' @param K the requested cluster numbers (vector of integers)
#' @param merge logical indicating whether cluster merging with
#' \pkg{flowMerge} should be attempted
#' @param B max. num. of EM iterations 
#' @param tol tolerance for EM convergence
#' @param lambda intial Box-Cox trafo
#' @param nu initial Box-Cox trafo, Inf for pure Gaussian
#' @param nu.est 0: no, 1: non-specific, 2: cluster-specific estimation of nu
#' @param trans 0: no, 1: non-specific, 2: cluster-specific estim. of lambda
#' @param ... further parameter to \code{flowClust}
#' @export
flowclusterTimeseries <- function(tset, ncpu=1, K=10, merge=FALSE,
                                  B=500, tol=1e-5, lambda=1,
                                  nu=4, nu.est=0, trans=1, ...) {

    if ( !requireNamespace("flowMerge", quietly = TRUE) )
      stop("`flowclusterTimeseries' requires the bioconductor package `flowMerge' (incl. its dependencies  `flowCore' and `flowClust'")
    
    dat <- tset$dat
    rm.vals <- tset$rm.vals
    clsDat <- dat[!rm.vals,]

    ##if ( ncpu>1 ) # TODO: how to avoid parallel mode?
    oldcpu <- unlist(options("cores"))
    options(cores=ncpu)

    fcls <- flowClust::flowClust(clsDat, K=K, B=B, tol=tol, lambda=lambda,
                                 nu=nu, nu.est=nu.est, trans=trans, ...)

    ## collect clusterings
    cluster.matrix <- matrix(0, nrow=nrow(dat), ncol=length(K))
    colnames(cluster.matrix) <- as.character(K)
    bic <- rep(NA, length(K))
    names(bic) <- as.character(K)
    icl <- bic
    for ( i in 1:length(fcls) ) {
      if ( length(fcls) > 1 ) fc <- fcls[[i]]
      else fc <- fcls
      cl.num <- as.character(fc@K)
      cluster <- flowClust::Map(fc,rm.outliers=F)
      cluster.matrix[!rm.vals, cl.num] <- cluster
      bic[cl.num] <- fc@BIC
      icl[cl.num] <- fc@ICL
    }
    ## max BIC and ICL
    max.bic <- max(bic, na.rm=T)
    max.clb <- K[which(bic==max.bic)]
    max.icl <- max(icl, na.rm=T)
    max.cli <- K[which(icl==max.icl)]
    ## best K selection
    ## use K with max BIC
    selected <- max.clb

    ## MERGE CLUSTERS, starting from best BIC by flowMerge
    mrg.cl <- mrg.id <- obj <- NULL
    if ( merge ) {
        best <- which(K==max.clb)
        if ( length(fcls) > 1 ) fc <- fcls[[best]]
        else fc <- fcls
        mcls <- rep(0, nrow(dat))
        mrg.id <- mrg.cl <- "NA"
        obj <- try(flowMerge::flowObj(fc, flowCore::flowFrame(clsDat)))
        if ( class(obj)!="try-error" ) {
            mrg <- try(flowMerge::merge(obj))
            if ( class(mrg)!="try-error" ) {
                mrg.cl <- flowMerge::fitPiecewiseLinreg(mrg)
                obj <- mrg[[mrg.cl]]
                mcls[!rm.vals] <- flowClust::Map(obj, rm.outliers=F)
                mrg.id <- paste(K[best],"m",mrg.cl,sep="")
            }
        }
        cluster.matrix <- cbind(cluster.matrix, mcls)
        colnames(cluster.matrix)[ncol(cluster.matrix)] <- mrg.id
        all <- append(fcls,obj)
    } else {
        all <- fcls
    }
    ## collect centers, Pci and Ccc corelation matrices (see clusterTimeseries)
    ## TODO: TEST FOR SEGMENTATION (currently only used for final
    ## segment time series)
    ## -> is `mu' really the same as centers and are Ccc and Pci
    ## correct? 
    centers <- Pci <- Ccc <- rep(list(NA),length(all))
    for ( i in 1:length(all) ) {
        if ( length(all) > 1 ) fc <- all[[i]]
        else fc <- all

        ## get cluster centers!
        x <- fc@mu
        rownames(x) <- 1:nrow(x)
        colnames(x) <- colnames(clsDat)
        
        centers[[i]] <- x
        ## C(c,c) - cluster X cluster cross-correlation matrix
        cr <- stats::cor(t(x))

        Ccc[[i]] <- cr
        
        ## P(c,i) - position X cluster correlation
        P <- matrix(NA,nrow=nrow(dat),ncol=nrow(fc@mu))
        P[!rm.vals,] <- clusterCor_c(clsDat, fc@mu)

        Pci[[i]] <- P
    }
    colnames(cluster.matrix) <- names(centers) <-
        names(Pci) <- names(Ccc) <- paste("K:",K,sep="")

    
    ## clustering data set for use in segmentCluster.batch 
    fcset <- list(clusters=cluster.matrix,
                  centers=centers, Pci=Pci, Ccc=Ccc,
                  K=K, usedk=K, selected=selected, warn=NULL,
                  bic=bic, icl=icl,
                  ids=colnames(cluster.matrix),
                  tsid=rep(tset$id,ncol(cluster.matrix)),
                  flowClust=fcls, flowMerge=obj, # flowClust/flowMerge results
                  max.clb=max.clb, max.cli=max.cli,
                  merged.K=mrg.cl, merged=mrg.id)
    class(fcset) <- "clustering" 

    ## add cluster colors
    fcset <- colorClusters(fcset)

    if ( ncpu>1 )
      options(cores=oldcpu)
    
    ## silent return
    tmp <- fcset
}

#' AIC/BIC for kmeans
#'
#' provides a log-likelihood method for \code{kmeans} results,
#' after Neal Fultz at \url{https://stackoverflow.com/a/33202188},
#' also available via the 
#' \href{https://rdrr.io/github/nfultz/stackoverflow/src/R/logLik_kmeans.R}{stackoverflow package}
#' @param object a \code{kmeans} object
#' @param ... unused
#' @export
logLik.kmeans <- function(object, ...)
    structure(object$tot.withinss,
              df = nrow(object$centers)*ncol(object$centers),
              nobs = length(object$cluster),
              class = 'logLik')


#' Cluster a processed time-series with k-means.
#' 
#' A wrapper for \code{\link[stats:kmeans]{kmeans}}, clustering
#' a time-series object \code{tset} provided by \code{\link{processTimeseries}},
#' where specifically the DFT of a time-series and requested data
#' transformation were calculated. The clustering is performed
#' on the \code{tset$dat} matrix.
#' @param tset a timeseries processed by \code{\link{processTimeseries}}
#' @param K selected cluster numbers, the argument \code{centers}
#' of \code{\link[stats:kmeans]{kmeans}} 
#' @param iter.max the maximum number of iterations allowed in
#' \code{\link[stats:kmeans]{kmeans}}, see there
#' @param nstart initialization \code{\link[stats:kmeans]{kmeans}}:
#' "how many random sets should be chosen?", see there
#' @param nui.thresh threshold correlation of a data point to a cluster
#' center; if below the data point will be added to nuissance cluster 0
#' @param verb level of verbosity, 0: no output, 1: progress messages
#'@export
clusterTimeseries <- function(tset, K=16, iter.max=100000, nstart=100,
                              nui.thresh=-Inf, verb=1) {


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

    ## BIC/AIC for kmeans
    bic <- rep(NA, length(K))
    names(bic) <- as.character(K)
    aic <- bic
    
    if ( verb>0 ) {
        cat(paste("Timeseries N\t",N,"\n",sep=""))
        cat(paste("Used datapoints\t",sum(!rm.vals),"\n",sep=""))
    }
    
    usedk <- K
    for ( k in 1:length(K) ) {
        
        ## get cluster number K
        Kused <- min(c(K[k],sum(!duplicated(dat[!rm.vals,]))))
        
        if ( verb>0 )
            cat(paste("Clusters K\t", Kused, "\n",sep=""))
        
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
        P <- matrix(NA,nrow=N,ncol=Kused)
        P[!rm.vals,] <- clusterCor_c(dat[!rm.vals,], km$centers)

        Pci[[k]] <- P

        ## calculate BIC/AIC
        bic[k] <- stats::BIC(km)
        aic[k] <- stats::AIC(km)
    }

    ## re-assign by correlation threshold
    ## NOTE: this only affects scoring function ccor
    for ( k in 1:ncol(clusters) ) {
        cls <- clusters[,k]
        for ( p in 1:nrow(Pci[[k]]) )
            if ( !any(Pci[[k]][p,] > nui.thresh, na.rm=TRUE) )
                cls[p] <- 0
        clusters[,k] <- cls
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
        names(Pci) <- names(Ccc) <- paste("K:",K,sep="") #paste(id,"_K:",K,sep="")

    ## max BIC and ICL
    max.bic <- max(bic, na.rm=T)
    max.clb <- K[which(bic==max.bic)[1]]
    max.aic <- max(aic, na.rm=T)
    max.cla <- K[which(aic==max.aic)[1]]
    ## best K selection
    ## use K with max BIC
    selected <- max.clb
    
    ## clustering data set for use in segmentCluster.batch 
    cset <- list(clusters=clusters, centers=centers, Pci=Pci, Ccc=Ccc,
                 K=K, usedk=usedk, selected=selected,
                 bic=bic, aic=aic, 
                 max.clb=max.clb, max.cla=max.cla,
                 warn=warn, ids=colnames(clusters),
                 tsid=rep(id,ncol(clusters)))
    class(cset) <- "clustering"

    ## add cluster colors
    cset <- colorClusters(cset)
    
    ## silent return
    tmp <- cset
}

#' Assign colors to clusters.
#' 
#' Takes a clustering set as returned by \code{\link{clusterTimeseries}} and
#' assigns colors to each cluster in each clustering along
#' the "hue" color wheel, as in \code{scale_colour_hue} in \code{ggplot2};
#' if \code{cset} contains a sorting (see
#' \code{clusterSort}), this sorting will be used to assing colors along
#' the color wheel, otherwise a sorting will be calculated first;
#' @param cset a clustering set as returned by \code{\link{clusterTimeseries}}
#'@export
colorClusters <- function(cset) {

    ## each column in the clustering matrix is one clustering
    cset$colors <- rep(list(NA), ncol(cset$clusters))

    ## sort if no sorting is present
    ## this requires the cluster-cluster similarity matrix Ccc
    if ( !"sorting" %in% names(cset) )
        cset <- sortClusters(cset, sort=TRUE)

    ## generate colors; use gray for nuissance
    for ( k in 1:ncol(cset$clusters) ) {
        cols <- c("#888888", color_hue(ncol(cset$Ccc[[k]])))
        names(cols) <- c("0", cset$sorting[[k]])
        cset$colors[[k]] <- cols
    }
    names(cset$colors) <- colnames(cset$clusters)
    cset
}

#' Sort clusters by similarity.
#' 
#' Takes a clustering set as returned by \code{\link{clusterTimeseries}} and
#' uses the cluster-cluster similarity matrix \code{Ccc} to sort
#' clusters by their similarity, starting with the first cluster `1'; the next
#' cluster is the first cluster (lowest cluster number) with the highest
#' similarity to cluster `1', and procedding from there. The final sorting is
#' added as \code{sorting} to \code{cset} and the annotated
#' \code{cset} is returned. This sorting is subsequently used to select
#' cluster colors. This simply allows for more informative plots
#' of the clustering underlying a segmentation but has no consequence
#' on segmentation itself.
#' @param cset a clustering set as returned by \code{\link{clusterTimeseries}}
#' @param sort if set to FALSE the clusters will be sorted merely numerically
#' @param verb level of verbosity, 0: no output, 1: progress messages
#'@export
sortClusters <- function(cset, sort=TRUE, verb=0) {

    ## each column in the clustering matrix is one clustering
    cset$sorting <- rep(list(NA), ncol(cset$clusters))

    ## merely generate numerical sorting if sort is FALS
    if ( !sort ) {
        sorting <- NULL
        for ( k in 1:ncol(cset$clusters) ) 
            sorting[[k]] <- sort(unique(cset$clusters[,k]))
        sorting <- lapply(sorting, function(x) x[x!=0])
        names(sorting) <- colnames(cset$clusters)
        cset$sorting <- sorting
        return(cset)
    }
    
    ## each clustering in \code{cset} comes with a cluster-cluster
    ## similarity matrix (used for scoring function \code{ccor});
    ## here we use it to get a rough sorting, simply starting with
    ## cluster '1'; the first cluster with the highest similarity
    ## to cluster '1' is taken as the next cluster
    for ( k in 1:ncol(cset$clusters) ) {
        Ccc <- cset$Ccc[[k]]
        #if ( is.null(colnames(Ccc)) )
        #  colnames(Ccc) <- rownames(Ccc) <- 1:ncol(Ccc)
        sorting <- colnames(Ccc)
        remaining <- sorting
        cl <- remaining[1]
        remaining <- remaining[remaining!=cl]
        cln.srt <- cl
        ## start at first cluster
        while ( length(remaining) > 1 ) {
            clcor <- Ccc[cl,remaining,drop=FALSE]
            ## get first cluster with highest correlation
            new <- colnames(clcor)[which.max(clcor)] 
            if ( verb>0 ) cat(paste("\t", cl, ">", new,
                                  round(max(clcor),2), "\n"))
            cl <- new
            remaining <- remaining[remaining!=cl]
            cln.srt <- c(cln.srt, new)
        }
        ## add last
        cln.srt <- c(cln.srt, remaining)
        cset$sorting[[k]] <- cln.srt
    }
    names(cset$sorting) <- colnames(cset$clusters)
    
    cset
}

#' Parameters for \code{\link{segmentCluster.batch}}.
#' 
#' Generates the parameter list (\code{varySettings}) for
#' \code{\link{segmentCluster.batch}}, using defaults
#' for all parameters not passed.
#' @param E exponent to scale similarity matrices, must be odd
#' to maintain negative correlations!
#' @param S the scoring function to be used: "ccor", "icor" or "cls"
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
#'@export
setVarySettings <- function(E=c(1,3),
                            S="ccor",
                            M=100,
                            Mn=100,
                            a=-2, nui=c(1,3),
                            nextmax=TRUE,
                            multi="max",
                            multib="max") {
    list(E=E, S=S, M=M, Mn=Mn, a=a, nui=nui, # scoring
         nextmax=nextmax, multi=multi, multib=multib) # backtracing
}

#' Batch wrapper for \code{\link{segmentClusters}}.
#' 
#' A high-level wrapper for multiple runs of segmentation by
#' \code{\link{segmentClusters}} for multiple clusterings and
#' multiple segmentation parameters. It additionally allows to
#' tag adjacent segments to be potentially fused due to similarity
#' of their clusters.
#' @param cset a clustering set as returned by \code{\link{clusterTimeseries}}
#' @param varySettings list of settings where each entry can be a vector;
#' the function will construct a matrix of all possible combinations of
#' parameter values in this list, call \code{\link{segmentClusters}} for
#' each, and report a matrix of segments where the segment `type' is
#' constructed from the varied parameters; see option \code{short.name}.
#' A varySettings list with all required (default) parameters can be
#' obtained via function \code{\link{setVarySettings}}.
#' @param fuse.threshold if adjacent segments are associated with clusters
#' the centers of which have a Pearson correlation \code{>fuse.threshold}
#' the field "fuse" will be set to 1 for the second segments (top-to-bottom
#' as reported)
#' @param rm.nui remove nuissance cluster segments from final results
#' @param type.name vector of strings selecting the parameters which will be
#' used as segment types. Note, that all parameters that are actually varied
#' will be automatically added (if missing). The list can include parameters
#' from time-series processing found in the "clustering" object \code{cset}
#' as \code{cset$tids}.
#' @param short.name default type name construction; if TRUE (default)
#' parameters that are not varied will not be part of the segment type and ID.
#' This argument has no effect if argument \code{type.name} is set.
#' @param id if set, the default segment IDs, constructed from numbered
#' segment types, are replaced by this
#' @param save.matrix store the total score matrix \code{S(i,c)} and the
#' backtracing matrix \code{K(i,c)}; useful in testing stage or for
#' debugging or illustration of the algorithm
#' TODO: save.matrix is currently not implemented, since batch function
#' returns a matrix only
#' @param verb level of verbosity, 0: no output, 1: progress messages
#' @details This is a high-level wrapper for \code{\link{segmentClusters}}
#' which allows segmentation over multiple clusterings as provided by the
#' function \code{\link{clusterTimeseries}} and over multiple segmentation
#' parameters. Each parameter in the list \code{varySettings} can be
#' a vector and ALL combinations of the passed parameter values will
#' be used for one run of \code{\link{segmentClusters}}.
#'@export
segmentCluster.batch <- function(cset, varySettings=setVarySettings(),
                                 fuse.threshold=0.2, rm.nui=TRUE, 
                                 type.name, short.name=TRUE, id,
                                 save.matrix=FALSE, verb=1) {

    ## TODO: allow defaults; getSettings to get full list!

    ## SETTING UP PARAMETER MATRIX
    ## 1) combine clusterings with segmentation parameters
    nk <- length(cset$K)
    vS <- append(list(K=colnames(cset$clusters)), varySettings)
    vL <- sapply(vS,length)
    names(vL) <- names(vS)
    rL <- c(1,vL)
    params <- as.data.frame(matrix(NA,ncol=length(vS),
                                 nrow=prod(sapply(vS,length))))
    colnames(params) <- names(vS)
    ## fill parameter matrix
    for ( j in 1:ncol(params) ) 
        params[,j] <- rep(rep(vS[[j]],prod(rL[1:(j)])),
                          each=prod(rL[(j+2):length(rL)]))

    ## TODO: add time-series processing info
    if ( !is.null(cset$tsid) ) {
        cllst <- strsplit(cset$tsid, "_")
        ## get class ids
        clid <- unique(unlist(lapply(cllst, function(x) sub(":.*","",x))))
        ## fill data.frame
        cltab <- data.frame(matrix(NA, nrow=length(cset$tsid),
                                   ncol=length(clid)))
        colnames(cltab) <- clid
        for ( i in 1:length(cllst) ) {
            tmp <- strsplit(cllst[[i]],":")
            class <- unlist(lapply(tmp, function(x) x[2]))
            names(class) <- unlist(lapply(tmp, function(x) x[1]))
            cltab[i, names(class)] <- class
        }
        ## copy each to all existing param values
        if ( any(rL>1) ){ # do only if necessary, because otherwise
                          # cltab will be converted to a vector!
            cltab <- apply(cltab, 2, function(x)
                           rep(rep(x,prod(rL[1])),
                               each=prod(rL[(1+2):length(rL)])))
        }
        params <- cbind(cltab, params)
        ## recalculate vL; number of types in each setting
        vL <- apply(params,2,function(x) length(unique(x)))
        names(vL) <- colnames(params)
    }

    ## TODO: store previous processing IDs explicitly, to be used in plots

    ## CONSTRUCT SEGMENT CLASS NAMES!
    if ( missing(type.name) ) { 
        type.name <- colnames(params)
        ## rm those with length==1 to keep short names
        ## UNLESS there is no variation
        if ( short.name )
          if ( sum(vL>1)>0 )
            type.name <- type.name[vL>1]
          else
            type.name <- "S" # DEFAULT ID: scoring function
    } else {
        ## add all varied parameters appear in type.name
        varied <- names(vL[vL>1])
        type.name <- unique(c(type.name, varied))
    }
    
    if ( verb>0 )
        cat(paste("SEGMENTATIONS\t",nrow(params),"\n",sep=""))

    allsegs <- sgtypes <- NULL
    ## internal matrices (S, K, S1)
    if ( save.matrix ) 
      SK <- rep(list(NA), nrow(params))
    else
      SK <- NULL
    ## colors & sorting
    seg.col <- rep(list(NA), nrow(params))
    seg.srt <- rep(list(NA), nrow(params))
    
    ## TODO: convert this loop to lapply and try parallel use!
    ## TODO: redirect messages to msgfile or store in results
    for ( i in 1:nrow(params) ) {

        ## construct segment type name  
        sgtype <- paste(paste(type.name,
                              unlist(lapply(params[i,type.name],as.character)),
                              sep=":"),
                        collapse="_")
        ## clustering K comes formatted already (K:id)
        sgtype <- sub("K:K:","K:", sgtype)
        sgtypes <- c(sgtypes, sgtype)

        ## clustering input
        K <- as.character(params[i,"K"])
        seq <- cset$clusters[,K]


        ## scoring params
        S <- as.character(params[i,"S"])
        E <-  as.numeric(params[i,"E"])
        M <-  as.numeric(params[i,"M"])
        Mn <- as.numeric(params[i,"Mn"])
        nui<- as.numeric(params[i,"nui"])
        a <-  as.numeric(params[i,"a"])

        ## back-tracing params
        multi   <- as.character(params[i,"multi"])
        multib  <- as.character(params[i,"multib"])
        nextmax <- as.logical(params[i,"nextmax"])

        if ( S=="ccor" ) csim <- cset$Ccc[[K]]
        if ( S=="icor" ) csim <- cset$Pci[[K]]
        if ( S=="ccls" ) csim <- NULL

        if ( verb>0 )
            cat(paste("SEGMENT TYPE\t",sgtype,
                      "\t", i,"of",nrow(params),"\n",sep=""))

        ## TODO: pass cset, inherit colors and set type ID there!
        seg <-segmentClusters(seq=seq,csim=csim,E=E,
                              S=S,M=M,Mn=Mn,nui=nui,a=a,
                              multi=multi,multib=multib,nextmax=nextmax,
                              save.matrix=save.matrix,rm.nui=rm.nui,verb=verb)

        ## retrieve algo-internal vectors and matrices
        ## S1: s(1,i,C) - scoring function from j=1 to i
        ## S: S(i,C) - score matrix
        ## K: K(i,C) - backtracing
        ## all can be plotted via dedicated functions
        if ( save.matrix ) 
            SK[[i]] <- seg$SK[[1]]

        ## pass on cluster coloring as segment coloring
        if ( "colors" %in% names(cset) )
          seg.col[[i]] <- cset$colors[[K]]
        if ( "sorting" %in% names(cset) )
          seg.srt[[i]] <- cset$sorting[[K]]

        ## tag adjacent segments from correlating clusters
        close <- fuseTagSegments(seg$segments, Ccc=cset$Ccc[[K]],
                                 fuse.threshold=fuse.threshold)
        if ( sum(close)>0 & verb>0 )
            cat(paste("Fused tags\t",sum(close), "\n",sep=""))

        ## collect results
        if ( nrow(seg$segments) > 0 ) {

            ## add colors as column
            ## TODO: redundant with cluster colors; but both are used
            colors <- seg$segments[,"CL"]
            if ( "colors" %in% names(cset) )
                colors <- cset$colors[[K]][as.character(seg$segments[,"CL"])]

            ##close <- rep(FALSE, nrow(seg$segments))
            sgids <- paste(sgtype, 1:nrow(seg$segments),sep="_")
            segs <- data.frame(ID=sgids,
                               type=rep(sgtype,length(sgids)),
                               seg$segments,
                               fuse=close,
                               color=colors)
            allsegs <- rbind(allsegs,segs)
            
        } 
        if ( verb>0 )
            cat(paste("Segments\t", nrow(seg$segments), "\n",sep=""))
    }
    if ( is.null(allsegs) & verb>0 )
        cat(paste("Total segments\t0\n",sep=""))
    else {
        cat(paste("Total segments\t",nrow(allsegs),"\n",sep=""))
        ## OVERRIDE ID
        if ( !missing(id) ) 
            allsegs[,"ID"] <- paste(id, 1:nrow(allsegs), sep="_")
    }

    ## name all stored data!
    if ( save.matrix )
      names(SK) <- sgtypes
    rownames(params) <- names(seg.col) <- names(seg.srt) <- sgtypes

    ## store length of sequence
    N <- nrow(cset$clusters)
    
    ## TODO: introduce and use classes for segment results
    sset <- list(segments=allsegs, N=N, colors=seg.col, sorting=seg.srt,
                 SK=SK, settings=params, ids=sgtypes)
    class(sset) <- "segments"
    return(sset)
}

## tags adjacent segments if they are from correlating (>\code{fuse.thresh})
## clusters
fuseTagSegments <- function(segs, Ccc, fuse.threshold=.2) {

    if ( nrow(segs)==0 ) return(NULL)
    if ( nrow(segs)==1 ) return(FALSE)
    
    fuse <- rep(-Inf,nrow(segs))
    for ( j in 2:nrow(segs) ) {
        previous <- segs[j-1,1]
        current <- segs[j,1]
        if ( !0 %in% c(previous,current) ) 
            fuse[j] <- Ccc[current,previous]
    }
            
    ## FUSE directly adjacent if clusters correlate?
    adj <- segs[2:nrow(segs),2] - segs[2:nrow(segs)-1,3] ==1
    close <- c(FALSE,adj) & fuse > fuse.threshold

    close
}

