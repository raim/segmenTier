### TIME-SERIES CLUSTERING PARAMETERS

prepareTimeseries <- function(ts,
                              smooth=FALSE, trafo="",
                              use.fft=TRUE, dft.range=2:7,
                              use.snr=TRUE, low.thresh=1) {
    tsd <-ts 

##tsd <- as.matrix(tsd[,4:ncol(tsd)])
    tsd[is.na(tsd)] <- 0
    zs <- apply(tsd,1,sum)==0

    ## smooth? TESTED, doesn't help to avoid fragmentation
    smooth <- FALSE
    if ( smooth ) {
        tsm <- apply(tsd,2,ma,150,FALSE)
        tsd[!zs,] <- tsm[!zs,]
    }
    
    ## transform raw data?
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

        ## filter low expression
        tot <- Re(fft[,1]) # NOTE: DC component = rowSums(tsd)
        low <- tot < low.thresh
        ## plot(tot,col=low+1)

        ## calculate SNR-scaling!
        if ( use.snr ) {
          amp <- abs(fft)
          snr <- fft
          for ( a in 2:ncol(snr) )
            snr[,a] <- fft[,a]/apply(amp[,-c(1,a)],1,mean)
          fft <- snr
        }
        
        ## filter selected components
        dat <- fft[,dft.range]

         
        ## get Real and Imaginary pars
        dat <- cbind(Re(dat),Im(dat))

    }else {
        dat <- tsd
        dat[zs,] <- NA
        ## sequence length
        N <- nrow(dat)

        ## filter by low expression
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

    list(ts=tsd, tot=tot, rm.vals=rm.vals, low.vals=low)
}
clusterTimeseries <- function(tset, selected=16, kiter=100000, nstart=100) {

    dat <- tset$ts
    rm.vals <- tset$rm.vals

    N <- nrow(dat)

    ## CLUSTERING
    ## stored data
    clusters <- matrix(NA, nrow=nrow(dat), ncol=length(selected))
    centers <- Pci <- Cca <- rep(list(NA), length(selected))
    
    allsegs <- NULL

    usedk <- selected
    for ( k in 1:length(selected) ) {

        ## get cluster number K
        K <- min(c(selected[k],sum(!duplicated(dat[!rm.vals,]))))
        cat(paste("clustering:", segid, ", N=",N,", K=",K, "\n"))
        
        ## cluster
        ## TODO: cluster SNR instead!
        km <- kmeans(dat[!rm.vals,],K,iter.max=kiter,nstart=nstart)
        ## use alternative algo if this error occured
        if (km$ifault==4) {
            km <- kmeans(dat[!rm.vals,],K,iter.max=kiter,nstart=nstart,
                         algorithm="MacQueen")
            cat(paste(segid,"\tWARNING: quick-transfer error\n"))
        }
        
        ## prepare cluster sequence
        seq <- rep(0, N) ## init. to nuissance cluster 0
        seq[!rm.vals] <- km$cluster
        
        ## store which K was used, the clustering and cluster centers
        usedk[k] <- K
        clusters[,k] <- seq
        centers[[k]] <- km$centers
        
        ## cluster center cross-correlation matrix
        cr <- cor(t(km$centers))

        Ccc[[k]] <- cr
        
        ## P(c,i) - position X cluster center correlation
        ## NOTE: we could also calculate P(i,c) for "low" values
        ##       but probably weaken the splits between adjacent ?
        P <- matrix(NA,nrow=N,ncol=K)
        P[!rm.vals,] <- clusterCor_c(dat[!rm.vals,], km$centers)

        Pci[[k]] <- P
    }
    colnames(clusters) <- names(centers) <- paste("K",usedk,sep="")

    list(clusters=clusters, Pci=P, Ccc=cr,
         usedk=usedk, centers=centers)
}

segmentTimeseries <- function(nui.cr=1, scores=c("ccor","icor"), M=175, Mn=15,
                              multi="max",multib= "max", nextmax=TRUE) {}


    

