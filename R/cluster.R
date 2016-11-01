### TIME-SERIES CLUSTERING PARAMETERS

## get Discrete Fourier Transformation
get.fft <- function(x) {
    n <- floor(ncol(x)/2) +1 ## Nyquist-freq
    fft <- t(mvfft(t(x)))[,1:n]
    colnames(fft) <- c("DC",as.character(1:(n-1)))
    fft
}


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

clusterTimeseries <- function(tset, selected=16, kiter=100000, nstart=100) {

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
        km <- kmeans(dat[!rm.vals,],K,iter.max=kiter,nstart=nstart)
        ## use alternative algo if this error occured
        if (km$ifault==4) {
            km <- kmeans(dat[!rm.vals,],K,iter.max=kiter,nstart=nstart,
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
        cr <- cor(t(km$centers))

        Ccc[[k]] <- cr
        
        ## P(c,i) - position X cluster correlation
        ## NOTE: we could also calculate P(i,c) for "low" values
        ## (see processTimeseries)
        ##       but probably weaken the splits between adjacent ?
        P <- matrix(NA,nrow=N,ncol=K)
        P[!rm.vals,] <- clusterCor_c(dat[!rm.vals,], km$centers)

        Pci[[k]] <- P
    }
    colnames(clusters) <- names(centers) <- paste("K",usedk,sep="")

    list(clusters=clusters, Pci=Pci, Ccc=Ccc,
         selected=selected, usedk=usedk, centers=centers)
}

segmentClusterset <- function(cset, csim.scale=1, scores="ccor",
                              M=175, Mn=20, a=2, nui=1,
                              nextmax=TRUE, multi="min",multib="min", 
                              ncpu=1, verb=1, save.mat="") {

        
    allsegs <- NULL

    ## TODO: generate param combinatoric matrix/list
    ## clusterings, ncol(cset$clusters),
    ## csim.scales
    ## scoring functions
    ## M
    ## Mn
    ## nui
    ## nextmax, multi, multib
    ## and loop through that instead

    for ( k in 1:ncol(cset$clusters) ) {

        seq <- cset$clusters[,k]
        selected <- cset$selected[k]

        ## segment type - clustering
        ktype <- paste("K",selected,"_", "k", k, sep="")

        segments <- NULL

        for ( score in scores  ) {

            for ( scale in csim.scale ) {
                
                if ( score=="ccor" ) csim <- cset$Ccc[[k]]
                if ( score=="icor" ) csim <- cset$Pci[[k]]
                if ( score=="xcor" ) csim <- cset$Ccc[[k]]

                ## segment type - scoring function
                sgtype <- paste(score, scale,sep="")
            
                seg <-segmentClusters(seq=seq,csim=csim,csim.scale=scale,
                                      score=score,M=M,Mn=Mn,nui=nui.cr,
                                      multi=multi,multib=multib,nextmax=nextmax,
                                      save.mat=save.mat,verb=2)
                
                ## store segments
                segs <- seg$segments
                
                ## FUSE correlating?
                ## CHECK HERE SINCE WE HAVE THE SIMILARITY MATRIX?
                ## TODO: move to extra function?
                if ( nrow(segs)>1 ) {
                    fuse <- rep(NA,nrow(segs))
                    for ( j in 2:nrow(segs) ) 
                        fuse[j] <- cset$Ccc[[k]][segs[j,1],segs[j-1,1]]
                    ## FUSE directly adjacent if clusters correlate?
                    adj <- segs[2:nrow(segs),2] - segs[2:nrow(segs)-1,3] ==1
                    close <- c(FALSE,adj) & fuse > fuse.thresh
                    if ( sum(close)>0 )
                        cat(paste(ktype, sgtype,
                                  "\t",sum(close), "segments could be fused\n"))
                }
                
                ## name segments
                if ( nrow(segs) > 0 ) {
                    if ( nrow(segs)==1 ) close <- FALSE 
                    segs <- cbind(segs,close) # bind fuse information
                    rownames(segs) <- paste(sgtype, 1:nrow(segs),sep="_")
                    segments <- rbind(segments,segs)
                }
            }
            if ( is.null(segments) ) {
                cat(paste("no segments for clustering", ktype,
                          "and scoring function",sgtype, "\n"))
                next
            }
        
        
            ## SEGMENT TYPE:
            ## K (cluster number), k (repeated  runs of same clustering),
            ## NOTE: naming by original K selected[k], the used[k] can be lower
            ## if not enough data was present
            ## TODO: instead of constructing a name
            ## just add all info as table cols here
            ## score, M, Mn, nui, (dyn.prog. settings), usedk, selectedk
        #sgtype <- paste("K",selected,"_", "k", k, sep="")
        
            ## storing results
            colnames(segments) <- c("cluster","start","end","fuse")
            segids <- paste(ktype, "_", rownames(segments),sep="")
            segtypes <- paste(ktype, "_", sgtype, sep="")
            segs <- data.frame(ID=segids,
                               type=segtypes,
                               CL=segments[,1],
                               segments[,2:3,drop=FALSE],
                               fuse=segments[,4,drop=FALSE])
            allsegs <- rbind(allsegs, segs)
        }
    }
    allsegs
}
    

