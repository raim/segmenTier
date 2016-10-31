### TIME-SERIES CLUSTERING PARAMETERS

processTimeseries <- function(ts,
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
    centers <- Pci <- Ccc <- rep(list(NA), length(selected))
    
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

    list(clusters=clusters, Pci=Pci, Ccc=Ccc,
         selected=selected, usedk=usedk, centers=centers)
}

segmentClusterset <- function(cset, csim.scale=1, scores="ccor",
                              M=175, Mn=20, a=2, nui=1,
                              nextmax=TRUE, multi="min",multib="min", 
                              ncpu=1, verb=1, save.mat="") {

        
    allsegs <- NULL

    for ( k in 1:ncol(cset$clusters) ) {

        seq <- cset$clusters[,k]
        selected <- cset$selected[k]

        scrR <- rep(list(NA),length(scores))
        names(scrR) <- scores
        segments <- NULL

        for ( score in scores  ) {
            multS <- rep(list(NA),length(multi) +1)
            names(multS) <- c("SM",multi)
            multS[[multi]] <- rep(list(NA),length(multib)+1)
            names(multS[[multi]]) <- c("SK",multib)
            
            if ( score=="ccor" ) csim <- cset$Ccc[[k]]
            if ( score=="icor" ) csim <- cset$Pci[[k]]
            if ( score=="xcor" ) csim <- cset$Ccc[[k]]
            
            ## scale csim!
            ## NOTE: should be odd number to maintain neg. values!
                                        #csim <- csim^scale
            
            seg <- segmentClusters(seq=seq,csim=csim,csim.scale=csim.scale,
                                   score=score,M=M,Mn=Mn,nui=nui.cr,
                                   multi=multi, multib=multib,nextmax=nextmax,
                                   save.mat=c("SK"),verb=2)
            ##multS$SM <- seg$SM
            multS[[multi]]$SK <- seg$SK
            ## store
            multS[[multi]][[multib]] <- seg$segments
            scrR[[score]] <- multS
            
            ## store segments
            segs <- seg$segments
            
            ## fuse by data
            #segs <- fuseSegments(segs, seq, dat)
            
            ## FUSE correlating?
            if ( nrow(segs)>1 ) {
                fuse <- rep(NA,nrow(segs))
                for ( j in 2:nrow(segs) ) 
                    fuse[j] <- cr[segs[j,1],segs[j-1,1]]
                ## FUSE directly adjacent if clusters correlate?
                adj <- segs[2:nrow(segs),2] - segs[2:nrow(segs)-1,3] ==1
                close <- c(FALSE,adj) & fuse > fuse.thresh
                if ( sum(close)>0 )
                    cat(paste("\t",sum(close), "segments could be fused\n"))
            }
            if ( nrow(segs) > 0 ) {
                if ( nrow(segs)==1 ) close <- FALSE 
                segs <- cbind(segs,close)
                rownames(segs) <- paste(paste(score,scale,sep=""),
                                        1:nrow(segs),sep="_")
                segments <- rbind(segments,segs)
            }
        }
        if ( is.null(segments) ) {
            cat(paste("no segments\n"))
            next
        }
        
        ## SEGMENT TYPE:
        ## K (cluster number), k (repeated  runs of same clustering),
        ## NOTE: naming by original K selected[k], the used[k] can be lower
        ## if not enough data was present
        ## TODO: add M, nui, (dyn.prog. settings)
        ## TODO: add usedk
        sgtype <- paste(#trafo,ifelse(trafo!="","_",""),
                        #ifelse(use.snr,"snr_",""),
                        "K",str_pad(selected,2,pad="0"),"_", "k", k, sep="")
        
        ## storing results
        colnames(segments) <- c("cluster","start","end","fuse")
        segids <- paste(segid, "_", sgtype, "_", rownames(segments),sep="")
        segtypes <- paste(sgtype, "_",sub("_.*", "",rownames(segments)),sep="")
        #segcoors <- cbind(chr = rep(NA, nrow(segments)),
        #                  segments[, 2:3,drop=FALSE] + primseg[i,"start"] - 1,
        #                  strand=rep(NA,nrow(segments)),
        #                  fuse=segments[,4,drop=FALSE])
        #segcoors <- index2coor(segcoors,chrS)
        ## store only the final 
        ##centers[[k]] <- centers[[k]][sort(unique(segments[,1])),]
        segs <- data.frame(ID=segids,
                           type=segtypes,
                           CL=segments[,1],
                           segments[,2:3,drop=FALSE],
                           fuse=segments[,4,drop=FALSE])
                           #segcoors)
        allsegs <- rbind(allsegs, segs)
    }
    allsegs
}
    

