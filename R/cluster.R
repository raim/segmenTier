### TIME-SERIES CLUSTERING PARAMETERS


### DATA TRANSFORMATION UTILS
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
ma <- function(x, n=5, circular=FALSE) {
    stats::filter(x,rep(1/n,n), sides=2, circular=circular)
}

### CHROMOSOME COORDINATE UTILS
## copied from genomeBrowser utils on 20161216
## https://gitlab.com/raim/genomeBrowser/blob/master/src/genomeBrowser_utils.R

#' convert chromosome coordinates to continuous index
#' @param features a table of chromosome features that must contain
#' the chromosome number (option \code{chrCol}), one or more chromosome
#' positions (option \code{cols}) and strand information (column
#'  \code{strandCol}).
#' @param chrS the chromosome index, indicating the start position
#' of each chromosome in the continuous index, derived from chromosome length
#' information
#' @param cols name of the columns giving coordinates that will be mapped
#' to continuous index
#' @param chrCol name of the column that gives the chromosome number
#' @param strandCol name of the column that gives forward/reverse strand
#' information
#' @param reverse a vector of possible reverse strand indicators
#' @export
coor2index <- function(features, chrS,
                       cols=c("start","end","coor"),
                       chrCol="chr", strandCol="strand",
                       reverse=c("-",-1)) {
  cols <- cols[cols%in%colnames(features)]
  for ( col in cols ) {
    features[,col] <- features[,col]+chrS[features[,chrCol]]
    minus <- features[,strandCol]%in%reverse
    features[minus,col] <- features[minus,col]+max(chrS)
  }
  features[,chrCol] <- 1
  features 
}

#' Simple version of \code{\link{index2coor}} for single values
#' @param pos the continuous index position that will be mapped to
#' chromosome coordinates
#' @param chrS the chromosome index, indicating the start position
#' of each chromosome in the continuous index, derived from chromosome length
#' information
#' @export
idx2coor <- function(pos, chrS) {
  coor <- cbind(chr=rep(1,length(pos)),coor=pos,strand=rep(NA,length(pos)))
  for ( i in 1:(length(chrS)-1) ) {
    ## frw strand
    current <- pos>chrS[i] & pos<=chrS[i+1]
    coor[current,"coor"] <- pos[current] - chrS[i]
    coor[current,"chr"] <- i
    coor[current,"strand"] <- 1
    ## rev strand
    current <- pos>(chrS[i]+max(chrS)) & pos<=(chrS[i+1]+max(chrS))
    coor[current] <- pos[current] - chrS[i] - max(chrS)
    coor[current,"chr"] <- i
    coor[current,"strand"] <- -1
  }
  coor  
}
#' get the chromosome from continuous index
#' @param idx index position for which chromosome information is reported
#' @param chrS the chromosome index, indicating the start position
#' of each chromosome in the continuous index, derived from chromosome length
#' information
#' @return returns the chromosome number
#' @export
idx2chr <- function(idx,chrS) {
    chr <- sapply(idx,function(x) which(chrS>x)[1]-1)
    chr[is.na(chr)] <- sapply(idx[is.na(chr)],function(x) # reverse strand
        which((chrS+max(chrS))>x)[1]-1)
    chr
}
#' get the strand from continuous index
#' @param idx index position for which strand information is reported
#' @param chrS the chromosome index, indicating the start position
#' of each chromosome in the continuous index, derived from chromosome length
#' information
#' @return returns the strand
#' @export
idx2str <- function(idx,chrS)
    ifelse(idx > max(chrS),-1,1)



#' convert continuous index to chromosome coordinates (reverse of
#' \code{\link{coor2index}})
#' @param features a table of chromosome features that must contain
#' the chromosome number (option \code{chrCol}), one or more chromosome
#' positions (option \code{coorCols}) and strand information (column
#'  \code{strandCol}).
#' @param chrS the chromosome index, indicating the start position
#' of each chromosome in the continuous index, derived from chromosome length
#' information
#' @param cols names of the columns giving coordinates that will be mapped
#' to continuous index
#' @param chrCol name of the column that gives the chromosome number
#' @param strandCol name of the column that gives forward/reverse strand
#' information
#' @param reverse a vector of possible reverse strand indicators
#' @export
index2coor <- function(features, chrS,
                       cols=c("start","end","coor"),
                       chrCol="chr", strandCol="strand",
                       strands=c(1,-1)) {

  cols <- cols[cols%in%colnames(features)]
  orig <- features[,cols,drop=FALSE]

  ## add chromosome and strand columns, if not present
  if ( !chrCol%in%colnames(features) )
    features <- cbind(chr=rep(NA,nrow(features)),features)
  if ( !strandCol%in%colnames(features) )
    features <- cbind(features,strand=rep(NA,nrow(features)))
  
  ## remap values back to original coordinates
  for ( i in 1:(length(chrS)-1) ) {
    ## forward strand
    current <- orig[,cols[1]]>chrS[i] & orig[,cols[1]]<=chrS[i+1]
    for ( col in cols )
      features[current,col] <- orig[current,col] - chrS[i]
    features[current,chrCol] <- i
    features[current,strandCol] <- strands[1]
    ## reverse strand
    current <- orig[,cols[1]]>(chrS[i]+max(chrS)) &
      orig[,cols[1]]<=(chrS[i+1]+max(chrS))
    for ( col in cols )
      features[current,col] <- orig[current,col] - chrS[i] - max(chrS)
    features[current,chrCol] <- i
    features[current,strandCol] <- strands[2]
  }
  features
}

## PLOT UTILITIES
#' Switch between plot devices
#' @param file.name file name without suffix (.png, etc)
#' @param type plot type: pdf, png or eps
#' @param width figure width in inches
#' @param height figure height in inches
#' @param res resolution in ppi (pixels per inch), only for 'png'
#' @export
plotdev <- function(file.name="test", type="png", width=5, height=5, res=100) {
  file.name <- paste(file.name, type, sep=".")
  if ( type == "png" )
    png(file.name, width=width, height=height, units="in", res=res)
  if ( type == "eps" )
    postscript(file.name, width=width, height=height, paper="special")
  if ( type == "pdf" )
    pdf(file.name, width=width, height=height)
}

#' pre-segmentation of whole-genome data into chunks that can
#' be handle by segmenTier.
#' @param ts the time-series of readcounts for the complete chromosome,
#' rows are chromosomal positions and columns are time-points; reverse
#' strand rows at the bottom of the matrix. Option \code{chrS} can be
#' used to handle chromosome ends and to optionally (\code{map2chrom})
#' map the resulting primary segment coordinates to chromosome coordinates.
#' @param chrS a chromosome index, indicating at wich positions
#' chromosomes start; this is required for handling chromosome ends
#' and forward and reverse strand values, but can be omitted.
#' @param avg the broad moving average of read-count presence
#' (number of time-points with >0 reads) for a first broad segmentation
#' @param minrd the minimal number of time-points with reads in the broad
#' moving average used as cutoff between segments
#' @param favg as \code{avg}, but a narrower moving average used in
#' end scanning that can result in fusing back segments w/o good separation
#' @param minds minimum distance between two segments (will be fused otherwise)
#' @param map2chrom if true, argument \code{chrS} is required to map
#' the segment coordinates to chromosomal coordinates
#' @param fig.path a directory path for plots of the segment end scanning;
#' no figures will be plotted if \code{fig.path} is not provided.
#' @param fig.type image type, "png" or "pdf"
#' @param seg.path a directory path where individual segments' data will
#' be written to as tab-delimited .csv files; no files will be written if
#' \code{seg.path} is not provided.
#' @export
presegment <- function(ts, chrS, avg=1000, favg=100, minrd=8, minds=250,
                       map2chrom=FALSE,
                       seg.path, fig.path, fig.type="png", verb=1) {

    if ( verb> 0 )
        cat(paste("Calculating total read-counts and moving averages.\n"))
    
    ## total time series
    numts <- rowSums(ts > 0) ## timepoints with reads
    
    ## moving averages of read-count presence 
    avgts <- ma(numts,n=avg,circular=TRUE) # long mov.avg
    avgfn <- ma(numts,n=favg,circular=TRUE) # short mov.avg

    ## main primary segment definition! 
    segs <- avgts > minrd # smoothed total read number is larger then threshold
    
    ## set chromosome ends to FALSE as well
    if ( !missing(chrS) )
        segs[c(chrS[2:length(chrS)],chrS+1)] <- FALSE

    ## areas below expression threshold
    empty <- which(!segs) 

    ## distance between empty areas
    emptycoor <- empty
    if ( !missing(chrS) ) # accounts for chromosome ends - TODO: could be faster
        emptycoor <- idx2coor(empty, chrS)[,"coor"]
    distn <- diff(emptycoor) 

    ## PRE-SEGMENTS: ALL WHERE DISTANCE BETWEEN EMPTY AREAS IS >1
    start <- empty[which(distn>1)]+1
    end <- start + distn[which(distn>1)]-2

    ## fuse close segments, < minds
    close <- start[2:length(end)] - end[2:length(end)-1] < minds

    if ( verb>0 )
        cat(paste("Fusing", sum(close), "segments with distance <",minds,"\n"))

    start <- start[c(TRUE,!close)]
    end <- end[c(!close,TRUE)]
    ## remove too small segments
    small <- end-start < minds
    start <- start[!small]
    end <- end[!small]
    primseg <- cbind(start,end)


    ## (2) expand ends in both directions until mov.avg. (n=10) of signal is 0
    ## TODO: analyze gradients and minima, and add to appropriate segments
    if ( verb>0 )
        cat(paste("Scanning borders.\n"))
    fused <- 0
    for ( sg in 2:nrow(primseg) ) {
        rng <- primseg[sg-1,2]:primseg[sg,1]
        ## scan from both sides, and if overlapping
        ## take minimum between
        k <- j <- min(rng)
        i <- max(rng)
        while ( k<(i-1) ) { # expand left segment to right
            if ( avgfn[k]==0 ) break
            k <- k+1
        }
        while ( i>(j+1) ) { # expand right segment to left
            if ( avgfn[i]==0 ) break
            i <- i-1
        }
        primseg[sg-1,2] <- k
        primseg[sg,1] <- i
        
        if ( i <= k )  {
            if (  verb > 1 )
                cat(paste("segment #",sg,i-k,
                          "overlap, will be fused with",sg-1,"\n"))
            fused <- fused +1
        }
        
        if ( missing(fig.path) ) next
    
        ## plot borders
        bord <- range(rng)
        rng<- max(rng[1]-1000,1):(rng[length(rng)]+1000)
        file.name <-file.path(fig.path,
                              ifelse(i<=k,
                                     paste("fused_",fused,sep=""),
                                     paste("border_",sg-1-fused,sep="")))
        plotdev(file.name,width=4,height=4,type=fig.type)
        plot(rng,numts[rng],type="l",ylim=c(-2,24),main=ifelse(i<=k,"fuse",""))
        lines(rng,avgts[rng],col=3)
        lines(rng,avgfn[rng],col=2);
        abline(h=8,col=3)
        if ( bord[1]==bord[2] ) {
            #cat(paste(sg, "borders equal\n"))
            abline(v=bord[1],col=3)
        }else
            arrows(x0=bord[1],x1=bord[2],y0=-2,y1=-2,col=3)
        if ( k==i ) {
            #cat(paste(sg, "k==i equal\n"))
            abline(v=k,col=3)
        } else
            arrows(x0=k,x1=i,y0=-1,y1=-1,col=2)
        dev.off()
    }

    ## (3) fuse primary segments with distance <=1, incl.
    ## those where ends where swapped in end extension
    start <- primseg[,1]
    end <- primseg[,2]
    close <- start[2:length(end)] - end[2:length(end)-1] < 2

    if ( verb>0 )
        cat(paste("Fusing", sum(close), "more segments\n"))
    
    start <- start[c(TRUE,!close)]
    end <- end[c(!close,TRUE)]
    
    ## (4) split chromosome ends!
    ## TODO: why is multiple chromosome end handling required?
    ## 
    ## get chromosomes of starts and ends via chrS
    if ( !missing(chrS) ) {
        
        if ( verb>0 )
            cat(paste("Splitting segments that still span chromosome ends.\n"))

        schr <- idx2chr(start,chrS) # forward strand
        echr <- idx2chr(end,chrS) # reverse strand
        splt <- which(echr!=schr) # which are spanning chromosome ends?
    
        ## split chromosome-spanning segments, and fuse with rest
        old <- cbind(start[-splt], end[-splt])
        
        str <- idx2str(start,chrS)[splt]
        ## (str==-1)*max(chrS) adds minus strand to end
        new<-rbind(cbind(start[splt],chrS[schr[splt]+1] + (str==-1)*max(chrS)),
                   cbind(chrS[schr[splt]+1]+1 + (str==-1)*max(chrS) ,end[splt]))
        seg <- rbind(old,new)
        start <- seg[,1] # re-assign start/end of segments
        end <- seg[,2]
        
        ## re-order
        end <- end[order(start)]
        start <- sort(start)

        ## remove too small segments again
        small <- end-start < minds
        start <- start[!small]
        end <- end[!small]
    }
    
    primseg <- cbind(start,end)  ## DONE - PRIMARY SEGMENTS v3 DEFINED!

    ## write out data for each segment, if requested
    if ( !missing(seg.path) ) {
        if ( verb>0 )
            cat(paste("Writing segment data to single files.\n"))
        writeSegments(data=ts, segments=primseg, name="primseg", path=seg.path)
    }
    
    ## map back to original chromosome coordinates
    if ( !missing(chrS) & map2chrom ) {
        #primseg <-cbind(start=primseg[,1], end=primseg[,2])
        primseg <- index2coor(primseg,chrS)
    }
    primseg
}

#' writing out data for segments to files, used for writing the primary
#' segments by \code{\link{presegment}} to single files that are
#' then further processed by segmenTier.
#' @param data a data matrix to which coordinates in \code{segments}
#' refer to
#' @param segments a matrix that must contain "start" and "end" (columns)
#' of segments; these coordinates will be extracted from \code{data}
#' and written to individual files, numbered by the row number in segments.
#' @param path optional output path where files will be written, if not supplied
#' files will end up in the current working directory (`getwd`)
#' @export
writeSegments <- function(data, segments, name, path) {

    if ( missing(name) ) name <- "segment"
    for ( i in 1:nrow(segments) ) {
        rng <- segments[i,"start"]:segments[i,"end"]
        if ( length(rng) < 100 )
            cat(paste("segment",i, length(rng),"\n"))
        tsd <- ts[rng,]
        id <- ifelse("ID"%in%colnames(segments), segments[i,"ID"], i)
        file.name <- paste(name, "_",id,".csv",sep="")
        if ( !missing(path) )
            file.name <- file.path(path, file.name)
        write.table(tsd,file.name,row.names=FALSE,sep="\t")
    }
 }

#' process a time-series apt for the \code{\link{segmenTier}}
#' clustering wrapper \code{\link{clusterTimeseries}}
#' @param ts the timeseries as a matrix, where columns are the timepoints
#' and rows individual measurements (e.g., genomic positions for transcriptome
#' data)
#' @param smooth.space integer, if set a moving average is calculated for
#' each time-point between adjacent data points using stats
#' package's \code{link[stats:smooth]{smooth}} with span \code{smooth.space}
#' @param smooth.time integer, if set the time-series will be smoothed
#' using stats package's \code{link[stats:filter]{filter}} to calculate a
#' moving average with span \code{smooth.time} and
#' \code{link[stats:smoothEnds]{smoothEnds}} to extrapolate smoothed first
#' and last time-points (again using span \code{smooth.time})
#' @param trafo prior data transformation, pass any function name, e.g.,
#' "log", or the package functions "ash" (\code{asinh = ln(x + sqrt(x^2+1))})
#' or "log_1" for (\code{ln(ts+1)}) 
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
#' "ash" (\code{asinh= ln(x + sqrt(x^2+1))}) or "log_1" for (\code{ln(ts+1)}) 
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
                              use.snr=TRUE, low.thresh=-Inf, 
                              smooth.space, smooth.time, keep.zeros=FALSE) {

    ## processing ID - this will be inherited to clusters
    ## and from there to segment ID and type
    processing <- paste("T:",trafo,sep="")
    if ( use.fft )
      processing <- paste(processing,"_",
                          paste("D:dft",paste(range(dft.range),collapse="-"),
                                sep=""),".",
                          paste("dc",dc.trafo,sep=""),".",
                          ifelse(use.snr,"snr","raw"),
                          sep="")
    
    tsd <-ts 
    tsd[is.na(tsd)] <- 0 # set NA to zero (will become nuissance cluster)
    zs <- apply(tsd,1,sum)==0 # remember all zeros

    ## smooth time-points between adjacent positions
    ## currently not used, doesn't help to avoid fragmentation!
    if ( !missing(smooth.space) ) {
        if ( smooth.space>1 ) {
            tsm <- apply(tsd[!zs,], 2 ,ma, smooth.space,FALSE)
            tsd[!zs,] <- tsm
        }
    }
    ## smooth time-series
    if ( !missing(smooth.time) ) {
        ## currently used only in clustering final segment time series
        if ( smooth.time>1 ) {
            tsm <- t(apply(tsd[!zs,], 1, ma, n=smooth.time, circular=FALSE))
            ends <- stats::smoothEnds(tsd[!zs,], k=smooth.time)
            tsd[!zs,] <- tsm
            tsd[!zs,c(1,ncol(tsd))] <- ends[,c(1,ncol(tsd))]
        }
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
        re <- Re(dat)
        im <- Im(dat)
        if ( 1 %in% dft.range ) # rm 0 DC component from Im
            im <- im[,-1]
        dat <- cbind(re,im)

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

## TODO: adapt to be used in segmentation as well, is fcls@mu
## equal/similar to kmeans' 'centers'? Are Ccc and Pci calculated correctly?
#' wrapper for \code{\link[flowClust]{flowClust}}, currently only used for
#' clustering of final segment time-series; it could in principle also
#' be used for segmentation, but that has not been tested.
#' @param tset processed time-series as provided by
#' \code{\link{processTimeseries}}
#' @param B max. num. of EM iterations 
#' @param tol tolerance for EM convergence
#' @param lambda intial Box-Cox trafo
#' @param nu initial Box-Cox trafo, Inf for pure Gaussian
#' @param nu.est 0: no, 1: non-specific, 2: cluster-specific estimation of nu
#' @param trans 0: no, 1: non-specific, 2: cluster-specific estim. of lambda
#' @param ... further parameter to \code{\link[flowClust]{flowClust}}
#' @export
flowclusterTimeseries <- function(tset, ncpu=1, K=10, B=500, tol=1e-5, lambda=1,
                                nu=4, nu.est=0, trans=1, ...) {

    require("flowClust")
    require("flowMerge")
    
    dat <- tset$dat
    rm.vals <- tset$rm.vals
    clsDat <- dat[!rm.vals,]

    if ( ncpu>1 )
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
   

    ## MERGE CLUSTERS, starting from best BIC by flowMerge
    best <- which(K==max.clb)
    if ( length(fcls) > 1 ) fc <- fcls[[best]]
    else fc <- fcls
    obj <- flowObj(fc, flowFrame(clsDat))
    mrg <- merge(obj)
    mrg.cl <- fitPiecewiseLinreg(mrg)
    obj <- mrg[[mrg.cl]]
    mcls <- rep(0, nrow(dat))
    mcls[!rm.vals] <- flowClust::Map(obj, rm.outliers=F)
    mrg.id <- paste(K[best],"m",mrg.cl,sep="")
    cluster.matrix <- cbind(cluster.matrix, mcls)
    colnames(cluster.matrix)[ncol(cluster.matrix)] <- mrg.id

    ## collect centers, Pci and Ccc corelation matrices (see clusterTimeseries)
    ## TODO: TEST FOR SEGMENTATION (currently only used for final
    ## segment time series)
    ## -> is `mu' really the same as centers and are Ccc and Pci
    ## correct? 
    all <- append(fcls,obj)
    centers <- Pci <- Ccc <- rep(list(NA),length(all))
    for ( i in 1:length(all) ) {
        if ( length(all) > 1 ) fc <- all[[i]]
        else fc <- all
        centers[[i]] <- fc@mu
        ## C(c,c) - cluster X cluster cross-correlation matrix
        cr <- stats::cor(t(fc@mu))

        Ccc[[i]] <- cr
        
        ## P(c,i) - position X cluster correlation
        P <- matrix(NA,nrow=nrow(dat),ncol=nrow(fc@mu))
        P[!rm.vals,] <- clusterCor_c(clsDat, fc@mu)

        Pci[[i]] <- P
    }
    
    ## clustering data set for use in segmentCluster.batch 
    fcset <- list(clusters=cluster.matrix,
                  centers=centers, Pci=Pci, Ccc=Ccc,
                  K=K, usedk=K, warn=NULL,
                  flowClust=fcls, flowMerge=obj, # flowClust/flowMerge results
                  max.clb=max.clb, max.cli=max.cli, merged=mrg.id,
                  bic=bic, icl=icl)
    class(fcset) <- "clustering" 

    if ( ncpu>1 )
        options(cores=1)
    ## silent return
    tmp <- fcset
}

#' wrapper for \code{\link[stats:kmeans]{kmeans}} clustering
#' of a time-series preprocessed by \code{\link{processTimeseries}}.
#' @param tset a timeseries processed by \code{\link{processTimeseries}}
#' @param K selected cluster numbers, the argument \code{centers}
#' of \code{\link[stats:kmeans]{kmeans}} 
#' @param iter.max the maximum number of iterations allowed in
#' \code{\link[stats:kmeans]{kmeans}}, see there
#' @param nstart initialization \code{\link[stats:kmeans]{kmeans}}:
#' "how many random sets should be chosen?", see there
#' @param nui.thresh threshold correlation of a data point to a cluster
#' center; if below the data point will be added to nuissance cluster 0
#'@export
clusterTimeseries <- function(tset, K=16, iter.max=100000, nstart=100, nui.thresh=-Inf) {


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
        P <- matrix(NA,nrow=N,ncol=Kused)
        P[!rm.vals,] <- clusterCor_c(dat[!rm.vals,], km$centers)

        Pci[[k]] <- P
    }

    ## re-assign by correlation threshold
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
        names(Pci) <- names(Ccc) <- paste(id,"_K",K,sep="")

    ## clustering data set for use in segmentCluster.batch 
    cset <- list(clusters=clusters, centers=centers, Pci=Pci, Ccc=Ccc,
                 K=K, usedk=usedk, warn=warn)
    class(cset) <- "clustering"

    ## silent return
    tmp <- cset
}

#' A high-level wrapper for multiple runs of segmentation by
#' \code{\link{segmentClusters}} for multiple clusterings and
#' multiple segmentation parameters. It additionally allows to
#' tag adjacent segments to be potentially fused due to similarity
#' of their clusters.
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
        close <- fuseTagSegments(seg$segments, Ccc=cset$Ccc[[K]],
                                 fuse.threshold=fuse.threshold)
        if ( sum(close)>0 & verb>0 )
          cat(paste("\t",sum(close), "segments could be fused\n"))

        ## collect results
        if ( nrow(seg$segments) > 0 ) {

            sgids <- paste(sgtype, 1:nrow(seg$segments),sep="_")
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
    else {
        ## OVERRIDE ID
        if ( !missing(id) ) 
            allsegs[,"ID"] <- paste(id, 1:nrow(allsegs), sep="_")
    }


    allsegs
}

## tags adjacent segments if they are from correlating (>\code{fuse.thresh})
## clusters
fuseTagSegments <- function(segs, Ccc, fuse.threshold=.2) {

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

