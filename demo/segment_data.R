
library("segmenTier")
library("Rcpp")
source("~/programs/segmenTier/R/plot.R")
source("~/programs/segmenTier/R/segment.R")
sourceCpp("~/programs/segmenTier/src/segment.cpp")
source("~/programs/segmenTier/R/cluster.R")
sourceCpp("~/programs/segmenTier/src/cluster.cpp")
load("~/programs/segmenTier/data/primseg436.rda")

## load time-series data
## contains tsd from primseg436 for
## a 7.6 kb genomic region
data(primseg436)

### TIME-SERIES CLUSTERING PARAMETERS
use.fft <- TRUE # cluster discrete Fourier transform of data?
use.snr <- TRUE # take SNR of DFT
trafo <- "raw" # "ash" # "log_1" #
dc.trafo <- "raw" # "ash" # NOTE: add component 1 (DC) to DFT range to use
low.thresh <- -Inf # -Inf/0 # minimal mean value (DC component of DFT if use.fft)
dft.range <- 2:7 # range of DFT to use for clustering
K <- c(16) # cluster number K
iter.max <- 100000 # max. iterations in kmeans
nstart <- 100   # number of initial configurations tested in kmeans

### SEGMENTATION PARAMETERS
vary <- setVarySettings(
  ## SCORING
  E=2, #c(1,3), # scale exponent of similarity matrices csim
  S=c("ccor","icor","ccls"), # SCORING FUNCTIONS
  M=100, #c(30,175), # scoring function minimal length penalty
  Mn=100, ## for nuissance clusters: allow smaller segments!?
  a=-2, 
  nui=2, #-/+ correlation of nuissance cluster with others and itself
    ## BACKTRACING
  nextmax=TRUE, # in back-tracing, search for the next non-decreasing S(i,c)
  multi="max", # "min" # handling of multiple max. score k in scoring
  multib="max" # "min" # multiple max. score clusters in back-tracing
  )


## PRE-PROCESS TIME SERIES FOR CLUSTERING
## take DFT and scale amplitudes, and
## select components of DFT
tset <- processTimeseries(ts=tsd, trafo=trafo, dc.trafo=dc.trafo,
                          use.fft=use.fft, dft.range=dft.range,
                          use.snr=use.snr, low.thresh=low.thresh)

## CLUSTER PRE-PROCESSED TIME SERIES
cset <- clusterTimeseries(tset, K=K, iter.max=iter.max, nstart=nstart)

## CALCULATE SEGMENTS FOR ALL CLUSTERINGS and
## FOR CHOSEN SEGMENTATION PARAMETERS
sset <- segmentCluster.batch(cset, varySettings=vary, 
                             verb=1, save.matrix=TRUE)
allsegs <- sset$segments
## PLOT RESULTS


## get time-series data

if ( !interactive() )
    png("segment_data.png",res=300,units="in", width=10,height=5)


nsg <- length(sset$ids)
nk <- length(cset$ids)
## number of plots
## two for time-series (total and heatmap)
## for each clustering: 1xclustering, 1x all segments, S/S1 for each segmentation type
nplots <- 2 + length(K) * (2 + 2*nsg/nk)
par(mfcol=c(nplots,1),mai=c(.3,1.5,.01,.01),mgp=c(1.3,.5,0),xaxs="i")

## TIME-SERIES PLOT UTILITY: plot both the total signal (optionally used
## for threshold) and a heatmap of the time-series
plot.tset(tset, plot=c("total","timeseries"))
## CLUSTERING PLOT UTILITY: 
## each clustering can have multiple segmentations; plot each
## NOTE that clusterings are sorted (by their similarity matrix `Ccc`)
## and colored along a color-wheel
for ( k in 1:ncol(cset$clusters) ) {
    plot.cset(cset, k)
    ## TODO:
    ## plot.sset(sset, c("segments")) - 1x
    ## plot.sset(sset, c("SV","SK")) -  2x each type for k 
    ## for each clustering, plot SV, SK and segments
    ##plot.sset(sset, c("segments"), types=tps) 
    columns <- c(name="ID", type="type", start="start", end="end",
                 color="color")
    ## filter allsegs by segments for the current clustering
    kid <- cset$ids[k]
    tps <- rownames(sset$settings)[sset$settings[,"K"] %in% kid]
    ypos <- segment.plotFeatures(allsegs[allsegs[,"type"]%in%tps,],
                                 coors=coors, typord=TRUE,cuttypes=TRUE,
                                 ylab="", names=FALSE,columns=columns,tcx=.5)
    axis(1)
    ## plot fuse tag
    fuse <- allsegs[allsegs[,"fuse"],]
    points(fuse[,"start"],ypos[fuse[,"type"]],col="black",pch=4, lwd=1,cex=1.5)
    ## plot S/S1
    plot.sset(sset, c("S","S1"), types=tps) 
    
}

if ( !interactive() )
    dev.off()
