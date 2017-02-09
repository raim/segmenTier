
library("segmenTier")

## load time-series data
## contains tsd from primseg436 for
## a 7.6 kb genomic region
data(primseg436)

### TIME-SERIES CLUSTERING PARAMETERS
use.fft <- TRUE # cluster discrete Fourier transform of data?
use.snr <- TRUE # take SNR of DFT
trafo <- "raw" # "ash" # "log_1" #
dc.trafo <- "ash" # "ash" # NOTE: add component 1 (DC) to DFT range to use
low.thresh <- -Inf # -Inf/0 # minimal mean value (DC component of DFT if use.fft)
dft.range <- 2:7 # range of DFT to use for clustering
K <- c(16) # cluster number K
iter.max <- 100000 # max. iterations in kmeans
nstart <- 100   # number of initial configurations tested in kmeans
nui.thresh <- 0.6 # threshold of position-cluster correlation below which the position will be assigned to the nuissance cluster

### SEGMENTATION PARAMETERS
vary <- setVarySettings(
    E=1:3,    # scale exponent of similarity matrices csim
    S="icor", # SCORING FUNCTIONS
    M=100,    # scoring function minimal length penalty
    Mn=100,   # M for nuissance clusters
    nui=1:3   #-/+ correlation of nuissance cluster with others and itself
)

## Here we vary the parameters with strongest effect
## on segmentation of our budding yeast data set.
## Paramters with little effect on similarity-based scoring functions
## icor and ccor are discussed in demo(segment_test)

## PRE-PROCESS TIME SERIES FOR CLUSTERING
## take DFT and scale amplitudes, and
## select components of DFT
tset <- processTimeseries(ts=tsd, trafo=trafo, dc.trafo=dc.trafo,
                          use.fft=use.fft, dft.range=dft.range,
                          use.snr=use.snr, low.thresh=low.thresh)

## CLUSTER PRE-PROCESSED TIME SERIES
cset <- clusterTimeseries(tset, K=K, iter.max=iter.max, nstart=nstart,
                          nui.thresh=nui.thresh)

## CALCULATE SEGMENTS FOR ALL CLUSTERINGS and
## FOR CHOSEN SEGMENTATION PARAMETERS
sset <- segmentCluster.batch(cset, varySettings=vary, 
                             verb=1, save.matrix=TRUE)

## NOTE: segments are in sset$segments
head(sset$segments)


## PLOT RESULTS

## plot segmentation
if ( !interactive() )
    png("segment_data.png",res=300,units="in", width=10,height=5)

# plot.matrix=TRUE will additionally plot the internal scoring matrices
plotSegmentation(tset, cset, sset, plot.matrix=FALSE, cex=.5, lwd=2) 


if ( !interactive() )
    dev.off()
