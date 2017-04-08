
debug <- FALSE
if ( debug ) {
    library("Rcpp")
    source("~/programs/segmenTier/R/plot.R")
    source("~/programs/segmenTier/R/cluster.R")
    source("~/programs/segmenTier/R/segment.R")
    sourceCpp("~/programs/segmenTier/src/segment.cpp")
    sourceCpp("~/programs/segmenTier/src/cluster.cpp")
    load("~/programs/segmenTier/data/primseg436.rda")
} else {
    library("segmenTier")


## load time-series data
## contains tsd from primseg436 for
## a 7.6 kb genomic region
    data(primseg436)
}

plot.pdf <- FALSE

## EXAMPLE DATASET FROM BUDDING YEAST
## NOTE Here we vary the parameters with strongest effect
## on segmentation of our budding yeast data set.
## Paramters with little effect on similarity-based scoring functions
## icor and ccor are discussed in demo(segment_test)


### TIME-SERIES PROCESSING PARAMETERS
## NOTE that the current pipe-line for batch processing
## allows only one configuration of time-series processing
## while all downstream steps (clustering, segmentation)
## can be varied over parameter ranges.
## NOTE that currently segmenTier is not tested for clustering
## of time-series directly (without Discrete Fourier Transform)
## but can be tested by setting use.fft to FALSE
trafo <- "raw"     # transformation of the raw data 
use.fft <- TRUE    # cluster discrete Fourier transform of data?
use.snr <- TRUE    # use DFT scaling (SNR is described as
                   # relative amplitude scaling in Machne&Murray, PLoS ONE 2012)
dft.range <- 1:7   # range of DFT to use for clustering
dc.trafo <- "ash"  # transformation of the first (DC) component of the DFT
                   # NOTE: add component 1 (DC) to DFT range to use
low.thresh <- -Inf # minimal total signal (DC component of DFT if use.fft)

### CLUSTERING PARAMETERS
K <- c(16)         # cluster number K; multiple allowed; specifically, note
                   # that k-means has a random effect at initialization
                   # and replicates of the same K can  yield different
                   # results for otherwise 
nui.thresh <- 0.6  # threshold of position-cluster correlation below which
                   # the position will be assigned to the nuissance cluster
## k-means initialization
iter.max <- 100000 # max. iterations in kmeans
nstart <- 100      # number of initial configurations tested in kmeans

### SEGMENTATION PARAMETERS
## segmenTier parameters are handled via the settings function,
## where all parameters can be passed as vectors.
vary <- setVarySettings(
    E=1:3,    # scale exponent of similarity matrices csim
    S="icor", # SCORING FUNCTIONS
    M=c(100,150), # scoring function minimal length penalty
    Mn=100,   # M for nuissance clusters
    nui=1:3   #-/+ correlation of nuissance cluster with others and itself
)

## PRE-PROCESS TIME SERIES FOR CLUSTERING
## take DFT and scale amplitudes, and
## select components of DFT
tset <- processTimeseries(ts=tsd, trafo=trafo, dc.trafo=dc.trafo,
                          use.fft=use.fft, dft.range=dft.range,
                          use.snr=use.snr, low.thresh=low.thresh)

## CLUSTER PRE-PROCESSED TIME SERIES
set.seed(15) # stable kmeans clustering
cset <- clusterTimeseries(tset, K=K, iter.max=iter.max, nstart=nstart,
                          nui.thresh=nui.thresh)

## CALCULATE SEGMENTS FOR ALL CLUSTERINGS and
## FOR CHOSEN SEGMENTATION PARAMETERS
sset <- segmentCluster.batch(cset, varySettings=vary,type.name=c("E","M","nui"), verb=1, save.matrix=FALSE)
## NOTE: segments are in sset$segments
head(sset$segments)

## PLOT RESULTS
## plot segmentation
if ( plot.pdf )
  plotdev("segment_data_exponents",res=300,width=10,height=5,type="pdf")

# plot.matrix=TRUE will additionally plot the internal scoring matrices
plotSegmentation(tset, cset, sset, plot.matrix=FALSE, cex=.5, lwd=2) 

if ( plot.pdf )
  dev.off()

## VARY SETTINGS

## vary M; E=nui=2
vary <- setVarySettings(
    E=2,    # scale exponent of similarity matrices csim
    S="icor", # SCORING FUNCTIONS
    M=c(100,150,200), # scoring function minimal length penalty
    Mn=100,   # M for nuissance clusters
    nui=2   #-/+ correlation of nuissance cluster with others and itself
)
varM <- segmentCluster.batch(cset, varySettings=vary,type.name=c("E","M","nui"), verb=1, save.matrix=FALSE)

## vary E; nui=2, M=150
vary <- setVarySettings(
    E=1:3,    # scale exponent of similarity matrices csim
    S="icor", # SCORING FUNCTIONS
    M=150,   # scoring function minimal length penalty
    Mn=100,   # M for nuissance clusters
    nui=2   #-/+ correlation of nuissance cluster with others and itself
)
varE <- segmentCluster.batch(cset, varySettings=vary,type.name=c("E","M","nui"), verb=1, save.matrix=FALSE)

## vary nui; E=2, M=150
vary <- setVarySettings(
    E=2,    # scale exponent of similarity matrices csim
    S="icor", # SCORING FUNCTIONS
    M=150,   # scoring function minimal length penalty
    Mn=100,   # M for nuissance clusters
    nui=1:3   #-/+ correlation of nuissance cluster with others and itself
)
varNui <- segmentCluster.batch(cset, varySettings=vary,type.name=c("E","M","nui"), verb=1, save.matrix=FALSE)



## worst at E=nui=3, M=75
## worst at E=nui=1, M=200


## BAD PARAMETER SETS:
## UNDER-FRAGMENTATION - red cluster in fig 2
vary <- setVarySettings(
    E=1,    # scale exponent of similarity matrices csim
    S="ccor", # SCORING FUNCTIONS
    M=200,   # scoring function minimal length penalty
    Mn=100,   # M for nuissance clusters
    nui=1   #-/+ correlation of nuissance cluster with others and itself
)
bad1 <- segmentCluster.batch(cset, varySettings=vary,type.name=c("E","M","nui"))
## OVER-FRAGMENATION - magenta
vary <- setVarySettings(
    E=3,    # scale exponent of similarity matrices csim
    S="ccor", # SCORING FUNCTIONS
    M=75,   # scoring function minimal length penalty
    Mn=100,   # M for nuissance clusters
    nui=3   #-/+ correlation of nuissance cluster with others and itself
)
bad2 <- segmentCluster.batch(cset, varySettings=vary,type.name=c("E","M","nui"))

## BEST FRAGMENATION - cyan
vary <- setVarySettings(
    E=2,    # scale exponent of similarity matrices csim
    S="icor", # SCORING FUNCTIONS
    M=150,   # scoring function minimal length penalty
    Mn=100,   # M for nuissance clusters
    nui=2   #-/+ correlation of nuissance cluster with others and itself
)
best <- segmentCluster.batch(cset, varySettings=vary,type.name=c("E","M","nui"))

## use layout to combine plots
if ( plot.pdf )
  plotdev("segment_data_examples",res=300,width=10,height=6,type="pdf")
layout(matrix(1:9,ncol=1),heights=c(.25,.5,.5,.3,.3,.3,.1,.1,.1))
par(mai=c(0.1,2,0.05,0.01),xaxs="i",yaxs="r")
par(cex=1) 
plot(tset)
par(cex=.6) 
plot(cset,axes=2,cex=.7)
par(cex=1.2) # increase axis labels
par(mai=c(0.01,2,0.01,0.01))
plot(varM,"segments",lwd=3)
plot(varE,"segments",lwd=3)
plot(varNui,"segments",lwd=3)
plot(bad1,"segments",lwd=3)
plot(best,"segments",lwd=3)
plot(bad2,"segments",lwd=3)
if ( plot.pdf )
  dev.off()


### MULTIPLE CLUSTERINGS
## Here we generate multiple clusterings, and segmentations
## (here with constant parameters) will be calculated for all of them.
## NOTE that nui.thresh acts as a noise filter for clustering, 
## based on a minimal position-cluster similarity in marix cset$Pci
## NOTE that random effects of k-means clustering could potentially
## be utilized to clean data.
## cluster
nui.thresh <- nui.thresh
vK <- rep(16,6) # c(16,16,16) #,20,20,20)
set.seed(10) # stable kmeans clustering
kset <- clusterTimeseries(tset, K=vK, iter.max=iter.max, nstart=nstart,
                          nui.thresh=nui.thresh)

## calculate segments
vary <- setVarySettings(
    E=3,    # scale exponent of similarity matrices csim
    S="icor", # SCORING FUNCTIONS
    M=100,   # scoring function minimal length penalty
    Mn=100,   # M for nuissance clusters
    nui=3   #-/+ correlation of nuissance cluster with others and itself
)
#vary$nui <- vary$E <- 3
vark <- segmentCluster.batch(kset, varySettings=vary, 
                             verb=1, save.matrix=TRUE)
## plot segmentation
if ( plot.pdf )
    plotdev("segment_data_clusterings",res=300,width=10,height=5,type="pdf")

# plot.matrix=TRUE will additionally plot the internal scoring matrices
#par(mai=c(0.1,2,0.05,0.01),xaxs="i",yaxs="r")
plotSegmentation(NULL, kset, vark, plot.matrix=FALSE, cex=.5, lwd=2, mai=c(0.1,2,0.05,0.01)) 
#plot(cset)

if ( plot.pdf )
    dev.off()


## CLUSTER-FREE - NOTE: this takes a lot of memory and time!
## each position x_i is treated as its own cluster
##N <- nrow(tset$dat)
##D <- sum(!tset$rm.vals)
##seq <- rep(0, N)
##seq[!tset$rm.vals] <- 1:D
##P <- matrix(NA,nrow=N,ncol=D)
#### position-position cross-correlation
##P[!tset$rm.vals,] <- clusterCor_c(tset$dat[!tset$rm.vals,],
##                                  tset$dat[!tset$rm.vals,])
##cset <- list()
##cset$clusters <- matrix(seq,ncol=1)
##cset$Pci <- list(P)
##colnames(cset$clusters) <- names(cset$Pci) <- "all"
##
##vary$S <- "icor" # only 'icor' is possible for this approach!
##vary$M <- 150
##vary$nui <- vary$E <- 3
##sset <- segmentCluster.batch(cset, varySettings=vary, 
##                             verb=1, save.matrix=TRUE)
