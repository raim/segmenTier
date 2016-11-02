
library("segmenTier")

## load time-series data
## contains tsd from primseg436 for
## a 7.6 kb genomic region
data(primseg436)

### TIME-SERIES CLUSTERING PARAMETERS
use.fft <- TRUE # cluster discrete Fourier transform of data?
use.snr <- TRUE # take SNR of DFT
trafo <- "" # "ash" # "log" # 
low.thresh <- 1 # -Inf/0 # minimal mean value (DC component of DFT if use.fft)
dft.range <- 2:7 # range of DFT to cluster to use for clustering
selected <- c(8,12,16) # cluster number K
iter.max <- 100000 # max. iterations in kmeans
nstart <- 100   # number of initial configurations tested in kmeans

### SEGMENTATION PARAMETERS
nui.cr <- 1 #  -/+ correlation of nuissance cluster with others and itself

## SCORING FUNCTION
## which scoring functions to use
scores <- c("ccor","icor")
csim.scale <- c(1,3) # 3 # 1 ## scale exponent of similarity matrices csim
## scoring function minimal length penalty
M <- c(175) # 30 # for empty set?
Mn <- 15 # for nuissance clusters: allow smaller segments!

## TOTAL SCORE and BACK-TRACING
## NOTE: these have little effect on real data sets, except perhaps nextmax?
## handling of multiple max. score k in scoring
multi <- "max" # c("max","min")
## handling of multiple max. score clusters in back-tracing
multib <- "max" # c("max","skip","min")
## in back-tracing, search for the next non-decreasing S(i,c)
nextmax <-TRUE


## PRE-PROCESS TIME SERIES FOR CLUSTERING
## take DFT and scale amplitudes, and
## select components of DFT

tset <- processTimeseries(ts=tsd, smooth=FALSE, trafo=trafo,
                          use.fft=use.fft, dft.range=dft.range,
                          use.snr=use.snr, low.thresh=low.thresh)

## CLUSTER PRE-PROCESS TIME SERIES
cset <- clusterTimeseries(tset, selected=selected,iter.max=iter.max, nstart=nstart)

## CALCULATE SEGMENTS FOR ALL CLUSTERINGS and
## FOR CHOSEN SEGMENTATION PARAMETERS
allsegs <- segmentCluster.batch(cset, csim.scale=csim.scale, score=scores,
                                M=M, Mn=Mn, a=2, nui=nui.cr,
                                nextmax=nextmax, multi=multi,multib=multib, 
                                ncpu=1, verb=1, save.mat="")

## PLOT RESULTS

N <- nrow(tsd)
coors <- c(chr=1,start=1,end=N) # "chromosome" coordinates

## get time-series data
tsd <- tset$ts # incl. all trafos and zeros set to NA
tot <- tset$tot # total of the time-series
low <- tset$low.vals # cut-off for non-clustered positions

colors0 <- gray.colors(100) ## heatmap colors for the timeseries 
colors0[1] <- "#FFFFFF" ## replace minimal by white

par(mfcol=c(3,1),mai=c(.01,1.5,.01,.01),mgp=c(1.7,.5,0),xaxs="i")
plot(1:N,tot,log=ifelse(trafo!="","","y"),type="l",lwd=2,axes=FALSE)
points((1:N)[low],tot[low],col=2,cex=.5)
axis(2);
segment.plotHeat(tsd,coors=coors,chrS=0,colors=colors0, colnorm=TRUE)
columns <- c(name="ID", type="type", start="start", end="end")
segment.plotFeatures(allsegs, coors=coors,
                     typord=TRUE,cuttypes=TRUE,
                     ylab="", names=FALSE,columns=columns,tcx=.5)
