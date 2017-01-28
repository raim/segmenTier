
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
nui.cr <- 2 #  -/+ correlation of nuissance cluster with others and itself
a <- -2

## SCORING FUNCTION
## which scoring functions to use
scores <- c("ccor","icor","ccls") #"icor" #
csim.scale <- c(1,3) # 3 # 1 ## scale exponent of similarity matrices csim
## scoring function minimal length penalty
M <- c(30,175) # 30 # for empty set?
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
tset <- processTimeseries(ts=tsd, trafo=trafo, dc.trafo=dc.trafo,
                          use.fft=use.fft, dft.range=dft.range,
                          use.snr=use.snr, low.thresh=low.thresh)

## CLUSTER PRE-PROCESSED TIME SERIES
cset <- clusterTimeseries(tset, K=K, iter.max=iter.max, nstart=nstart)

## CALCULATE SEGMENTS FOR ALL CLUSTERINGS and
## FOR CHOSEN SEGMENTATION PARAMETERS
allsegs <- segmentCluster.batch(cset, csim.scale=csim.scale, score=scores,
                                M=M, Mn=Mn, a=a, nui=nui.cr,
                                nextmax=nextmax, multi=multi,multib=multib, 
                                ncpu=1, verb=1, save.mat=TRUE)

## PLOT RESULTS


## get time-series data
ts <- tset$ts # incl. all trafos and zeros set to NA
ts[tset$zero.vals,] <- NA
tot <- tset$tot # total of the time-series

N <- nrow(ts)
coors <- c(chr=1,start=1,end=N) # "chromosome" coordinates

colors0 <- rev(gray.colors(100)) ## heatmap colors for the timeseries 
colors0[1] <- "#FFFFFF" ## replace minimal by white

if ( !interactive() )
    png("segment_data.png",res=300,units="in", width=10,height=5)

par(mfcol=c(3,1),mai=c(.3,1.5,.01,.01),mgp=c(1.3,.5,0),xaxs="i")
plot(1:N,tot,log=ifelse(trafo!="","","y"),type="l",lwd=2,axes=FALSE,ylab=NA)
polygon(x=c(1,1,N,N),y=c(min(tot,na.rm=TRUE),rep(low.thresh,2),min(tot,na.rm=TRUE)),col="#00000055",border=NA)
abline(h=low.thresh,col="#000000BB")
lines(1:N,tot)
axis(2);
axis(1)
mtext("total signal", 2, 2)
segment.plotHeat(ts,coors=coors,chrS=0,colors=colors0, colnorm=TRUE)
axis(2,at=1:ncol(ts))
axis(1)
mtext("time points", 2, 2)
## TODO: plot clustering
columns <- c(name="ID", type="type", start="start", end="end")
ypos <- segment.plotFeatures(allsegs, coors=coors,
                             typord=TRUE,cuttypes=TRUE,
                             ylab="", names=FALSE,columns=columns,tcx=.5)
axis(1)
## plot fuse tag
fuse <- allsegs[allsegs[,"fuse"],]
points(fuse[,"start"],ypos[fuse[,"type"]],col="red",pch=3, lwd=2,cex=2)

if ( !interactive() )
    dev.off()
