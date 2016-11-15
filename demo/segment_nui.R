
library("segmenTier")

## load time-series data
## contains tsd from primseg436 for
## a 7.6 kb genomic region
data(primseg436)

### TIME-SERIES CLUSTERING PARAMETERS
use.fft <- TRUE # cluster discrete Fourier transform of data?
use.snr <- TRUE # take SNR of DFT
trafo <- "" # "ash" # "log" #
## nuissance assignment test: 
low.thresh <- -Inf #1/0 # minimal mean value (DC component of DFT if use.fft)
keep.zeros <- FALSE
dft.range <- 1:7 # range of DFT to cluster to use for clustering
dc.trafo <- "ash"
selected <- c(16) # cluster number K
iter.max <- 100000 # max. iterations in kmeans
nstart <- 100   # number of initial configurations tested in kmeans

## re-assign clusters
nui.thresh <- 0.6

### SEGMENTATION PARAMETERS
nui.cr <- 2 #  -/+ correlation of nuissance cluster with others and itself
    
## SCORING FUNCTION
## which scoring functions to use
scores <- c("ccor","icor")
csim.scale <- c(1,3,5) # 3 # 1 ## scale exponent of similarity matrices csim
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
tset <- processTimeseries(ts=tsd,
                          smooth=FALSE, trafo=trafo, keep.zeros=keep.zeros,
                          use.fft=use.fft, dft.range=dft.range,
                          use.snr=use.snr, dc.trafo=dc.trafo,
                          low.thresh=low.thresh)

## CLUSTER PRE-PROCESSED TIME SERIES
cset <- clusterTimeseries(tset, selected=selected,iter.max=iter.max, nstart=nstart)

## filter low-correlating by theta
for ( k in 1:ncol(cset$clusters) ) {
    cls <- cset$clusters[,k]
    Pci <- cset$Pci[[k]]
    for ( i in 1:nrow(Pci) )
        if ( !any(Pci[i,] > nui.thresh, na.rm=TRUE) )
            cls[i] <- -1
    cls[cls==-1] <- 0
    cset$clusters[,k] <- cls
}

## CALCULATE SEGMENTS FOR ALL CLUSTERINGS and
## FOR CHOSEN SEGMENTATION PARAMETERS
allsegs <- segmentCluster.batch(cset, csim.scale=csim.scale, score=scores,
                                M=M, Mn=Mn, a=2, nui=nui.cr,
                                nextmax=nextmax, multi=multi,multib=multib, 
                                ncpu=1, verb=1, save.mat="")



## break-count
starts <- ends <- rep(0,N)
tab <- table(allsegs[,"start"])
starts[as.numeric(names(tab))] <- tab
tab <- table(allsegs[,"end"])
ends[as.numeric(names(tab))] <- tab

## TODO bin starts


## type-count
ntype <- length(unique(allsegs[,"type"]))
## consensus segments: in more then half of types
n.thresh <- ntype/2
constarts <- which(starts>=n.thresh)
conends <- which(ends>=n.thresh)

consegs <- data.frame(matrix(NA,ncol=ncol(allsegs),nrow=length(constarts)))
colnames(consegs) <- colnames(allsegs)

consegs[,"start"] <- constarts
consegs[,"end"] <- conends
consegs[,"fuse"] <- FALSE
consegs[,"type"] <- "consensus"


consegs <- consegs[consegs[,"end"]-consegs[,"start"] > 1,]

segs <- rbind(allsegs,consegs)

## search closest end:

## moing average of hits for plots
starts <- ma(starts,5)*5
ends <- ma(ends,5)*5
starts[starts==0] <- NA
ends[ends==0] <- NA




## PLOT RESULTS


## get time-series data
ts <- tset$ts # incl. all trafos and zeros set to NA
ts[tset$zero.vals,] <- NA
tot <- tset$tot # total of the time-series

N <- nrow(ts)
coors <- c(chr=1,start=1,end=N) # "chromosome" coordinates

colors0 <- rev(gray.colors(100)) ## heatmap colors for the timeseries 
colors0[1] <- "#FFFFFF" ## replace minimal by white

png(paste("~/programs/segmenTier/test/testset_nui",nui.thresh,"_low",low.thresh,ifelse(keep.zeros,"_with0","_no0"),".png",sep=""), width=8,height=6, res=300, units="in")

par(mfcol=c(4,1),mai=c(.3,1.5,.01,.01),mgp=c(1.3,.5,0),xaxs="i")

#1
plot(1:N,tot,log=ifelse(trafo!="","","y"),type="l",lwd=2,axes=FALSE,ylab=NA)
polygon(x=c(1,1,N,N),y=c(min(tot,na.rm=TRUE),rep(low.thresh,2),min(tot,na.rm=TRUE)),col="#00000055",border=NA)
abline(h=low.thresh,col="#000000BB")
lines(1:N,tot)
axis(2);
axis(1)
mtext("total signal", 2, 2)
#2
segment.plotHeat(ts,coors=coors,chrS=0,colors=colors0, colnorm=TRUE)
axis(2,at=1:ncol(ts))
axis(1)
mtext("time points", 2, 2)
## TODO: plot clustering
columns <- c(name="ID", type="type", start="start", end="end")
#3
ypos <- segment.plotFeatures(segs, coors=coors,
                             typord=TRUE,cuttypes=TRUE,
                             ylab="", names=FALSE,columns=columns,tcx=.5)
axis(1)
## plot fuse tag
fuse <- segs[segs[,"fuse"],]
points(fuse[,"start"],ypos[fuse[,"type"]],col="red",pch=3, lwd=2,cex=2)

#4 hit number
plot(1:N,starts,type="h",lwd=2,xlim=coors[2:3],ylim=c(0,max(c(starts,ends)+1,na.rm=TRUE)),axes=FALSE,ylab=NA,col=NA)
abline(h=ntype/2)
lines(1:N +0.5,starts,type="h",col="#0000FFAA",lwd=2)
lines(1:N -0.5,ends,  type="h",col="#FF0000AA",lwd=2)
axis(1)
axis(2)
mtext("hit count", 2, 2)
dev.off()
