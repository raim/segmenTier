
library("segmenTier")
source("~/programs/segmenTier/R/cluster.R")
source("~/programs/segmenTier/R/segment.R")

## load time-series data
## contains tsd from primseg436 from unpublished yeast transcriptome timeseries
data(primseg436)

### TIME-SERIES CLUSTERING PARAMETERS
use.fft <- TRUE # cluster discrete Fourier transform of data?
use.snr <- TRUE # take SNR of DFT
trafo <- "" # "ash" # "log" # 
low.thresh <- 1 # -Inf/0 # filter by log(DC-component=total signal over time)
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


tset <- processTimeseries(ts=tsd, smooth=FALSE, trafo=trafo,
                          use.fft=use.fft, dft.range=dft.range,
                          use.snr=use.snr, low.thresh=low.thresh)


cset <- clusterTimeseries(tset, selected=selected,iter.max=iter.max, nstart=nstart)

allsegs <- segmentCluster.batch(cset, csim.scale=csim.scale, score=scores,
                                M=M, Mn=Mn, a=2, nui=nui.cr,
                                nextmax=nextmax, multi=multi,multib=multib, 
                                ncpu=1, verb=1, save.mat="")


    ## plot all
#out <- file.path(paste(outname,"_",segid,sep=""))
#coors <- c(chr=1,unlist(primseg[i,1:2])) 
#plotdev(out,width=6,height=3.5,type=fig.type)

browser.path <- sub("GENBRO=","",system("env|grep GENBRO",intern=TRUE))
source(file.path(browser.path,"src/genomeBrowser_utils.R")) # for coor2index
source(file.path(browser.path,"src/genomeBrowser.R")) ## for chrS

N <- nrow(tsd)
coors <- c(chr=1,start=1,end=N)

tsd <- tset$ts # incl. all trafos and zeros set to NA
#tsd[tset$zero.vals,] <- NA
tot <- tset$tot
low <- tset$low.vals

library("colorRamps")
colors0 <- matlab.like(100)  ## read count
colors0[1] <- "#FFFFFF" ## replace minimal by white


png("testset_srg1_scales.png",width=10,height=5,units="in",res=300)
par(mfcol=c(3,1),mai=c(.01,1.5,.01,.01),mgp=c(1.7,.5,0),xaxs="i")
plot(1:N,tot,log=ifelse(trafo!="","","y"),type="l",lwd=2,axes=FALSE)
points((1:N)[low],tot[low],col=2,cex=.5)
axis(2);
plotHeat(tsd,coors=coors,chrS=0,colors=colors0, colnorm=TRUE)
columns <- c(name="ID", type="type", start="start", end="end")
tmp <- plotFeatures(allsegs, coors=coors,
                    typord=TRUE,cuttypes=TRUE,
                    ylab="", names=FALSE,columns=columns,tcx=.5)
dev.off()
