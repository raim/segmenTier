
### TIME-SERIES CLUSTERING PARAMETERS
use.fft <- TRUE # cluster discrete Fourier transform of data?
use.snr <- TRUE # take SNR of DFT
trafo <- "" # "ash" # "log" # 
low.thresh <- 1 # -Inf/0 # filter by log(DC-component=total signal over time)
dft.range <- 2:7 # range of DFT to cluster to use for clustering
selected <- c(8,12,16) # cluster number K
kiter <- 100000 # max. iterations in kmeans
nstart <- 100   # number of initial configurations tested in kmeans

### SEGMENTATION PARAMETERS
nui.cr <- 1 #  -/+ correlation of nuissance cluster with others and itself

## SCORING FUNCTION
## which scoring functions to use
scores <- c("ccor","icor")
scale <- 3 # 3 # 1 ## scale exponent of similarity matrices csim
## scoring function minimal length penalty
M <- 175 # 30 # for empty set?
Mn <- 15 # for nuissance clusters: allow smaller segments!

## TOTAL SCORE and BACK-TRACING
## NOTE: these have little effect on real data sets, except perhaps nextmax?
## handling of multiple max. score k in scoring
multi <- "max" # c("max","min")
## handling of multiple max. score clusters in back-tracing
multib <- "max" # c("max","skip","min")
## in back-tracing, search for the next non-decreasing S(i,c)
nextmax <-TRUE

## TODO: load time-series data

tset <- processTimeseries(ts=tsd,
                          smooth=FALSE, trafo=trafo,
                          use.fft=use.fft, dft.range=dft.range,
                          use.snr=use.snr, low.thresh=low.thresh)


cset <- clusterTimeseries(tset, selected=selected,kiter=kiter, nstart=nstart)

allsegs <- segmentClusterset(cset, csim.scale=scale,
                             scores=scores,
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

par(mfcol=c(5,1),mai=c(.01,.8,.01,.01),mgp=c(1.7,.5,0),xaxs="i")
plot(1:N,tot,log=ifelse(trafo!="","","y"),type="l",lwd=2,axes=FALSE)
points((1:N)[low],tot[low],col=2,cex=.5)
axis(2);
plotHeat(tsd,coors=c(chr=1,start=1,end=N),chrS=0,colors=colors0,
         colnorm=TRUE)
segcols <- columns; segcols["color"] <- "CL"
typs <- sort(unique(allsegs[,"type"])); sgtypes <- typs;#sgtypes <- NULL;
                                        #for ( score in scores ) sgtypes <- c(sgtypes,grep(score,typs,value=TRUE))
tmp <- plotFeatures(allsegs,
                    coors=index2coor(t(as.matrix(c(coors,strand=NA))),chrS),
                    types=sgtypes,typord=TRUE,cuttypes=TRUE,
                    ylab="", names=FALSE,columns=segcols,tcx=.5)
#tmp <- plotFeatures(feats, coors=coors,
#                    types=types,typord=TRUE,cuttypes=TRUE,
#                    ylab="", names=FALSE,columns=columns,tcx=.5)
#axis(1)
#tmp <- plotFeatures(transcripts, coors=coors,
#                    typord=TRUE,cuttypes=TRUE,ylab="", names=TRUE,
#                    columns=columns)
