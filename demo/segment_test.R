
library("segmenTier")
library("Rcpp")
source("~/programs/segmenTier/R/plot.R")
source("~/programs/segmenTier/R/cluster.R")
source("~/programs/segmenTier/R/segment.R")
sourceCpp("~/programs/segmenTier/src/segment.cpp")

## a sequence of clusters - note that `0' will be treated
## as the nuissance cluster, which will not result in segments
seq <- c(4,4,1,2,4,4,4,2,4,4,3,4,4,2,4,1,4,4,0,0,0,0,0,0,1,2,2,2,4,1,1,1,1,
         0,0,1,1,1,1,1,3,3,3,0,3,3,3,3,0,1,3,4,3,2,4,4,1,3)

### clustering can also be letters
## TODO: use this to test allowing character clusters
#seq <- letters[seq]
## BUT: nuissance cluster must be 0 or "0"!
#seq[seq=="a"]  <- "0"

## list of non-nuissance clusters
C <- unique(seq)
C <- C[as.character(C)!="0"]

## SCORING FUNCTIONS to be tested
scores <- c("ccor","icor", "ccls") # ,"cls"

## SCORING FUNCTION PARAMETERS
## segment size and penalty parameters
M <- 3 # minimal segment size; note: this is not a strict cutoff
a <- -2 # penalty for non-matching clusters in scoring function "ccls"

## CLUSTER SIMILARITIES
## scoring function "ccor":
## cluster X cluster correlation matrix csim for "ccor"
Ccc <- matrix(0,ncol=length(C),nrow=length(C))
colnames(Ccc) <- rownames(Ccc) <- as.character(C)
diag(Ccc)<- 1 # set diagonal to 1
Ccc[1,2] <- Ccc[2,1] <- -.6 # just enough to avoid merging of 2 with 1
Ccc[2,4] <- Ccc[4,2] <- -.4 # just enought to avoid merging of 2 with 4
Ccc[1,3] <- Ccc[3,1] <- -.5

## scoring function "icor":
## position X cluster correlation matrix csim for "icor"
Pci <- matrix(0,ncol=length(C),nrow=length(seq))
for ( i in 1:length(seq) ) {
    for ( j in 1:ncol(Pci) )
        set.seed(42) # to keep results constant
        Pci[i,j] <- sample(seq(-1,.1,.01))[3] # set bad correlation
    if ( seq[i]!=0 ) { # set good correlation to its own cluster
        set.seed(42) # to keep results constant
        Pci[i,seq[i]] <- sample(seq(.3,1,.01))[1]
    }
}

## construct cset
cset <- list()
class(cset) <- "clustering"
cset$clusters <- matrix(seq,ncol=1)
colnames(cset$clusters) <- paste("K",length(C),sep="")
cset$Ccc <- cset$Pci <- list()
cset$Ccc[[1]] <- Ccc
cset$Pci[[1]] <- Pci

## cluster sorting & coloring
cset <- colorClusters(cset)

## TOTAL SCORING MATRIX S(i,c)
## handling of multiple max. score k 
## "max" for highest k (favors shorter segments) or
## "min" for lowest k (favors longer segments)
multis <- c("max","min") # total scoring matrix, take

## BACK-TRACING:
## when back-tracing search the next highest score before starting a segment
nextmax <-TRUE
## handling of multiple max. score clusters in back-tracing
## "skip" or 
## "max" for highest k (favors shorter segments) or
## "min" for lowest k (favors longer segments)
multibs <- c("max","skip","min") # back-tracing


scrR <- rep(list(NA),length(scores)) # result list for scoring function
names(scrR) <- scores
for ( score in scores ) {

    ## set csim according to scoring function
    if ( score == "ccor" ) csim <- Ccc
    else if ( score == "icor") csim <- Pci
    else csim <- NULL # will be constructed for 'ccls'

    ## result list for total scoring, will also hold the used
    ## scoring function matrices
    multS <- rep(list(NA),length(multis) +1) 
    names(multS) <- c("SM",multis)
    
    for ( multi in multis ) {

        ## result list for back-tracing, will also hold
        ## the total scoring and back-tracing matrices
        multS[[multi]] <- rep(list(NA),length(multibs)+1) 
        names(multS[[multi]]) <- c("SK",multibs)
        
        for ( multib in multibs ) {

            ## run the algorithm with current parameters, and
            ## store SM and SK matrices for later plots
            seg <- segmentClusters(seq=seq,csim=csim,score=score,M=M,Mn=M,
                                   a=a,nui=1,
                                   multi=multi, multib=multib,nextmax=nextmax,
                                   save.mat=TRUE)
            ##multS$SM <- seg$SM # TODO: test whether they are the same?
            multS[[multi]]$SK <- seg$SK
            ## store segments!
            multS[[multi]][[multib]] <- seg$segments
        }
    }
    ## store all results for the current scoring function
    scrR[[score]] <- multS
}

## plotting segments and score traces (from total scoring matrix)
plotSegments(scrR, seq=seq) #,out="segment_test")

## plot scoring function matrices
##for ( score in scores ) 
##    plotScoring(scrR[[score]]$SM, seq=seq, score=score) 
