
library("segmenTier")
##source("~/programs/segmenTier/R/segment.R")

## a sequence of clusters - note that `0' will be treated
## as the nuissance cluster, which will not result in segments
seq <- c(4,4,1,2,4,4,4,2,4,4,3,4,4,2,4,1,4,4,0,0,0,0,0,0,1,2,2,2,4,1,1,1,1,
         0,0,1,1,1,1,1,3,3,3,0,3,3,3,3,0,1,3,4,3,2,4,4,1,3)

## list of non-nuissance clusters
C <- 1:4

## scoring function "ccor":
## cluster X cluster correlation matrix, for "ccor"
cr <- matrix(0,ncol=length(C),nrow=length(C))
diag(cr)<- 1 # set diagonal to 1
cr[1,2] <- cr[2,1] <- -.4 # just enough to avoid merging of 2 with 1
cr[1,3] <- cr[3,1] <- -.5
cr[2,4] <- cr[4,2] <- -.4 # just enought to avoid merging of 2 with 4

## scoring function "icor":
## position X cluster correlation matrix, for "icor"
ic <- matrix(0,ncol=length(C),nrow=length(seq))
for ( i in 1:length(seq) ) {
    for ( j in 1:ncol(ic) )
        set.seed(42) # to keep results constant
        ic[i,j] <- sample(seq(-1,.1,.01))[3] # set bad correlation
    if ( seq[i]!=0 ) { # set good correlation to its own cluster
        set.seed(42) # to keep results constant
        ic[i,seq[i]] <- sample(seq(.3,1,.01))[1]
    }
}

## SCORING FUNCTIONS to be tested
scores <- c("ccor","cls","icor")

## SCORING FUNCTION PARAMETERS
## segment size and penalty parameters
M <- 3 # minimal segment size; note: this is not a strict cutoff
a <- 2 # penalty for non-matching clusters in scoring function "cls"

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
    if ( score == "ccor" ) csim <- cr
    else if ( score == "icor") csim <- ic
    else csim <- NULL

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
            seg <- segmentClusters(seq=seq,csim=csim,score=score,M=M,a=a,nui=1,
                                   multi=multi, multib=multib,nextmax=nextmax,
                                   save.mat=c("SM","SK"))
            multS$SM <- seg$SM # TODO: test whether they are the same?
            multS[[multi]]$SK <- seg$SK
            ## store segments!
            multS[[multi]][[multib]] <- seg$segments

        }
    }
    ## store all results for the current scoring function
    scrR[[score]] <- multS
}

## plot scoring function matrices
for ( score in scores ) 
    plotScoring(scrR[[score]]$SM, seq=seq) 
## plotting segments and total scoring matrices
plotSegments(scrR, seq=seq) #,out="segment_test")
