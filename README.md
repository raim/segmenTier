![segmenTier](doc/logo.png) 

# Similarity-Based Segmentation of Multi-Dimensional Signals

## Theory

The theory behind the package is outlined in
[Machne, Murray & Stadler 2017](http://www.nature.com/articles/s41598-017-12401-8).

## Installation

```R
library(devtools)
install_github("raim/segmenTier")
```

## Usage

### Quick Guide

```R
library(segmenTier)

data(primseg436) # RNA-seq time-series data

## cluster timeseries:
tset <- processTimeseries(ts=tsd, dc.trafo="ash", dft.range=1:7)
cset <- clusterTimeseries(tset, K=16, nui.thresh=0.6)
## and segment it:
segments <- segmentClusters(seq=cset, M=150, E=2, nui=3, S="icor")
## and inspect results:
plotSegmentation(tset,cset,segments
```

### Demos

Usage of the package is further demonstrated in two R demos:

#### Direct Interface to Algorithm

The main low level interface to the algorithm, function `segmentClusters`,
is demonstrated in the file [demo/segment_test.R](demo/segment_test.R). 
It produces Supplemental Figure S1 of the preprint 
manuscript.

To run it as a demo in R simply type:
```
library(segmenTier)
demo("segment_test", package = "segmenTier")
```

#### Interfaces to Time-Series Processing, Clustering and Batch Segmentation 

A real-life data set is processed, clustered and 
segmented with varying parameters in 
[demo/segment_data.R](demo/segment_data.R).

This demo runs quite long, since it calculates many 
segmentations. It produces Figure 3 and Supplemental Figures
S4a and S4b of the preprint manuscript.

```
demo("segment_data", package = "segmenTier")
```

