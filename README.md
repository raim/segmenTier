![segmenTier](doc/logo.png) 

# Similarity-Based Segmentation of Multi-Dimensional Signals

`segmenTier` is a dynamic programming solution to segmentation based
 on maximization of arbitrary similarity measures within segments.
 The general idea, theory and this implementation are described in
 [Machne, Murray & Stadler
 (2017)](http://www.nature.com/articles/s41598-017-12401-8).  In
 addition to the core algorithm, the package provides time-series
 processing and clustering functions as described in the
 publication. These are generally applicable where a `kmeans`
 clustering yields meaningful results, and have been specifically
 developed for clustering of the Discrete Fourier Transform of
 periodic gene expression data (`circadian' or `yeast metabolic
 oscillations'). This clustering approach is outlined in the
 supplemental material of [Machne & Murray (2012)]
 (https://doi.org/10.1371/journal.pone.0037906), and here is used as a
 basis of segment similarity measures.

## Theory


The theory behind the package is outlined in detail in
[Machne, Murray & Stadler 2017](http://www.nature.com/articles/s41598-017-12401-8) and summarized in the package vignette.

## Installation

```R
library(devtools)
install_github("raim/segmenTier", subdir = "pkg")
```

## Usage

### Quick Guide

```R
library(segmenTier)

data(primseg436) # RNA-seq time-series data

## cluster timeseries:
tset <- processTimeseries(ts=tsd, dft.range=1:7, dc.trafo="ash")
cset <- clusterTimeseries(tset)
## and segment it:
segments <- segmentClusters(seq=cset, M=100, E=2, nui=3, S="icor")
## and inspect results:
plotSegmentation(tset,cset,segments)
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

