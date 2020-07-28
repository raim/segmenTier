[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/segmenTier)](https://cran.r-project.org/package=segmenTier)
[![Downloads](https://cranlogs.r-pkg.org/badges/segmenTier)](https://cran.r-project.org/package=segmenTier)


![segmenTier](pkg/vignettes/logo.png) 

# Similarity-Based Segmentation of Multi-Dimensional Signals

`segmenTier` is a dynamic programming solution to segmentation based
 on maximization of arbitrary similarity measures within segments.
 The general idea, theory and this implementation are described in
 [Machne, Murray & Stadler (2017)](http://www.nature.com/articles/s41598-017-12401-8).
 In addition to the core algorithm, the package provides time-series
 processing and clustering functions as described in the
 publication. These are generally applicable where a `kmeans`
 clustering yields meaningful results, and have been specifically
 developed for clustering of the Discrete Fourier Transform of
 periodic gene expression data (circadian or yeast metabolic
 oscillations). This clustering approach is outlined in the
 supplemental material of [Machne & Murray (2012)]
 (https://doi.org/10.1371/journal.pone.0037906), and here is used as a
 basis of segment similarity measures.

## News

* Version 0.1.2: 
    - more general defaults in `processTimeseries` 
      (`use.fft=FALSE`, `na2zero=FALSE`) allow to set-up time-series 
      without any transformations for clustering
    - Doc and vignette have been substantially re-worked

## Theory


The theory behind the package is outlined in detail in
[Machne, Murray & Stadler 2017](http://www.nature.com/articles/s41598-017-12401-8) and summarized in the package vignette.

## Installation

The development version can be installed from github using
[`devtools`](https://cran.r-project.org/package=devtools):

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
tset <- processTimeseries(ts=tsd, na2zero=TRUE, use.fft=TRUE,
                          dft.range=1:7, dc.trafo="ash", use.snr=TRUE)
cset <- clusterTimeseries(tset, K=12)
## and segment it:
segments <- segmentClusters(seq=cset, M=100, E=2, nui=3, S="icor")
## inspect results:
print(segments)
plotSegmentation(tset,cset,segments)
## and get segment border table for further processing
segments$segments
```

### Demos

Usage of the package is further demonstrated in two R demos:

#### Demo I: Direct Interface to Algorithm

The main low level interface to the algorithm, function
`segmentClusters`, is demonstrated in the file
[demo/segment_test.R](demo/segment_test.R).  It produces Supplemental
Figure S1 of [Machne, Murray & Stadler
2017](http://www.nature.com/articles/s41598-017-12401-8).

To run it as a demo in R simply type:
```
library(segmenTier)
demo("segment_test", package = "segmenTier")
```

#### Demo II: Clustering, Batch Segmentation & Parameter Scans

A real-life data set is processed, clustered and 
segmented with varying parameters in 
[demo/segment_data.R](demo/segment_data.R).

This demo runs quite long, since it calculates many segmentations. It
provides a comprehensive overview of the effects of segmentation
parameters `E`, `M` and `nui`, and produces (among others) Figure 3
and Supplemental Figures S4a and S4b of [Machne, Murray & Stadler
2017](http://www.nature.com/articles/s41598-017-12401-8).

```
demo("segment_data", package = "segmenTier")
```

# Karl, the segmenTier

![](pkg/vignettes/anadenobolus_arboreus_credit.png)
