# ![segmenTier](doc/logo.png) : Similarity-Based Segmentation of Multi-Dimensional Signals

## Quick Guide

```R
library(devtools)
install_github("raim/segmenTier")
library(segmenTier)
demo("segment_test") # direct algorithm interface and plotting
demo("segment_data") # time-series processing, clustering and batch segmentation
```

## Theory

The theory behind the package is outlined in 
the manuscript at the preprint server of 
the [Santa Fe Institute](https://www.santafe.edu/).

## Installation

```R
library(devtools)
install_github("raim/segmenTier")
```

## Usage

Usage of the package is demonstrated in two R demos.

### Direct Interface to Algorithm

The main low level interface to the algorithm, function `segmentClusters`,
is demonstrated in the file [demo/segment_test.R](demo/segment_test.R). 
It produces Supplemental Figure S1 of the preprint 
manuscript.

To run it as a demo in R grep simply type:
```
library(segmenTier)
demo("segment_test", package = "segmenTier")
```

### Interfaces to Time-Series Processing, Clustering and Batch Segmentation 

A real-life data set is processed, clustered and 
segmented with varying parameters in 
[demo/segment_data.R](demo/segment_data.R).

This demo runs quite long, since it calculates many 
segmentations. It produces Figure 3 and Supplemental Figures
S4a and S4b of the preprint manuscript.

```
demo("segment_data", package = "segmenTier")
```

