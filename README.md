# ![segmenTier](doc/logo.png)

Similarity-Based Segmentation of
Multi-Dimensional Signals

## Theory

The theory behind the package is outlined in 
the manuscript at the preprint server of the
the [Santa Fe Institute](https://www.santafe.edu/).

## Installation

```R
library(devtools)
install_github("raim/segmenTier")
```

## Usage

Usage of the package is demonstrated in two R demos.

Artificial data to demonstrate the low level functions
are in file [demo/segment_test.R]([demo/segment_test.R]). 
It produces Supplemental Figure S1 of the preprint 
manuscript.

To run it in R grep simply type:
```
library(segmenTier)
demo("segment_test", package = "segmenTier")
```

A real-life data set is segmented with varying
parameters in [demo/segment_data.R](demo/segment_data.R).

This demo runs quite long, since it calculates many 
segmentations. It produces Figure 3 and Supplemental Figures
S4a and S4b of the preprint manuscript.

```
demo("segment_data", package = "segmenTier")
```

