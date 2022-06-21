
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TimeSeriesAnalysis

<!-- badges: start -->
<!-- badges: end -->

TimeSeriesAnalysis is a transcriptomic analysis tool for both RNA
sequencing and microarray data

## Overview

TimeSeries (TS) is an analysis and visualization package for RNAseq and
microarray data. TS extracts significant genes from time course
transcriptomic data by performing differential gene expression on both
the conditional and temporal axes. It then employs partitioning
algorithm based on recursive thresholding (PART) clustering to identify
small genomic clusters of relevance, followed by running the clusters
through gprofiler to reveal the biological relevance of each cluster.

TS performs: data normalization and processing PCA plots Differential
gene expression (conditional and temporal) PART clustering Heatmaps
Trajectory of identified clusters Gprofiler (functional enrichment)
analysis of clusters Dotplots and MDS plots of Gprofiler results Nearest
ancestor clustering of GOs GO ancestor queries

## Installation

You can install the development version of TimeSeriesAnalysis like so:

``` r
install.packages("devtools")
devtools::install_github("Ylefol/TimeSeriesAnalysis")
```

Certain bioconductor packages may have to be installed before
installation of TimeSeriesAnalysis. This can be done with the following
code snippet:

``` r
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
bio_pkgs <- c(PACKAGES HERE)
BiocManager::install(bio_pkgs)
```

## Rmarkdown format

TimeSeriesAnalysis was developed with user-friendliness in mind, the
core script of the package is a Rmarkdown file with explanations of what
each code block performs. Computationally intensive code blocks save
their work in order to avoid loss of computation time.
TimeSeriesAnalysis comes with a test dataset. Users are recommended to
first run the package with the test data to then modify the Rmarkdown
script with the necessary information for their purposes (input files,
organism, genes of interest etc…). The rmarkdown script and the expected
result from an example run can be found and downloaded from this
repository, they are located in the ‘rmarkdown_method’ folder. If users
prefer a script approach, the equivalent of the rmarkdown report is
provided in two scripts, one for the computation tasks and a second for
the analyses. These two scripts are found in the ‘script_method’ folder.

## Data

The data used as an example is data generated from the LEMONAID
consortium More info about the data?
