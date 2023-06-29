% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{murine_TS_data}
\alias{murine_TS_data}
\title{RNAseq counts for AID activation cocktail in a murine timeseries exeriment time series experiment}
\format{
A SummarizedExperiment containing the count matrix as the assay and the sample data as the colData
\describe{
  \item{sample_dta}{Dataframe indicating which count files belong to which group, timepoint, and replicate}
  \item{counts}{A matrix of all the count files}
}
}
\source{
<https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html>

<https://www.github.com/Ylefol/TimeSeriesAnalysis>
}
\usage{
murine_TS_data
}
\description{
A data set with RNAseq results for a AID activation cocktail composed of LPS, CD40L, and TGFb-1.
Two groups exist: Control (no cocktail) and Treated (with activation cocktail).
Initial experiment harvested timepoints at 0, 15min, 30min, 1h, 2h, 3h, 6h, 12h, 24h, and 48h.
Only one replicate exists per group per timepoint, so subsequent timepoint of the same groups have been merged
as a singular timepoint. For example, timepoints 0 and 15m became timepoint 7.5 minutes.
}
\keyword{datasets}