% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{Celegans_TS_data}
\alias{Celegans_TS_data}
\title{RNAseq counts for Celegans and krill oil timeseries exeriment time series experiment}
\format{
An object of class \code{list} of length 2.
}
\source{
<https://www.github.com/Ylefol/TimeSeriesAnalysis>
}
\usage{
Celegans_TS_data
}
\description{
A dataset with RNAseq results for a timeseries experiment with three timepoints, at days 1, 3. and 6.
The dataset contains two Celegans strain, the BY273 strain (a parkinsons model) and the N2 strain (WT). Each strain
is split into two groups, one treated with krill oil and the other is the control. In total the dataset contains
four groups: BY273.ctrl, BY273.krill, N2.ctrl, N2.krill

#' @format A list containing a dataframe and another list of count dataframes
\describe{
  \item{sample_dta}{File indicating which count files belong to which group, timepoint, and replicate}
  \item{counts}{A list of all the count files}
}
}
\keyword{datasets}
