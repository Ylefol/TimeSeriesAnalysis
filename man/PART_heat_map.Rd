% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DE_PART_results_functions.R
\name{PART_heat_map}
\alias{PART_heat_map}
\title{Create and plot PART heatmap}
\usage{
PART_heat_map(object, heat_name = "custom_heat_map")
}
\arguments{
\item{object}{A time series object}

\item{heat_name}{The file name given to the saved heatmap}
}
\value{
none
}
\description{
The heatmap shows selected genes as rows and replicates of the time series experiment
as columns.
Gene clusters are identified by colored row annotation. Replicates are split
by their grouping (conditional) and ordered (within each other) via time points
The legend for groupings, clusters, timepoints, and heatmap values is given on
the right hand side of the heatmap.

the heatmap is saved twice, once
in png format, and the other in svg format
}
