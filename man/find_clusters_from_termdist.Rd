% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MDS_GO_results_functions.R
\name{find_clusters_from_termdist}
\alias{find_clusters_from_termdist}
\title{Cluster GO terms using Wang measurement}
\usage{
find_clusters_from_termdist(GO_df, semantic_data)
}
\arguments{
\item{GO_df}{The dataframe of GOs as returned by \code{gprofiler_cluster_analysis}}

\item{semantic_data}{semantic similarity data as created by the godata function}
}
\value{
clust_df dataframe showing each GO with it's cluster
}
\description{
Function which clusters GO terms and returns the result in dataframe format

This function uses the Wang measurement method with the ward.D2 aggregation method
clustering parameters are pre-set within the function.
}
