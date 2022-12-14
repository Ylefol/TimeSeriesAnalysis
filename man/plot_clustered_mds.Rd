% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MDS_GO_results_functions.R
\name{plot_clustered_mds}
\alias{plot_clustered_mds}
\title{plot clustered MDS}
\usage{
plot_clustered_mds(main_matrix, cluster_info)
}
\arguments{
\item{main_matrix}{Matrix containing the data to be plotted as created by  \code{SS_GO_clusters}}

\item{cluster_info}{A dataframe with cluster names and the number of GOs in each cluster}
}
\value{
the plotly object for the MDS plot
}
\description{
This function was adapted from the MDSplot() function from the visEAGO package
}
