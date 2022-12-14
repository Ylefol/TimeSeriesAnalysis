% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DE_PART_results_functions.R
\name{wrapper_cluster_trajectory}
\alias{wrapper_cluster_trajectory}
\title{Plot cluster trajectory for all PART clusters}
\usage{
wrapper_cluster_trajectory(
  object,
  cluster_traj_dta,
  mean_cluster_traj_dta,
  log_TP = F,
  plot_name = "Ctraj"
)
}
\arguments{
\item{cluster_traj_dta}{The trajectory data for all clusters being calculated
The data is calculated/obtained from \code{calculate_cluster_traj_data} function}

\item{plot_name}{The name given to the plot file as it is saved}

\item{mean_cluster_traj_data}{The mean value for each clusters trajectory,
the data is calculated/obtained from \code{calculated_mean_cluster_traj} function}
}
\description{
Function which plots the trajectory of all the clusters given
The function will split the clusters into groups of 8 to keep
the figures clean. If multiple figures are created, they will be names
1_of_x until it reaches x_of_x.
}
