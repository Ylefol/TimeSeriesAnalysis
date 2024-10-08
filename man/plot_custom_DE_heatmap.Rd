% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DE_PART_results_functions.R
\name{plot_custom_DE_heatmap}
\alias{plot_custom_DE_heatmap}
\title{Plots DE heatmaps}
\usage{
plot_custom_DE_heatmap(
  heat_mat,
  col_split,
  row_splits,
  l2fc_col,
  log_transform,
  color_vect,
  legend_value = "counts",
  plot_file_name = "custom_heatmap",
  custom_width = 15,
  custom_height = 5,
  do_SVG = TRUE
)
}
\arguments{
\item{heat_mat}{The main count matrix of the heatmap}

\item{col_split}{The named vector used to split the columns into it's segments}

\item{row_splits}{The named vector used to split the rows into the two groups}

\item{l2fc_col}{The log2foldchange vector used to create the histogram}

\item{log_transform}{boolean indicating if the results were log transformed or not}

\item{color_vect}{A vector for the colors to be used in the heatmap}

\item{legend_value}{string indicating the type of value used
(counts for RNAseq, intensity value for microarray)}

\item{plot_file_name}{The name given to the saved heatmap}

\item{custom_width}{The width of the heatmap}

\item{custom_height}{The height of the heatmap}

\item{do_SVG}{Boolean to check if the SVG file should be created or not}
}
\value{
none or heatmap plot
}
\description{
The function is intended to be called within the \code{custom_heatmap_wrapper}

The function creates a segmented heatmap where rows are samples/patients and columns
are genes. The heatmap is column segmented based on the number of differential expression
analyses contained within it. A log2foldChange histogram is attached to the bottom
of the heatmap to indicate significance of the genes being visualized.

The column/segment legend is given at the bottom of the heatmap with the number
of genes in each group being also provided.

The gradient legend of the values within the heatmap is given on the right hand side.

The heatmap is saved in both png and svg format.
}
\examples{
TS_object<-create_example_object_for_R()
TS_object <- normalize_timeSeries_with_deseq2(time_object=TS_object)
#Perform conditional differential gene expression analysis
TS_object<-conditional_DE_wrapper(TS_object,vignette_run=TRUE)
heat_dta<-create_conditional_heatmap_matrix(TS_object)
heat_dta<-prepare_heat_data(heat_dta,log_transform=TRUE)

log_l2fc<-log_transform_l2fc_vect(heat_dta[['l2fc_vector']])

heat_plot<- plot_custom_DE_heatmap(heat_dta[['heat_matrix']],heat_dta[['region_split']],
                                   heat_dta[['group_split']],log_l2fc,log_transform = TRUE,
                                   legend_value='counts', plot_file_name = NULL,
                                   do_SVG=FALSE)

}
