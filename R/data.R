#' @title RNAseq counts for AID activation cocktail in PBMC time series experiment
#'
#' @description A data set with RNAseq results for various AID activation cocktails on a PBMC cell line.
#' Four groups exist: Control (no cocktail), LPS, IgM activation, and TGFb activation. Every group
#' has two replicates per timepoint. Timepoints are days 1, 3, and 9. The data was subsetted to only contain
#' the genes found to be differentially expressed in the PBMC analysis. Subsetting was performed
#' in order to alleviate the vignettes and examples.
#'
#' @format A SummarizedExperiment containing the count matrix as the assay and the sample data as the colData
#' \describe{
#'   \item{sample_dta}{Dataframe indicating which count files belong to which group, timepoint, and replicate}
#'   \item{counts}{A matrix of all the count files}
#' }
#' @source <https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html>
#' @source <https://www.github.com/Ylefol/TimeSeriesAnalysis>
"PBMC_TS_data"

#' @title Pre-computed PART, gprofiler,semantic data, and differential expression results for the PBMC dataset
#'
#' @description A list containing the gprofiler results for the PBMC example.
#' This data exists to allow for the creation of Vignettes in the event of a network failure.
#' Also contains a list of lists with various PART results, from the inputed matrix to the
#' computed clusters and the various parameters used for the clustering.
#'
#' @format A list of lists
#' @source <https://www.github.com/Ylefol/TimeSeriesAnalysis>
"PBMC_pre_loaded"
