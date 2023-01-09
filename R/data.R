#' @title RNAseq counts for AID activation cocktail in PBMC time series experiment
#'
#' @description A data set with RNAseq results for various AID activation cocktails on a PBMC cell line.
#' Four groups exist: Control (no cocktail), LPS, IgM activation, and TGFb activation. Every group
#' has two replicates per timepoint. Timepoints are days 1, 3, and 9.
#'
#' @format A SummarizedExperiment containing the count matrix as the assay and the sample data as the colData
#' \describe{
#'   \item{sample_dta}{Dataframe indicating which count files belong to which group, timepoint, and replicate}
#'   \item{counts}{A matrix of all the count files}
#' }
#' @source <https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html>
#' @source <https://www.github.com/Ylefol/TimeSeriesAnalysis>
"AID_TS_data"


#' @title RNAseq counts for AID activation cocktail in a murine timeseries exeriment time series experiment
#'
#' @description A data set with RNAseq results for a AID activation cocktail composed of LPS, CD40L, and TGFb-1.
#' Two groups exist: Control (no cocktail) and Treated (with activation cocktail).
#' Initial experiment harvested timepoints at 0, 15min, 30min, 1h, 2h, 3h, 6h, 12h, 24h, and 48h.
#' Only one replicate exists per group per timepoint, so subsequent timepoint of the same groups have been merged
#' as a singular timepoint. For example, timepoints 0 and 15m became timepoint 7.5 minutes.
#'
#' @format A SummarizedExperiment containing the count matrix as the assay and the sample data as the colData
#' \describe{
#'   \item{sample_dta}{Dataframe indicating which count files belong to which group, timepoint, and replicate}
#'   \item{counts}{A matrix of all the count files}
#' }
#' @source <https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html>
#' @source <https://www.github.com/Ylefol/TimeSeriesAnalysis>
"murine_TS_data"


#' @title RNAseq counts for Celegans and krill oil timeseries exeriment time series experiment
#'
#' @description A dataset with RNAseq results for a timeseries experiment with three timepoints, at days 1, 3. and 6.
#' The dataset contains two Celegans strain, the BY273 strain (a parkinsons model) and the N2 strain (WT). Each strain
#' is split into two groups, one treated with krill oil and the other is the control. In total the dataset contains
#' four groups: BY273.ctrl, BY273.krill, N2.ctrl, N2.krill
#'
#' @format A SummarizedExperiment containing the count matrix as the assay and the sample data as the colData
#' \describe{
#'   \item{sample_dta}{Dataframe indicating which count files belong to which group, timepoint, and replicate}
#'   \item{counts}{A matrix of all the count files}
#' }
#' @source <https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html>
#' @source <https://www.github.com/Ylefol/TimeSeriesAnalysis>
"Celegans_TS_data"
