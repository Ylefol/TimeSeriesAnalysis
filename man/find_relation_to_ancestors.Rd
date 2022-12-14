% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MDS_GO_results_functions.R
\name{find_relation_to_ancestors}
\alias{find_relation_to_ancestors}
\title{Ancestor querying for GOs}
\usage{
find_relation_to_ancestors(target_ancestors, GOs_to_check, ontology = "BP")
}
\arguments{
\item{target_ancestors}{A vector of GO IDs to use as ancestors}

\item{GOs_to_check}{A dataframe created by either \code{gprofiler_for_conditional_DE}
or \code{gprofiler_cluster_analysis}}

\item{ontology}{The ontology to be used to search for ancestry}
}
\value{
A updated dataframe where all non-ancestor related GOs were removed
GOs with a found ancestor are associated to it via a new column called 'ancestor'
}
\description{
Function which finds if a term is associated with a target ancestor

The function will retrieve ancestry information from the annotationDbi of the
loaded library as well as the inputted ontology (BP by default). It will then
isolate and preserve the GO terms found to be related to one or several of the target
ancestors. If the GO is related to several ancestors, only the one highest in the
tree is linked to that GO.

The function also associates a color to each ancestor, due to this, there is a
maximum of 9 ancestors that can be selected
}
