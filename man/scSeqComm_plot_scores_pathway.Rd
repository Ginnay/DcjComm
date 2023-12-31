% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{scSeqComm_plot_scores_pathway}
\alias{scSeqComm_plot_scores_pathway}
\title{Plot intercellular and intracellular scores with pathway information (heatmap)}
\usage{
scSeqComm_plot_scores_pathway(
  data,
  title,
  facet_grid_x = "cluster_R",
  annotation_GO = NULL,
  cutoff = 0.05,
  topGO = 5,
  GO_ncol = 1
)
}
\arguments{
\item{data}{A data frame containing results of signaling analysis.}

\item{title}{A character string to be used as title of plot.}

\item{facet_grid_x}{A variable defining faceting group on columns dimension.}

\item{annotation_GO}{A data.frame or a list of data.frame with results of enrichment analysis.}

\item{cutoff}{A numeric value defining the significance value for GO terms.}

\item{topGO}{A numeric value defining the maximum number of terms to be visualized.}

\item{GO_ncol}{A numeric value defining the number of columns on which show the GO terms table.}
}
\value{
Plot showing intracellular and intercellular scores and pathway
}
\description{
Function to plot the intercellular (mapped in point size scale) and intracellular scores (mapped in color scale) of the given input data with pathway information on y axis and cluster pairs on x axis.
Grey color is associated to NA values of intracellular score.
}
