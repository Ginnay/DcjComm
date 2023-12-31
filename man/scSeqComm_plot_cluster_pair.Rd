% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{scSeqComm_plot_cluster_pair}
\alias{scSeqComm_plot_cluster_pair}
\title{Plot intercellular and intracellular scores of ligand-receptor pairs between a given cell clusters pair (heatmap)}
\usage{
scSeqComm_plot_cluster_pair(data, title, selected_cluster_pair = NULL)
}
\arguments{
\item{data}{A data frame containing results of signaling analysis.}

\item{title}{A character string to be used as title of plot.}

\item{selected_cluster_pair}{String or array of strings containg cell clusters pairs to plot.}
}
\value{
A plot (or a list of plots)
}
\description{
Function to plot intercellular (mapped in point size scale) and intracellular scores (mapped in color scale) of ligand-receptor pairs between a given cell clusters pair.
If input data contains more than one cell clusters pair, a plot will be generated for each cell clusters pair.
Grey color is associated to NA values of intracellular score.
}
