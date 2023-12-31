% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Tfactors_score_all_clusters.R
\name{Tfactors_score_all_clusters}
\alias{Tfactors_score_all_clusters}
\title{Compute the tfactors score in all cluster}
\usage{
Tfactors_score_all_clusters(regulatoryDB, genes_score, cell_cluster, cutoff)
}
\arguments{
\item{regulatoryDB}{The input regulatoryDB.}

\item{genes_score}{The TGpvalues obtained above.}

\item{cell_cluster}{The input cell_group.}

\item{cutoff}{Select genes after Bonferroni correction according to cutoff, default "0.05.}
}
\description{
Compute the tfactors score in all cluster
}
