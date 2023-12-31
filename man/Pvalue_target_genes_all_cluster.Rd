% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Pvalue_target_genes_all_cluster.R
\name{Pvalue_target_genes_all_cluster}
\alias{Pvalue_target_genes_all_cluster}
\title{Compute the pvalue of target genes in all cluster}
\usage{
Pvalue_target_genes_all_cluster(
  exprMatr,
  cell_cluster,
  target_genes,
  DEmethod = "MAST",
  H0_cells = NULL,
  unreliable_genes_all_clusters,
  N_cores = 1,
  backend = "doParallel"
)
}
\arguments{
\item{exprMatr}{The input expression matrix.}

\item{target_genes}{The input target genes.}

\item{cell_clusters}{The clusters obtained by cell cluster.}
}
\description{
Compute the pvalue of target genes in all cluster
}
