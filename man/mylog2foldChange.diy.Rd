% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mylog2foldChange.R
\name{mylog2foldChange.diy}
\alias{mylog2foldChange.diy}
\title{calculate foldchange for each celltype}
\usage{
mylog2foldChange.diy(inData, cell.type, method = "median", probs = 0.75)
}
\arguments{
\item{inData}{a dataframe of gene expression}

\item{cell.type}{the cell type which you want to calculate the foldchange}

\item{method}{default is "median",eg "median" or "mean"}

\item{probs}{the quantile for median calculate}
}
\value{
the foldchange value of each celltype in \code{cellwave object}
}
\description{
calculate foldchange for each celltype
}
