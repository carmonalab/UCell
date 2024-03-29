% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{sample.matrix}
\alias{sample.matrix}
\title{Sample dataset to test UCell installation}
\format{
A sparse matrix of 30 cells and 20729 genes.
}
\source{
\url{https://doi.org/10.1016/j.cell.2021.04.048}
}
\usage{
sample.matrix
}
\description{
A sparse matrix (class "dgCMatrix") of single-cell transcriptomes
(scRNA-seq) for 30 cells and 20729 genes.
Single-cell UMI counts were normalized using a standard log-normalization:
counts for each cell were divided by the total counts for that cell and
multiplied by 10,000, then natural-log transformed using \code{log1p}. \cr\cr
This a subsample of T cells from the
large scRNA-seq PBMC dataset published
by \href{https://doi.org/10.1016/j.cell.2021.04.048}{Hao et al.} and
available as UMI counts at
\url{https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat}
}
\keyword{datasets}
