\name{NEWS}
\title{News for package \pkg{UCell}}

\section{Changes in version 1.3.1}{
\itemize{
    \item Restructure code to conform to BioC standards.
    \item Switch from \code{future} to \code{BiocParallel} to parallelize jobs.
    \item Add support for SingleCellExperiment - new function \code{ScoreSignatures_UCell_sce()} interacts directly with sce objects.
    \item Signatures cannot be larger than \code{maxRank} parameter.
    \item Do not rank more genes (\code{maxRank}) than there are in the input matrix.
}}