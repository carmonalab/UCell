\name{NEWS}
\title{News for package \pkg{UCell}}

\section{New in version 2.0.0}{
\itemize{
    \item Update code to pass all BioC checks.
    \item The function \code{ScoreSignatures_UCell()} and \code{StoreRankings_UCell()} accept directly \code{sce} objects.
    \item Takes custom BiocParallel::bpparam() object as input to specify parallelisation.
}}


\section{New in version 1.3.1}{
\itemize{
    \item Restructure code to conform to BioC standards.
    \item Switch from \code{future} to \code{BiocParallel} to parallelize jobs.
    \item Add support for SingleCellExperiment - new function \code{ScoreSignatures_UCell_sce()} interacts directly with sce objects.
    \item Signatures cannot be larger than \code{maxRank} parameter.
    \item Do not rank more genes (\code{maxRank}) than there are in the input matrix.
}}