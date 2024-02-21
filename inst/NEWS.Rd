\name{NEWS}
\title{News for package \pkg{UCell}}

\section{New in version 2.7.4}{
\itemize{
    \item Add support for Seurat v5 assay and datasets in multiple layers. Added dependency on Seurat >= 5.0.
    \item Change default for \code{chunk.size} to 100. On parallelized jobs, this can improve up to 2-fold execution time.
    \item Smaller test dataset (30 cells), to speed up package function checks.  
}}

\section{New in version 2.1.2}{
\itemize{
    \item New function \code{SmoothKNN()} for k-nearest neighbor smoothing of UCell scores. It can be applied both on SingleCellExperiment and Seurat objects (S3 method).
    \item Add two new vignettes: along with basic usage (vignette 1), there are now dedicated vignettes for running UCell with SingleCellExperiment objects (vignette 2) and Seurat objects (vignette 3). kNN smoothing is illustrated for both object types.
    \item Fixing a bug that prevented storing of feature ranks.
   
}}

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