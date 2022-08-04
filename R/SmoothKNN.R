#' Smooth signature scores by KNN
#'
#' This function performs smoothing of single-cell scores by weighted
#' average of the k-nearest neighbors. It can be useful to 'impute' scores by
#' neighboring cells and partially correct data sparsity. While this function
#' has been designed to smooth UCell scores, it can be applied to any numerical
#' metadata contained in SingleCellExperiment or Seurat objects
#'
#' @param obj Input object - either a [SingleCellExperiment] object
#'     or a [Seurat] object.
#' @param reduction Which dimensionality reduction to use for kNN smoothing.
#'     It must be already present in the input object.
#' @param k Number of neighbors for kNN smoothing
#' @param BNPARAM A [BiocNeighborParam] object specifying the algorithm to use
#'     for kNN calculation.
#' @param signature.names The names of the signatures (or any numeric metadata
#'     column) for which to calculate kNN-smoothed scores
#' @param suffix The suffix to append to metadata columns for the new
#'     knn-smoothed scores      
#' @examples
#' # Run UCell
#' gene.sets <- list(Tcell = c("CD2","CD3E","CD3D"),
#'                 Myeloid = c("SPI1","FCER1G","CSF1R"))
#' data(sample.matrix)
#' obj <- Seurat::CreateSeuratObject(sample.matrix)                
#' 
#' obj <- AddModuleScore_UCell(obj,features = gene.sets, name=NULL)
#' # Run PCA
#' obj <- RunPCA(obj)
#' # Smooth signatures
#' obj <- SmoothKNN(obj, reduction="pca", signature.names=names(gene.sets))
#' head(obj[[]])
#'
#' @importFrom methods setMethod setGeneric
#' @import BiocNeighbors
#' @export SmoothKNN
SmoothKNN <- function(
    obj=NULL,
    reduction="pca",
    k=10,
    BNPARAM=AnnoyParam(),
    signature.names=NULL,
    suffix="_KNN") {
  
  setGeneric("SmoothKNN")
}

setMethod("SmoothKNN", signature(obj="ANY"), SmoothKNN_unsupported)
setMethod("SmoothKNN", signature(obj="SingleCellExperiment"), SmoothKNN_sce)
setMethod("SmoothKNN", signature(obj="Seurat"), SmoothKNN_Seurat)
