#' Smooth signature scores by kNN
#'
#' This function performs smoothing of single-cell scores by weighted
#' average of the k-nearest neighbors. It can be useful to 'impute' scores by
#' neighboring cells and partially correct data sparsity. While this function
#' has been designed to smooth UCell scores, it can be applied to any numerical
#' metadata contained in SingleCellExperiment or Seurat objects
#'
#' @param obj Input object - either a [SingleCellExperiment] object
#'     or a Seurat object.
#' @param signature.names The names of the signatures (or any numeric metadata
#'     column) for which to calculate kNN-smoothed scores
#' @param reduction Which dimensionality reduction to use for kNN smoothing.
#'     It must be already present in the input object.
#' @param k Number of neighbors for kNN smoothing
#' @param BNPARAM A [BiocNeighborParam] object specifying the algorithm to use
#'     for kNN calculation.
#' @param suffix The suffix to append to metadata columns for the new
#'     knn-smoothed scores  
#' @param sce.expname For sce objects only - which experiment stores the
#'    signature       
#' @examples
#' # Run UCell
#' library(Seurat)
#' gene.sets <- list(Tcell = c("CD2","CD3E","CD3D"),
#'                 Myeloid = c("SPI1","FCER1G","CSF1R"))
#' data(sample.matrix)
#' obj <- Seurat::CreateSeuratObject(sample.matrix)                
#' 
#' obj <- AddModuleScore_UCell(obj,features = gene.sets, name=NULL)
#' # Run PCA
#' obj <- FindVariableFeatures(obj) |> ScaleData() |> RunPCA()
#' # Smooth signatures
#' obj <- SmoothKNN(obj, reduction="pca", signature.names=names(gene.sets))
#' head(obj[[]])
#'
#' @importFrom methods setMethod setGeneric is
#' @import BiocNeighbors
#' @export SmoothKNN
SmoothKNN <- function(
    obj=NULL,
    signature.names=NULL,
    reduction="pca",
    k=10,
    BNPARAM=AnnoyParam(),
    suffix="_kNN",
    sce.expname="UCell") {
  
  if (methods::is(obj, "SingleCellExperiment")) {
    SmoothKNN_sce(obj=obj, signature.names = signature.names,
                  reduction=reduction, k=k, BNPARAM=BNPARAM,
                  suffix=suffix, sce.expname=sce.expname)
  }
  else if (methods::is(obj, "Seurat")) {
    SmoothKNN_Seurat(obj=obj, signature.names = signature.names,
                  reduction=reduction, k=k, BNPARAM=BNPARAM,
                  suffix=suffix)
  } else {
    stop("Unsupported format. Please provide a 'sce' or 'Seurat' object")
  }
}

