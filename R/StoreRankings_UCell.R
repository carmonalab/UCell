#' Calculate and store gene rankings for a single-cell dataset
#'
#' Given a gene vs. cell matrix, calculates the rankings of expression for all genes in each cell. 
#' 
#' While \code{\link{ScoreSignatures_UCell}} can be used 'on the fly' to evaluate signatures in a query dataset, it requires recalculating gene
#' ranks at every execution. If you have a large dataset and plan to experiment with multiple signatures, evaluating the same dataset multiple times,
#' this function allows you to store pre-calculated ranks so they do not have to be recomputed every time. Pre-calculated ranks can then be applied to the
#' function \code{\link{ScoreSignatures_UCell}} to evaluate gene signatures in a significantly faster way on successive iterations.
#'
#' @param matrix Input matrix, either stored in a [SingleCellExperiment] object or as a raw matrix. \code{dgCMatrix} format supported.
#' @param maxRank Maximum number of genes to rank per cell; above this rank, a given gene is considered as not expressed
#' @param assay Assay where the data is to be found (for input in 'sce' format)
#' @param chunk.size Number of cells to be processed simultaneously (lower size requires slightly more computation but reduces memory demands)
#' @param ncores Number of processors to parallelize computation
#' @param ties.method How ranking ties should be resolved (passed on to [data.table::frank])
#' @param force.gc Explicitly call garbage collector to reduce memory footprint
#' @param seed Integer seed
#' @return Returns a sparse matrix of pre-calculated ranks that can be used multiple times to evaluate different signatures
#' @examples
#' ## Not run:
#' data(sample.matrix)
#' ranks <- StoreRankings_UCell(sample.matrix)
#' ranks[1:5,1:5]
#' gene.sets <- list( Tcell_signature = c("CD2","CD3E","CD3D"),
#'                  Myeloid_signature = c("SPI1","FCER1G","CSF1R"))
#' scores <- ScoreSignatures_UCell(features=gene.sets, precalc.ranks=ranks)
#' ## End (Not run)
#' @importFrom methods is 
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment assay
#' @import Matrix
#' @export
StoreRankings_UCell <- function(matrix, maxRank=1500, chunk.size=1000,
                                ncores=1, assay='counts', ties.method="average",
                                force.gc=FALSE, seed=123) {
    
    #Check type of input
    if (methods::is(matrix, "SingleCellExperiment")) { # sce object
        if (!assay %in% names(matrix@assays)) {
            stop(sprintf("Assay %s not found in sce object.", assay))
        }
        m <- SummarizedExperiment::assay(matrix, assay)
    } else if (methods::is(matrix, "matrix") | #matrix or DF
               methods::is(matrix, "dgCMatrix") |
               methods::is(matrix, "data.frame")) { 
        m <- matrix
    } else {
        stop("Unrecognized input format.")
    }
    
    features <- rownames(m)[1]  #placeholder signature
    meta.list <- calculate_Uscore(m, features=features, maxRank=maxRank,
                                  chunk.size=chunk.size, ncores=ncores,
                                  ties.method=ties.method, storeRanks=TRUE,
                                  force.gc=force.gc)
    
    ranks.all <- lapply(meta.list,function(x) rbind(x[["cells_rankings"]]))
    ranks.all <- Reduce(cbind, ranks.all)
    
    return(ranks.all)
    
}