#' Calculate module enrichment scores from single-cell data
#'
#' Given a gene vs. cell matrix, calculates module/signature enrichment scores
#' on single-cell level using Mann-Whitney U statistic.
#' UCell scores are normalized U statistics (between 0 and 1), and they are
#' mathematically related to the Area under the ROC curve (see
#' [Mason and Graham](https://doi.org/10.1256/003590002320603584))
#' These scores only depend on the gene expression ranks of individual cell, 
#' and therefore they are robust across datasets regardless of dataset
#' composition.
#'
#' @param matrix Input matrix, either stored in a [SingleCellExperiment] object
#'     or as a raw matrix. \code{dgCMatrix} format supported.
#' @param features A list of signatures, for example:
#'     \code{list(Tcell_signature = c("CD2","CD3E","CD3D"),
#'     Myeloid_signature = c("SPI1","FCER1G","CSF1R"))}
#'     You can also specify positive and negative gene sets by adding a + or - 
#'     sign to genes in the signature; see an example below
#' @param precalc.ranks If you have pre-calculated ranks using
#'     \code{\link{StoreRankings_UCell}}, you can specify the pre-calculated
#'     ranks instead of the gene vs. cell matrix.
#' @param maxRank Maximum number of genes to rank per cell; above this rank, a
#'     given gene is considered as not expressed. Note: this parameter is 
#'     ignored if \code{precalc.ranks} are specified
#' @param assay The sce object assay where the data is to be found
#' @param chunk.size Number of cells to be processed simultaneously (lower size
#'     requires slightly more computation but reduces memory demands)
#' @param w_neg Weight on negative genes in signature. e.g. `w_neg=1` weighs
#'     equally up- and down-regulated genes, `w_neg=0.5` gives 50% less
#'     importance to negative genes
#' @param BPPARAM A [BiocParallel::bpparam()] object that tells UCell
#'     how to parallelize. If provided, it overrides the `ncores` parameter.     
#' @param ncores Number of processors to parallelize computation. If
#'     \code{BPPARAM = NULL}, the function uses
#'     \code{BiocParallel::MulticoreParam(workers=ncores)}
#' @param ties.method How ranking ties should be resolved - passed on to
#'     [data.table::frank]
#' @param missing_genes How to handle missing genes in matrix:
#'     "impute": impute expression to zero; "skip": remove missing
#'     genes from signature 
#' @param name Name suffix appended to signature names
#' @param force.gc Explicitly call garbage collector to reduce memory footprint
#' @return Returns input SingleCellExperiment object with UCell scores
#'     added to altExp
#' @examples
#' library(UCell)
#' # Using sparse matrix
#' data(sample.matrix)
#' gene.sets <- list( Tcell_signature = c("CD2","CD3E","CD3D"),
#'                  Myeloid_signature = c("SPI1","FCER1G","CSF1R"))
#' scores <- ScoreSignatures_UCell(sample.matrix, features=gene.sets)
#' head(scores)
#' 
#' # Using sce object
#' library(SingleCellExperiment)
#' data(sample.matrix)
#' my.sce <- SingleCellExperiment(list(counts=sample.matrix))
#' gene.sets <- list( Tcell_signature = c("CD2","CD3E","CD3D"),
#'                  Myeloid_signature = c("SPI1","FCER1G","CSF1R"))
#' my.sce <- ScoreSignatures_UCell(my.sce, features=gene.sets)
#' altExp(my.sce, 'UCell')
#'
#' @importFrom methods is 
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment assay assays SummarizedExperiment
#' @import Matrix
#' @export
ScoreSignatures_UCell <- function(
        matrix=NULL,features, precalc.ranks=NULL,
        maxRank=1500, w_neg=1, name="_UCell",
        assay="counts", chunk.size=100,
        missing_genes = c("impute","skip"),
        BPPARAM=NULL, ncores=1,
        ties.method="average", force.gc=FALSE) {
    
    features <- check_signature_names(features)
    missing_genes <- match.arg(missing_genes)
    
    #Check type of input
    if (methods::is(matrix, "SingleCellExperiment")) { # sce object
        if (!assay %in% names(SummarizedExperiment::assays(matrix))) {
            stop(sprintf("Assay %s not found in sce object.", assay))
        }
        m <- SummarizedExperiment::assay(matrix, assay) 
    } else if (methods::is(matrix, "matrix") | #matrix or DF
        methods::is(matrix, "dgCMatrix") |
        methods::is(matrix, "data.frame")) { 
            m <- matrix
    } else {
        m <- NULL
    }
    
    if (is.null(m) & is.null(precalc.ranks)) {
        stop("Unrecognized input format.")
    }
    
    #Run on pre-calculated ranks ('m' can be NULL)
    if (!is.null(precalc.ranks)) {
        u.list <- rankings2Uscore(precalc.ranks, features=features,
            chunk.size=chunk.size,w_neg=w_neg,
            ncores=ncores, BPPARAM=BPPARAM,
            missing_genes=missing_genes,
            force.gc=force.gc, name=name)
    } else {
        u.list <- calculate_Uscore(m, features=features, maxRank=maxRank,
            chunk.size=chunk.size, w_neg=w_neg,
            ties.method=ties.method, ncores=ncores,
            missing_genes=missing_genes,
            BPPARAM=BPPARAM, force.gc=force.gc, name=name)
    }
    u.merge <- lapply(u.list,function(x) rbind(x[["cells_U"]]))
    u.merge <- Reduce(rbind, u.merge)
    
    if (methods::is(matrix, "SingleCellExperiment")) {
        SingleCellExperiment::altExp(matrix, "UCell") <- 
            SummarizedExperiment(assays = list("UCell" = t(u.merge)))
        return(matrix)
    } else {
        return(u.merge)
    }
}
