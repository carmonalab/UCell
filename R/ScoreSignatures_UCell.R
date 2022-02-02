#' Calculate module enrichment scores from single-cell data
#'
#' Given a gene vs. cell matrix, calculates module/signature enrichment scores on single-cell level using Mann-Whitney U statistic.
#' UCell scores are normalized U statistics (between 0 and 1), and they are mathematically related to the Area under the ROC curve (see [Mason and Graham](https://doi.org/10.1256/003590002320603584))
#'
#' These scores only depend on the gene expression ranks of individual cell, and therefore they are robust across datasets regardless of dataset composition.
#'
#' @param matrix A gene vs. cell data matrix, either in sparse or dense format. Leave empty if providing a rank matrix with \code{precalc.ranks}
#' @param features A list of signatures, for example: \code{list( Tcell_signature = c("CD2","CD3E","CD3D"), Myeloid_signature = c("SPI1","FCER1G","CSF1R"))}
#'     You can also specify positive and negative gene sets by adding a + or - sign to genes in the signature; see an example below
#' @param precalc.ranks A sparse matrix of pre-calculated ranks, obtained with \code{\link{StoreRankings_UCell}}
#' @param maxRank Maximum number of genes to rank per cell; above this rank, a given gene is considered as not expressed.
#'     Note: this parameter is ignored if \code{precalc.ranks} are specified
#' @param chunk.size Number of cells to be processed simultaneously (lower size requires slightly more computation but reduces memory demands)
#' @param w_neg Weight on negative genes in signature. e.g. `w_neg=1` weighs equally up- and down-regulated genes, `w_neg=0.5` gives 50% less importance to negative genes
#' @param ncores Number of processors to parallelize computation. Requires package \code{future}
#' @param ties.method How ranking ties should be resolved (passed on to [data.table::frank])
#' @param name Name suffix appended to signature names
#' @param force.gc Explicitly call garbage collector to reduce memory footprint
#' @param seed Integer seed for [future.apply::future_lapply] parallel execution
#' @return Returns a dataframe of signature scores for each cell
#' @examples
#' ## Not run:
#' library(UCell)
#' my.matrix <- UCell::sample.matrix
#' gene.sets <- list( Tcell_signature = c("CD2","CD3E","CD3D"),
#'                  Myeloid_signature = c("SPI1","FCER1G","CSF1R"))
#' scores <- ScoreSignatures_UCell(my.matrix, features=gene.sets)
#' scores[1:5,]
#' ## End (Not run)
#' @import Matrix
#' @export
ScoreSignatures_UCell <- function(matrix=NULL, features, precalc.ranks=NULL, maxRank=1500, w_neg=1, name="_UCell",
                                  chunk.size=1000, ncores=1, ties.method="average", force.gc=FALSE, seed=123) {
  
  features <- check_signature_names(features)
  
  if (!is.null(precalc.ranks)) {
     meta.list <- rankings2Uscore(precalc.ranks, features=features, chunk.size=chunk.size, w_neg=w_neg,
                                  ncores=ncores, force.gc=force.gc, name=name)
  } else {
     meta.list <- calculate_Uscore(matrix, features=features, maxRank=maxRank, chunk.size=chunk.size, w_neg=w_neg,
                                   ties.method=ties.method, ncores=ncores, force.gc=force.gc, name=name)
  }
  meta.merge <- lapply(meta.list,function(x) rbind(x[["cells_AUC"]]))
  meta.merge <- Reduce(rbind, meta.merge)
  
  return(meta.merge)

}
