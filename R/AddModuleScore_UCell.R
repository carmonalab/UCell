#' Calculate module enrichment scores from single-cell data (Seurat interface)
#'
#' Given a Seurat object, calculates module/signature enrichment scores at single-cell level using the Mann-Whitney U statistic.
#' UCell scores are normalized U statistics (between 0 and 1), and they are mathematically related to the Area under the ROC curve (see [Mason and Graham](https://doi.org/10.1256/003590002320603584))
#' 
#' In contrast to Seurat's AddModuleScore, which is normalized by binning genes of similar expression at the population level, UCell scores depend 
#' only on the gene expression ranks of individual cell, and therefore they are robust across datasets regardless of dataset composition.
#'
#' @param object Seurat object
#' @param features A list of signatures, for example: \code{list( Tcell_signature = c("CD2","CD3E","CD3D"), Myeloid_signature = c("SPI1","FCER1G","CSF1R"))}
#'     You can also specify positive and negative gene sets by adding a + or - sign to genes in the signature; see an example below
#' @param chunk.size Number of cells to be processed simultaneously (lower size requires slightly more computation but reduces memory demands)
#' @param maxRank Maximum number of genes to rank per cell; above this rank, a given gene is considered as not expressed.
#' @param ncores Number of processors to parallelize computation. Requires package \code{future}
#' @param storeRanks Store ranks matrix in Seurat object (@misc slot) for fast subsequent computations. This option will demand large amounts of RAM.
#' @param w_neg Weight on negative genes in signature. e.g. `w_neg=1` weighs equally up- and down-regulated genes, `w_neg=0.5` gives 50% less importance to negative genes
#' @param assay Pull out data from this assay of the Seurat object (if NULL, use \code{DefaultAssay(obj)})
#' @param slot Pull out data from this slot of the Seurat object
#' @param ties.method How ranking ties should be resolved (passed on to [data.table::frank])
#' @param force.gc Explicitly call garbage collector to reduce memory footprint
#' @param seed Integer seed for [future.apply::future_lapply] parallel execution
#' @param name Name tag that will be appended at the end of each signature name, "_UCell" by default (e.g. signature score in meta data will be named: Myeloid_signature_UCell)
#' @return Returns a Seurat object with module/signature enrichment scores added to object meta data; each score is stored as the corresponding signature name provided in \code{features} followed by the tag given in \code{name} (or "_UCell" by default )
#' @examples
#' ## Not run:
#' library(UCell)
#' gene.sets <- list(Tcell_signature = c("CD2","CD3E","CD3D"),
#'                 Myeloid_signature = c("SPI1","FCER1G","CSF1R"))
#' SeuratObject <- AddModuleScore_UCell(SeuratObject,features = gene.sets)
#' SeuratObject$Tcell_signature_UCell
#' head(SeuratObject@meta.data)
#' 
#' ## Using positive and negative gene sets
#' markers <- list()
#' markers$Tcell_gd <- c("TRDC+", "TRGC1+", "TRGC2+", "TRDV1+","TRAC-","TRBC1-","TRBC2-")
#' markers$Tcell_NK <- c("FGFBP2+", "SPON2+", "KLRF1+", "FCGR3A+", "CD3E-","CD3G-")
#' SeuratObject <- AddModuleScore_UCell(SeuratObject, features = markers)
#' FeaturePlot(SeuratObject, features=c("Tcell_gd_UCell","Tcell_NK_UCell"))
#' ## End (Not run)
#' @export
AddModuleScore_UCell <- function(obj, features, maxRank=1500, chunk.size=1000, ncores=1, storeRanks=F, w_neg=1,
                                 assay=NULL, slot="data", ties.method="average", force.gc=FALSE, seed=123, name="_UCell") {
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Function 'AddModuleScore_UCell' requires the Seurat package. Please install it.", call. = FALSE)
  }  
  
  features <- check_signature_names(features)
  
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(obj)
  }
  precomputedRanks <- obj@misc[["UCell"]][[assay]][["cells_rankings"]]
  
  #If rank matrix was pre-computed, evaluate the new signatures from these ranks
  #Else, calculate new ranks to score signatures (optionally storing ranks, takes up memory but become very fast to evaluate further signatures)
  if (!is.null(precomputedRanks)) {
    meta.list <- rankings2Uscore(precomputedRanks, features=features, chunk.size=chunk.size, w_neg=w_neg,
                                 ncores=ncores, force.gc=force.gc, name=name)
    
  } else {
    meta.list <- calculate_Uscore(Seurat::GetAssayData(obj, slot, assay=assay), features=features, maxRank=maxRank, chunk.size=chunk.size, w_neg=w_neg,
                                  ncores=ncores, ties.method=ties.method, force.gc=force.gc, storeRanks=storeRanks, name=name)
    
    #store ranks matrix?
    if (storeRanks==T){
      cells_rankings.merge <- lapply(meta.list,function(x) rbind(x[["cells_rankings"]]))
      cells_rankings.merge <- Reduce(cbind, cells_rankings.merge)
      
      obj@misc[["UCell"]][[assay]] <- list(cells_rankings=cells_rankings.merge)
    }
  }
  
  meta.merge <- lapply(meta.list,function(x) rbind(x[["cells_AUC"]]))
  meta.merge <- Reduce(rbind, meta.merge)
  
  obj <- Seurat::AddMetaData(obj, as.data.frame(meta.merge))
  
  return(obj)
}