#' Calculate Mann Whitney U from a vector of ranks
#' 
#' @param rank_value A vector of ranks
#' @param maxRank Max number of features to include in ranking
#' @param sparse Whether the vector of ranks is in sparse format
#' 
#' @return Normalized AUC (as U statistic) for the vector
u_stat <- function(rank_value, maxRank=1000, sparse=FALSE){
    
    if (sparse==TRUE){
        rank_value[rank_value==0] <- maxRank+1
    }
    
    insig <- rank_value > maxRank
    if (all(insig)) {
        return(0L)
    } else {
        rank_value[insig] <- maxRank+1
        rank_sum <- sum(rank_value)
        len_sig <- length(rank_value)
        
        u_value <- rank_sum - (len_sig * (len_sig + 1))/2
        auc <- 1 - u_value/(len_sig * maxRank)
        return(auc)
    }
}

#' Calculate U scores for a list of signatures, given a rank matrix
#' 
#' @param   sig_list      A list of signatures
#' @param   ranks_matrix  Matrix of pre-computed ranks
#' @param   maxRank       Max number of features to include in ranking,
#'     for u_stat function
#' @param   sparse        Whether the vector of ranks is in sparse format
#' @param   w_neg         Weight on negative signatures
#' 
#' @return A matrix of U scores
#' @import data.table
#' 
u_stat_signature_list <- function(sig_list, ranks_matrix, maxRank=1000,
    sparse=FALSE, w_neg=1) {
    
    dim <- ncol(ranks_matrix)-1
    u_matrix <- vapply(sig_list, FUN.VALUE = numeric(dim), FUN=function(sig) {
        sig_neg <- grep('-$', unlist(sig), perl=TRUE, value=TRUE)
        sig_pos <- setdiff(unlist(sig), sig_neg)
        
        if (length(sig_pos)>0) {
            sig_pos <- gsub('\\+$','',sig_pos,perl=TRUE)
            u_p <- as.numeric(ranks_matrix[
                sig_pos,
                lapply(.SD, function(x)
                    u_stat(x,maxRank = maxRank,sparse=sparse)),
                .SDcols=-1, on="rn"])
        } else {
            u_p <- rep(0, dim(ranks_matrix)[2]-1)
        }
        if (length(sig_neg)>0) {
            sig_neg <- gsub('-$','',sig_neg,perl=TRUE)
            u_n <- as.numeric(ranks_matrix[
                sig_neg,
                lapply(.SD, function(x)
                    u_stat(x,maxRank = maxRank,sparse=sparse)),
                .SDcols=-1, on="rn"])
        } else {
            u_n <- rep(0, dim(ranks_matrix)[2]-1)
        }
        
        diff <- u_p - w_neg*u_n   #Subtract negative sets, if any
        diff[diff<0] <- 0
        return(diff)
    })
    
    rownames(u_matrix) <- colnames(ranks_matrix)[-1]
    return (u_matrix)
}

#' Calculate rankings and scores for query data and given signature set
#' 
#' @param   matrix        Input data matrix 
#' @param   features      List of signatures
#' @param   maxRank       Rank cutoff (1500) 
#' @param   chunk.size    Cells per sub-matrix (1000)
#' @param   BPPARAM       A BioParallel object to instruct UCell how to
#'    parallelize  
#' @param   ncores        Number of cores to use for parallelization
#' @param   w_neg         Weight on negative signatures
#' @param   ties.method   How to break ties, for data.table::frankv
#'     method ("average")
#' @param   storeRanks    Store ranks? (FALSE) 
#' @param   force.gc      Force garbage collection? (FALSE) 
#' @param   name          Suffix for metadata columns ("_UCell") 
#' 
#' @return  A list of signature scores
#' @importFrom methods is 
#' @import  Matrix
#' @import  BiocParallel
calculate_Uscore <- function(
        matrix, features,  maxRank=1500, chunk.size=1000,
        BPPARAM = NULL, ncores=1, w_neg=1, ties.method="average",
        storeRanks=FALSE, force.gc=FALSE, name="_UCell"){
    
    #Make sure we have a sparse matrix
    if (!methods::is(matrix, "dgCMatrix")) {
        matrix <- Matrix::Matrix(as.matrix(matrix),sparse = TRUE)
    }
    
    #Check if all genes in signatures are present in the data matrix
    matrix <- check_genes(matrix, features)
    
    #Do not evaluate more genes than there are
    if (!is.numeric(maxRank)) {
        stop("Rank cutoff (maxRank) must be a number")
    }
    if (maxRank > nrow(matrix)) {
        maxRank <- nrow(matrix)
    }
    
    #Weight on neg signatures must be >=0
    if (is.null(w_neg)) {w_neg <- 1}
    if (w_neg<0) {stop("Weight on negative signatures (w_neg) must be >=0")}
    
    #Signatures cannot be larger than maxRank parameter
    sign.lgt <- lapply(features, length)
    if (any(sign.lgt > maxRank)) {
        stop("One or more signatures contain more genes than maxRank parameter.
            Increase maxRank parameter or make shorter signatures")
    }
    
    #Split into manageable chunks
    split.data <- split_data.matrix(matrix=matrix, chunk.size=chunk.size)
    
    #Either take a BPPARAM object, or make one on the spot using 'ncores'
    if (is.null(BPPARAM)) {
        BPPARAM <- BiocParallel::MulticoreParam(workers=ncores)
    }
    meta.list <- BiocParallel::bplapply(
        X = split.data, 
        BPPARAM =  BPPARAM,
        FUN = function(x) {
            cells_rankings <- data_to_ranks_data_table(x,
                ties.method=ties.method)
            cells_AUC <- u_stat_signature_list(features, cells_rankings, 
                maxRank=maxRank, sparse=FALSE,
                w_neg=w_neg)
            colnames(cells_AUC) <- paste0(colnames(cells_AUC),name)
            if (storeRanks==TRUE){
                gene.names <- as.character(as.matrix(cells_rankings[,1]))
                #make sparse
                cells_rankings[cells_rankings>maxRank] <- 0
                ranks.sparse <- Matrix::Matrix(as.matrix(
                    cells_rankings[,-1]),sparse = TRUE)
                dimnames(ranks.sparse)[[1]] <- gene.names
                if (force.gc) {
                    cells_rankings <- NULL
                    gc()
                }
                return(list(cells_rankings=ranks.sparse, cells_AUC=cells_AUC))
            } else {
                if (force.gc) {
                    cells_rankings <- NULL
                    gc()
                }
                return(list(cells_AUC=cells_AUC))
            }
            
        })
    return(meta.list)
}

#' Get signature scores from pre-computed rank matrix
#' 
#' @param ranks_matrix  A rank matrix
#' @param features      List of signatures
#' @param chunk.size    How many cells per matrix chunk
#' @param w_neg         Weight on negative signatures
#' @param BPPARAM       A BioParallel object to instruct UCell how to
#'     parallelize  
#' @param ncores        How many cores to use for parallelization?
#' @param force.gc      Force garbage collection to recover RAM? (FALSE)
#' @param name          Name suffix for metadata columns ("_UCell")
#' 
#' @return                    A list of signature scores
#' @import    data.table
rankings2Uscore <- function(ranks_matrix, features, chunk.size=1000, w_neg=1,
                            BPPARAM = NULL,ncores=1, force.gc=FALSE,
                            name="_UCell") {
    
    #Check if all genes in signatures are present in the stored signatures
    ranks_matrix <- check_genes(ranks_matrix, features)
    
    #Weight on neg signatures must be >=0
    if (is.null(w_neg)) {w_neg <- 1}
    if (!is.numeric(w_neg) | w_neg<0) {
        stop("Weight on negative signatures (w_neg) must be >=0")}
    
    maxRank <- max(ranks_matrix)
    split.data <- split_data.matrix(matrix=ranks_matrix, chunk.size=chunk.size)
    rm(ranks_matrix)
    
    #Either take a BPPARAM object, or make one on the spot using 'ncores'
    if (is.null(BPPARAM)) {
        BPPARAM <- BiocParallel::MulticoreParam(workers=ncores)
    }
    meta.list <- BiocParallel::bplapply(
        X = split.data, 
        BPPARAM =  BPPARAM,
        FUN = function(x) {
            
            dense <- as.matrix(x)
            dense <- as.data.table(dense, keep.rownames=TRUE)
            setkey(dense, "rn", physical=FALSE)
            
            cells_AUC <- u_stat_signature_list(features, dense,
                maxRank=maxRank, sparse=TRUE,
                w_neg=w_neg)
            colnames(cells_AUC) <- paste0(colnames(cells_AUC),name)
            
            if (force.gc) {
                dense <- NULL
                gc()
            }
            return(list(cells_AUC=cells_AUC))
        }
    )
    return(meta.list)
}

#' Check genes
#'
#' Check if all genes in signatures are found in data matrix - otherwise
#' add zero counts in data-matrix to complete it
#' 
#' @param matrix Input data matrix
#' @param features List of genes that must be present
#'     (otherwise they are added)
#' 
#' @return Same input matrix, extended to comprise any missing genes
check_genes <- function(matrix, features) {
    features <- unlist(features)
    features <- gsub("[-+]$","",features,perl=TRUE)
    missing <- setdiff(features, rownames(matrix))
    ll <- length(missing)
    
    if (ll/length(features) > 0.5) {
        n <- round(100*ll/length(features))
        mess <- sprintf("Over half of genes (%s%%)", n)
        mess <- paste(mess, "in specified signatures are missing from data.",
            "Check the integrity of your dataset.") 
        warning(mess, immediate.=TRUE, call.=FALSE, noBreaks.=TRUE)
    }
    
    if (ll>0) {
        dim1 <- length(missing)
        dim2 <- ncol(matrix)
        add.mat <-  Matrix::Matrix(data=min(matrix),
            nrow=dim1, ncol=dim2, sparse = TRUE)

        rownames(add.mat) <- missing
        matrix <- rbind(matrix, add.mat)
        
        missing.concatenate <- paste(missing, collapse=",")
        mess <- sprintf("The following genes were not found and will be
                        imputed to exp=0:\n* %s",missing.concatenate)
        warning(mess, immediate.=TRUE, call.=FALSE, noBreaks.=TRUE)
    }
    return(matrix)
}

#' Check signature names and add standard names is missing
#' 
#' @param features List of signatures for scoring
#' 
#' @return The input list of signatures, with standard names if
#'     provided un-named

check_signature_names <- function(features) {
    defaultSigName <- paste0(rep("signature_",length(features)),
        seq_along(features))
    if(is.null(names(features))){
        names(features) <- defaultSigName
    } else {
        invalidNames <- names(features) == "" | duplicated(names(features))
        names(features)[invalidNames] <- defaultSigName[invalidNames]
    }
    return(features)
}


#' Calculate per-cell feature rankings
#' 
#' @param data              Expression data matrix 
#' @param ties.method       How to break ties (passed on to data.table::frankv)
#' 
#' @return                  A data.table of ranks 
#' @import data.table
data_to_ranks_data_table <- function(data, ties.method="average") {
    dt <- as.data.table(as.matrix(data))
    rnaDT.ranks.dt <- dt[, lapply(.SD, function(x)
        frankv(x,ties.method=ties.method,order=c(-1L)))]
    rnaDT.ranks.rownames <- rownames(data)
    rnaDT.ranks.dt.rn <- cbind(rn=rnaDT.ranks.rownames, rnaDT.ranks.dt)
    setkey(rnaDT.ranks.dt.rn, "rn", physical = FALSE)
    return(rnaDT.ranks.dt.rn)
}

#' Split data matrix into smaller sub-matrices ('chunks')
#'
#' @param   matrix      Input data matrix 
#' @param   chunk.size  How many cells to include in each sub-matrix
#' 
#' @return  A list of sub-matrices, each with size {n_features x chunk_size}
split_data.matrix <- function(matrix, chunk.size=1000) {
    ncols <- dim(matrix)[2]
    nchunks <- (ncols-1) %/% chunk.size + 1
    
    split.data <- list()
    min <- 1
    for (i in seq_len(nchunks)) {
        if (i == nchunks-1) {  #make last two chunks of equal size
            left <- ncols-(i-1)*chunk.size
            max <- min+round(left/2)-1
        } else {
            max <- min(i*chunk.size, ncols)
        }
        split.data[[i]] <- matrix[,min:max]
        min <- max+1    #for next chunk
    }
    return(split.data)
}

#' Smoothing scores by KNN
#'
#' @param   matrix  Input data matrix 
#' @param   nn      A nearest neighbor object returned by
#'   [BiocNeighbors::findKNN]
#' 
#' @return  A dataframe of knn-smoothed scores
knn_smooth_scores <- function(
    matrix=NULL,
    nn=NULL
) {
  
  sig.cols <- colnames(matrix)
  w.df <- vapply(sig.cols, FUN.VALUE=numeric(nrow(matrix)), FUN=function(s) {
    
    ss.scores <- matrix[,s]
    weighted.scores <- vapply(X = 1:nrow(nn$index),
                              FUN.VALUE = numeric(1),
                              FUN = function(x) {
                                r <- nn$index[x,]
                                r <- c(x,r)
                                
                                d <- nn$distance[x,]
                                d <- c(d[1],d)
                                
                                w <- 1/(0.01+d)
                                
                                sum(w * ss.scores[r])/sum(w)
                              })
  })
  rownames(w.df) <- rownames(matrix)
  as.data.frame(w.df)
}  


#' Smooth signature scores by kNN (Seurat)
#'
#' This function performs smoothing of single-cell scores by weighted
#' average of the k-nearest neighbors. It can be useful to 'impute' scores by
#' neighboring cells and partially correct data sparsity. While this function
#' has been designed to smooth UCell scores, it can be applied to any numerical
#' metadata contained in Seurat objects
#'
#' @param obj Input Seurat object.
#' @param signature.names The names of the signatures (or any numeric metadata
#'     column) for which to calculate kNN-smoothed scores
#' @param reduction Which dimensionality reduction to use for kNN smoothing.
#'     It must be already present in the input object.
#' @param k Number of neighbors for kNN smoothing
#' @param BNPARAM A [BiocNeighborParam] object specifying the algorithm to use
#'     for kNN calculation.
#' @param suffix The suffix to append to metadata columns for the new
#'     knn-smoothed scores  
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
#' @import BiocNeighbors
SmoothKNN_Seurat <- function(
    obj=NULL,
    signature.names=NULL,
    reduction="pca",
    k=10,
    BNPARAM=AnnoyParam(),
    suffix="_kNN"
) {
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Function 'SmoothKNN_UCell' requires the Seurat package.
            Please install it.", call. = FALSE)
  } 
  
  if (!reduction %in% Seurat::Reductions(obj)) {
    stop(sprintf("Could not find reduction %s in this object", reduction))
  }
  
  if (is.null(signature.names)) {
    stop("Please provide the metadata column names that you want to smooth")
  }
  found <- intersect(signature.names, colnames(obj[[]]))
  notfound <- setdiff(signature.names, found)
  
  if (length(found)==0) {
    stop("Could not find any of the given signatures in this object")
  }
  if (length(notfound)>0) {
    nf <- paste(notfound, collapse=",")
    mess <- sprintf("The following signature were found in metadata:\n* %s",nf)
    warning(mess, immediate.=TRUE, call.=FALSE, noBreaks.=TRUE)
  }
  
  m <- obj[[found]]
  
  # Find kNNs
  space <- Seurat::Embeddings(obj, reduction=reduction)
  nn <- findKNN(space, k=k, BNPARAM=BNPARAM)
  
  # Do smoothing
  smooth.df <- knn_smooth_scores(matrix=m, nn=nn)  
  
  colnames(smooth.df) <- paste0(colnames(smooth.df), suffix)
  obj <- Seurat::AddMetaData(obj, metadata = smooth.df)
  return(obj)
}

#' Smooth signature scores by kNN (SingleCellExperiment)
#'
#' This function performs smoothing of single-cell scores by weighted
#' average of the k-nearest neighbors. It can be useful to 'impute' scores by
#' neighboring cells and partially correct data sparsity. While this function
#' has been designed to smooth UCell scores, it can be applied to any numerical
#' metadata contained in a SingleCellExperiment object
#'
#' @param obj Input object stored in [SingleCellExperiment]
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
#' library(SingleCellExperiment)
#' library(scater)
#' data(sample.matrix)
#' sce <- SingleCellExperiment(list(counts=sample.matrix))
#' gene.sets <- list( Tcell_signature = c("CD2","CD3E","CD3D"),
#'                   Myeloid_signature = c("SPI1","FCER1G","CSF1R"))
#' sce <- ScoreSignatures_UCell(sce, features=gene.sets, name=NULL)
#' altExp(sce, 'UCell')
#' # Run PCA
#' sce <- logNormCounts(sce)
#' sce <- runPCA(sce, scale=TRUE, ncomponents=20)
#' # Smooth signatures
#' sce <- SmoothKNN(sce, reduction="PCA", signature.names=names(gene.sets))
#' altExp(sce, 'UCell_kNN')
#' 
#' @import BiocNeighbors
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment assay assays SummarizedExperiment
SmoothKNN_sce <- function(
    obj=NULL,
    signature.names=NULL,
    reduction="PCA",
    k=10,
    BNPARAM=AnnoyParam(),
    suffix="_kNN",
    sce.expname="UCell"
) {
  
  if (! reduction %in% reducedDimNames(obj)) {
    stop(sprintf("Could not find reduction %s in this object", reduction))
  }
  
  if (is.null(signature.names)) {
    stop("Please provide the metadata column names that you want to smooth")
  }
  
  if (!sce.expname %in% altExpNames(obj)) {
    stop(sprintf("Cannot find summarized experiment name: %s", sce.expname))
  } 
  
  exp <- altExp(obj, sce.expname)
  
  found <- intersect(signature.names, rownames(exp))
  notfound <- setdiff(signature.names, found)
  
  if (length(found)==0) {
    stop("Could not find any of the given signatures in this object")
  }
  if (length(notfound)>0) {
    nf <- paste(notfound, collapse=",")
    mess <- sprintf("The following signature were found:\n* %s",nf)
    warning(mess, immediate.=TRUE, call.=FALSE, noBreaks.=TRUE)
  }
  
  m <- as.data.frame(SummarizedExperiment::assay(exp))
  m <- t(m[found, ,drop=FALSE])
  
  # Find kNNs
  space <- reducedDim(obj, reduction)
  nn <- findKNN(space, k=k, BNPARAM=BNPARAM)
  
  # Do smoothing
  smooth.df <- knn_smooth_scores(matrix=m, nn=nn)  
  
  new.altExp <- paste0(sce.expname, suffix)
  
  altExp(obj, new.altExp) <-
    SummarizedExperiment(assays = list(new.altExp = t(smooth.df)))
  
  return(obj)
}
