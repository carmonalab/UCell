#' Calculate Mann Whitney U from a vector of ranks
#' 
#' @param rank_value A vector of ranks
#' @param maxRank Max number of features to include in ranking
#' @param sparse Whether the vector of ranks is in sparse format
#' 
#' @return Normalized AUC (as U statistic) for the vector
u_stat <- function(rank_value, maxRank=1000, sparse=FALSE){
    
    if (sparse==TRUE){
        rank_value[rank_value==0] <- maxRank
    }
    
    insig <- rank_value > maxRank
    if (all(insig)) {
        return(0L)
    } else {
        rank_value[insig] <- maxRank
        rank_sum <- sum(rank_value)
        len_sig <- length(rank_value)
        lfac <- len_sig*(len_sig + 1)/2
        u_value <- rank_sum - lfac
        auc <- 1 - u_value/(len_sig*maxRank - lfac)
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
    if (is.vector(u_matrix)) {  # Case of ncells=1
      u_matrix <- t(as.matrix(u_matrix))
    }
    
    rownames(u_matrix) <- colnames(ranks_matrix)[-1]
    return (u_matrix)
}

#' Calculate rankings and scores for query data and given signature set
#' 
#' @param   matrix        Input data matrix 
#' @param   features      List of signatures
#' @param   maxRank       Rank cutoff (1500) 
#' @param   chunk.size    Cells per sub-matrix (100)
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
        matrix, features,  maxRank=1500, chunk.size=100,
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
                cells_rankings[cells_rankings>=maxRank] <- 0
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
rankings2Uscore <- function(ranks_matrix, features, chunk.size=100, w_neg=1,
                            BPPARAM = NULL,ncores=1, force.gc=FALSE,
                            name="_UCell") {
    
    #Check if all genes in signatures are present in the stored signatures
    ranks_matrix <- check_genes(ranks_matrix, features)
    
    #Weight on neg signatures must be >=0
    if (is.null(w_neg)) {w_neg <- 1}
    if (!is.numeric(w_neg) | w_neg<0) {
        stop("Weight on negative signatures (w_neg) must be >=0")}
    
    maxRank <- max(ranks_matrix)+1
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
split_data.matrix <- function(matrix, chunk.size=100) {
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
        split.data[[i]] <- matrix[,min:max,drop=FALSE]
        min <- max+1    #for next chunk
    }
    return(split.data)
}

#' Smoothing scores by KNN
#'
#' @param   matrix Input data matrix 
#' @param   nn    A nearest neighbor object returned by
#'   [BiocNeighbors::findKNN]
#' @param   decay Exponential decay for nearest neighbor weight: (1-decay)^n
#' @param   up.only If set to TRUE, smoothed scores will only be
#'     allowed to increase by smoothing
#' 
#' @return  A dataframe of knn-smoothed scores
knn_smooth_scores <- function(
    matrix=NULL,
    nn=NULL,
    decay=0.1,   #decay must be bound between 0 and 1
    up.only=FALSE #scores can only increase
) {
  
  sig.cols <- colnames(matrix)
  
  w.df <- vapply(sig.cols, FUN.VALUE=numeric(nrow(matrix)), FUN=function(s) {
    ss.scores <- matrix[,s]
    weighted.scores <- vapply(X = seq_len(nrow(nn$index)),
                              FUN.VALUE = numeric(1),
                              FUN = function(x) {
                                r <- nn$index[x,]
                                r <- c(x,r)
                                i <- seq(0, length(r)-1)
                                w <- (1-decay)**i
                                sum(w * ss.scores[r])/sum(w)
                              })
    if (up.only) {
      pmax(weighted.scores, ss.scores)
    } else {
      weighted.scores
    }
  })
  rownames(w.df) <- rownames(matrix)
  as.data.frame(w.df)
}  




#' @rdname SmoothKNN
#' @method SmoothKNN Seurat
#' @export
SmoothKNN.Seurat <- function(
    obj=NULL,
    signature.names=NULL,
    reduction="pca",
    k=10,
    decay=0.1,
    up.only=FALSE,
    BNPARAM=AnnoyParam(),
    BPPARAM=SerialParam(),
    suffix="_kNN",
    assay=NULL,
    slot="data",
    sce.expname=NULL,
    sce.assay=NULL
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
  
  if (is.null(assay)) {  # Work on metadata
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
  } else {  # Work directly on features
    exp <- Seurat::GetAssayData(obj, layer=slot, assay=assay)
    feats <- rownames(exp)
    found <- intersect(signature.names, feats)
    notfound <- setdiff(signature.names, found)
    
    if (length(found)==0) {
      stop("Could not find any of the given features in this object")
    }
    if (length(notfound)>0) {
      nf <- paste(notfound, collapse=",")
      mess <- sprintf("The following features were not found in assay %s:\n* %s",
                      assay, nf)
      warning(mess, immediate.=TRUE, call.=FALSE, noBreaks.=TRUE)
    }
    m <- t(exp[found, , drop=FALSE])
  }
  ncells <- ncol(obj)
  
  if (decay<0 | decay>1) {
    stop("decay parameter must be a number between 0 and 1")
  }
    
  if (k<=0) {  #this behavior disables kNN smoothing
    k=1
    decay=1
  }
  
  if (ncells <= k) {
    k <- ncells-1
    warning("'k' capped at the number of observations minus 1")
  }
  
  if (ncells>1) {
    # Find kNNs
    space <- Seurat::Embeddings(obj, reduction=reduction)
    nn <- findKNN(space, k=k, BNPARAM=BNPARAM, BPPARAM=BPPARAM)
    
    # Do smoothing
    smooth.df <- knn_smooth_scores(matrix=m, nn=nn,
                                   decay=decay, up.only=up.only)  
  } else {
    smooth.df <- m
  }
  
  if (is.null(assay)) {  #metadata
    colnames(smooth.df) <- paste0(colnames(smooth.df), suffix)
    obj <- Seurat::AddMetaData(obj, metadata = smooth.df)
  } else {  #new assay
    nas <- paste0(assay, suffix)
    obj[[nas]] <- Seurat::CreateAssayObject(data=t(smooth.df))
  }
  return(obj)
}

#' @rdname SmoothKNN
#' @method SmoothKNN SingleCellExperiment
#' @importFrom stats setNames
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment assay assays SummarizedExperiment assayNames
#' @export
SmoothKNN.SingleCellExperiment <- function(
    obj=NULL,
    signature.names=NULL,
    reduction="PCA",
    k=10,
    decay=0.1,
    up.only=FALSE,
    BNPARAM=AnnoyParam(),
    BPPARAM=SerialParam(),
    suffix="_kNN",
    assay=NULL,
    slot="data",
    sce.expname=c("UCell","main"),
    sce.assay=NULL
) {
  
  if (!reduction %in% reducedDimNames(obj)) {
    stop(sprintf("Could not find reduction %s in this object", reduction))
  }
  
  if (is.null(signature.names)) {
    stop("Please provide the metadata column names that you want to smooth")
  }
  
  sce.expname <- sce.expname[1]
  if (!sce.expname %in% c(altExpNames(obj), "main")) {
    stop(sprintf("Cannot find summarized experiment name: %s", sce.expname))
  } 
 
  if (sce.expname == "main") {
     exp <- obj
  } else {
     exp <- altExp(obj, sce.expname)
  }
  found <- intersect(signature.names, rownames(exp))
  notfound <- setdiff(signature.names, found)
  
  if (length(found)==0) {
    stop("Could not find any of the given signatures in this object")
  }
  if (length(notfound)>0) {
    nf <- paste(notfound, collapse=",")
    mess <- sprintf("The following signatures were not found:\n* %s",nf)
    warning(mess, immediate.=TRUE, call.=FALSE, noBreaks.=TRUE)
  }
  
  if (is.null(sce.assay)) {
    sce.assay <- 1
  } else if (!sce.assay %in% assayNames(obj)) {
    stop(sprintf("Cannot find assay: %s", sce.assay))
  }
  m <- SummarizedExperiment::assay(exp, sce.assay)
  m <- t(m[found, ,drop=FALSE])
  ncells <- nrow(m)
  
  if (decay<0 | decay>1) {
    stop("decay parameter must be a number between 0 and 1")
  }
  
  if (k<=0) {  #this behavior disables kNN smoothing
    k=1
    decay=1
  }
  
  if (ncells <= k) {
    k <- ncells-1
    warning("'k' capped at the number of observations minus 1")
  }
  
  if (ncells>1) {
    # Find kNNs
    space <- reducedDim(obj, reduction)
    nn <- findKNN(space, k=k, BNPARAM=BNPARAM, BPPARAM=BPPARAM)
    # Do smoothing
    m.smooth <- knn_smooth_scores(matrix=m, nn=nn,
                                  decay=decay, up.only=up.only) 
  } else {
    m.smooth <- m
  }
  
  #New experiment with smoothed scores
  colnames(m.smooth) <- paste0(colnames(m.smooth), suffix)
  sce.newexp <- paste0(sce.expname, suffix)
  
  l <- list("sce" = t(m.smooth))
  SingleCellExperiment::altExp(obj, sce.newexp) <- 
    SummarizedExperiment(assays = setNames(l, sce.newexp))

  return(obj)
}
