
#Calculate AUC as Mann–Whitney U statistic from a vector of ranks
u_stat <- function(rank_value, maxRank=1000, sparse=F){

  if(sparse==T){
    rank_value[rank_value==0] <- maxRank+1
  }

  insig <- rank_value > maxRank
  if(all(insig)) {
    return(0L)
  } else {
    rank_value[insig] <- maxRank+1
    rank_sum = sum(rank_value)
    len_sig <- length(rank_value)
    u_value = rank_sum - (len_sig * (len_sig + 1))/2
    auc = u_value/(len_sig * maxRank)
    auc = 1 - auc
    return(auc)
  }
}

#Check if all genes in signatures are found in data matrix - otherwise add zero counts in data-matrix to complete it
check_genes <- function(matrix, features) {
   features <- unlist(features)
   missing <- setdiff(features, rownames(matrix))
   ll <- length(missing)
   
   if (ll/length(features) > 0.5) {
      warning(sprintf("Over half of genes (%s%%) in specified signatures are missing from data. Check the integrity of your dataset\n", round(100*ll/length(features))))
   }
   
   if (ll>0) {
      add.mat <- Matrix::sparseMatrix(length(missing), ncol(matrix))
      rownames(add.mat) <- missing
      matrix <- rbind(matrix, add.mat)
      
      missing.concatenate <- paste(missing, collapse=",")
      warning(sprintf("The following genes were not found and will be imputed to exp=0:\n* %s",missing.concatenate))
   }
   return(matrix)
}

#Check signature names and add standard names is missing

check_signature_names <- function(features) {
  defaultSigName <- paste0(rep("signature_",length(features)),seq_along(features))
  if(is.null(names(features))){
    names(features) <- defaultSigName
  } else {
    invalidNames <- names(features) == "" | duplicated(names(features))
    names(features)[invalidNames] <- defaultSigName[invalidNames]
  }
  return(features)
}

#Calculate AUC for a list of signatures, from a ranks matrix
u_stat_signature_list <- function(sig_list, ranks_matrix, maxRank=1000, sparse=F) {

  u_matrix <- sapply(sig_list, function(sig) {
    as.numeric(ranks_matrix[sig, lapply(.SD, function(x) u_stat(x,maxRank = maxRank,sparse=sparse)),.SDcols=-1, on="rn"])
    })

  rownames(u_matrix) <- colnames(ranks_matrix)[-1]
  return (u_matrix)
}

# Calculate features' ranks from expression data matrices
data_to_ranks_data_table = function(data, ties.method="average") {
  dt <- as.data.table(as.matrix(data))
  rnaDT.ranks.dt <- dt[, lapply(.SD, function(x) frankv(x,ties.method=ties.method,order=c(-1L)))]
  rnaDT.ranks.rownames <- rownames(data)
  rnaDT.ranks.dt.rn <- cbind(rn=rnaDT.ranks.rownames, rnaDT.ranks.dt)
  setkey(rnaDT.ranks.dt.rn, rn, physical = F)
  return(rnaDT.ranks.dt.rn)
}

#split data matrix into cell chunks
split_data.matrix <- function(matrix, chunk.size=1000) {
  ncols <- dim(matrix)[2]
  nchunks <- (ncols-1) %/% chunk.size + 1

  split.data <- list()
  for (i in 1:nchunks) {
    min <- 1 + (i-1)*chunk.size
    max <- min(i*chunk.size, ncols)
    split.data[[i]] <- matrix[,min:max]
  }
  return(split.data)
}

#Get signature scores from precomputed rank matrix

rankings2Uscore <- function(ranks_matrix, features, chunk.size=1000, 
                            ncores=1, force.gc=FALSE, name="_UCell") {

  #Check if all genes in signatures are present in the stored signatures
  ranks_matrix <- check_genes(ranks_matrix, features)

  maxRank <- max(ranks_matrix)
  split.data <- split_data.matrix(matrix=ranks_matrix, chunk.size=chunk.size)
  rm(ranks_matrix)

  if (ncores>1) {
    plan(future::multisession(workers=future_param_ncores))

    meta.list <- future_lapply(
      X = split.data,
      FUN = function(x) {

        dense <- as.matrix(x)
        dense <- as.data.table(dense, keep.rownames=TRUE)
        setkey(dense, rn, physical=F)

        cells_AUC <- u_stat_signature_list(features, dense, maxRank=maxRank, sparse=T)
        colnames(cells_AUC) <- paste0(colnames(cells_AUC),name)
        
        if (force.gc) {
          dense <- NULL
          gc()
        }
        return(list(cells_AUC=cells_AUC))
      },
      future.seed = future_param_seed
    )
    plan(strategy = "sequential")

  } else {

    meta.list <- lapply(
      X = split.data,
      FUN = function(x) {

        dense <- as.matrix(x)
        dense <- as.data.table(dense, keep.rownames=TRUE)
        setkey(dense, rn, physical=F)

        cells_AUC <- u_stat_signature_list(features, dense, maxRank=maxRank, sparse=T)
        colnames(cells_AUC) <- paste0(colnames(cells_AUC),name)
        
        if (force.gc) {
          dense <- NULL
          gc()
        }

        return(list(cells_AUC=cells_AUC))
      } )
  }
  return(meta.list)
}

#Calculate rankings and scores for query data and given signature set
calculate_Uscore <- function(matrix, features,  maxRank=1500, chunk.size=1000, ncores=1, 
                             ties.method="average", storeRanks=FALSE, force.gc=FALSE, name="_UCell") {

  #Make sure we have a sparse matrix
  require(Matrix)
  if (class(matrix) != "dgCMatrix") {
    matrix <- Matrix::Matrix(as.matrix(matrix),sparse = T)
  }
  #Check if all genes in signatures are present in the data matrix
  matrix <- check_genes(matrix, features)

  #Split into manageable chunks
  split.data <- split_data.matrix(matrix=matrix, chunk.size=chunk.size)

  #Parallelize?
  if (ncores>1) {
    plan(future::multisession(workers=future_param_ncores))

    meta.list <- future_lapply(
      X = split.data,
      FUN = function(x) {

        cells_rankings <- data_to_ranks_data_table(x, ties.method = ties.method)
        cells_AUC <- u_stat_signature_list(features, cells_rankings, maxRank=maxRank, sparse=F)

        colnames(cells_AUC) <- paste0(colnames(cells_AUC),name)

        if (storeRanks==T){
          gene.names <- as.character(as.matrix(cells_rankings[,1]))
          #make sparse
          cells_rankings[cells_rankings>maxRank] <- 0
          ranks.sparse <- Matrix::Matrix(as.matrix(cells_rankings[,-1]),sparse = T)
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
      },
      future.seed = future_param_seed
    )
    plan(strategy = "sequential")

  } else {
    meta.list <- lapply(
      X = split.data,
      FUN = function(x) {
        cells_rankings <- data_to_ranks_data_table(x, ties.method = ties.method)
        cells_AUC <- u_stat_signature_list(features, cells_rankings, maxRank=maxRank, sparse=F)
        colnames(cells_AUC) <- paste0(colnames(cells_AUC),name)

        if (storeRanks==T){
          gene.names <- as.character(as.matrix(cells_rankings[,1]))
          #make sparse
          cells_rankings[cells_rankings>maxRank] <- 0
          ranks.sparse <- Matrix::Matrix(as.matrix(cells_rankings[,-1]),sparse = T)
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

      } )
  }
  return(meta.list)
}
