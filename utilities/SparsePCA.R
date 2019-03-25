# Author: Shristi Pandey
# The following set of functions prepare the data and run sparse PCA 

#========================================
  # Parameters:
  #   pc.compute: The number of sparse PCs to compute 
  #   num.nonneg: Number of non-zero components in each PC. This is essentially the sparsity parameter of SparsePCA
  # ======================================
# Returns:
#   A Seurat Object with the PCA slot updated with the results of running sparse PCA
# ======================================

set.ifnull = function(x,y) {
  if(is.null(x)) x=y
  return(x)
}

# Finds the set of sparse components that can optimally reconstruct the gene expression data. 
# The amount of sparsity is controlled the number of non-zero components
RunSparsePCA <- function (object, pc.genes = NULL, pcs.compute = 20, use.imputed = FALSE, 
                          rev.pca = FALSE, weight.by.var = TRUE, do.print = TRUE, pcs.print = 1:5, 
                          genes.print = 30, reduction.name = "pca", reduction.key = "PC", 
                          assay.type = "RNA", num.nonneg = 200, ...) 
{
  data.use <- PrepDR2(object = object, genes.use = pc.genes, 
                     use.imputed = use.imputed, assay.type = assay.type)

    pcs.compute <- min(pcs.compute, nrow(x = data.use) - 1)
    pca.results <- irlba::ssvd(x = t(x = data.use), k = pcs.compute, n=num.nonneg, maxit = 200)
    gene.loadings <- pca.results$v
    sdev <- diag(pca.results$d)/sqrt(max(1, ncol(data.use) - 1))
    if (weight.by.var) {
      cell.embeddings <- pca.results$u %*% pca.results$d
    }
    else {
      cell.embeddings <- pca.results$u
    }

  rownames(x = gene.loadings) <- rownames(x = data.use)
  colnames(x = gene.loadings) <- paste0(reduction.key, 1:pcs.compute)
  rownames(x = cell.embeddings) <- colnames(x = data.use)
  colnames(x = cell.embeddings) <- colnames(x = gene.loadings)
  pca.obj <- new(Class = "dim.reduction", gene.loadings = gene.loadings, 
                 cell.embeddings = cell.embeddings, sdev = sdev, key = reduction.key)
  eval(expr = parse(text = paste0("object@dr$", reduction.name, 
                                  "<- pca.obj")))
  parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("RunPCA"))]
  object <- SetCalcParams2(object = object, calculation = "RunPCA", 
                          ... = parameters.to.store)
  if (is.null(object@calc.params$RunPCA$pc.genes)) {
    object@calc.params$RunPCA$pc.genes <- rownames(data.use)
  }
  if (do.print) {
    PrintDim(object = object, dims.print = pcs.print, genes.print = genes.print, 
             reduction.type = reduction.name)
  }
  return(object)
}

#Takes a seurat object and returns a gene expression matrix of only the most variable genes across the dataset. 
#The output is then used to run Sparse PCA
PrepDR2 <- function (object, genes.use = NULL, use.imputed = FALSE, assay.type = "RNA") 
{
  if (length(object@var.genes) == 0 && is.null(x = genes.use)) {
    stop("Variable genes haven't been set. Run MeanVarPlot() or provide a vector\n          of genes names in genes.use and retry.")
  }
  if (use.imputed) {
    data.use <- t(x = scale(x = t(x = object@imputed)))
  }
  else {
    data.use <- GetAssayData(object, assay.type = assay.type, 
                             slot = "scale.data")
  }
  genes.use <- set.ifnull(x = genes.use,  object@var.genes)
  genes.use <- unique(x = genes.use[genes.use %in% rownames(x = data.use)])
  genes.var <- apply(X = data.use[genes.use, ], MARGIN = 1, 
                     FUN = var)
  genes.use <- genes.use[genes.var > 0]
  genes.use <- genes.use[!is.na(x = genes.use)]
  data.use <- data.use[genes.use, ]
  return(data.use)
}

SetCalcParams2 <- function (object, calculation, time = TRUE, ...) 
{
  object@calc.params[calculation] <- list(...)
  object@calc.params[[calculation]]$object <- NULL
  object@calc.params[[calculation]]$object2 <- NULL
  if (time) {
    object@calc.params[[calculation]]$time <- Sys.time()
  }
  return(object)
}

