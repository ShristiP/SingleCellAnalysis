# Author: Shristi Pandey
# This function explores all sparsity parameters for PCA and all clustering resolution to generate the best clustering results out of a single cell transcriptomes
#
# ======================================
# Input: A seurat Object
# Parameters:
#   num.nonneg: Number of non-zero components in each PC. This is essentially the sparsity parameter of SparsePCA
#   nonZeroFactors: a vector of all the sparsity paramaters to try. 
#   ClusterRes = a vector of all the clustering resolutions to try. 

# ======================================
# Returns:
#   a data.frame of all the clustering results. 
#   
# ======================================

load("/Users/ShristiPandey/Dropbox/SingleCellAnalysis/ForPaper/FinalObjects/LarvaSs2With10XIdent.RObj")
load("/Users/ShristiPandey/Dropbox/SingleCellAnalysis/ForPaper/FinalObjects/FinalLarvaSs2_justgng8L.RObj")

# Calculate sparse PCA across many non-negative component numbers
non.zero.use <- c(500, 250, 100, 80, 60, 40, 30, 20, 10, 5)

sparse.pca.results <- lapply(object, non.zero.use, function(nn) {
  print(paste0(Sys.time(), ": Starting ", nn, " non-negative components."))
  LarvaSs2_final <- RunSparsePCA(object, pc.genes = object@var.genes, pcs.compute = 30, use.imputed = FALSE, num.nonneg = nn, do.print = F)
  return(object@dr$pca)
})
names(sparse.pca.results) <- non.zero.use

#Perform Clustering across all parameters
res = c(0.8, 1.0, 1.2, 1.4, 1.6, 1.8)
all.clustering.results. = ClusterAcrossAllParameters(nonZeroFactors = non.zero.use, ClusterRes = res)

#Perform clustering with all parameters and store the results in a dataframe 
#=============================================================================

ClusterAcrossAllParameters <- function(object, nonZeroFactors, ClusterRes){
  clustering.results <- matrix(0, nrow = length(nonZeroFactors), ncol = length(ClusterRes))
  for (i in 1:length(nonZeroFactors)){
    currNonZero <- nonZeroFactors[i]
    object@dr$pca <- sparse.pca.results[[as.character(nonZeroFactors[i])]]
    for (j in 1: length(ClusterRes)){
      currRes = ClusterRes[j]
      print(currRes)
      object_final_1 <- FindClusters(object, dims.use = 1:30, reduction.type = "pca", force.recalc = T, resolution = currRes)
      numClusters <- length(levels(object@ident))
      print(numClusters)
      #ari <- adjustedRandIndex(a, LarvaSs2_final_1@ident)
      #print(ari)
      clustering.results[i,j] <- numClusters
      #clustering.results[i,j] <- ari
      print(clustering.results)
    }
  }
  clustering.results <- as.data.frame(clustering.results, row.names = as.character(nonZeroFactors))
  return(clustering.results)
}

