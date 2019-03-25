# Author: Shristi Pandey
#================================
# This will run a basic version of Seurat taking user input along the way for parameters
# User input parameters:
#     y.cutoff = threshold for variable gene detection (Seurat Method)
#     diffCV.cutoff = threshold for selection of variable genes (Pandey et al., method)
#     num.pc = The number of PCs to be used for clustering. 

# Returns: A seurat object where the cells have been clustered 
#         A matrix of containing a ranked list of markers and the associated statistics for all classes
#================================



RunBasicSeurat <- function(raw.data, project.name , num.pc, clustering.res){
  source("C:/Users/ShristiPandey/Dropbox/SourceCodeRepos/VariableGeneSelection/VarGenes.R")
  seuratObject <- CreateSeuratObject(raw.data = raw.data, project = project.name, min.cells = 3, min.genes = 300)
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = seuratObject@data), value = TRUE)
  percent.mito <- Matrix::colSums(seuratObject@raw.data[mito.genes, ])/Matrix::colSums(seuratObject@raw.data)
  seuratObject <- AddMetaData(object = seuratObject, metadata = percent.mito, col.name = "percent.mito")
  VlnPlot(object = seuratObject, features.plot = c("nUMI", "percent.mito"), nCol = 3)
  seuratObject <- ScaleData(object = seuratObject, vars.to.regress = c("nUMI", "percent.mito"))
  seuratObject <- FilterCells(object = seuratObject, subset.names = c("nGene"), low.thresholds = c(200), high.thresholds = c(2500))
  seuratObject <- NormalizeData(object = seuratObject, normalization.method = "LogNormalize", scale.factor = 10000)
  
  #select variable genes 
  varGenes <- new.meanCVfit(count.data = seuratObject@raw.data,  diffCV.cutoff = 0.4, mean.min = 0.005, mean.max = 100, main.use = "", do.plot = T)
  Cutoff.val  <- readline(prompt = "Enter revised value for diffcv.cutoff: ")
  varGenes <- new.meanCVfit(count.data = seuratObject@raw.data,  diffCV.cutoff = Cutoff.val, mean.min = 0.005, mean.max = 100, main.use = "", do.plot = T)
  seuratObject <- FindVariableGenes(seuratObject, y.cutoff = 0.5, do.plot = T)
  yCutoff <- readline(prompt = "Enter a revised value for y.cutoff based on the plot: ")
  seuratObject <- FindVariableGenes(seuratObject, y.cutoff = yCutoff, do.plot = T)
  
  totalVarGenes <- union(seuratObject@varGenes, varGenes)
  seuratObject@var.genes <- totalVarGenes
  
  seuratObject <- RunPCA(object = seuratObject, pc.genes = seuratObject@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 40)
  PCAPlot(seuratObject, 1, 2)
  PCElbowPlot(seuratObject, 40)
  num.pc  <- readline(prompt = "Enter the number of PCs you want to use for clustering: ")
  num.pc <- as.integer(num.pc)
  seuratObject <- RunTSNE(seuratObject, dims.use = 1:num.pc, do.fast = TRUE)
  seuratObject <- FindClusters(seuratObject, reduction.type = "pca", dims.use = 1:num.pc, resolution = 2, print.output = T, save.SNN = T, force.recalc = TRUE)
  markers = FindAllMarkers(seuratObject, test.use = 'bimod', only.pos = T, min.pct = 0.25, do.print = T)
  return(seuratObject, markers)
}
