# Author: Shristi Pandey
#===========================
# This function finds markers that are enriched in a clusters by the area 
# under the precision-recall curve (AUCPR). AUCPR is a quantitative measure of “cluster-specificity” of a marker 
# by balancing both recall (the sensitivity of marker gene detection within the cluster of interest) 
# and precision (accuracy of the quantitative levels of gene expression as a predictor of the correct cell type). 
#========================================
# Input: Seurat Object 
# Parameters:
#   thresh.use: 
#   min.pct: only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations
#   min.diff.pct: only test genes that show a minimum difference in the fraction of detection between the two groups. Set to -Inf by default
#   only.pos: Only return positive markers (FALSE by default)
#   max.cells.per.ident Down sample each identity class to a max number.
# ======================================
# Returns:
#   Matrix containing a ranked list of putative markers, and associated statistics (p-values, AUCPR)
# ======================================


FindMarkersAUCPR = function (object, ident.1, ident.2 = NULL, genes.use = NULL, 
          thresh.use = 0.25, min.pct = 0.1, min.diff.pct = -Inf, 
          print.bar = TRUE, only.pos = FALSE, max.cells.per.ident = Inf, 
          random.seed = 1, latent.vars = "nUMI", min.cells = 3) 
  {
  library(pbapply)
  genes.use <- set.ifnull(x = genes.use, y = rownames(x = object@data))
  if (max.cells.per.ident < Inf) {
    object <- SubsetData(object = object, max.cells.per.ident = max.cells.per.ident, 
                         random.seed = random.seed)
  }
  if (length(x = as.vector(x = ident.1) > 1) && any(as.character(x = ident.1) %in% 
                                                    object@cell.names)) {
    cells.1 <- intersect(x = ident.1, y = object@cell.names)
  }
  else {
    cells.1 <- WhichCells(object = object, ident = ident.1)
  }
  if (length(x = as.vector(x = ident.2) > 1) && any(as.character(x = ident.2) %in% 
                                                    object@cell.names)) {
    cells.2 <- intersect(x = ident.2, y = object@cell.names)
  }
  else {
    if (is.null(x = ident.2)) {
      cells.2 <- object@cell.names
    }
    else {
      cells.2 <- WhichCells(object = object, ident = ident.2)
    }
  }
  cells.2 <- setdiff(x = cells.2, y = cells.1)
  if (length(x = cells.1) == 0) {
    print(paste("Cell group 1 is empty - no cells with identity class", 
                ident.1))
    return(NULL)
  }
  if (length(x = cells.2) == 0) {
    print(paste("Cell group 2 is empty - no cells with identity class", 
                ident.2))
    return(NULL)
  }
  thresh.min <- 0
  data.temp1 <- round(x = apply(X = object@data[genes.use, 
                                                cells.1, drop = F], MARGIN = 1, FUN = function(x) {
                                                  return(sum(x > thresh.min)/length(x = x))
                                                }), digits = 3)
  data.temp2 <- round(x = apply(X = object@data[genes.use, 
                                                cells.2, drop = F], MARGIN = 1, FUN = function(x) {
                                                  return(sum(x > thresh.min)/length(x = x))
                                                }), digits = 3)
  data.alpha <- cbind(data.temp1, data.temp2)
  colnames(x = data.alpha) <- c("pct.1", "pct.2")
  alpha.min <- apply(X = data.alpha, MARGIN = 1, FUN = max)
  names(x = alpha.min) <- rownames(x = data.alpha)
  genes.use <- names(x = which(x = alpha.min > min.pct))
  alpha.diff <- alpha.min - apply(X = data.alpha, MARGIN = 1, 
                                  FUN = min)
  genes.use <- names(x = which(x = alpha.min > min.pct & alpha.diff > 
                                 min.diff.pct))
  if (length(cells.1) < min.cells) {
    stop(paste("Cell group 1 has fewer than", as.character(min.cells), 
               "cells in identity class", ident.1))
  }
  if (length(cells.2) < min.cells) {
    stop(paste("Cell group 2 has fewer than", as.character(min.cells), 
               " cells in identity class", ident.2))
  }
  data.1 <- apply(X = object@data[genes.use, cells.1, drop = F], 
                  MARGIN = 1, FUN = ExpMean)
  data.2 <- apply(X = object@data[genes.use, cells.2, drop = F], 
                  MARGIN = 1, FUN = ExpMean)
  total.diff <- (data.1 - data.2)
  genes.diff <- names(x = which(x = abs(x = total.diff) > thresh.use))
  genes.use <- intersect(x = genes.use, y = genes.diff)
 
  
  to.return <- MarkerTestPR(object = object, cells.1 = cells.1, 
                            cells.2 = cells.2, genes.use = genes.use, print.bar = print.bar)
  
  to.return[, "avg_diff"] <- total.diff[rownames(x = to.return)]
  to.return <- cbind(to.return, data.alpha[rownames(x = to.return), 
                                           ])
  
  to.return <- to.return[order(-to.return$myAUC, -to.return$avg_diff), 
                           ]

  if (only.pos) {
    to.return <- subset(x = to.return, subset = avg_diff > 
                          0)
  }
  return(to.return)
}



MarkerTestPR = function (object, cells.1, cells.2, genes.use = NULL, print.bar = TRUE) 
{
  genes.use <- set.ifnull(x = genes.use, y = rownames(x = object@data))
  to.return <- AUCPRMarkerTest(data1 = object@data[, cells.1], 
                             data2 = object@data[, cells.2], mygenes = genes.use, 
                             print.bar = print.bar)
  return(to.return)
}

AUCPRMarkerTest = function (data1, data2, mygenes, print.bar = TRUE) 
{
  myAUC <- unlist(x = lapply(X = mygenes, FUN = function(x) {
    return(DifferentialAUCPR(x = as.numeric(x = data1[x, ]), 
                           y = as.numeric(x = data2[x, ])))
  }))
  myAUC[is.na(x = myAUC)] <- 0
  if (print.bar) {
    iterate.fxn <- pblapply
  }
  else {
    iterate.fxn <- lapply
  }
  avg_diff <- unlist(x = iterate.fxn(X = mygenes, FUN = function(x) {
    return(ExpMean(x = as.numeric(x = data1[x, ])) - ExpMean(x = as.numeric(x = data2[x, 
                                                                                      ])))
  }))
  toRet <- data.frame(cbind(myAUC, avg_diff), row.names = mygenes)
  toRet <- toRet[rev(x = order(toRet$myAUC)), ]
  return(toRet)
}

DifferentialAUCPR = function (x, y) 
{
  prediction.use <- ROCR::prediction(predictions = c(x, y), labels = c(rep(x = 1, 
                                                                     length(x = x)), rep(x = 0, length(x = y))), label.ordering = 0:1)
  
  perf1 <- ROCR::performance(prediction.use, "prec", "rec")
  
  is.nanvals = is.nan(perf1@x.values[[1]]) | is.nan(perf1@y.values[[1]])
  require(caTools)
  auc = trapz(perf1@x.values[[1]][!is.nanvals], 
        perf1@y.values[[1]][!is.nanvals])
  #plot(perf1@x.values[[1]], perf1@y.values[[2]])
  auc.use <- round(x = auc, digits = 3)
  return(auc.use)
}