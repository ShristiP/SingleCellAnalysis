# Function to plot the relative gene expression patterns among single cell clusters in the form of a dot plot 

set.ifnull=function(x,y) {
  if(is.null(x)) x=y
  return(x)
}

nogrid = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# max.val.exp 
# max.val.perc


Perc.pos.by.ident = function(object,features.use=NULL, Count.mat = NULL, ident.use=NULL, thresh.use=NULL,do.plot=TRUE,do.transpose=FALSE,max.val.perc=NULL, max.val.exp=NULL,cols.use=gray.colors(10),max.size=10,min.perc=0,norm.exp=NULL,...) {
  
  features.use=set.ifnull(features.use, object@var.genes[1:20])
  features.use=features.use[features.use %in% rownames(object@data)]
  ident.use=set.ifnull(ident.use, levels(object@ident))
  thresh.use = set.ifnull(thresh.use, object@is.expr)
  
  #Matrix of percent expressing cells
  PercMat = matrix(0, nrow=length(features.use), ncol = 0)
  rownames(PercMat) = features.use; 
  
  #Matrix of average transcript levels
  ExpMat = PercMat;
  
  #Count mat
  Count.mat = set.ifnull(Count.mat, exp(object@data[features.use, colnames(object@data)]) - 1)
  
  
  for (i in ident.use){
    cells.in.cluster = names(object@ident)[which(object@ident== i)]
    vec.exp = apply(Count.mat[features.use, cells.in.cluster], 1, function(x) sum(x>thresh.use)/length(x)) 
    PercMat = cbind(PercMat,vec.exp)
    
    # Median in only expressing cells
    #vec.exp = apply(Count.mat[features.use, cells.in.cluster], 1, function(x) if (sum(x>0) > 1){ median(x[x>0]) } else {sum(x)})
    
    # Median in all cells
    vec.exp = apply(Count.mat[features.use, cells.in.cluster], 1, function(x) if (sum(x>0) > 1){ mean(x) } else {sum(x)})
    ExpMat = cbind(ExpMat, vec.exp)
  }
  colnames(ExpMat) = ident.use
  colnames(PercMat) = ident.use
  print(dim(ExpMat))
  
  if (!is.null(norm.exp)){
    if (norm.exp < 0 | norm.exp > 1){
      print("Warning : norm.exp should be a value betwen (0,1). Skipping normalization")
      next
    } else{
      quant.vals = apply(ExpMat, 1, function(x) quantile(x, norm.exp))
      ExpMat = t(scale(t(ExpMat), center=FALSE, scale=quant.vals))
      
    }
    
  }
  
  rows.use = rownames(PercMat)[apply(PercMat, 1, function(x) max(x) >= min.perc)];
  PercMat = PercMat[rows.use,]
  ExpMat = ExpMat[rows.use,]
  features.use = rows.use
  if (!is.null(max.val.perc)) PercMat[PercMat > max.val.perc] = max.val.perc
  if (!is.null(max.val.exp)) ExpMat[ExpMat > max.val.exp] = max.val.exp
  
  if (do.plot){
    
    ExpVal = melt(ExpMat)
    PercVal = melt(PercMat)
    colnames(ExpVal) = c("gene","cluster","nTrans")
    ExpVal$percExp = PercVal$value*100
    
    if (!do.transpose){
      ExpVal$gene = factor(ExpVal$gene, levels=rev(features.use))
      ExpVal$cluster = factor(ExpVal$cluster, levels= ident.use)
      p=ggplot(ExpVal, aes(y = factor(gene),  x = factor(cluster))) + geom_point(aes(colour = nTrans,  size =percExp)) + 
        scale_color_gradient(low ="blue",   high = "red", limits=c(min(ExpVal$nTrans), max(ExpVal$nTrans) ))+scale_size(range = c(1, max.size))+   theme_bw() +nogrid
      p = p + xlab("Cluster") + ylab("Gene") + theme(axis.text.x=element_text(size=12, face="italic", hjust=1)) + 
        theme(axis.text.y=element_text(size=12, face="italic"))  
      print(p)
    } else {
      ExpVal$gene = factor(ExpVal$gene, levels=features.use)
      ExpVal$cluster = factor(ExpVal$cluster, levels= rev(ident.use))
      p=ggplot(ExpVal, aes(y = factor(cluster),  x = factor(gene))) + geom_point(aes(colour = nTrans,  size =percExp)) + 
        scale_color_gradient(low ="blue",   high = "red", limits=c( min(ExpVal$nTrans), max(ExpVal$nTrans) ))+scale_size(range = c(1, max.size))+   theme_bw() +nogrid
      p = p + ylab("Cluster") + xlab("Gene") + theme(axis.text.x=element_text(size=12, face="italic", angle=45, hjust=1)) + 
        theme(axis.text.y=element_text(size=12, face="italic"))
      print(p)
      
    }
    
  }else {
    to.return=list()
    to.return$ExpMat = ExpMat;
    to.return$PercMat = PercMat;
    return(to.return)
  }
}