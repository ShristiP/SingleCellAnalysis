# Author: Shristi Pandey
#============================================================================================
# This script is used to train a random forest model that is used to map clusters across the two developmental time points. 


# load the relevant seurat Objects to compare the gene expression patterns in cell types across larval and adult forebrains
load("/Users/ShristiPandey/Dropbox/SingleCellAnalysis/Forebrain/Objects/AdultForebrainWTAll_08212018.RObj")
load("/Users/ShristiPandey/Dropbox/SingleCellAnalysis/Forebrain/Objects/Znf_wt_forebrainonlyFinal.RObj")

##: ++++++Function to plot the confusion matrix++++++++++++++++++++++++++++++++++++++++++++++++
library(reshape)
plotConfusionMatrix = function(X,row.scale=TRUE, col.scale=FALSE, cols.use=gray.colors(10), max.size=5, ylab.use="Known", xlab.use="Predicted"){
  
  if (!col.scale & row.scale){ X = t(scale(t(X), center=FALSE, scale=rowSums(X)));  X=X*100 }
  if (col.scale & !row.scale){ X = scale(X, center=FALSE, scale=colSums(X)); X = X*100 }
  if(col.scale & row.scale){
    print("Only one of row.scale or col.scale should be true. performing row scaling by default")
    X = t(scale(t(X), center=FALSE, scale=rowSums(X)))
    X=X*100
  }
  X = melt(X)
  colnames(X) = c("Known", "Predicted", "Percentage")
  X$Known = factor(X$Known, levels=rev(unique(X$Known)));
  X$Predicted = as.factor(X$Predicted)
  p = ggplot(X, aes(y = Known,  x = Predicted)) + geom_point(aes(colour = Percentage,  size =Percentage)) + 
    scale_color_gradient(low ="blue",   high = "red", limits=c(0, 100 ))+scale_size(range = c(1, max.size))+theme_bw() #+nogrid
  p = p + xlab(xlab.use) + ylab(ylab.use) + theme(axis.text.x=element_text(size=12, face="italic", hjust=1)) + 
    theme(axis.text.y=element_text(size=12, face="italic"))  
  print(p)
}

#predicting Fish from Mouse

# Random Forest
library(randomForest)
AdultForebrain = AdultFB
LarvalBrain = Znf_wt_forebrain

# Random forest will be performed using the union of the genes that are variable aross both the larval and adult forebrains
genes.use = union(AdultForebrain@var.genes, LarvalBrain@var.genes)
genes.use = genes.use[genes.use %in% rownames(LarvalBrain@data)]
genes.use = genes.use[genes.use %in% rownames(AdultForebrain@data)]
#genes.use  = intersect(AdultForebrain@var.genes, LarvalBrain@var.genes)

#++Train on adult and predict larva++++++++++++++++++++++++++++++++++++++++++ 
training.set = c(); test.set=c() 
training.label = c(); test.label=c();
#Larva10x = set.all.ident(Larva10x, id = 'm')

for (i in as.numeric(levels(AdultForebrain@ident))){
  #print(i)
  cells.in.clust = WhichCells(AdultForebrain,i);
  n = min(500, round(length(cells.in.clust)*0.7))
  train.temp = cells.in.clust[sample(length(cells.in.clust))][1:n]
  test.temp = setdiff(cells.in.clust, train.temp)
  training.set = c(training.set,train.temp); test.set=c(test.set,test.temp)
  training.label = c(training.label, rep(i,length(train.temp))); test.label = c(test.label, rep(i, length(test.temp)));
  #print(test.label)
}

predictor_Data = t(scale(t(as.matrix(AdultForebrain@data[genes.use,])),center=TRUE, scale=TRUE))
# Option 2 (No z-scoring):
# predictor_Data = as.matrix(AdultForebrain@data[genes.use,])

tmp = as.vector(table(training.label))
sampsizes = rep(min(tmp),length(tmp))
 
rf_output=randomForest(x=t(predictor_Data[,training.set]), y=factor(training.label), importance = TRUE, ntree = 501, proximity=TRUE, sampsize=sampsizes, keep.inbag=TRUE, replace=FALSE) 
Conf_OOB0 = rf_output$confusion

test.predict.prob = predict(rf_output,t(predictor_Data[,test.set]), type="prob")
thresh = 0.08 # The class with the maximum probability needs to have at least this margin
test.predict = apply(test.predict.prob,1,function(x) if (max(x) > thresh){ which.max(x)-1 } else {40}) # 1: indicates that the functiona applies over rows: which.max(x) returns index but have to subtract 1 because indexing starts from 1 in R 
Conf_test = table(test.label,test.predict)

postscript("ConfusionMatrixMouseToMouse.eps", width = 7, height = 6)
plotConfusionMatrix(Conf_test,row.scale=TRUE, max.size = 12, xlab.use="Predicted RF Test", ylab.use="Mouse Habenula")
dev.off()

# Now use the trained RF model to predict adult labels for larval cells. 
Larva.rf = as.matrix(Znf_wt_forebrain@data[genes.use,])

# Scaling the data together with Adult
Larva.rf = t(scale(t(Larva.rf), center=rowMeans(as.matrix(LarvalBrain@data[genes.use,])), scale=TRUE))

# No scaling
Larva.rf[is.na(Larva.rf)] = 0
Larva.ident = factor(LarvalBrain@ident)
Larva.predict.prob = predict(rf_output,t(Larva.rf), type="prob")
thresh = 0.08 # The class with the maximum probability needs to have at least this margin
Larva.predict = apply(Larva.predict.prob,1,function(x) if (max(x) > thresh){ which.max(x) -1 } else {40})
Conf_test = table(Larva.ident, Larva.predict) 

postscript("ConfusionMatrixMouseToFIsh.eps", width = 9, height = 12)
plotConfusionMatrix(Conf_test,row.scale=TRUE, max.size = 12, xlab.use="predicted id", ylab.use="actual larval id")
dev.off()