library(Hmisc)
# A function to merge cells of two clusters into a single cluster
force.merge = function(object,merge.clust=NULL,new.name=NULL,char=FALSE,set.data.info=NULL,...){
            
            if (is.null(merge.clust)) return(object)
            if (is.null(new.name)) new.name = merge.clust[1]
            if (!char){
              if (is.null(set.data.info)){
                set.data.info="res.0.6"
                
              }
              
              ident.use= as.numeric(object@ident)
              names(ident.use) = names(object@ident)
              ident.use[ident.use %in% merge.clust] = new.name
              ident.use = factor(ident.use)
              if ((min(as.numeric(ident.use))) == 0){
                levels(ident.use) = c(0:(length(levels(ident.use))-1))
              } else {
                levels(ident.use) = c(1:length(levels(ident.use)))
              }
              
              object@ident = ident.use
              object@data.info[,set.data.info] = ident.use
            }
            
            if (char){
              if (is.null(set.data.info)){
                set.data.info="res.0.6"
                
              }
              ident.use= as.character(object@ident)
              names(ident.use) = names(object@ident)
              ident.use[ident.use %in% merge.clust] = new.name
              ident.use = factor(ident.use)
              object@data.info[,set.data.info] = ident.use
            }
            print("Returning Object")
            return(object)
          
}


cluster.dist.by.ident = function(object,ident1.use=NULL, ident2.use="orig.ident", ident2.counts=NULL, ident1.names=NULL,  ident2.names=NULL, cells.use=NULL, ylab.use="%", xlab.use="Sample", legend.title="cluster") 
{
            
            if (is.null(cells.use)){
             cells.use = colnames(object@data)
    
             }
            if (is.null(ident1.use)){
              ident1.use = "res.0.6"
              
            }
            ident1.use = set.ifnull(ident1.use,"res.0.6")
            if (!(ident1.use %in% colnames(object@data.info)) | !(ident2.use %in% colnames(object@data.info))){
              stop("One of ident1.use or ident2.use is invalid")
            }
            
            
            ident2.counts = set.ifnull(ident2.counts, table(object@data.info[cells.use,c(ident2.use)]))
            
            if (length(ident2.counts) != length(table(object@data.info[cells.use,ident2.use])) ){
              stop("The number of entries in ident2.counts must match the number of unique labels in ident2")
            }
            
            if (is.null(names(ident2.counts))){
              names(ident2.counts) = names(table(object@data.info[cells.use,ident2.use]))
            }
            
            
            print("Using the following base counts for ident.2")
            print(ident2.counts)
            
            # Counts
            object@data.info[,ident1.use] = as.factor(object@data.info[,ident1.use])
            object@data.info[,ident2.use] = as.factor(object@data.info[,ident2.use])
            df = table(object@data.info[cells.use,c(ident1.use, ident2.use)])
            df = scale(df, center=FALSE, scale=ident2.counts)
            df=df*100# 100 to convert to percentages
            df.melt = melt(df)
            df.melt[,ident1.use] = factor(df.melt[,ident1.use])
            df.melt[,ident2.use] = factor(df.melt[,ident2.use])
            if ((ylab.use) == "%") ylab.use = paste0("Percent", capitalize(ident2.use))
            ylab.use = gsub(" ","", ylab.use)
            colnames(df.melt)[3] = ylab.use
            require(RColorBrewer)
            nColors = length(unique(df.melt[,ident1.use]))
            if (nColors <=8 ){
              col.use= brewer.pal(nColors,"Set2")
            } else {
              require(randomcoloR)
              col.use = unname(distinctColorPalette(nColors))
            }
            p = ggplot() + geom_bar(aes_string(x=ident2.use,y=ylab.use,fill=ident1.use), data=df.melt, stat="identity") + xlab(xlab.use) +
              ylab(ylab.use) + scale_fill_manual(values = col.use) 
            #rainbow(length(unique(df.melt[,ident1.use])))
            #+  scale_fill_brewer(palette="Set2", guide=guide_legend(title = legend.title))  
            print(p)
}

set.ifnull = function(a,b){
  if (is.null(a)){
    a=b
  }
  return(a)
}



subset.data = function(object,cells.use=NULL,subset.name=NULL,ident.use=NULL,accept.low=-Inf, accept.high=Inf,do.center=F,do.scale=F,max.cells.per.ident=Inf, random.seed = 1,...) {
  data.use=NULL
  cells.use=set.ifnull(cells.use,object@cell.names)
  if (!is.null(ident.use)) {
    cells.use=WhichCells(object,ident.use)
  }
  if (!is.null(subset.name)) {
    data.use=FetchData(object,subset.name,...)
    if (length(data.use)==0) return(object)
    subset.data=data.use[,subset.name]
    pass.inds=which((subset.data>accept.low) & (subset.data<accept.high))
    cells.use=rownames(data.use)[pass.inds]
  }
  cells.use=ainb(cells.use,object@cell.names)
  cells.use = WhichCells(object, cells.use = cells.use, max.cells.per.ident = max.cells.per.ident, random.seed = random.seed)
  object@data=object@data[,cells.use]
  if(!(is.null(object@scale.data))) {
    if (length(colnames(object@scale.data)>0)) {
      object@scale.data[,cells.use]
      object@scale.data=object@scale.data[complete.cases(object@scale.data),cells.use]
    }
  }
  if (do.scale) {
    object=ScaleData(object,do.scale = do.scale,do.center = do.center)
    object@scale.data=object@scale.data[complete.cases(object@scale.data),cells.use]
  }
  object@ident=drop.levels(object@ident[cells.use])
  object@tsne.rot=object@tsne.rot[cells.use,]
  object@pca.rot=object@pca.rot[cells.use,]
  object@cell.names=cells.use
  
  object@gene.scores=data.frame(object@gene.scores[cells.use,]); colnames(object@gene.scores)[1]="nGene"; rownames(object@gene.scores)=colnames(object@data)
  object@data.info=data.frame(object@data.info[cells.use,])
  #object@mix.probs=data.frame(object@mix.probs[cells.use,]); colnames(object@mix.probs)[1]="nGene"; rownames(object@mix.probs)=colnames(object@data)
  
  return(object)
}