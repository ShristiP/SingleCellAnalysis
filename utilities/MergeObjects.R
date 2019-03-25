# Merges to Seurat objects and returns a new object with all the cells both objects. 
MergeObjects = function(object1, object2, id1_use = "res.0.5", 
                        id2_use="res.0.5", new_id_name = "res", 
                        object1_name = NULL, object2_name = NULL, norm.val=1e4, use.regressed=FALSE){
  
    # Assumes that object1 and object2 have been aligned to the same reference
    newobject = object1
    if (is.null(object1_name)){
      object1_name = "S1"
    }
    
    if (is.null(object2_name)){
      object2_name = "S2"
    }
    
    
    col_names_raw = c(paste0(object1_name, "_", colnames(object1@raw.data)), paste0(object2_name,"_", colnames(object2@raw.data))) 
    col_names = c(paste0(object1_name,"_", colnames(object1@data)), paste0(object2_name,"_", colnames(object2@data))) 
    newobject@raw.data = cbind(object1@raw.data, object2@raw.data)
    colnames(newobject@raw.data) = col_names_raw
    
    genes.use = union(rownames(object1@data), rownames(object2@data))
    
    #Normalize data
    col_sums = Matrix::colSums(newobject@raw.data[,col_names])
    newobject@data = norm.val * scale(newobject@raw.data[,col_names],center=FALSE, scale=col_sums)
    colnames(newobject@data) = col_names
    newobject@data = newobject@data[genes.use,]
    
    if (use.regressed){
      object1_data = object1@data
      colnames(object1_data) = paste0(object1_name,"_", colnames(object1@data))
      newobject@data[rownames(object1@data), colnames(object1_data)] = as.matrix(object1_data);
      
      object2_data = object2@data
      colnames(object2_data) = paste0(object2_name,"_", colnames(object2@data))
      newobject@data[rownames(object2@data), colnames(object2_data)] = as.matrix(object2_data);
    }
    
    
    # scale data
    newobject@scale.data = t(scale(t(newobject@data), center=TRUE, scale=TRUE))
    newobject@scale.data[is.na(newobject@scale.data)] = 0
    
    # var genes
    newobject@var.genes = union(object1@var.genes, object2@var.genes)
    
    # metaData
    metaData1 = object1@meta.data[,c(1:4)]
    metaData2 = object2@meta.data[,c(1:4)]
    metaData1$orig.ident = paste0(object1_name,"_",as.character(metaData1$orig.ident))
    metaData2$orig.ident = paste0(object2_name,"_", as.character(metaData2$orig.ident))
    newobject@meta.data = rbind(metaData1, metaData2)
    rownames(newobject@meta.data) = col_names
    newobject@meta.data$orig.ident = factor(newobject@meta.data$orig.ident)
    
    res1 = paste0(object1_name,"_",as.numeric(object1@meta.data[,id1_use]))
    res2 = paste0(object2_name,"_",as.numeric(object2@meta.data[,id2_use]))
    
    newobject@meta.data[,new_id_name] = factor(c(res1, res2))
    newobject@ident = newobject@meta.data[,new_id_name]
    return(newobject)

}

