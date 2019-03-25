# Single cell RNA-seq Analysis 
This repository contains scripts and functions used to perform single cell analysis in:

Pandey S., et al. Comprehensive Identification and Spatial Mapping of Habenular Neuronal Types Using Single-Cell RNA-Seq. 2018. Current Biology 28(7):1052-1065. DOI: https://doi.org/10.1016/j.cub.2018.02.040



## Selection of variable Genes 
Single cell dataset is noisy and sparse. Therefore, downstream PCA and clustering analysis is only done on those genes that exhibit the greater variability in expression patterns among the singe cells. These genes have variability in expression caused by both biological and technical reasons. To detect only biologically meaningful variation in expression, a null mathematical model is built to model the relationship between average UMI counts and  coefficient of variation (CV) across all genes, based on a negative binomial distribution that incorporates sampling noise and relative library size. Those genes that have a CV greater than the null model (threshold determined by diff.cutoff) are chosen as variable. Those represent the genes with the greatest biological variabilities. 


## Marker Identification using AUCPR

Markers for single cell clusters are nominated by performing a differential expression analysis between the cells in the cluster of interest and the rest of the cells in the dataset. Markers’ specificity and precision can be quantified using a statistical test based on the area under theprecision-recall curve (AUCPR). AUCPR is a quantitative measure of the balance between recall (the sensitivity of marker gene detection within the cluster of interest) and precision (accuracy of the quantitative levels of gene as a predictor of the correct cell type).Markers can be ‘digital’ (expressed only in the marked cluster) typically with AUCPR values > 0.8, or analog(expressed at a higher level in the marked cluster, but also detectable in other clusters) with AUCPR values between 0.6-0.8. 

## Random Forest for mapping clusters from larva to adult 

To evaluate the correspondence between the clusters in larva and adult, a multi class random forest classifier is trained on 70% of the larval data. This trained classifier is then used to assign a cluster label for the remaining 30% of the data. We assigned a class label to each cell, but only if a minimum of 15% of trees in the forest converged onto a decision(given that there are 16 classes, 6.25% vote would constitute a majority). Otherwise, the cells were labeled unassigned. 

This classifier was then used to predict the larval labels for adult cells.




