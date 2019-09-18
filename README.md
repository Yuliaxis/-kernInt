# R package 'kernInt'

-----------------------------

- Current version: 0.1.0
- Author: Elies Ramon
- e-mail: elies.ramon@cragenomica.es

## Purpose

**kernInt** aim is to integrate microbiome data of heterogeneous origins (metagenomics, metabolomics...) and, also, to use the kernel matrices to unify the diverse approaches (visualization, clustering, prediction...) used in the analysis of microbiome.


## Installation

In R console:  
															
`install.packages("devtools")`

`devtools::install_bitbucket("elies_ramon/kernint")`

## Package Overview

### Main features

- Implementation of quantitative Jaccard (Ruzicka similarity) with or without weights
- Implementation of Aitchison-RBF kernel
- Outliers / Novelty detection via SVMs.
- Integration of kernel matrices (the only approach available right now is doing the mean)
- Training/test splitting, k-Cross Validation and SVM classification and regression

### Example data

[**Soil data**](https://qiita.ucsd.edu/study/description/103)

- `soilDataRaw`: Bacterial abundances (raw counts) in 89 soils from across North and South America. 

- `soilMetaData`: Soil pH, annual season precipitation and temperature, country, elevation, etc.
 
[**IBDMDB** ](https://ibdmdb.org/tunnel/public/HMP2/WGS/1818/products)

- `speMGX`: Species abundances of intestinal microbiota of 84 individuals, 41 with Chron's Disease, 24 with Ulcerative Colitis, and 19 without inflamatory inflammatory bowel disease.

- `genMGX`: Genus abundances of intestinal microbiota of 84 individuals, 41 with Chron's Disease, 24 with Ulcerative Colitis, and 19 without inflamatory inflammatory bowel disease.
 
## Usage

### Visualization / Ordination

`library(catkern)`

`cRBFmatrix <- clrRBF(soilDataRaw[-89,])`

Plotting the kernel PCA: 

`cplot(cRBFmatrix,y=soilMetaData$ph[-89], col=c("aquamarine2","orchid3"),title = "Soil - Aitchison RBF kPCA",legend = TRUE)`

A gradient from more acid (greenish) soils to more basic (violet) soils can be observed.

We can also try the non-compositional kernels like quantitative Jaccard:

We first perform a Cumulative Sum Scaling normalisation (CSS):

`soilData <- CSSnorm(data=t(soilDataRaw[-89,]))`

`Jmatrix <- qJacc(soilData)`

`cplot(Jmatrix,y=soilMetaData$ph[-89], col=c("aquamarine2","orchid3"),title = "Soil - Ruzicka kPCA",legend = TRUE)`

or with weights:

`wJmatrix <- wqJacc(soilData,y=soilMetaData$ph[-89])`

`cplot(wJmatrix,y=soilMetaData$ph[-89], col=c("aquamarine2","orchid3"),title = "Soil - weighted Ruzicka kPCA",legend = TRUE)`


### Clustering

To perform a hierarchical clustering:

`rownames(cRBFmatrix) <- soilMetaData$ph[-89]`

`distances <- as.dist(1-cRBFmatrix)`

`clustering <- hclust(distances,method = "ward.D2")`

Plotting the dendogram: 

`plot(clustering)`

`clusters <- cutree(clustering, k = 3)`

Three different colors to represent acid pH, "intermediate" pH and basic pH: 

`rect.hclust(clustering , k = 3, border = c("aquamarine2", "orchid3","azure3"))`

### Outlier detection

`outliers(data=soilDataRaw,kernel="cRBF",nu=0.2)`

### SVM regression 

We want to predict the pH of soil (`y`) from the abundances (`data`). 
Being `cRBF` the Aitchison RBF kernel and `p` the proportion of data instances for the training set:

`regress(data=soilDataRaw[-89,], y=soilMetaData$ph[-89], kernel="cRBF",p=0.8)`

We can perform k-Cross Validation to train the hyperparameters (Cost):

`regress(data=soilDataRaw[-89,], y=soilMetaData$ph[-89], kernel="cRBF", C=c(0.1,1,10), k=10, p=0.8)`

### SVM classification

We want to predict if a certain individual has IBD or not.

`diag <- as.numeric(speMGX[,1])`
`diag[diag == 3] <- 1`
`classify(data=speMGX[,7:ncol(speMGX)],y=diag,kernel="qJac",C=c(0.1,1,10), k=10)`

`kernInt` supports several methods to deal with imbalanced data:

-Class weighting: 

`classify(data=speMGX[,7:ncol(speMGX)],diag,kernel="qJac",classimb="weights",C=c(0.001,0.01),k=10)`

-Undersampling:

`classify(data=speMGX[,7:ncol(speMGX)],diag,kernel="qJac",classimb="data",type="ubUnder",C=c(0.001,0.01),k=10`

-Oversampling: 

`classify(data=speMGX[,7:ncol(speMGX)],diag,kernel="qJac",classimb="data",type="ubOver",C=c(0.001,0.01),k=10`

-SMOTE: 

`classify(data=speMGX[,7:ncol(speMGX)],diag,kernel="qJac",classimb="data",type="ubSMOTE",C=c(0.001,0.01),k=10`

-One-class SVM:

`outliers(data=speMGX[,7:ncol(speMGX)],y=diag,kernel="qJac",p=0.8,k=10)`

-Probabilistic SVM:

`classify(data=speMGX[,7:ncol(speMGX)],diag,kernel="wqJac",C=c(0.1,1),CUT=c(0.3,0.4,0.5),k=10,prob = TRUE)`

-Unbalanced SVM:

(falta)

