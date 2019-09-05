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

`cplot(cRBFmatrix,y=soilMetaData$ph[-89], col=c("aquamarine2","orchid3"),title = "Soil - AitchisonRBF PCA",legend = TRUE)`

A gradient from more acid (greenish) soils to more basic (violet) soils can be observed.

### Clustering

To perform a hierarchical clustering:

`rownames(cRBFmatrix) <- soilMetaData$ph[-89]`

`distances <- as.dist(1-cRBFmatrix)`

`clustering <- hclust(distances,method = "ward.D2")`

Plotting the dendogram: 

`plot(clustering)`

`clusters <- cutree(clustering, k = 3)`

Three different colors to represent acid pH, "medium" pH and basic pH: 

`rect.hclust(clustering , k = 3, border = c("aquamarine2", "azure3","orchid3"))`

### Outlier detection

`outliers(data=soilDataRaw,kernel="cRBF",nu=0.9)`

### SVM regression 

We want to predict the pH of soil (`y`) from the abundances (`data`). 
Being `cRBF` the Aitchison RBF kernel and `p` the proportion of data instances for the training set:

`regress(data=soilDataRaw[-89,], y=soilMetaData$ph[-89], kernel="cRBF",p=0.8)`

We can perform k-Cross Validation to train the hyperparameters (Cost):

`regress(data=soilDataRaw[-89,], y=soilMetaData$ph[-89], kernel="cRBF", C=c(0.1,1,10), k=10, p=0.8)`

