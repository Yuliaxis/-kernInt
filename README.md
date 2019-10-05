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

- Implementation of quantitative Jaccard (Ruzicka similarity) with or without weights.
- Implementation of Aitchison-RBF kernel.
- Outliers / Novelty detection via SVMs.
- Integration of kernel matrices via Multiple Kernel Learning.
- Training/test splitting, k-Cross Validation and SVM classification and regression.

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


### SVM regression 

SVM regression is performed via the `regress()` function.

The most basic call to this function needs only three arguments: `data` (predictor variables), `y` (target variable) and `kernel` (the kernel function used).

For example, if we want to predict the pH of soil (`y`) from the abundances (`data`): 

`regress(data=soilDataRaw[-89,], y=soilMetaData$ph[-89], kernel="cRBF")`

If we have a pre-computed kernel matrix at hand, it can be passed as input to `data`. `kernel` should then be turned to `kernel="matrix"`.

The SVM hyperparameters Cost (`C`) and Epsilon (`E`) can be specified, and also the proportion of data instances for the training set (`p`).

`regress(data=soilDataRaw[-89,], y=soilMetaData$ph[-89], kernel="cRBF", p=0.6, C=5, E=0.001)`

In addition, a generic kernel hyperparameter (`H`) can be specified. For example, if the chosen kernel is RBF, `H` will be interpreted as *gamma*: (*RBF(x,y) = exp(-gamma * ||x-y||^2*)

`regress(data=soilDataRaw[-89,], y=soilMetaData$ph[-89], kernel="cRBF", C=5, H=0.1)`

We can perform k-Cross Validation to train the hyperparameters. This is done providing an argument to `k`:

`regress(data=soilDataRaw[-89,], y=soilMetaData$ph[-89], kernel="cRBF", C=c(1,5,10), E=c(0.001,0.1), k=10)`

If the input data has repeated rownames, `regress()` will consider that the row names that share id are repeated measures 
coming from the same individual. The function will ensure that all repeated measures are used either to train
or to test the model, but not for both, thus preserving the independence between the training and tets sets.

### SVM classification

The classical SVM classification is performed via the `classify()` function. One-class classification is available in the `outliers()` function.

The usage of `classify()` is for the most part similar to that of `regress()`. For example, if we want to predict if a certain individual has IBD or not:

`diag <- as.numeric(speMGX[,1])`

`diag[diag == 3] <- 1`

`classify(data=speMGX[,7:ncol(speMGX)],y=diag,kernel="qJac",C=c(0.1,1,10), k=10)`

Probabilistic classification is available setting `prob=TRUE`:

`classify(data=speMGX[,7:ncol(speMGX)],y=diag,kernel="qJac",prob=TRUE,C=c(0.1,1,10), k=10)`

Both classify() and outliers() have the same treatment regarding repeated row names in the input data than `regress()`. 
Also, `classify()` supports several methods to deal with imbalanced data:

-Class weighting: 

`classify(data=speMGX[,7:ncol(speMGX)],diag,kernel="qJac",classimb="weights",C=c(0.001,0.01),k=10)`

-Undersampling:

`classify(data=speMGX[,7:ncol(speMGX)],diag,kernel="qJac",classimb="data",type="ubUnder",C=c(0.001,0.01),k=10)`

-Oversampling: 

`classify(data=speMGX[,7:ncol(speMGX)],diag,kernel="qJac",classimb="data",type="ubOver",C=c(0.001,0.01),k=10)`

-SMOTE: 

`classify(data=speMGX[,7:ncol(speMGX)],diag,kernel="qJac",classimb="data",type="ubSMOTE",C=c(0.001,0.01),k=10)`

-One-class SVM:

`outliers(data=speMGX[,7:ncol(speMGX)],y=diag,kernel="qJac",p=0.8,k=10)`

-Probabilistic SVM with a cutoff different to 0.5:

`classify(data=speMGX[,7:ncol(speMGX)],diag,kernel="wqJac",C=c(0.1,1),CUT=c(0.3,0.4,0.5),k=10,prob = TRUE)`

-Unbalanced SVM:

(falta)

### Outlier detection

The `outliers()` function can be used either in a supervised or in an unsupervised way.

In the latter approach, the most basic call to this function needs only two arguments: `data` (predictor variables) and `kernel` (the kernel function used). 
Then, the function will return the data outliers:

`outliers(data=soilDataRaw,kernel="cRBF")`

The nu hyperparameter (`nu`) and a generic kernel hyperparamete `H` can be entered:

`outliers(data=soilDataRaw,kernel="cRBF",nu=0.2,H=05)`

If an argument for `y` is provided, `outliers()` functions as an one-class SVM. In that case, cross-validation will be performed if `k` has an argument. Also,
`p` stands for the proportion of data instances reserved for the training set 

`outliers(data=soilDataRaw,kernel="cRBF",y=soilMetaData$ph,nu=0.2)`


## MKL

MKL (Multiple Kernel Learning) is available to both `classify()` and `regress()`. All features of these two functions ara available when performing MKL. 

To do MKL, the `data` argument must be a list of length > 1 or a tridimensional array. Each element of the list should be a data.frame or matrix:

`falta`

A vector of kernel names can be passed to the `kernel` argument. That way a different kernel will be applied to each data type: 

`falta`

The `coeff` argument is for the weight of each data type in the kernel combination. When absent, the mean across all kernel matrices is performed.

`falta`

Kernel(s)' generic hyperparameter `H` can be a number (same hyperparameter for all data types), a vector (different hyperparameter for each data type 
and not cross-validation) or a matrix (different set of hyperparametes for each kernel to chose through cross-validation).



## Data fusion

MKL (Multiple Kernel Learning) is available through `KInt()` and `fuseData()`. The two return a fused kernel matrix, but the use kernel matrices as input while in the latter a list with the different data sources is needed.

`d <- list()`

`d[[1]] <- matrix(abs(rnorm(20)),nrow=4,ncol=5)`

`d[[2]] <- matrix(abs(rnorm(20)),nrow=4,ncol=5)`

We can use different kernel functions for each type of data; for example, Ruzicka for the first data type and Aitchison-RBF for the second:

`fuseData(DATA=d,kernel=c("qJac","cRBF"))`

The former command consider the two sources equally important. If not, we can state the weights:

`fuseData(DATA=d,kernel=c("qJac","cRBF"),coeff=c(0.9,0.1))`
