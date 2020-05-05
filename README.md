# R package 'kernInt'

-----------------------------

- Current version: 0.1.0
- Author: Elies Ramon
- e-mail: elies.ramon@cragenomica.es

## Purpose

**kernInt** uses the kernel framework to unify supervised and unsupervised microbiome analyses, while paying special attention to spatial and temporal integration. 

## Installation

In R console:  
															
`install.packages("devtools")`

`devtools::install_bitbucket("elies_ramon/kernint")`

## Package Overview

### Main features

- Implementation of compositional kernels like Aitchison-RBF kernel.
- Implementation of kernels suitable for functional data.
- Kernels derived from classical ecology distances, as Jaccard and Jensen-Shannon, are also available.
- Previously unpublished longitudinal pig gut microbiome dataset
- Automatic training/test splitting of the input data, k-Cross Validation and SVM classification and regression.
- Microbial signatures of the classification and regression models.
- Outliers / Novelty detection via SVMs.
- Integration of data from different sources via Multiple Kernel Learning (MKL)

### Example data

[**Soil data**](https://qiita.ucsd.edu/study/description/103): Bacterial abundances (raw counts) in 88 soils from across North and South America. Metadata as soil pH, annual season precipitation and temperature, country, elevation, etc. is available.

[**Smokers**](https://qiita.ucsd.edu/study/description/524): Microorganism abundances of oro and nasopharynx in 29 smokers and 33 nonsmokers.
 
[**Growth**](): Berkeley longitudinal height data of 54 girls and 39 boys (93 individuals in total) from ages 11 to 18.

**Pig data**:  Previously unpublished longitudinal gut microbiome dataset of 153 piglets during their first week of life.

### Vignette

An in-depth vignette covering kernel background, with step-to-step usage and detailed examples can be found *here*.
 
## Usage

### kernel PCA

As the standard PCA, kernel PCA can be used to summarize a dataset, to visualize a dataset or to extract features of a dataset. Data can be projected in a linear or nonlinear way, depending on the kernel used. When the kernel is the standard linear, kernel PCA is equivalent to standard PCA.

The `kernPCA()` function have two mandatory arguments: `data` and `kernel`:

`kernPCA(data=soil$abund,kernel="clin")`

The rest of arguments customize the plot. For example, the dots can be coloured according to a desired variable, which can be continuous or discrete. Here we show a kernel PCA of soil samples in which each sample is coloured accoirding to its pH: 

`kernPCA(data=soil$abund,kernel="clin",  y=soil$metadata$ph, col=c("aquamarine2","orchid3"),title = "Soil kernel PCA",legend = TRUE)`

The projected data can be retrieved setting `plot=FALSE`.

### Clustering

A dendogram plot presenting a hierarchical clustering can be obtained with:

`hklust(data=soil$abund,kernel="jac")`

Additional arguments allow changing the agglomeration method, draw the clusters or customizing the plot:

`hklust(data=soil$abund,kernel="jac", title = "Soil data cluster dendogram",cut=3,colors=2:4)`

The dendogram object can be retrieved setting `plot=FALSE`.


### SVM regression 

SVM regression is performed with the `regress()` function. `regress()` performs automatic training/test splitting of the input data, k-cross validation if requested, and regression with optimal hyperparameters. To perform MKL, go to *MKL section*.

The most basic call only needs three arguments: `data` (predictor variables; e.g. taxonomic abundances), `y` (target variable; e.g. a phenotype) and `kernel`.

For example, if we want to predict the pH of soil (`y`) from the abundances (`data`) using the compositional linear kernel:

`regress(data=soil$abund, kernel="clin", y=soil$metadata$ph)`

If the user has a pre-computed kernel matrix at hand, it can be passed as input to `data`. `kernel` should then be turned to `kernel="matrix"`.

The SVM hyperparameters Cost (`C`) and Epsilon (`E`) can be specified, thus tuning how the model will adjust to the data.

`regress(data=soil$abund, kernel="clin", y=soil$metadata$ph, C=5, E=0.001)`

In addition, a generic kernel hyperparameter (`H`) can be specified. For example, if the chosen kernel is RBF, `H` will be interpreted as *gamma*: (*RBF(x,y) = exp(-gamma * ||x-y||^2*)

`regress(data=soil$abund, y=soil$metadata$ph, kernel="crbf", C=5, H=0.1)`

k-Cross Validation can be performed to train the hyperparameters. This is done providing an argument to `k`:

`regress(data=soil$abund, y=soil$metadata$ph, kernel="clin", C=c(1,5,10), E=c(0.001,0.1), k=10)`

Training/test splitting is controlled with the `p` argument and the rownames of `data`. If given a numeric value between 0 and 1, `regress()` will consider it the proportion of data instances for the test set, and will do a random splitting. Default is `p=0.2`. If the input data has repeated rownames, `regress()` will consider that the rows that share id are repeated measures coming from the same individual. The function will ensure that all repeated measures are used either to trainor to test the model, but not for both, thus preserving the independence of the training and tets sets. However, users can enter the test partition, setting `p` to be a numeric (row indexes) or character (rownames) vector. The remainig data will be used as training.

####   Output 

A list containing:

- `$nmse`: Normalized mean squared error over the test data. This permits evaluating how good the model is at predicting.

- `$hyperparam`: Hyperparameter values used to build the model and cross-validation errors, if applicable.

- `$prediction`: Predicted and true values (test set). Rownames correspond to the indexes in the original data.

- `$var.imp`: The variable importance (e.g. microbial signature) if a linear or linear-like kernel is used.


### SVM classification

SVM classification is performed via the `classify()` function. Both binary or multiclass classification (one-vs-one) are supported. One-class classification is available in the `outliers()` function.

The usage of `classify()` is for the most part identical to that of `regress()`. For example, to predict if a certain soil came from forest, tropical, shrubland or grassland:

`classify(data=soil$abund ,y=soil$metadata[ ,"env_feature"],kernel="clin")`

Probabilistic classification is available setting `prob=TRUE`:

`classify(data=soil$abund ,y=soil$metadata[ ,"env_feature"],kernel="clin",prob=TRUE)`

Also, `classify()` supports several methods to deal with imbalanced data:

-Class weighting: 

`classify(data=soil$abund ,y=soil$metadata[ ,"env_feature"],kernel="clin", classimb="weights",C=c(0.001,0.01),k=10)`

-Undersampling:

`classify(data=soil$abund ,y=soil$metadata[ ,"env_feature"],kernel="clin", classimb="ubUnder",C=c(0.001,0.01),k=10)`

-Oversampling: 

`classify(data=soil$abund ,y=soil$metadata[ ,"env_feature"],kernel="clin", classimb="ubOver",C=c(0.001,0.01),k=10)`
 
####   Output 

A list containing:

- `$conf.matrix`: Confusion matrix (true versus predicted) for the test data. This permits evaluating how good the model is at predicting.

- `$hyperparam`: Hyperparameter values used to build the model and cross-validation errors, if applicable.

- `$prediction`: Predicted and true values (test set). Rownames correspond to the indexes in the original data. If `prob=TRUE`, the probability of each observation to belong to a given class.

- `$var.imp`: The variable importance (e.g. microbial signature) if a linear or linear-like kernel is used.

### Outlier detection

The `outliers()` function can be used either in a supervised or in an unsupervised way.

In the latter approach, the most basic call to this function needs two arguments: `data` (predictor variables) and `kernel` (the kernel function used). Then, the function will return the data outliers:

`outliers(data=soil$abund ,kernel="clin")`

The nu hyperparameter (`nu`) and a gamma hyperparamete `H` can be entered:

`outliers(data=soil$abund,kernel="crbf",nu=0.2,H=0.05)`

If an argument for `y` is provided, `outliers()` functions as an one-class SVM. In that case, cross-validation will be performed if `k` has an argument. Also, `p` stands for the proportion of data instances reserved for the training set 

`outliers(data=soil$abund,y=soil$metadata[ ,"env_feature"],kernel="clin",nu=c(0.1,0.2),k=5)`


## MKL

MKL (Multiple Kernel Learning) is available to both `classify()` and `regress()`. All features of these two functions are available when performing MKL. 

To do MKL, the `data` argument must be a list of length > 1. Each element of the list should be a data.frame or matrix, and rows should coincide. If `kernel="matrix"` data may be a tridimensional array. `kernel` argument may contain only one kernel name (thus implying that the kernel is the same for all datasets) or a vector of kernel names. That way a different kernel will be applied to each data type. For example, if we have a list with two pre-computed kernel matrices.

`regress(data=grKern[c(1,3)], y=growth[,2], kernel=c("matrix","matrix"))` 

The `coeff` argument is for the weight of each data type in the kernel combination. When absent, the mean across all kernel matrices is performed.

`regress(data=grKern[c(1,3)], y=growth[,2], coeff=c(0.6,0.4), kernel=c("matrix","matrix"),C=0.1)` 

The use of additional arguments as C, E, p... remains the same. Kernel(s)' generic hyperparameter In the MKL usage, `H` must be NULL or, either, a list so each element is the hyperparameter applied of each kernel function:

`classify(data=grKern[c(2,4)],y=growth[,1],kernel=c("linear","cov"), H=list(a=NULL,b=0.25))`

In the case of k-Cross-Validation:

`classify(data=grKern[c(2,4)],y=growth[,1],kernel=c("linear","cov"), H=list(dataset1=NULL,dataset2=c(0.25,0.5,0.75)),k=5)`


#### Side functions for data fusion

`KInt()` and `fuseData()` return a fused kernel matrix, but the use kernel matrices as input while in the latter a list with the different data sources is needed.

`d <- list()`

`d[[1]] <- matrix(abs(rnorm(20)),nrow=4,ncol=5)`

`d[[2]] <- matrix(abs(rnorm(20)),nrow=4,ncol=5)`

We can use different kernel functions for each type of data; for example, the Jaccard kernel for the first data type and the linear kernel for the second:

`fuseData(DATA=d,kernel=c("jac","lin"))`

The former command consider the two sources equally important. If not, we can state the weights:

`fuseData(DATA=d,kernel=c("jac","lin"),coeff=c(0.9,0.1))`

## Functional and Longitudinal data

Two kernels for longitudinal kernels are implemented: "flin" and "frbf". When used, `domain` argument has to be provided.
