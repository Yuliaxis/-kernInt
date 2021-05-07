# R package 'kernInt'

-----------------------------

- Current version: 0.1.0
- Author: Elies Ramon
- e-mail: eramongurrea@gmail.com

## Purpose

**kernInt** uses the kernel framework to unify supervised and unsupervised microbiome analyses, while paying special attention to spatial and temporal integration.  If you find our package useful, please cite:

Ramon E, Belanche-Muñoz L, Molist F, Quintanilla R, Perez-Enciso M, and Ramayo-Caldas Y (2021) kernInt: A Kernel Framework for Integrating Supervised and Unsupervised Analyses in Spatio-Temporal Metagenomic Datasets. Front. Microbiol. 12:609048.doi: 10.3389/fmicb.2021.609048

## Installation and usage 

In R console:  

```														
if (!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("elies-ramon/kernInt")
```

If metagenomeSeq was not installed previously:

```
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install("metagenomeSeq")
```

Once the package is installed, it can be used anytime typing:

```
library(kernInt)
```

## Package Overview

### Main features

- Integration of supervised (classification, regression) and unsupervised (kernel PCA, hierarchical clustering, outlier detection) analyses.
- Microbial signatures of the classification and regression models.
- Automatic training/test splitting of the input data, k-Cross Validation and SVM classification and regression.
- Integration of data from different sources via Multiple Kernel Learning (MKL)
- Previously unpublished longitudinal pig gut microbiome dataset
- Implementation of kernels for compositional data (Aitchison-RBF kernel, compositional linear)
- Implementation of kernels suitable for functional data (functional RBF, functional linear)
- Kernels derived from classical ecology distances, as Jaccard and Jensen-Shannon, are also available.


### Example datasets

We offer three metagenomic datasets with the package: a single point soil dataset, a human health dataset with an spatial component, and a novel longitudinal dataset concerning pig production. Also, to better illustrate the longitudinal treatment of data, we include the classical Berkeley Growth Dataset.

[**Soil data**](https://qiita.ucsd.edu/study/description/103): Bacterial abundances in 88 soils from across North and South America. Metadata as soil pH, annual season precipitation and temperature, country, elevation, etc. is available.

[**Smokers**](https://qiita.ucsd.edu/study/description/524): Microorganism abundances of right and left oro- and nasopharynx in 29 smokers and 33 nonsmokers.

**Pig data**: Previously unpublished longitudinal gut microbiome dataset of 153 piglets during their first week of life.
 
[**Growth**](https://europepmc.org/article/med/13217130): Berkeley longitudinal height data of 54 girls and 39 boys (93 individuals in total) from ages 0 to 18.


### Vignette

An online in-depth vignette covering the kernel framework, with step-to-step usage and detailed examples can be found [here.](https://elies-ramon.github.io/)

The same vignette can be also accessed offline when the package is loaded, typing:

``` 
browseVignettes("kernInt")
```

If no vignette is found, try again after doing this:

``` 
devtools::install_github("elies-ramon/kernInt", build_vignettes = TRUE)
``` 

