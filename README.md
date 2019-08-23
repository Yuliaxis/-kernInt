# R package 'kernInt'

-----------------------------

- Current version: 0.1.0
- Author: Elies Ramon
- e-mail: elies.ramon@cragenomica.es

## Purpose

*kernInt* aim is to integrate microbiome data of heterogeneous origins (metagenomics, metabolomics...) and, also, to use the kernel matrices to unify the diverse approaches (visualization, clustering, prediction...) used in the analysis of microbiome.


## Installation

In R console:  
															
`install.packages("devtools")`  
`devtools::install_bitbucket("elies_ramon/kernint")`

## Package Overview

### Main features

- Implementation of quantitative Jaccard (Ruzicka similarity) with or without weights
- Outliers / Novelty detection via SVMs.
- Integration of kernel matrices (the only approach available right now is doing the mean)
- Training/test splitting, k-Cross Validation and SVM classification
