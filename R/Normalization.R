## DATA NORMALIZATION

#' Cumulative Sum Scaling Normalization
#' @param data Input data
#' @param pcount Pseudocount value (default value 0)
#' @return CSS normalized data
#' @examples
#' soilData <- CSSnorm(data=soil$abund)
#' @importFrom metagenomeSeq newMRexperiment cumNormStat cumNorm MRcounts
#' @export


CSSnorm <- function(data,pcount=0) {
  data <- t(data+pcount)
  metSeqData <-  newMRexperiment(data, #phenoData<- phenotypeData,
                                       featureData=NULL, libSize=NULL, normFactors=NULL) #using filtered data
  p <-  cumNormStat(metSeqData) #default is 0.5
  cumNormData <-  cumNorm(metSeqData, p=p)
  #cumNormData
  cssData <-  t(MRcounts(cumNormData, norm=TRUE, log=TRUE))
  return(cssData)
}


#' clr transform
#'
#' @param data A matrix or data.frame with compositional data (raw counts)
#' @return clr transformation of data
#' @param pcount Pseudocount value (default value 0.1)
#' @examples
#' clr(data=soil$abund)
#' @importFrom robCompositions cenLR
#' @export
#'

clr <- function(data,pcount=0.1) {
  return(cenLR(data+pcount)$x.clr)
}

