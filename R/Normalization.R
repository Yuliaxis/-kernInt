## DATA NORMALIZATION


#' Cumulative Sum Scaling Normalization
#' @param data Input data
#' @return Normalized data
#' @examples
#' soilData <- CSSnorm(data=soil$abund)
#' @import metagenomeSeq
#' @export


CSSnorm <- function(data) {
  data <- t(data)
  data.metagenomeSeq <-  newMRexperiment(data, #phenoData<- phenotypeData,
                                       featureData=NULL, libSize=NULL, normFactors=NULL) #using filtered data
  p <-  cumNormStat(data.metagenomeSeq) #default is 0.5
  data.cumnorm <-  cumNorm(data.metagenomeSeq, p=p)
  #data.cumnorm
  data.CSS <-  t(MRcounts(data.cumnorm, norm=TRUE, log=TRUE))
  return(data.CSS)
}


#' clr transform
#'
#' @param data A matrix or data.frame with compositional data (raw counts)
#' @return clr transformation of data
#' @examples
#' clr(data=soil$abund)
#' @importFrom robCompositions cenLR
#' @export
#'

clr <- function(data) {
  minv <- min(  data[data!=0] )
  minv <- minv/10
  return(cenLR(data+minv)$x.clr)
}

