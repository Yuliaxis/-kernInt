## DATA NORMALIZATION ("not the compositional way")

#' DATA NORMALIZATION
#' @param data Input data
#' @return Normalized data
#' @examples
#' soilData <- CSSnorm(data=t(soilDataRaw))
#' @import metagenomeSeq
#' @export


CSSnorm <- function(data) {
  data.metagenomeSeq <-  newMRexperiment(data, #phenoData<- phenotypeData,
                                       featureData=NULL, libSize=NULL, normFactors=NULL) #using filtered data
  p <-  cumNormStat(data.metagenomeSeq) #default is 0.5
  data.cumnorm <-  cumNorm(data.metagenomeSeq, p=p)
  #data.cumnorm
  data.CSS <-  t(MRcounts(data.cumnorm, norm=TRUE, log=TRUE))
  return(data.CSS)
}

