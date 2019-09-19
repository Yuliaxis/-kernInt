# CLASSIFICATION PERFORMANCE MEASURES

#' Accuracy
#' @param ct Confusion Matrix
#' @return Accuracy
#' @examples
#' diag <- as.numeric(speMGX[,1])
#' diag[diag == 3] <- 1  # De 3 a 2 classes: No Malalt /  malalt
#' CM <- classify(data=speMGX[,7:ncol(speMGX)],y=diag,kernel="qJac")
#' Acc(CM)
#' @export

Acc <- function(ct) sum(diag(ct))/sum(ct)


#' F1
#' @param ct Confusion Matrix
#' @return F1
#' @examples
#' diag <- as.numeric(speMGX[,1])
#' diag[diag == 3] <- 1  # De 3 a 2 classes: No Malalt /  malalt
#' CM <- classify(data=speMGX[,7:ncol(speMGX)],y=diag,kernel="qJac")
#' F1(CM)
#' @export

F1 <-  function(ct) {
  REC <- Rec(ct)
  PREC <- Prec(ct)
  f1 <- (2*PREC*REC)/(PREC+REC)
  if(is.nan(f1)) f1 <- 0
  return(f1)
}


#' Precision
#' @param ct Confusion Matrix
#' @return Precision
#' @examples
#' diag <- as.numeric(speMGX[,1])
#' diag[diag == 3] <- 1  # De 3 a 2 classes: No Malalt /  malalt
#' CM <- classify(data=speMGX[,7:ncol(speMGX)],y=diag,kernel="qJac")
#' Prec(CM)
#' @export

Prec <- function(ct) {
  pr <- ct[2,2]/sum(ct[,2])
  if(is.nan(pr)) pr <- 0
  return(pr)
}

#' Recall
#' @param ct Confusion Matrix
#' @return Recall
#' @examples
#' diag <- as.numeric(speMGX[,1])
#' diag[diag == 3] <- 1  # De 3 a 2 classes: No Malalt /  malalt
#' CM <- classify(data=speMGX[,7:ncol(speMGX)],y=diag,kernel="qJac")
#' Rec(CM)
#' @export

Rec <-  function(ct) {
  rc <- ct[2,2]/sum(ct[2,])
  if(is.nan(rc)) rc <- 0
  return(rc)
}

