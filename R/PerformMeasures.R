# PERFORMANCE MEASURES

##  NMSE (regression)
#' @keywords internal
#' @importFrom stats var

error.norm <- function(target,prediction) {
  N <- length(target)
  error <- sum((target-prediction)^2)/((N-1)*var(target))
  return(error)
}

#' Accuracy
#' @param ct Confusion Matrix
#' @return Accuracy
#' @examples
#' diag <- as.numeric(speMGX[,1])
#' diag[diag == 3] <- 1  # De 3 a 2 classes: No Malalt /  malalt
#' CM <- classify(data=speMGX[,7:ncol(speMGX)],y=diag,kernel="qJac")
#' Acc(CM$"conf.matrix")
#' @export

Acc <- function(ct) sum(diag(ct))/sum(ct)


#' F1
#' @param ct Confusion Matrix
#' @param min.class Minority class
#' @return F1
#' @examples
#' diag <- as.numeric(speMGX[,1])
#' diag[diag == 3] <- 1  # De 3 a 2 classes: No Malalt /  malalt
#' CM <- classify(data=speMGX[,7:ncol(speMGX)],y=diag,kernel="qJac")
#' F1(CM$"conf.matrix")
#' @export

F1 <-  function(ct,min.class=2) {
  REC <- Rec(ct,min.class)
  PREC <- Prec(ct,min.class)
  f1 <- (2*PREC*REC)/(PREC+REC)
  if(is.nan(f1)) f1 <- 0
  return(f1)
}


#' Precision
#' @param ct Confusion Matrix
#' @param min.class Minority class
#' @return Precision
#' @examples
#' diag <- as.numeric(speMGX[,1])
#' diag[diag == 3] <- 1  # De 3 a 2 classes: No Malalt /  malalt
#' CM <- classify(data=speMGX[,7:ncol(speMGX)],y=diag,kernel="qJac")
#' Prec(CM$"conf.matrix")
#' @export

Prec <- function(ct,min.class=2) {
  pr <- ct[min.class,min.class]/sum(ct[,min.class])
  if(is.nan(pr)) pr <- 0
  return(pr)
}

#' Recall
#' @param ct Confusion Matrix
#' @param min.class Minority class
#' @return Recall
#' @examples
#' diag <- as.numeric(speMGX[,1])
#' diag[diag == 3] <- 1  # De 3 a 2 classes: No Malalt /  malalt
#' CM <- classify(data=speMGX[,7:ncol(speMGX)],y=diag,kernel="qJac")
#' Rec(CM$"conf.matrix")
#' @export

Rec <-  function(ct,min.class=2) {
  rc <- ct[min.class,min.class]/sum(ct[min.class,])
  if(is.nan(rc)) rc <- 0
  return(rc)
}

