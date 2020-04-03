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
#' model <- classify(data=soil$abund,y=soil$metadata[,"env_feature"],kernel="clin")
#' Acc(model$"conf.matrix")
#' @export

Acc <- function(ct) sum(diag(ct))/sum(ct)


#' F1
#' @param ct Confusion Matrix
#' @param min.class Minority class
#' @return F1
#' @examples
#' model <- classify(data=soil$abund,y=soil$metadata[,"env_feature"],kernel="clin")
#' F1(model$"conf.matrix")
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
#' model <- classify(data=soil$abund,y=soil$metadata[,"env_feature"],kernel="clin")
#' Prec(model$"conf.matrix")
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
#' model <- classify(data=soil$abund,y=soil$metadata[,"env_feature"],kernel="clin")
#' Rec(model$"conf.matrix")
#' @export

Rec <-  function(ct,min.class=2) {
  rc <- ct[min.class,min.class]/sum(ct[min.class,])
  if(is.nan(rc)) rc <- 0
  return(rc)
}

