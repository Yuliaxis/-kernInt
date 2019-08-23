# CLASSIFIER

#' SVM classifier
#'
#'
#' @param data Input data
#' @param classes Number of classes (2 or 3)
#' @param kernel "qJac" for quantitative Jaccard and "wqJacc" for quant Jaccard with weights
#' @param p Proportion of total data instances in the training set
#' @param C A vector with the possible costs to evaluate via k-Cross-Val. If no argument is provided cross-validation is not performed.
#' @param k The k for the k-Cross Validation. Minimum k = 2.
#' @param classimb "weights" to introduce class weights in the SVM algorithm and "data" to oversampling. If no argument provided, nothing is done.
#' @return The test error (accuracy)
#' @examples
#' classify(data=MGXdata,kernel="qJac",C=c(0.0001,0.001),k=10)
#' @importFrom kernlab SVindex as.kernelMatrix predict
#' @export




classify <- function(data,classes=2, kernel, p=0.8, C, k=10,  classimb="no") {

  # 1. Classes
  diagn <- as.numeric(data[,1])
  if(classes == 2)   diagn[diagn == 3] <- 1 # De 3 a 2 classes: No Malalt /  malalt
  diagn[diagn == 2] <- 0
  diagn <- as.factor(diagn)

  # 1. TR/TE
  N <- nrow(data)
  all.indexes <- 1:N

  learn.indexes <- trainIndx(n=N,ptrain=p)
  test.indexes <- all.indexes[-learn.indexes]

  nlearn <- length(learn.indexes)
  ntest <- N - nlearn


  # 2. Compute kernel matrix
  if(kernel == "qJac") {
    Jmatrix <- qJacc(data[,7:ncol(data)]) ##ATENCIÓ CANVIAR!!!!
  } else if(kernel == "wqJac") {
    Jmatrix <- wqJacc(data[,7:ncol(data)],y=data[,1])
  } else {

  }
  trMatrix <- Jmatrix[learn.indexes,learn.indexes]

  #3. Class imbalance


  # 4. Do R x k-Cross Validation
  if(hasArg(C)) {
    if(k<2) stop("k must be equal to or higher than 2")
    bh <- kCV(COST = C, K=trMatrix, Yresp=diagn[learn.indexes], k=k, R=k)
    cost <- bh$cost
  } else {
    cost <- 1
  }

  if(classimb=="data")  {

  } else if(classimb == "weights") {
    model <- ksvm(trMatrix, diagn[learn.indexes], kernel="matrix", type="C-svc",
                  C=cost,class.weights=c("1"=summary(diagn[learn.indexes])[1],"0"=summary(diagn[learn.indexes])[2]))
  } else {
    model <- ksvm(trMatrix, diagn[learn.indexes], kernel="matrix", type="C-svc",C=cost )
  }

  # 5. Prediction
  teMatrix <- Jmatrix[test.indexes,learn.indexes]

  teMatrix <- teMatrix[,SVindex(model),drop=FALSE]
  teMatrix <- as.kernelMatrix(teMatrix)
  pred <- kernlab::predict(model,teMatrix)

  ### Confusion matrix
  print(ct <- table(Truth=diagn[test.indexes], Pred=pred))

  # Test error
  te.error <- round(1-sum(diag(ct))/sum(ct),4)


  return(te.error)
}
