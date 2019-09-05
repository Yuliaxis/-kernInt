# REGRESSION

#' SVM regression
#'
#' @param data Input data
#' @param y Reponse variable (continuous)
#' @param kernel "cRBF" for clrRBF, "qJac" for quantitative Jaccard and  "wqJacc" for quantitative Jaccard with weights
#' @param g Gamma hyperparameter
#' @param p Proportion of total data instances in the training set
#' @param C A vector with the possible costs to evaluate via k-Cross-Val. If no argument is provided cross-validation is not performed.
#' @param k The k for the k-Cross Validation. Minimum k = 2.
#' @return NMSE (normalized mean squared error)
#' @examples
#' regress(data=soilDataRaw[-89,],soilMetaData$ph[-89],kernel="cRBF",C=c(1,10,50),k=10)
#' @importFrom kernlab as.kernelMatrix kernelMatrix predict rbfdot SVindex
#' @export

regress <- function(data, y, kernel, g=1, p=0.8, C=1, k) {

  # 1. TR/TE
  N <- nrow(data)
  all.indexes <- 1:N
  learn.indexes <- trainIndx(n=N,ptrain=p)
  test.indexes <- all.indexes[-learn.indexes]
  nlearn <- length(learn.indexes)
  ntest <- N - nlearn

  # 2. Compute kernel matrix
  if(kernel == "qJac") {
    cat("quantJaccard kernel \n")
    Jmatrix <- qJacc(data)
  } else if(kernel == "wqJac") {
    Jmatrix <- wqJacc(data,y=y)
    cat("quantJaccard kernel + weights \n")
  }  else if(kernel == "cRBF") {
    Jmatrix <- clrRBF(data)
    cat("clr + RBF \n")
  }   else {
    cat("standard RBF \n")

    # Jmatrix <-  kernelMatrix(rbfdot(sigma = g),data)

  }

  trMatrix <- Jmatrix[learn.indexes,learn.indexes]

  # 4. Do R x k-Cross Validation
  if(hasArg(k)) {
    if(k<2) stop("k must be equal to or higher than 2")
    bh <- kCV.reg(COST = C, K=trMatrix, Yresp=y[learn.indexes], k=k, R=k)
    cost <- bh$cost
  } else {
    if(length(C)>1) paste("C > 1 - Only the first element will be used")
    cost <- C[1]
  }

  model <- ksvm(trMatrix, y[learn.indexes],type="eps-svr", kernel="matrix", C=cost )

  # 5. Prediction
  teMatrix <- Jmatrix[test.indexes,learn.indexes]
  teMatrix <- teMatrix[,SVindex(model),drop=FALSE]
  teMatrix <- as.kernelMatrix(teMatrix)
  pred <- kernlab::predict(model,teMatrix)
  print(pred)

  return( error.norm(y[test.indexes],pred))
}

