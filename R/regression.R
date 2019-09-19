# REGRESSION

#' SVM regression
#'
#'regress() automatically trains a Support Vector Regrssion model, tests it and returns the Normalized Mean Squared Error.
#'
# Cross-validation is available to choose the best hyperparameters (e.g. Cost, Epsilon) during the training step.
#'
#' @param data Input data
#' @param y Reponse variable (continuous)
#' @param kernel "cRBF" for clrRBF, "qJac" for quantitative Jaccard and  "wqJacc" for quantitative Jaccard with weights.
#' "matrix" if a pre-calculated kernel matrix is given as input.
#' @param p Proportion of total data instances in the training set
#' @param C A vector with the possible costs to evaluate via k-Cross-Val. If no argument is provided cross-validation is not performed.
#' @param G Gamma hyperparameter
#' @param E Epsilon hyperparameter
#' @param k The k for the k-Cross Validation. Minimum k = 2.
#' @return NMSE (normalized mean squared error)
#' @examples
#' # Data normalization
#' soilData <- CSSnorm(data=t(soilDataRaw[-89,]))
#' # Simple regression without tuning the hyperparameters
#' regress(data=soilData,soilMetaData$ph[-89],kernel="qJac")
#' # The percentage of data for training can be changed (default: 0.8 Training / 0.2 Test):
#' regress(data=soilData,soilMetaData$ph[-89],kernel="qJac",p=0.6)
#' # Regression with 10-Cross-Validation to choose the best Cost and Epsilon:
#' regress(data=soilData, y=soilMetaData$ph[-89], kernel="qJac", C=c(0.1,1,10), E = c(0.01,0.1), k=10)
#' @importFrom kernlab as.kernelMatrix kernelMatrix predict rbfdot SVindex
#' @export

regress <- function(data, y, kernel, p=0.8, C=1, G=0, E=0.1, k) {

  # 1. TR/TE
  N <- nrow(data)
  all.indexes <- 1:N
  learn.indexes <- trainIndx(n=N,ptrain=p)
  test.indexes <- all.indexes[-learn.indexes]
  nlearn <- length(learn.indexes)
  ntest <- N - nlearn

  # 2. Compute kernel matrix
  Jmatrix <- kernelSelect(kernel,data,y)

  trMatrix <- Jmatrix[learn.indexes,learn.indexes]

  # 4. Do R x k-Cross Validation
  if(hasArg(k)) {
    if(k<2) stop("k must be equal to or higher than 2")
    bh <- kCV.reg(GAMMA = G, EPS = E, COST = C, K=trMatrix, Yresp=y[learn.indexes], k=k, R=k)
    cost <- bh$cost
  } else {
    if(length(C)>1) paste("C > 1 - Only the first element will be used")
    cost <- C[1]
  }

  if(G != 0)  {
    trMatrix <- exp(G * trMatrix)/exp(G)
    teMatrix <- exp(G * teMatrix)/exp(G)
  }

  model <- ksvm(trMatrix, y[learn.indexes],type="eps-svr", kernel="matrix", C=cost )

  # 5. Prediction
  teMatrix <- Jmatrix[test.indexes,learn.indexes]
  teMatrix <- teMatrix[,SVindex(model),drop=FALSE]
  teMatrix <- as.kernelMatrix(teMatrix)
  pred <- kernlab::predict(model,teMatrix)
  # print(pred)

  return( error.norm(y[test.indexes],pred))
}

