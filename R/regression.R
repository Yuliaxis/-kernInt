# REGRESSION

#' SVM regression
#'
#'regress() automatically trains a Support Vector Regrssion model, tests it and returns the Normalized Mean Squared Error.
#'
# Cross-validation is available to choose the best hyperparameters (e.g. Cost, Epsilon) during the training step.
#' If the input data has repeated rownames, classify() will consider that the row names that share id are repeated
#' measures from the same individual. The function will ensure that all repeated measures are used either to train
#' or to test the model, but not for both, thus preserving the independence between the training and tets sets.
#'
#' @param data Input data
#' @param y Reponse variable (continuous)
#' @param kernel "cRBF" for clrRBF, "qJac" for quantitative Jaccard and  "wqJacc" for quantitative Jaccard with weights.
#' "matrix" if a pre-calculated kernel matrix is given as input.
#' @param p Proportion of total data instances in the training set
#' @param C A vector with the possible costs to evaluate via k-Cross-Val. If no argument is provided cross-validation is not performed.
#' @param H Gamma hyperparameter
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

regress <- function(data, y, kernel, p=0.8, C=1, H=0, E=0.1, k) {

  if(class(data) == "list") {
    m <- length(data)
    if(m < 2) data <- unlist(data)
  } else if(class(data) == "array") {
    m <- dim(data)[3]
    if(m < 2) data <- unlist(data)
    data <- matrix(data[,,1],ncol=dim(data)[2],nrow=dim(data)[1])
  } else if(class(data) == "data.frame" | class(data) == "matrix") {
    m <- 1
  } else {
    stop("Wrong input data class.")
  }

  # 1. TR/TE
  index <- finalTRTE(data,p)
  learn.indexes <- index$li
  test.indexes <- index$ti

  # 2. Compute kernel matrix
  Jmatrix <- kernelSelect(kernel,data,y)

  trMatrix <- Jmatrix[learn.indexes,learn.indexes]
  teMatrix <- Jmatrix[test.indexes,learn.indexes]

  # 3. Do R x k-Cross Validation
  if(hasArg(k)) {
    if(k<2) stop("k must be equal to or higher than 2")
    bh <- kCV.core(H = H, method="svr", kernel=kernel,EPS = E, COST = C, K=trMatrix, Y=y[learn.indexes], k=k, R=k)
    cost <- bh$cost
    eps <- bh$epsilon
    H <- bh$h
  } else {
    if(length(C)>1)  paste("C > 1 - Only the first element will be used")
    if(length(E)>1) paste("E > 1 - Only the first element will be used")
    if(length(H)>1) paste("H > 1 - Only the first element will be used")
    cost <- C[1]
    eps <- E[1]
    H <- H[1]
  }

  trMatrix <- hyperkSelection(trMatrix,h=H,kernel=kernel)
  teMatrix <- hyperkSelection(teMatrix,h=H,kernel=kernel)

  model <- ksvm(trMatrix, y[learn.indexes],type="eps-svr", kernel="matrix", C=cost, epsilon = eps)

  # 5. Prediction
  teMatrix <- teMatrix[,SVindex(model),drop=FALSE]
  teMatrix <- as.kernelMatrix(teMatrix)
  pred <- kernlab::predict(model,teMatrix)

  return( error.norm(y[test.indexes],pred))
}

