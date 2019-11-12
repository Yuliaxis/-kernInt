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
#' @param data Input data: a matrix or data.frame with predictor variables. To perform MKL: a list of the *m* types of data to combine.
#' @param y Reponse variable (continuous)
#' @param kernel "linear" for linear kernel, cRBF" for clrRBF, "qJac" for quantitative Jaccard and  "wqJacc" for quantitative Jaccard with weights.
#' "matrix" if a pre-calculated kernel matrix is given as input. To perform MKL: Vector of *m* kernels to apply to each data type.
#' @param coeff ONLY IN MKL CASE: A *t·m* matrix of the coefficients, where *m* are the number of different data types and *t* the number of
#' different coefficient combinations to evaluate via k-CV. If absent, the same weight is given to all data sources.
#' @param p Proportion of total data instances in the training set
#' @param C A cost, or a vector with the possible costs to evaluate via k-Cross-Val.
#' @param H Gamma hyperparameter. A vector with the possible values to chose the best one via k-Cross-Val can be entered.
#' For the MKL, a list with *m* entries can be entered, being' *m* is the number of different data types. Each element on the list
#' must be a number or, if k-Cross-Validation is needed, a vector with the hyperparameters to evaluate for each data type.
#' @param E Epsilon hyperparameter, or a vector with the possible epsilons to evaluate via k-Cross-Val.
#' @param k The k for the k-Cross Validation. Minimum k = 2. If no argument is provided cross-validation is not performed.
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

regress <- function(data, y, coeff,  kernel, p=0.8, C=1, H=NULL, E=0.1, k) {

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
  if("time2" %in% kernel || "time" %in% kernel ) {
    print("Longitudinal")
    index <- longTRTE(data,p)
  } else {
    index <- finalTRTE(data,p) ## data és una matriu en aquest cas. passar-ho a MKL.
  }
  learn.indexes <- index$li
  test.indexes <- index$ti
  print(test.indexes)

  # 2. Compute kernel matrix

  if(m>1) {
    Jmatrix  <- seqEval(DATA=data, y=y, kernels=kernel,h=NULL) ## Sense especificar hiperparàmetre.
    trMatrix <- Jmatrix[learn.indexes,learn.indexes,]
    teMatrix <- Jmatrix[test.indexes,learn.indexes,]
  } else {
    print(dim(data))
    print(length(y))
    Jmatrix <- kernelSelect(kernel=kernel,data=data,y=y,h=NULL)
    trMatrix <- Jmatrix[learn.indexes,learn.indexes]
    teMatrix <- Jmatrix[test.indexes,learn.indexes]
  }

  # 3. Do R x k-Cross Validation
  if(hasArg(k)) {
    if(k<2) stop("k should be equal to or higher than 2")
    if(m>1)  {
      if(!hasArg(coeff)) coeff <- rep(1/m,m)
      bh <- kCV.MKL(ARRAY=trMatrix, COEFF=coeff, KERNH=H, kernels=kernel, method="svr", COST = C,EPS = E,
                     Y=y[learn.indexes], k=k,  R=1)
      coeff <- bh$coeff ##indexs
      print(coeff)

    } else {
      bh <- kCV.core(H = H, method="svr", kernel=kernel,EPS = E, COST = C, K=trMatrix, Y=y[learn.indexes], k=k, R=1)
    }
    eps <- bh$epsilon
    cost <- bh$cost
    H <- bh$h

  } else {
    if(length(C)>1)  paste("C > 1 - Only the first element will be used")
    if(length(E)>1) paste("E > 1 - Only the first element will be used")
    # if(length(H)>1) paste("H > 1 - Only the first element will be used")
    cost <- C[1]
    eps <- E[1]
    if(!is.null(H)) H <- kernHelp(H)$hyp

    # H <- H[1]
  }

  if(m>1) {
    for(j in 1:m) trMatrix[,,j] <- hyperkSelection(K=trMatrix[,,j], h=H[j],  kernel=kernel[j])
    for(j in 1:m) teMatrix[,,j] <- hyperkSelection(K=teMatrix[,,j], h=H[j],  kernel=kernel[j])
    print(trMatrix[1:10,1:10,4])

    trMatrix <- KInt(data=trMatrix,coeff=coeff)
    teMatrix <- KInt(data=teMatrix,coeff=coeff)
  }  else {
    trMatrix <- hyperkSelection(trMatrix,h=H,kernel=kernel)
    teMatrix <- hyperkSelection(teMatrix,h=H,kernel=kernel)
  }
  print(trMatrix[1:10,1:10])


  model <- ksvm(trMatrix, y[learn.indexes],type="eps-svr", kernel="matrix", C=cost, epsilon = eps)

  # 5. Prediction
  teMatrix <- teMatrix[,SVindex(model),drop=FALSE]
  teMatrix <- as.kernelMatrix(teMatrix)
  pred <- kernlab::predict(model,teMatrix)
  print(pred)
  print(y[test.indexes])
  return( error.norm(y[test.indexes],pred))
}

