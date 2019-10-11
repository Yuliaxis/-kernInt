# OUTLIER / NOVELTY DETECTION

#' SVM outlier detection / one class SVM
#'
#' outliers() has two principal usages: unsupervised detection of outliers in data, o supervised one-class SVC.
#'
#' If outliers() is used in a supervised way and the input data has repeated rownames, classify() will consider that the row names that share id are repeated
#' measures from the same individual. The function will ensure that all repeated measures are used either to train
#' or to test the model, but not for both, thus preserving the independence between the training and tets sets.
#'
#' @param data Input data
#' @param y Reponse variable. If a value is provided, outliers() functions as an one-class SVM.
#' @param kernel "qJac" for quantitative Jaccard and "wqJacc" for quant Jaccard with weights.
#' "matrix" if a pre-calculated kernel matrix is given as input.
#' @param nu Hyperparameter nu
#' @param p If a value for y is provided, p is the proportion of total data instances in the training set
#' @param k The k for the k-Cross Validation. Minimum k = 2.
#' @param H Hyperparameter gamma
#' @return The indexes of the outliers (outlier detection) or, if a value is provided for y, the confusion matrix (one-class SVM)
#' @examples
#' # Outlier detection
#' outliers(data=soilDataRaw,kernel="cRBF",nu=0.3)
#' ## One-class SVM:
#' ## Preparing the y
#' diag <- as.numeric(speMGX[,1])
#' diag[diag == 3] <- 1  # De 3 a 2 classes: No Malalt /  malalt
#' diag[diag == 2] <- 0
#' diag <- as.numeric(!diag) # No malalt classe 1, malalt classe 0
#' ## One-class SVM changing the percentage of data for training (70%) and the hyperparameter nu:
#' outliers(data=speMGX[,7:ncol(speMGX)],y=diag,kernel="qJac",nu=0.2,p=0.7)
#' ## One-class SVM with 10-Cross-Validation:
#' outliers(data=speMGX[,7:ncol(speMGX)],y=diag,kernel="qJac",nu=c(0.45,0.5),H=c(0.1,1),k=10)
#' @importFrom kernlab ksvm predict
#' @export



outliers <- function(data,y,kernel,nu,p=0.8,k,H=0) {

  if(hasArg(y)) {
    y <- as.factor(y)
    levels(y) <- c("0","1")
  }

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

  Jmatrix <- kernelSelect(kernel,data,y)

  if(hasArg(y)){
    index <- finalTRTE(data,p) ## data és una matriu en aquest cas. passar-ho a MKL.
    learn.indexes <- index$li
    test.indexes <- index$ti

    trMatrix <- Jmatrix[learn.indexes,learn.indexes]
    teMatrix <- Jmatrix[test.indexes,learn.indexes]

    if(hasArg(k)) {
      if(k<2) stop("k must be equal to or higher than 2")
      bh <- kCV.core(method="one",K=trMatrix,  kernel=kernel,Y=y[learn.indexes], NU=nu, H=H, k=k, R=k)
      nu <- bh$nu
      H <- bh$h
      print(c(nu,H))
    } else {
      if(length(nu)>1) paste("Nu length > 1 but not k provided - Only the first element will be used")
      nu <- nu[1]
      if(length(H)>1) paste("Gamma length > 1 but not k provided - Only the first element will be used")
      H <- H[1]
    }

    trMatrix <- hyperkSelection(trMatrix,h=H,kernel=kernel)
    teMatrix <- hyperkSelection(teMatrix,h=H,kernel=kernel)

    model <- ksvm(trMatrix,nu=nu, type="one-svc", kernel="matrix")

    teMatrix <- teMatrix[,SVindex(model),drop=FALSE]
    teMatrix <- as.kernelMatrix(teMatrix)
    pred <- kernlab::predict(model,teMatrix)

    pred <- as.factor(as.numeric(pred))
    levels(pred) <- c("0","1")

    ct <- table(Truth=y[test.indexes], Pred=pred)
    return(ct)

  } else {
    model <- ksvm(Jmatrix,nu=nu, type="one-svc", kernel="matrix")
    get_index <- kernlab::predict(model)
    return(which(get_index[,1]==TRUE))
  }
}

