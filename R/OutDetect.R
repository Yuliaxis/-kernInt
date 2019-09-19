# OUTLIER / NOVELTY DETECTION

#' SVM outlier detection / one class SVM
#'
#' outliers() has two principal usages: unsupervised detection of outliers in data, o supervised one-class SVC.
#'
#' @param data Input data
#' @param y Reponse variable. If a value is provided, outliers() functions as an one-class SVM.
#' @param kernel "qJac" for quantitative Jaccard and "wqJacc" for quant Jaccard with weights.
#' "matrix" if a pre-calculated kernel matrix is given as input.
#' @param nu Hyperparameter nu
#' @param p If a value for y is provided, p is the proportion of total data instances in the training set
#' @param k The k for the k-Cross Validation. Minimum k = 2.
#' @param G Hyperparameter gamma
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
#' outliers(data=speMGX[,7:ncol(speMGX)],y=diag,kernel="qJac",nu=c(0.45,0.5),G=c(0.1,1),k=10)
#' @importFrom kernlab ksvm predict
#' @export



outliers <- function(data,y,kernel,nu,p=0.8,k,G=0) {

  if(hasArg(y)) y <- as.factor(y)

  Jmatrix <- kernelSelect(kernel,data,y)

  if(hasArg(y)){

    N <- nrow(data)
    all.indexes <- 1:N

    learn.indexes <- trainIndx(n=N,ptrain=p)
    test.indexes <- all.indexes[-learn.indexes]

    nlearn <- length(learn.indexes)
    ntest <- N - nlearn

    trMatrix <- Jmatrix[learn.indexes,learn.indexes]
    teMatrix <- Jmatrix[test.indexes,learn.indexes]

    if(hasArg(k)) {
      if(k<2) stop("k must be equal to or higher than 2")
      bh <- kCV.one(K=trMatrix, Yresp=y[learn.indexes], NU=nu, GAMMA=G, k=k, R=k)
      nu <- bh$nu
      g <- bh$g
      print(c(nu,g))
    } else {
      if(length(nu)>1) paste("Nu length > 1 but not k provided - Only the first element will be used")
      nu <- nu[1]
      if(length(G)>1) paste("Gamma length > 1 but not k provided - Only the first element will be used")
      g <- G[1]
    }

    if(g != 0)  {
      trMatrix <- exp(g * trMatrix)/exp(g)
      teMatrix <- exp(g * teMatrix)/exp(g)
    }

    model <- ksvm(trMatrix,nu=nu, type="one-svc", kernel="matrix")

    teMatrix <- teMatrix[,SVindex(model),drop=FALSE]
    teMatrix <- as.kernelMatrix(teMatrix)
    pred <- kernlab::predict(model,teMatrix)

    pred <- as.factor(as.numeric(pred))

    levels(y) <- c("0","1")
    levels(pred) <- c("0","1")

    ct <- table(Truth=y[test.indexes], Pred=pred)
    return(ct)

  } else {
    model <- ksvm(Jmatrix,nu=nu, type="one-svc", kernel="matrix")
    get_index <- kernlab::predict(model)
    return(which(get_index[,1]==TRUE))
  }

}

