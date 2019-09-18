# OUTLIER / NOVELTY DETECTION

#' SVM outlier detection
#'
#'
#' @param data Input data
#' @param y Reponse variable
#' @param kernel "qJac" for quantitative Jaccard and "wqJacc" for quant Jaccard with weights
#' @param nu Hyperparameter nu
#' @param G Hyperparameter gamma
#' @param p If a value is provided, outliers() functions as an one-class SVM. p is the proportion of total data instances in the training set
#' @param k The k for the k-Cross Validation. Minimum k = 2.
#' @return The indexes of the outliers (outlier detection) or, if a value is provided for p, the confusion matrix (one-class SVM)
#' @examples
#' outliers(data=soilDataRaw,kernel="cRBF",nu=0.3)
#'
#' diag <- as.numeric(speMGX[,1])
#' diag[diag == 3] <- 1  # De 3 a 2 classes: No Malalt /  malalt
#' diag[diag == 2] <- 0
#' diag <- as.numeric(!diag) # No malalt classe 1, malalt classe 0
#' outliers(data=speMGX[,7:ncol(speMGX)],y=diag,kernel="qJac",nu=0.2,p=0.8,k=10)
#' outliers(data=speMGX[,7:ncol(speMGX)],y=diag,kernel="qJac",nu=c(0.1,0.2,0.3),G=c(0.01,0.1,1,5,10),p=0.8,k=10)
#' @importFrom kernlab ksvm predict
#' @export


outliers <- function(data,y,kernel,nu,p,k,G=0) {

  y <- as.factor(y)

  if(kernel == "qJac") {
    Jmatrix <- qJacc(data)
  } else if(kernel == "wqJac") {
    Jmatrix <- wqJacc(data,y=y)
    cat("quantJaccard kernel + weights \n")
  }  else if(kernel == "cRBF") {
    Jmatrix <- aitch.dist(data)
    cat("clr + RBF \n")
  }   else {
    cat("standard RBF \n")

    # Jmatrix <-  kernelMatrix(rbfdot(sigma = g),data)

  }

  if(hasArg(p)){

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
      if(length(nu)>1) paste("Hyperparameters length > 1 but not k provided - Only the first element will be used")
      nu <- nu[1]
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

