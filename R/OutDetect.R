# OUTLIER / NOVELTY DETECTION

#' SVM outlier detection
#'
#'
#' @param data Input data
#' @param y Reponse variable
#' @param kernel "qJac" for quantitative Jaccard and "wqJacc" for quant Jaccard with weights
#' @param nu Hyperparameter nu
#' @param p If a value is provided, outliers() functions as an one-class SVM. p is the proportion of total data instances in the training set
#' @return The indexes of the outliers
#' @examples
#' outliers(data=soilDataRaw,kernel="cRBF",nu=0.3)
#' diag <- as.numeric(speMGX[,1])
#' diag[diag == 3] <- 1  # De 3 a 2 classes: No Malalt /  malalt
#' diag[diag == 2] <- 0
#' diag <- as.numeric(!diag) # No malalt classe 1, malalt classe 0
#' outliers(data=speMGX[,7:ncol(speMGX)],y=diag,kernel="cRBF",nu=0.2,p=0.8)
#' @importFrom kernlab ksvm predict
#' @export


outliers <- function(data,y,kernel,nu,p) {

  if(kernel == "qJac") {
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

  if(hasArg(p)){

    N <- nrow(data)
    all.indexes <- 1:N

    learn.indexes <- trainIndx(n=N,ptrain=p)
    test.indexes <- all.indexes[-learn.indexes]

    nlearn <- length(learn.indexes)
    ntest <- N - nlearn

    trMatrix <- Jmatrix[learn.indexes,learn.indexes]

    model <- ksvm(trMatrix,nu=nu, type="one-svc", kernel="matrix")
    teMatrix <- Jmatrix[test.indexes,learn.indexes]

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

