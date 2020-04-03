# OUTLIER / NOVELTY DETECTION

#' SVM outlier detection / one class SVM
#'
#' outliers() has two principal usages: unsupervised detection of outliers in data, o supervised one-class SVC.
#'
#' If outliers() is used in a supervised way and the input data has repeated rownames, classify() will consider that the row names that share id are repeated
#' measures from the same individual. The function will ensure that all repeated measures are used either to train
#' or to test the model, but not for both, thus preserving the independence between the training and tets sets.
#'
#' @param data Input data: a matrix or data.frame with predictor variables/features as columns.
#' To perform MKL: a list of *m* datasets. All datasets should have the same number of rows
#' @param y Reponse variable (factor)
#' @param kernel "lin" or rbf" to standard Linear and RBF kernels. "clin" for compositional linear and "crbf" for Aitchison-RBF
#' kernels. "jac" for quantitative Jaccard / Ruzicka kernel. "jsk" for Jensen-Shannon Kernel. "flin" and "frbf" for functional linear
#' and functional RBF kernels. "matrix" if a pre-computed kernel matrix is given as input.
#' To perform MKL: Vector of *m* kernels to apply to each dataset.
#' @param nu Hyperparameter nu
#' @param p If a value for y is provided, p is the proportion of total data instances in the test set
#' @param k The k for the k-Cross Validation. Minimum k = 2.
#' @param H Hyperparameter gamma
#' @param domain Only used in "frbf" or "flin".
#' @return The indexes of the outliers (outlier detection) or, if a value is provided for y, the confusion matrix (one-class SVM)
#' @examples
#' # Outlier detection
#' outliers(data=soil$abund,kernel="clin",nu=0.2)
#' ## One-class SVM:
#' outliers(data=soil$abund ,y=soil$metadata[ ,"env_feature"],kernel="clin")
#' ## One-class SVM with 10-Cross-Validation:
#' outliers(data=soil$abund ,y=soil$metadata[ ,"env_feature"],kernel="clin",nu=c(0.45,0.5),k=10)
#' @importFrom kernlab ksvm predict
#' @importFrom methods hasArg
#' @export



outliers <- function(data,y,kernel,nu=0.2,p=0.2,k,domain=NULL,H=NULL) {

  if(hasArg(y)) y <- as.factor(y)

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

  Jmatrix <- kernelSelect(kernel,data,domain,h=NULL)

  if(hasArg(y)){
    index <- finalTRTE(data,1-p) ## data és una matriu en aquest cas. passar-ho a MKL.
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

