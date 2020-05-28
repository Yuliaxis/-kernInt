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
#' @param y Reponse variable (continuous)
#' @param kernel "lin" or rbf" to standard Linear and RBF kernels. "clin" for compositional linear and "crbf" for Aitchison-RBF
#' kernels. "jac" for quantitative Jaccard / Ruzicka kernel. "jsk" for Jensen-Shannon Kernel. "flin" and "frbf" for functional linear
#' and functional RBF kernels. "matrix" if a pre-computed kernel matrix is given as input.
#' To perform MKL: Vector of *m* kernels to apply to each dataset.
#' @param coeff ONLY IN MKL CASE: A *t·m* matrix of the coefficients, where *m* are the number of different data types and *t* the number of
#' different coefficient combinations to evaluate via k-CV. If absent, the same weight is given to all data sources.
#' @param p The proportion of data reserved for the test set. Otherwise, a vector containing the indexes or the names of the rows for testing.
#' @param H Gamma hyperparameter (only in RBF-like functions). A vector with the possible values to chose the best one via k-Cross-Val can be entered.
#' For the MKL, a list with *m* entries can be entered, being' *m* is the number of different data types. Each element on the list
#' must be a number or, if k-Cross-Validation is needed, a vector with the hyperparameters to evaluate for each data type.
#' @param nu Hyperparameter nu
#' @param k The k for the k-Cross Validation. Minimum k = 2. If no argument is provided cross-validation is not performed.
#' @param domain Only used in "frbf" or "flin".
#' @return The indexes of the outliers (outlier detection) or, if a value is provided for y, the confusion matrix (one-class SVM)
#' @examples
#' # Outlier detection
#' outliers(data=soil$abund,kernel="clin",nu=0.2)
#' ## One-class SVM:
#' outliers(data=soil$abund ,y=soil$metadata[ ,"env_feature"],kernel="clin")
#' ## One-class SVM with 10-Cross-Validation:
#' outliers(data=soil$abund ,y=soil$metadata[ ,"env_feature"],kernel="clin",nu=c(0.45,0.5),k=10)
#' ## With data of multiple sources
#'outliers(data=smoker$abund,kernel="crbf",H=list(nasL=0.01,nasR=0.01,oroL=0.1,oroR=0.1))
#' @importFrom kernlab as.kernelMatrix ksvm predict
#' @importFrom methods hasArg
#' @export



outliers <- function(data,y=NULL, kernel, coeff, nu=0.2,p=0.2,k,domain=NULL,H=NULL) {

  ### Checking data
  check <- checkinput(data,kernel)
  m <- check$m
  data <- check$data
  kernel <- check$kernel

  # 2. Compute kernel matrix

   Jmatrix<- seqEval(DATA=data,domain=domain, kernels=kernel,h=NULL) ## Sense especificar hiperparàmetre.
   if(!hasArg(coeff)) coeff <- rep(1/m,m)

  if(hasArg(y)){
    if(length(y) != check$n) stop("Length of the target variable do not match with the row number of predictors")
    y <- as.factor(y)

    inds <- checkp(p=p,data=data)
    learn.indexes <- inds$learn.indexes
    test.indexes <- inds$test.indexes

    try <- y[learn.indexes]
    tey <- y[test.indexes]

    if(m>1) {
      trMatrix <- Jmatrix[learn.indexes,learn.indexes,]
      teMatrix <- Jmatrix[test.indexes,learn.indexes,]
    } else {
      trMatrix <- Jmatrix[learn.indexes,learn.indexes]
      teMatrix <- Jmatrix[test.indexes,learn.indexes]
    }

    if(hasArg(k)) {
      if(k<2) stop("k should be equal to or higher than 2")
      if(m>1)  {
        bh <- kCV.MKL(ARRAY=trMatrix, COEFF=coeff, KERNH=H, kernels=kernel, method="one-svc", NU = nu,
                      Y=try, k=k, R=k)
        coeff <- bh$coeff ##indexs
        conserv <- c("coeff","nu","error")
        if(!is.null(H)) conserv <- c(conserv,"h")
        bh <- bh[conserv]
      } else {
        bh <- kCV.core(method="one-svc",NU = nu, H = H, kernel=kernel, K=trMatrix,
                       Y=try, k=k, R=k)
        bh <- bh[-which(is.na(bh))]
      }
      nu <- bh$nu
      H <- bh$h

    } else {
      if(!is.null(H))   {
        H <- kernHelp(H)$hyp
        bh <- data.frame(H)
      } else {
        bh <- NULL
      }
      if(length(nu)>1) warning("Multiple nu and no k provided - Only the first element will be used")
      nu <- nu[1]
      if(m>1)   {
        bh <- list(coeff=coeff,h=H,nu=nu)
      }  else {
        bh <- cbind(bh,nu)
      }
    }

    if(m>1) {
      for(j in 1:m) trMatrix[,,j] <- hyperkSelection(K=trMatrix[,,j], h=H[j],  kernel=kernel[j])
      for(j in 1:m) teMatrix[,,j] <- hyperkSelection(K=teMatrix[,,j], h=H[j],  kernel=kernel[j])
      trMatrix <- KInt(data=trMatrix,coeff=coeff)
      teMatrix <- KInt(data=teMatrix,coeff=coeff)
    }  else {
      trMatrix <- hyperkSelection(trMatrix,h=H,kernel=kernel)
      teMatrix <- hyperkSelection(teMatrix,h=H,kernel=kernel)
    }

    model <- ksvm(trMatrix,nu=nu, type="one-svc", kernel="matrix")

    teMatrix <- teMatrix[,SVindex(model),drop=FALSE]
    teMatrix <- as.kernelMatrix(teMatrix)
    pred <- kernlab::predict(model,teMatrix)
    pred <- as.factor(pred)
    levels(pred) <- levels(y)
    ct <- table(Truth=tey, Pred=pred)  ### Confusion matrix
    test <- data.frame(true=tey,predicted=pred)
    return(list("conf.matrix"=ct,"hyperparam"=bh,"prediction"=test))
  } else {
    if(m>1) {
      H <- unlist(H)
      for(j in 1:m)Jmatrix[,,j] <- hyperkSelection(K=Jmatrix[,,j], h=H[j],  kernel=kernel[j])
      Jmatrix <- KInt(data=Jmatrix,coeff=coeff)
    }
    Jmatrix <- as.kernelMatrix(Jmatrix)
    model <- ksvm(Jmatrix,nu=nu, type="one-svc", kernel="matrix")
    get_index <- kernlab::predict(model)
    return(which(get_index[,1]))
  }
}

