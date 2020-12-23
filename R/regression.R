# REGRESSION

#' SVM regression
#'
#'regress() automatically trains a Support Vector Regression model, tests it and returns the Normalized Mean Squared Error.
#'
# Cross-validation is available to choose the best hyperparameters (e.g. Cost, Epsilon) during the training step.
#' If the input data has repeated rownames, classify() will consider that the row names that share id are repeated
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
#' @param C A cost, or a vector with the possible costs to evaluate via k-Cross-Val.
#' @param H Gamma hyperparameter (only in RBF-like functions). A vector with the possible values to chose the best one via k-Cross-Val can be entered.
#' For the MKL, a list with *m* entries can be entered, being' *m* is the number of different data types. Each element on the list
#' must be a number or, if k-Cross-Validation is needed, a vector with the hyperparameters to evaluate for each data type.
#' @param E Epsilon hyperparameter, or a vector with the possible epsilons to evaluate via k-Cross-Val.
#' @param k The k for the k-Cross Validation. Minimum k = 2. If no argument is provided cross-validation is not performed.
#' @param domain Only used in "frbf" or "flin".
#' @return NMSE (normalized mean squared error), chosen hyperparameters, test set predicted and observed values, and variable importances (only with linear-like kernels)
#' @examples
#' # Simple regression without tuning the hyperparameters
#' regress(data=soil$abun,soil$metadata$ph,kernel="clin")
#' # The percentage of data for training can be changed (default: 0.8 Training / 0.2 Test):
#' regress(data=soil$abun,soil$metadata$ph,kernel="clin",p=0.6)
#' # Regression with 10-Cross-Validation to choose the best Cost and Epsilon:
#' regress(data=soil$abun,soil$metadata$ph,kernel="clin", C=c(0.1,1,10), E = c(0.01,0.1), k=10)
#' # Regression with MKL:
#' Nose <- list()
#' Nose$left <- CSSnorm(smoker$abund$nasL)
#' Nose$right <- CSSnorm(smoker$abund$nasR)
#' age <- smoker$metadata$age[seq(from=1,to=62*4,by=4)]
#' w <- matrix(c(0.5,0.1,0.9,0.5,0.9,0.1),nrow=3,ncol=2)
#' regress(data=Nose,kernel="jac",y=age,C=c(1,10,100), coeff = w, k=10)
#' @importFrom kernlab alpha alphaindex as.kernelMatrix kernelMatrix predict rbfdot SVindex
#' @importFrom methods hasArg is
#' @export

regress <- function(data, y,  coeff,  kernel, p=0.2,  C=1, H=NULL, E=0.01, domain=NULL, k) {

  ### Checking data
  check <- checkinput(data,kernel)
  m <- check$m
  data <- check$data
  kernel <- check$kernel

  if(length(y) != check$n) stop("Length of the target variable do not match with the row number of predictors")

  # 1. TR/TE
  inds <- checkp(p=p,data=data)
  learn.indexes <- inds$learn.indexes
  test.indexes <- inds$test.indexes

  try <- y[learn.indexes]
  tey <- y[test.indexes]


  # 2. Compute kernel matrix

  if(m>1) {
    Jmatrix<- seqEval(DATA=data,domain=domain, kernels=kernel,h=NULL) ## Sense especificar hiperparàmetre.
    trMatrix <- Jmatrix[learn.indexes,learn.indexes,]
    teMatrix <- Jmatrix[test.indexes,learn.indexes,]
    if(!hasArg(coeff)) coeff <- rep(1/m,m)
  } else {
    Jmatrix <- kernelSelect(kernel=kernel,domain=domain,data=data,h=NULL)
    trMatrix <- Jmatrix[learn.indexes,learn.indexes]
    teMatrix <- Jmatrix[test.indexes,learn.indexes]
  }


  # 3. Do R x k-Cross Validation
  if(hasArg(k)) {
    if(k<2) stop("k should be equal to or higher than 2")
    if(m>1)  {
      bh <- kCV.MKL(ARRAY=trMatrix, COEFF=coeff, KERNH=H, kernels=kernel, method="svr", COST = C,EPS = E,
                     Y=try, k=k,  R=k)
      coeff <- bh$coeff ##indexs
      conserv <- c("coeff","cost","error")
      if(!is.null(H)) conserv <- c(conserv,"h")
      bh <- bh[conserv]
    } else {
      bh <- kCV.core(H = H, method="svr", kernel=kernel,EPS = E, COST = C, K=trMatrix, Y=try, k=k, R=k)
    }
    eps <- bh$eps
    cost <- bh$cost
    H <- bh$h
  } else {
    if(length(C)>1)  warning("Multiple C and no k - Only the first element will be used")
    if(length(E)>1) warning("Multiple E and no k - Only the first element will be used")
    if(!is.null(H))   {
      H <- kernHelp(H)$hyp
      bh <- data.frame(H)
    } else {
      bh <- NULL
    }
    cost <- C[1]
    eps <- E[1]
    if(m>1)   {
      bh <- list(coeff=coeff,h=H,cost=cost,eps=eps)
    }  else {
      bh <- cbind(bh,cost,eps)
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

  model <- ksvm(trMatrix, try,type="eps-svr", kernel="matrix", C=cost, epsilon = eps)

  ##### Importances (only linear)

  alphaids <- alphaindex(model) # Indices of SVs in original data
  alphaids <-  learn.indexes[unlist(alphaids)]
  alphas <- alpha(model)

  importances <- impCore(kernel=kernel,alphas=alphas,alphaids=alphaids,data=data,ys=as.numeric(y[alphaids]), m=m, coeff=coeff)

  # 5. Prediction
  teMatrix <- teMatrix[,SVindex(model),drop=FALSE]
  teMatrix <- as.kernelMatrix(teMatrix)
  pred <- kernlab::predict(model,teMatrix)
  err <- error.norm(tey,pred)
  test <- data.frame(true=tey,predicted=pred)
  rownames(test) <- test.indexes
  return(list("nmse"=err,"hyperparam"=bh,"prediction"=test,"var.imp"=importances))
}

