# CLASSIFIER

#' SVM classifier
#'
#' classify() automatically trains a Support Vector Classification model, tests it and returns the confusion matrix.
#'
#' Cross-validation is available to choose the best hyperparameters (e.g. Cost) during the training step.
#'
#' The classification can be hard (predicting the class) or soft (predicting the probability of belonging to a given class)
#'
#' Another feature is the possibility to deal with imbalanced data in the target variable with several techniques:
#' \describe{
#'   \item{Data Resampling}{Oversampling techniques (oversample the minority class, generating synthetic data with SMOTE)
#'   or undersampling the majority class.}
#'   \item{Class weighting}{Giving more weight to the minority class}
#' }
#' To use one-class SVM to deal with imbalanced data, see: outliers()
#'
#' If the input data has repeated rownames, classify() will consider that the row names that share id are repeated
#' measures from the same individual. The function will ensure that all repeated measures are used either to train
#' or to test the model, but not for both, thus preserving the independance between the training and tets sets.
#'
#' Currently, the classification can be only performed if the target variable is binary (two classes).
#'
#' @param data Input data: a matrix or data.frame with predictor variables. To perform MKL: a list of the *m* types of data to combine.
#' @param y Reponse variable (binary)
#' @param kernel "cRBF" for clrRBF, "qJac" for quantitative Jaccard,  "wqJac" for quantitative Jaccard with weights.
#' "matrix" if a pre-calculated kernel matrix is given as input. To perform MKL: Vector of *m* kernels to apply to each data type.
#' @param coeff ONLY IN MKL CASE: A *t·m* matrix of the coefficients, where *m* are the number of different data types and *t* the number of
#' different coefficient combinations to evaluate via k-CV. If absent, the same weight is given to all data sources.
#' @param prob if TRUE class probabilities (soft-classifier) are computed instead of a True-or-false assignation (hard-classifier)
#' @param classimb "weights" to introduce class weights in the SVM algorithm and "data" to oversampling. If other arguments are provided nothing is done.
#' @param type If classimb = "data", the procedure to data oversampling or undersampling ("ubOver","ubUnder" or "ubSMOTE")
#' @param p Proportion of total data instances in the training set
#' @param k The k for the k-Cross Validation. Minimum k = 2. If no argument is provided cross-validation is not performed.
#' @param C The cost. A vector with the possible costs (SVM hyperparameter) to evaluate via k-Cross-Val can be entered too.
#' @param H Gamma hyperparameter. A vector with the possible values to chose the best one via k-Cross-Val can be entered.
#' For the MKL, a list with *m* entries can be entered, being' *m* is the number of different data types. Each element on the list
#' must be a number or, if k-Cross-Validation is needed, a vector with the hyperparameters to evaluate for each data type.
#' @param CUT Cut-off if prob = TRUE. If CUT is a vector, the best cut-off can be obtained via cross-validation.
#' @return Confusion matrix or, if prob = TRUE and not cutoff is set, a data.frame with the class probability and the actual class.
#' @examples
#' #Preparing the target variable:
#' diag <- as.numeric(speMGX[,1])
#' diag[diag == 3] <- 1  # 3 classes to 2: Disease (1) /  no Disease (2)
#' # Classification with 10-Cross-Validation
#' classify(data=speMGX[,7:ncol(speMGX)],diag,kernel="qJac",C=c(0.1,1),k=10)
#' # Probabilistic Classification with no cross-validation:
#' classify(data=speMGX[,7:ncol(speMGX)],diag,kernel="wqJac",prob = TRUE)
#' # Classification (with 10-CV) accounting for the imbalanced data: an example of data resampling
#' classify(data=speMGX[,7:ncol(speMGX)],diag,kernel="qJac",classimb="data",type="ubUnder",C=c(0.001,0.01),k=10)
#' # Classification (with 10-CV) accounting for the imbalanced data: class weighting
#' classify(data=speMGX[,7:ncol(speMGX)],diag,kernel="qJac",classimb="weights",C=c(0.001,0.01),k=10)
#' @importFrom kernlab as.kernelMatrix kernelMatrix predict rbfdot SVindex
#' @importFrom unbalanced ubBalance
#' @importFrom ROSE roc.curve
#' @export



classify <- function(data, y, coeff, kernel,  prob=FALSE, classimb="no", type="ubOver", p=0.8, k, C=1, H=NULL, CUT=NULL) {

  # y class
  diagn <- as.factor(y)
  levels(diagn) <- c("1","2")
  # data class
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
  if(classimb == "weights") {
    wei <- c("1"=as.numeric(summary(diagn[learn.indexes])[2]),"2"=as.numeric(summary(diagn[learn.indexes])[1]))
  } else {
    wei <- NULL
  }
print(test.indexes)
  if(classimb=="data")  {
    s <- sampl(data=data,diagn=diagn,learn.indexes=learn.indexes,test.indexes=test.indexes,kernel=kernel,type=type)
    data <- s$data
    diagn <- s$y
    learn.indexes <- s$li
    test.indexes <- s$ti
  }

    # 2. Compute kernel matrix

  if(m>1) {
    Jmatrix  <- seqEval(DATA=data, y=diagn, kernels=kernel,h=NULL) ## Sense especificar hiperparàmetre.
    trMatrix <- Jmatrix[learn.indexes,learn.indexes,]
    teMatrix <- Jmatrix[test.indexes,learn.indexes,]
  } else {
    Jmatrix <- kernelSelect(kernel=kernel,data=data,y=diagn,h=NULL)
    trMatrix <- Jmatrix[learn.indexes,learn.indexes]
    teMatrix <- Jmatrix[test.indexes,learn.indexes]
  }

  # 3. Do R x k-Cross Validation

  if(hasArg(k)) {
    if(k<2) stop("k should be equal to or higher than 2")
    if(m>1)  {
      if(!hasArg(coeff)) coeff <- rep(1/m,m)
      bh <- kCV.MKL(ARRAY=trMatrix, COEFF=coeff, KERNH=H, kernels=kernel, method="svc", COST = C,
                    CUT=CUT, Y=diagn[learn.indexes], k=k,  prob=prob, R=1,classimb=wei)
      coeff <- bh$coeff ##indexs
      print(coeff)

    } else {
    bh <- kCV.core(method="svc",COST = C, H = H, kernel=kernel, CUT=CUT, K=trMatrix, prob=prob,
                   Y=diagn[learn.indexes], k=k, R=1,classimb=wei)
    }
    CUT <- bh$cut
    cost <- bh$cost
    H <- bh$h

  } else {
    if(length(C)>1) paste("C > 1 and no k provided - Only the first element will be used")
    cost <- C[1]
    # if(length(H)>1) paste("H > 1 and no k provided- Only the first element will be used")
    # H <- H[1]
    if(!is.null(CUT) && length(CUT)>1) {
      paste("CUT > 1 and no k provided - Only the first element will be used")
      CUT <- CUT[1]
    }
    if(!is.null(H)) H <- kernHelp(H)$hyp

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

  # 4. Model
  print(H)
  print(dim(trMatrix))
  print(trMatrix[1:10,1:10])
  print(length(diagn[learn.indexes]))
  print(prob)
  print(cost)
  print(wei)
  model <- ksvm(trMatrix, diagn[learn.indexes], kernel="matrix", type="C-svc", prob.model = prob, C=cost, class.weights=wei)

  # 5. Prediction

  teMatrix <- teMatrix[,SVindex(model),drop=FALSE]
  teMatrix <- as.kernelMatrix(teMatrix)

  if(prob)  {
    pred <- predict(model,teMatrix,type = "probabilities")
    if(is.null(CUT)) {
      return(data.frame(Actual=diagn[test.indexes],Predicted = as.factor(pred)))
    }
    print(paste("Best cut is", cut))
    pred <- (pred[,2] < cut)
    pred[pred] <- 1
    pred[pred==0] <- 2
  }   else    {
      pred <- kernlab::predict(model,teMatrix)
    }
   pred <- as.factor(pred)
   levels(pred) <- c("1","2")
   print(pred)

  ### Confusion matrix
  ct <- table(Truth=diagn[test.indexes], Pred=pred)
  return(ct)
}

