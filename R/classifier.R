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
#' @param data Input data
#' @param y Reponse variable (binary)
#' @param kernel "cRBF" for clrRBF, "qJac" for quantitative Jaccard,  "wqJacc" for quantitative Jaccard with weights.
#' "matrix" if a pre-calculated kernel matrix is given as input.
#' @param prob if TRUE class probabilities (soft-classifier) are computed instead of a True-or-false assignation (hard-classifier)
#' @param classimb "weights" to introduce class weights in the SVM algorithm and "data" to oversampling. If other arguments are provided nothing is done.
#' @param type If classimb = "data", the procedure to data oversampling or undersampling ("ubOver","ubUnder" or "ubSMOTE")
#' @param p Proportion of total data instances in the training set
#' @param k The k for the k-Cross Validation. Minimum k = 2. If no argument is provided cross-validation is not performed.
#' @param C The cost. A vector with the possible costs (SVM hyperparameter) to evaluate via k-Cross-Val can be entered too.
#' @param G Gamma hyperparameter. A vector with the possible gammas to evaluate via k-Cross-Val can be entered too.
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



classify <- function(data, y, kernel,  prob=FALSE, classimb="no", type="ubOver", p=0.8, k, C=1, G=0, CUT) {

  # 1. Classes
  diagn <- as.factor(y)
  levels(diagn) <- c("1","2")

  # 1. TR/TE

  ids <- as.factor(rownames(data))
  N <-  nlevels(ids)
  # N <- nrow(data)
  all.indexes <- 1:N

  learn.indexes <- trainIndx(n=N,ptrain=p)
  test.indexes <- all.indexes[-learn.indexes]

  ##Mostres vinculades
  if(length(ids) > nlevels(ids)) {
    trNames <- levels(ids)[learn.indexes]
    teNames <-  levels(ids)[test.indexes]
    learn.indexes <- which(ids %in% trNames)
    test.indexes <- which(ids %in% teNames)
    #Unique test - si es vol llevar, comentar aquestes 4 línies.
    names(test.indexes) <- ids[which(ids %in% teNames)]
    test.indexes <- sample(test.indexes)[teNames]
    names(test.indexes) <- NULL
    test.indexes <- sort(test.indexes)
  }

  if(classimb=="data")  {

    nlearn <- length(learn.indexes)
    ntest <- length(test.indexes)
    N <- nlearn + ntest

    diagn <- diagn[c(learn.indexes,test.indexes)]
    if(kernel == "matrix") {
      if(type == "ubSMOTE") stop("Kernel matrix as input is not compatible with SMOTE. Original dataset is required.")

      dades <- data[c(learn.indexes,test.indexes),c(learn.indexes,test.indexes)]
      rownames(dades) <- 1:N

      if(type == "ubOver")  SobrDadesTr <- ubBalance(dades[1:nlearn,], diagn[1:nlearn], type=type, positive=2, k=0)
      if(type == "ubUnder")  SobrDadesTr <- ubBalance(dades[1:nlearn,], diagn[1:nlearn], type=type, positive=2)

      ii <- c(as.numeric(rownames(SobrDadesTr$X)),(nlearn+1):N)
      data <- data[ii,ii]
      diagn <- diagn[ii]
      N <- nrow(data)
      nlearn <- length(SobrDadesTr$Y)
      learn.indexes <- 1:nlearn
      test.indexes <- (nlearn+1):N
    } else {
    dades <- data[c(learn.indexes,test.indexes),]

    if(type == "ubUnder") SobrDadesTr <- ubBalance(dades[1:nlearn,], diagn[1:nlearn], type=type, positive=2)
    if(type == "ubOver") SobrDadesTr <- ubBalance(dades[1:nlearn,], diagn[1:nlearn], type=type, positive=2,  k=0)
    if(type == "ubSMOTE") SobrDadesTr <- ubBalance(dades[1:nlearn,], diagn[1:nlearn], type=type, positive=2)
    data <- rbind(SobrDadesTr$X,dades[(nlearn+1):N,])
    nlearn <- length(SobrDadesTr$Y)
    N <- nrow(data)
    diagn <- c(SobrDadesTr$Y, diagn[test.indexes])
    diagn <- as.factor(diagn)
    learn.indexes <- 1:nlearn
    test.indexes <- (nlearn+1):N
    }

  }

  if(classimb == "weights")  {
    wei <- TRUE
  } else {
    wei <- FALSE
  }

  # 2. Compute kernel matrix
  Jmatrix <- kernelSelect(kernel,data,y)

  trMatrix <- Jmatrix[learn.indexes,learn.indexes]
  teMatrix <- Jmatrix[test.indexes,learn.indexes]


  # 4. Do R x k-Cross Validation
  if(hasArg(k)) {
    if(k<2) stop("k must be equal to or higher than 2")
    if(hasArg(CUT)) {
       bh <- kCV(COST = C, GAMMA = G, CUT=CUT, K=trMatrix, prob=prob, Yresp=diagn[learn.indexes], k=k, R=k,classimb=wei)
       cut <- bh$cut
       } else {
       bh <- kCV(COST = C, GAMMA = G, K=trMatrix, prob=prob, Yresp=diagn[learn.indexes], k=k, R=k,classimb=wei)
      }
    cost <- bh$cost
    G <- bh$gamma

  } else {
    if(length(C)>1) paste("C > 1 and no k provided - Only the first element will be used")
    cost <- C[1]
    if(length(G)>1) paste("G > 1 and no k provided- Only the first element will be used")
    G <- G[1]
    if(hasArg(CUT) && length(CUT)>1) {
      paste("CUT > 1 and no k provided - Only the first element will be used")
      CUT <- CUT[1]
    }
  }

  if(G != 0)  {
    trMatrix <- exp(G * trMatrix)/exp(G)
    teMatrix <- exp(G * teMatrix)/exp(G)
  }

  if(wei) {
    model <- ksvm(trMatrix, diagn[learn.indexes], kernel="matrix", type="C-svc",
                  C=cost,GAMMA = G, class.weights=c("1"=as.numeric(summary(diagn[learn.indexes])[2]),"2"=as.numeric(summary(diagn[learn.indexes])[1])))
  } else {
    model <- ksvm(trMatrix, diagn[learn.indexes], prob.model = prob, kernel="matrix", type="C-svc",C=cost,GAMMA = G )
  }

  # 5. Prediction

  teMatrix <- teMatrix[,SVindex(model),drop=FALSE]
  teMatrix <- as.kernelMatrix(teMatrix)

  if(prob)  {
    pred <- predict(model,teMatrix,type = "probabilities")
    if(!hasArg(CUT)) {
      return(data.frame(Actual=diagn[test.indexes],Predicted = as.factor(pred)))
    }
    print(paste("Best cut is", cut))
    pred <- (pred[,2] < cut)
    pred[pred] <- 1
    pred[pred==0] <- 2
    }
  else    { pred <- kernlab::predict(model,teMatrix) }
   pred <- as.factor(pred)
   levels(pred) <- c("1","2")
   print(pred)

  ### Confusion matrix
  ct <- table(Truth=diagn[test.indexes], Pred=pred)
  return(ct)
}
