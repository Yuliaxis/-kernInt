# CLASSIFIER

#' SVM classifier
#'
#' @param data Input data
#' @param y Reponse variable (categoric)
#' @param classes Number of classes (2 or 3)
#' @param kernel "qJac" for quantitative Jaccard and "wqJacc" for quant Jaccard with weights
#' @param p Proportion of total data instances in the training set
#' @param C A vector with the possible costs to evaluate via k-Cross-Val. If no argument is provided cross-validation is not performed.
#' @param k The k for the k-Cross Validation. Minimum k = 2.
#' @param classimb "weights" to introduce class weights in the SVM algorithm and "data" to oversampling. If other arguments are provided nothing is done.
#' @param type Procedure to data oversampling ("ubOver" or "ubSMOTE")
#' @return Confusion matrix
#' @examples
#' classify(data=speMGX[,7:ncol(speMGX)],speMGX[,1],kernel="qJac",C=c(0.1,1),k=10)
#' classify(data=speMGX[,7:ncol(speMGX)],speMGX[,1],kernel="qJac",C=1,classimb="data", type="ubOver")
#' @importFrom kernlab SVindex as.kernelMatrix predict
#' @importFrom unbalanced ubBalance
#' @importFrom ROSE roc.curve
#' @export




classify <- function(data, y, classes=2, kernel, p=0.8, C, k=10,  classimb="no", type="ubOver") {

  # 1. Classes
  diagn <- as.numeric(y)
  if(classes == 2)   diagn[diagn == 3] <- 1 # De 3 a 2 classes: No Malalt /  malalt
  diagn <- as.factor(diagn)

  # 1. TR/TE
  N <- nrow(data)
  all.indexes <- 1:N

  learn.indexes <- trainIndx(n=N,ptrain=p)
  test.indexes <- all.indexes[-learn.indexes]

  nlearn <- length(learn.indexes)
  ntest <- N - nlearn

  if(classimb=="data")  {
    dades <- data[c(learn.indexes,test.indexes),]
    diagn <- diagn[c(learn.indexes,test.indexes)]
    print(summary(diagn))

    if(type == "ubOver")  SobrDadesTr <- ubBalance(dades[1:nlearn,], diagn[1:nlearn], type=type, positive=2,  k=0)
    if(type == "ubSMOTE")  SobrDadesTr <- ubBalance(dades[1:nlearn,], diagn[1:nlearn], type=type, positive=2)

    data <- rbind(SobrDadesTr$X,dades[(nlearn+1):N,])
    nlearn <- length(SobrDadesTr$Y)
    N <- nrow(data)
    diagn <- c(SobrDadesTr$Y, diagn[test.indexes])
    print(diagn)
    learn.indexes <- 1:nlearn
    test.indexes <- (nlearn+1):N
  }


  # 2. Compute kernel matrix
  if(kernel == "qJac") {
    Jmatrix <- qJacc(data)
  } else if(kernel == "wqJac") {
    Jmatrix <- wqJacc(data,y=diagn)
  } else {

  }
  trMatrix <- Jmatrix[learn.indexes,learn.indexes]

  # 4. Do R x k-Cross Validation
  if(hasArg(C)) {
    if(k<2) stop("k must be equal to or higher than 2")
    bh <- kCV(COST = C, K=trMatrix, Yresp=diagn[learn.indexes], k=k, R=k)
    cost <- bh$cost
  } else {
    cost <- 1
  }

  if(classimb == "weights") {
    model <- ksvm(trMatrix, diagn[learn.indexes], kernel="matrix", type="C-svc",
                  C=cost,class.weights=c("1"=summary(diagn[learn.indexes])[2],"2"=summary(diagn[learn.indexes])[1]))
  } else {
    model <- ksvm(trMatrix, diagn[learn.indexes], kernel="matrix", type="C-svc",C=cost )
  }

  # 5. Prediction
  teMatrix <- Jmatrix[test.indexes,learn.indexes]

  teMatrix <- teMatrix[,SVindex(model),drop=FALSE]
  teMatrix <- as.kernelMatrix(teMatrix)
  pred <- kernlab::predict(model,teMatrix)

  ### Confusion matrix
  ct <- table(Truth=diagn[test.indexes], Pred=pred)
  cat(paste("Accuracy:",round(Acc(ct),digits=4),"\n"))
  pr <- Prec(ct)
  cat(paste("Precision:",round(pr,digits=4),"\n"))
  rc <- Rec(ct)
  cat(paste("Recall:",round(rc,digits=4),"\n"))
  cat(paste("F1:",round(F1(Prec=pr,Rec=rc),digits=4),"\n"))

  print(roc.curve(diagn[test.indexes], pred))

  return(ct)
}
