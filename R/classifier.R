# CLASSIFIER

#' SVM classifier
#'
#' @param data Input data
#' @param y Reponse variable (categoric)
#' @param kernel "cRBF" for clrRBF, "qJac" for quantitative Jaccard and  "wqJacc" for quantitative Jaccard with weights
#' @param p Proportion of total data instances in the training set
#' @param C A vector with the possible costs to evaluate via k-Cross-Val. If no argument is provided cross-validation is not performed.
#' @param G Gamma hyperparameter
#' @param CUT if prob == TRUE the possible cut-offs to consider one class or another.
#' @param k The k for the k-Cross Validation. Minimum k = 2.
#' @param prob if TRUE builds a model for calculating class probabilities
#' @param classimb "weights" to introduce class weights in the SVM algorithm and "data" to oversampling. If other arguments are provided nothing is done.
#' @param type Procedure to data oversampling ("ubOver","ubUnder" or "ubSMOTE")
#' @return Confusion matrix
#' @examples
#' diag <- as.numeric(speMGX[,1])
#' diag[diag == 3] <- 1  # De 3 a 2 classes: No Malalt /  malalt#'
#' classify(data=speMGX[,7:ncol(speMGX)],diag,kernel="qJac",C=c(0.1,1),k=10)
#' classify(data=speMGX[,7:ncol(speMGX)],diag,kernel="qJac",C=1,classimb="data", type="ubOver")
#' @importFrom kernlab as.kernelMatrix kernelMatrix predict rbfdot SVindex
#' @importFrom unbalanced ubBalance
#' @importFrom ROSE roc.curve
#' @export



classify <- function(data, y, kernel,  p=0.8, C=1, G=0, CUT=0.5, k, prob=FALSE, classimb="no", type="ubOver") {

  # 1. Classes
  diagn <- as.factor(y)

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
    if(type == "ubUnder")  SobrDadesTr <- ubBalance(dades[1:nlearn,], diagn[1:nlearn], type=type, positive=2)
    if(type == "ubSMOTE")  SobrDadesTr <- ubBalance(dades[1:nlearn,], diagn[1:nlearn], type=type, positive=2)
    data <- rbind(SobrDadesTr$X,dades[(nlearn+1):N,])
    nlearn <- length(SobrDadesTr$Y)
    N <- nrow(data)
    diagn <- c(SobrDadesTr$Y, diagn[test.indexes])
    diagn <- as.factor(diagn)
    print(diagn)
    learn.indexes <- 1:nlearn
    test.indexes <- (nlearn+1):N
   }


  # 2. Compute kernel matrix
  if(kernel == "qJac") {
    cat("quantJaccard kernel \n")
    Jmatrix <- qJacc(data)
  } else if(kernel == "wqJac") {
    Jmatrix <- wqJacc(data,y=diagn)
    cat("quantJaccard kernel + weights \n")

  }  else if(kernel == "cRBF") {
    Jmatrix <- clrRBF(data)
    cat("clr + RBF \n")

  }   else {
    cat("standard RBF \n")

    # Jmatrix <-  kernelMatrix(rbfdot(sigma = g),data)

  }

  trMatrix <- Jmatrix[learn.indexes,learn.indexes]
  teMatrix <- Jmatrix[test.indexes,learn.indexes]


  # 4. Do R x k-Cross Validation
  if(hasArg(k)) {
    if(k<2) stop("k must be equal to or higher than 2")
    bh <- kCV(COST = C, GAMMA = G, CUT=CUT, K=trMatrix, prob=prob, Yresp=diagn[learn.indexes], k=k, R=k)
    if(classimb == "weights") {
      # diagn <- as.numeric(diagn)   ############# culpa weights
      bh <- kCV(COST = C,GAMMA = G, K=trMatrix, Yresp=diagn[learn.indexes], k=k, R=k,classimb=TRUE)
    }
    cost <- bh$cost
    G <- bh$gamma
    cut <- bh$cut
  } else {
    if(length(C)>1) paste("C > 1 - Only the first element will be used")
    cost <- C[1]
    G <- G[1]
    CUT <- CUT[1]
  }

  if(G != 0)  {
    trMatrix <- exp(G * trMatrix)/exp(G)
    teMatrix <- exp(G * teMatrix)/exp(G)
  }


  if(classimb == "weights") {
    model <- ksvm(trMatrix, diagn[learn.indexes], kernel="matrix", type="C-svc",
                  C=cost,GAMMA = G, class.weights=c("1"=as.numeric(summary(diagn[learn.indexes])[2]),"2"=as.numeric(summary(diagn[learn.indexes])[1])))
  } else {
    model <- ksvm(trMatrix, diagn[learn.indexes], prob.model = prob, kernel="matrix", type="C-svc",C=cost,GAMMA = G )
  }


  # 5. Prediction

  teMatrix <- teMatrix[,SVindex(model),drop=FALSE]
  teMatrix <- as.kernelMatrix(teMatrix)

  levels(diagn) <- c("1","2")
  if(prob)  {
    pred <- predict(model,teMatrix,type = "probabilities")
    print(cut)
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
  # print(ct)
  # cat(paste("Accuracy:",round(Acc(ct),digits=4),"\n"))
  # pr <- Prec(ct)
  # cat(paste("Precision:",round(pr,digits=4),"\n"))
  # rc <- ct[2,2]/sum(ct[2,])
  # cat(paste("Recall:",round(rc,digits=4),"\n"))
  # f1 <- F1(Prec=pr,Rec=rc)
  # cat(paste("F1:",round(f1,digits=4),"\n"))
   # AUC <- roc.curve(diagn[test.indexes], pred)
  # print(AUC$auc)
}
