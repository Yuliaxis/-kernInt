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
#' @param data Input data: a matrix or data.frame with predictor variables/features as columns.
#' To perform MKL: a list of *m* datasets. All datasets should have the same number of rows.
#' @param y Reponse variable (factor)
#' @param kernel "lin" or rbf" to standard Linear and RBF kernels. "clin" for compositional linear and "crbf" for Aitchison-RBF
#' kernels. "jac" for quantitative Jaccard / Ruzicka kernel. "jsk" for Jensen-Shannon Kernel. "flin" and "frbf" for functional linear
#' and functional RBF kernels. "matrix" if a pre-computed kernel matrix is given as input.
#' To perform MKL: Vector of *m* kernels to apply to each dataset.
#' @param coeff ONLY IN MKL CASE: A *t·m* matrix of the coefficients, where *m* are the number of different data types and *t* the number of
#' different coefficient combinations to evaluate via k-CV. If absent, the same weight is given to all data sources.
#' @param prob if TRUE class probabilities (soft-classifier) are computed instead of a True-or-false assignation (hard-classifier)
#' @param classimb "weights" to introduce class weights in the SVM algorithm, "ubOver" to oversampling and "ubUnder" to undersample.
#' @param p The proportion of data reserved for the test set. Otherwise, a vector containing the indexes or the names of the rows for testing.
#' @param k The k for the k-Cross Validation. Minimum k = 2. If no argument is provided cross-validation is not performed.
#' @param C The cost. A vector with the possible costs (SVM hyperparameter) to evaluate via k-Cross-Val can be entered too.
#' @param H Gamma hyperparameter (only in RBF-like functions). A vector with the possible values to chose the best one via k-Cross-Val can be entered.
#' For the MKL, a list with *m* entries can be entered, being' *m* is the number of different data types. Each element on the list
#' must be a number or, if k-Cross-Validation is needed, a vector with the hyperparameters to evaluate for each data type.
#' @param domain Only used in "frbf" or "flin".
#' @return Confusion matrix or, if prob = TRUE and not cutoff is set, a data.frame with the class probability and the actual class.
#' @examples
#' classify(data=soil$abund ,y=soil$metadata[ ,"env_feature"],kernel="clin")
#' @importFrom kernlab alpha alphaindex as.kernelMatrix kernelMatrix predict rbfdot SVindex
#' @importFrom unbalanced ubBalance
#' @importFrom methods hasArg
#' @export



classify <- function(data, y,  coeff="mean", kernel,  prob=FALSE, classimb, p=0.2, k, C=1, H=NULL, domain=NULL) {

  ### Checking data
  check <- checkinput(data,y,kernel)
  m <- check$m
  data <- check$data
  kernel <- check$kernel

  # y class
  diagn <- as.factor(y)

  # 1. TR/TE
  if((length(p) == 1) && (p < 1)) { ### p és sa proporció de test.
    if(p<=0) stop("A test partition is mandatory")
    index <- finalTRTE(data,1-p)
    learn.indexes <- index$li
    test.indexes <- index$ti
  } else {                #### els índexs de test són entrats de forma manual
    if(class(p)=="character") {
      test.indexes <- which(rownames(data) %in% p)
    } else {
      test.indexes <- p
    }
    learn.indexes <- (1:nrow(data))[-test.indexes]
  }


  # 2. Compute kernel matrix

  try <- diagn[learn.indexes]
  tey <- diagn[test.indexes]

  if(m>1) {
    Jmatrix<- seqEval(DATA=data,domain=domain, kernels=kernel,h=NULL) ## Sense especificar hiperparàmetre.
    trMatrix <- Jmatrix[learn.indexes,learn.indexes,]
    teMatrix <- Jmatrix[test.indexes,learn.indexes,]
  } else {
    Jmatrix <- kernelSelect(kernel=kernel,domain=domain,data=data,h=NULL)
    trMatrix <- Jmatrix[learn.indexes,learn.indexes]
    teMatrix <- Jmatrix[test.indexes,learn.indexes]
  }

  # 3. Data imbalance

  wei <- NULL
  if(hasArg(classimb)) {
    if(classimb == "weights") {
      wei <- as.numeric(summary(try))
      names(wei) <- levels(try)
    } else if(classimb=="ubOver" || classimb=="ubUnder")  {
      s <- dataSampl(data=trMatrix,tedata=teMatrix, diagn=try,kernel=kernel,type=classimb)
      trMatrix <- s$data
      teMatrix <- s$tedata
      try <- s$diagn
    }  else {
      stop("Class balancing method misspelled or incorrect")
    }
  }


  # 4.  R x k-Cross Validation

  if(hasArg(k)) {
    if(k<2) stop("k should be equal to or higher than 2")
    if(m>1)  {
      if(class(coeff) == "character") {
        if(coeff == "mean") {
          coeff <- rep(1/m,m)
        } else {
          d <- aperm(trMatrix,c(3,1,2))
          x <- lapply(seq_len(nrow(d)), function(i) d[i,,]) ## transformar en llista
          coeff <- umkl(X=x,method=coeff)
        }
      }
      bh <- kCV.MKL(ARRAY=trMatrix, COEFF=coeff, KERNH=H, kernels=kernel, method="svc", COST = C,
                     Y=try, k=k,  prob=prob, R=k,classimb=wei)
      coeff <- bh$coeff ##indexs
      conserv <- c("coeff","cost","error")
      if(!is.null(H)) conserv <- c(conserv,"h")
      bh <- bh[conserv]
    } else {
    bh <- kCV.core(method="svc",COST = C, H = H, kernel=kernel, K=trMatrix, prob=prob,
                   Y=try, k=k, R=k,classimb=wei)
    bh <- bh[-which(is.na(bh))]
    }
    cost <- bh$cost
    H <- bh$h

  } else {
    if(!is.null(H))   {
      H <- kernHelp(H)$hyp
      bh <- data.frame(H)
    } else {
      bh <- NULL
    }
    if(length(C)>1) warning("Multiple C and no k provided - Only the first element will be used")
    cost <- C[1]
    if(m>1)   {
      bh <- list(coeff=coeff,h=H,cost=cost)
    }  else {
      bh <- cbind(bh,cost)
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

  # 5. Model

  model <- ksvm(trMatrix, try, kernel="matrix", type="C-svc", prob.model = prob, C=cost, class.weights=wei)

  ##### Importances (only linear)
  importances <- NULL

    if(identical(unique(kernel), "lin")) {
      if(nlevels(y) == 2) {
        alphaids <- alphaindex(model) # Indices of SVs in original data
        alphaids <-  learn.indexes[unlist(alphaids)]
        alphas <- alpha(model)
        alphas <- unlist(alphas)
        ys <-  as.numeric(y[alphaids])
        ys[ys==2] <- -1

        if(m>1) {
          coeff <- array(rep(coeff,each=length(data)/m),dim=dim(data))
          cosn <-  apply(data^2,3L,rowSums) ## cosine normalization
          cosn <- array(rep(cosn,each=dim(data)[2]),dim=c(dim(data)[2],dim(data)[1],dim(data)[3]))
          cosn <- aperm(cosn,c(2,1,3))
          svmatrix <- coeff * data * 1/(sqrt(cosn))
          svmatrix <- t(apply(svmatrix, 1L, c))
          svmatrix <- as.matrix(svmatrix[alphaids, ])

        } else {
          svmatrix <- as.matrix(data[alphaids, ])
          svmatrix <- svmatrix /sqrt(rowSums(svmatrix^2))  ### cosine normalization
        }
        importances  <- (colSums( matrix((ys * alphas),ncol=ncol(svmatrix),nrow=length(ys)) * svmatrix))^2
      }  else {
      print("not implemented")
      }
    }


  # 6. Prediction

  teMatrix <- teMatrix[,SVindex(model),drop=FALSE]
  teMatrix <- as.kernelMatrix(teMatrix)

  if(prob)  {
    pred <- predict(model,teMatrix,type = "probabilities")
    ct <- NULL
    test <- cbind(tey,pred)
    colnames(test)[1] <- "true"
  }  else    {
    pred <- kernlab::predict(model,teMatrix)
    pred <- as.factor(pred)
    levels(pred) <- levels(diagn)
    ct <- table(Truth=tey, Pred=pred)  ### Confusion matrix
    test <- data.frame(true=tey,predicted=pred)
  }
  rownames(test) <- test.indexes
  return(list("conf.matrix"=ct,"hyperparam"=bh,"prediction"=test,"var.imp"=importances))
}

