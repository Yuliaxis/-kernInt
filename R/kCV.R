## K-fold cross- validation (classification)
#' @keywords internal
#' @importFrom kernlab ksvm cross

kCV.class <- function(CUT, COST, K, Yresp, k, R, prob, classimb=FALSE) {

  # on Y és el vector resposta, i K.train la submatriu amb els individus de training
  min.error <- Inf
  if(is.null(CUT)) cut <- 0.5

  for (c in COST){

    outer.error <- vector(mode="numeric",length=R)
    for (o in 1:R) {
      unordered <- sample.int(nrow(K))
      Kmatrix <- K[unordered,unordered]
      Y <- Yresp[unordered]
      if(prob & hasArg(CUT)) {
        N <- trunc(nrow(Kmatrix)/k,digits=0)
        PRED <- matrix(NA,ncol=1,nrow=nrow(Kmatrix))
        # rownames(PRED) <- as.character(Y)
        for(p in 0:(k-1)) {
          if(p < (k-1)) {
             indexTE <- (1+(p*N)):((p+1)*N)
          } else {
             indexTE <- (1+(p*N)):nrow(Kmatrix)
          }
          TEST <- Kmatrix[indexTE,-indexTE]
          K.model <- ksvm(Kmatrix[-indexTE,-indexTE], Y[-indexTE], type="C-svc",kernel="matrix",prob.model=prob,C=c) # Rular el mètode
          TEST <- TEST[,SVindex(K.model),drop=FALSE]
          TEST <- as.kernelMatrix(TEST)
          pred <- predict(K.model,TEST,type = "probabilities")
          PRED[indexTE,1] <- pred[,2] #minority class
        }
        PRED <- na.omit(PRED)
        acu <- vector(mode="numeric",length=length(cut))
        for(cut in 1:length(CUT)) {
          pr <- (PRED < CUT[cut])
          pr[pr] <- 1
          pr[pr==0] <- 2
          acu[cut] <- sum(as.numeric(Y) == pr)/nrow(pr) # És ok aquesta aproximació??
        }
        cut <- CUT[which.max(acu)]
        outer.error[o] <- max(acu)
      } else {
        K.model <- ksvm(Kmatrix, Y, type="C-svc",kernel="matrix",prob.model=prob,class.weights=classimb,C=c,cross=k) # Rular el mètode
        }
      if(!prob) outer.error[o] <- cross(K.model) # La mitjana dels errors és l'error de CV
    }
    v.error <- mean(outer.error)
    print(v.error)
    if (min.error > v.error) {
        min.error <- v.error
        best.cost <- c
        best.cut <- cut
      }
    }
  print(best.cut)
  best.hyp <- data.frame(cost=best.cost,epsilon=NA,cut=best.cut,nu=NA, error= min.error)
  return(best.hyp)
}



## K-fold cross- validation (one-class SVM)
#' @keywords internal
#' @importFrom kernlab ksvm cross

kCV.one <- function(K, Yresp, NU, k=k, R=k) {
  min.error <- Inf
  for (nu in NU) {
    outer.error <- vector(mode="numeric",length=R)
    for (o in 1:R) {
      unordered <- sample.int(nrow(K))
      Kmatrix <- K[unordered,unordered]
      Y <- Yresp[unordered]
      K.model <- ksvm(Kmatrix, Y, type="one-svc", kernel="matrix",nu=nu,cross=k) # Rular el mètode
      outer.error[o] <- cross(K.model) # La mitjana dels errors és l'error de CV
    }
    v.error <- mean(outer.error)
    print(v.error)
    if (min.error > v.error) {   # < o <= ???
      min.error <- v.error
      best.h1 <- nu
    }
  }
  best.hyp <- data.frame(cost=NA,epsilon=NA,cut=NA,nu=best.h1, error= min.error)
  return(best.hyp)
}


## K-fold cross-validation (regression SVM)
#' @keywords internal
#' @importFrom kernlab ksvm cross

kCV.reg <- function(EPS, COST, K, Yresp, k, R) {
  min.error <- Inf
  for (e in EPS) {
    for (c in COST){
      outer.error <- vector(mode="numeric",length=R)
      for (o in 1:R) {
        unordered <- sample.int(nrow(K))
        Kmatrix <- K[unordered,unordered]
        Y <- Yresp[unordered]
        K.model <- ksvm(Kmatrix,Y,type="eps-svr",kernel="matrix",C=c,epsilon=e, cross=k) # Rular el mètode
        outer.error[o] <- cross(K.model) # La mitjana dels errors és l'error de CV
      }
      v.error <- mean(outer.error)
      print(v.error)
      if (min.error > v.error) {
        min.error <- v.error
        best.cost <- c
        best.e <- e
      }
    }
  }
  best.hyp <- data.frame(cost=best.cost,epsilon=best.e,cut=NA,nu=NA,error= min.error)
  return(best.hyp)
}


## Gamma hyperparameter - general kCV procedure if MKL is not being performed
#' @keywords internal
kCV.core <- function(H, kernel, method, K, ...) {
  min.error <- Inf
  for (h in H) {
    Kmatrix <- hyperkSelection(K=K,h=h,kernel=kernel)
    # Y <- Yresp
    if(method == "svc") {
      bh <- kCV.class(K=Kmatrix, ...) ## calls svm classification
    } else if(method == "svr") { ## smv calls regression
      bh <- kCV.reg(K=Kmatrix, ...)
    } else { ## calls one-class svm
      bh <- kCV.one(K=Kmatrix, ...)
    }
    if (min.error > bh$error) {
      best.h <- h
      best.cost <- bh$cost
      best.e <- bh$epsilon
      best.cut <- bh$cut
      best.nu <- bh$nu
      min.error <- bh$error
    }
  }

  best.hyp <- data.frame(h=best.h,cost=best.cost,epsilon=best.e,cut=best.cut,nu=best.nu,error= min.error)
  return(best.hyp)
}

## kCV procedure if MKL is being performed
#' @keywords internal
kCV.MKL <- function(ARRAY, COEFF, KERNH, kernels, method, ...) {
  min.error <- Inf
  nhyp <- length(KERNH)
  ARRAY2 <- matrix(0,dim=c(dim(ARRAY)[1],dim(ARRAY)[2],nhyp))
  for(k in 1:nhyp) {
    j <- ceiling(k/nrow(KERNH))
    ARRAY2[,,k] <- hyperkSelection(ARRAY[,,j],h=KERNH[[k]],kernel=kernels[j])
  }
  indexes <- expand.grid(rep(1:nrow(KERNH),ncol(KERNH)))
  m <- nrow(COEFF)
  for(i in 1:m) {
    for(k in 1:nrow(indexes)) {
      Kmatrix <- KInt(ARRAY2[,,indexes[k,]],coeff=COEFF[i,])
      if(method == "svc") {
        bh <- kCV.class(Kmatrix,...) ## calls svm classification
      } else if(method == "svr") { ## smv calls regression
        bh <- kCV.reg(Kmatrix,...)
      } else { ## calls one-class svm
        bh <- kCV.one(Kmatrix,...)
      }
      if (min.error > bh$error) {
        ii <- indexes[k,]
        best.cost <- bh$cost
        best.e <- bh$epsilon
        best.cut <- bh$cut
        best.nu <- bh$nu
        best.coeff <- COEFF[i,]
        min.error <- bh$error
      }
    }
  }
  ii <- matrix(c(ii,1:length(ii)),ncol=2)
  best.h <- rep(0,m)
  for(j in 1:m) best.h[j] <- KERNH[ii[j],j]
  best.hyp <- list(coeff=best.coeff,h=best.h,cost=best.cost,epsilon=best.e,cut=best.cut,nu=best.nu,error= min.error)
  return(best.hyp)
}
