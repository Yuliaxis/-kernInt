# HELPER FUNCTIONS

# Pairwise kernel evaluation indexes

#' @keywords internal
expand.grid.mod <- function(x, rep) { # x is a vector

  g <- function(i) {
    z <- setdiff(x, x[seq_len(i-rep)])
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  do.call(rbind, lapply(seq_along(x), g))
}

#Training indexes
#' @keywords internal
trainIndx <- function(n, ptrain = 0.8) {
  unord <- sample.int(n)
  return(sort(sample(unord,round(ptrain*n))))
}

#Kernel selection
#' @keywords internal
kernelSelect <- function(kernel,data,y) {
  if(kernel == "qJac") {
    cat("quantJaccard kernel\n")
    return(qJacc(data))
  } else if(kernel == "wqJac") {
    cat("quantJaccard kernel + weights \n")
    return(wqJacc(data,y=y))
  }  else if(kernel == "cRBF") {
    cat("clr + RBF \n")
    return(clrRBF(data)) ## s'ha d'arreglar això, perquè ara mateix la gamma no es pot tocar.
  } else if(kernel == "time") {
    cat("Time matrix \n")
    return(TimeK(data))
  } else if(kernel == "matrix") {
    cat("Pre-computed kernel matrix given \n")
    return(data)
  }   else {
    cat("standard RBF \n")

    # Jmatrix <-  kernelMatrix(rbfdot(sigma = G),data)
    return(qJacc(data)) ##temporalment

  }
}



## K-fold cross- validation
#' @keywords internal
#' @importFrom kernlab ksvm cross

kCV <- function(GAMMA, CUT, COST, K, Yresp, k, R, prob, classimb=FALSE) {

  # on Y és el vector resposta, i K.train la submatriu amb els individus de training
  min.error <- Inf
  cut <- 0.5
  for (g in GAMMA){
    if (g == 0) {
      Kmatrix <- K
    } else { Kmatrix <- exp(g*K)/exp(g) #Standardized Kernel Matrix. Otherwise exp(g*K)
    }

    for (c in COST){
        Kmatrix <- K
        Y <- Yresp
        outer.error <- vector(mode="numeric",length=R)
          for (o in 1:R) {
            unordered <- sample.int(nrow(Kmatrix))
            Kmatrix <- Kmatrix[unordered,unordered]
            Y <- Y[unordered]
            if(classimb)  {
              K.model <- ksvm(Kmatrix, Y, type="C-svc",kernel="matrix",C=c,cross=k,
                                         class.weights=c("1"=as.numeric(summary(Yresp)[2]),"2"=as.numeric(summary(Yresp)[2]))) # Rular el mètode
            } else if (prob & hasArg(CUT)) {
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
              K.model <- ksvm(Kmatrix, Y, type="C-svc",kernel="matrix",prob.model=prob,C=c,cross=k) # Rular el mètode

            }
            if(!prob) outer.error[o] <- cross(K.model) # La mitjana dels errors és l'error de CV
          }
          v.error <- mean(outer.error)
          print(v.error)
          if (min.error > v.error) {   # < o <= ???
             min.error <- v.error
             best.cost <- c
             best.g <- g
             best.cut <- cut
          }
    }
  }
  print(best.cut)
  best.hyp <- data.frame(cost=best.cost,gamma=best.g,cut=best.cut,error= min.error)
  return(best.hyp)
}


## K-fold cross- validation (one-class SVM)
#' @keywords internal
#' @importFrom kernlab ksvm cross

kCV.one <- function(GAMMA, K, Yresp, NU, k=k, R=k) {
  min.error <- Inf
  for (g in GAMMA){
    if (g == 0) {
      Kmatrix <- K
    } else {
      Kmatrix <- exp(g*K)/exp(g) #Standardized Kernel Matrix. Otherwise exp(g*K)
    }

    for (nu in NU) {
      Kmatrix <- K
      Y <- Yresp
      outer.error <- vector(mode="numeric",length=R)
      for (o in 1:R) {
        unordered <- sample.int(nrow(Kmatrix))
        Kmatrix <- Kmatrix[unordered,unordered]
        Y <- Y[unordered]
        K.model <- ksvm(Kmatrix, Y, type="one-svc", kernel="matrix",nu=nu,cross=k) # Rular el mètode
        outer.error[o] <- cross(K.model) # La mitjana dels errors és l'error de CV
      }
      v.error <- mean(outer.error)
      print(v.error)
      if (min.error > v.error) {   # < o <= ???
        min.error <- v.error
        best.h1 <- nu
        best.h2 <- g
      }
    }
  }
best.hyp <- data.frame(nu=best.h1,g=best.h2, error= min.error)
return(best.hyp)
}

## K-fold cross- validation (regression)
#' @keywords internal
#' @importFrom kernlab ksvm cross

kCV.reg <-function(EPS, COST, GAMMA, K, Yresp, k, R) {

  # on Y és el vector resposta, i K.train la submatriu amb els individus de training
  min.error <- Inf
  for (g in GAMMA){
    if (g == 0) {
      Kmatrix <- K
    } else { Kmatrix <- exp(g*K)/exp(g) #Standardized Kernel Matrix. Otherwise exp(g*K)
    }
    Y <- Yresp
    for (e in EPS) {
      for (c in COST){
        outer.error <- vector(mode="numeric",length=R)
        for (o in 1:R) {
          unordered <- sample.int(nrow(Kmatrix))
          Kmatrix <- Kmatrix[unordered,unordered]
          Y <- Y[unordered]
          K.model <- ksvm(Kmatrix,Y,type="eps-svr",kernel="matrix",C=c,epsilon=e, cross=k) # Rular el mètode
          outer.error[o] <- cross(K.model) # La mitjana dels errors és l'error de CV
        }
        v.error <- mean(outer.error)
        # print(v.error)
        if (min.error > v.error) {   # < o <= ???
          min.error <- v.error
          best.cost <- c
          best.e <- e
          best.g <- g
        }
      }
    }
  }
  best.hyp <- data.frame(cost=best.cost,epsilon=best.e,gamma=best.g,error= min.error)
  return(best.hyp)
}
##  NMSE (regression)
#' @keywords internal
error.norm <- function(target,prediction) {
  N <- length(target)
  error <- sum((target-prediction)^2)/((N-1)*var(target)) ##(norm.mse <- model$deviance/((N-1)*var(target)))
  return(error)
}


## Harmonic mean, to compute the F1 accuracy:
#' @keywords internal
harm <- function (a,b) { 2/(1/a+1/b) }
