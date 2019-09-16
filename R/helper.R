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


## K-fold cross- validation
#' @keywords internal
#' @importFrom kernlab ksvm cross

kCV <- function(GAMMA,  COST, K, Yresp, k, R,prob, classimb=FALSE) {

  # on Y és el vector resposta, i K.train la submatriu amb els individus de training
  min.error <- Inf
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
                                         class.weights=c("1"=19,"2"=65)) # Rular el mètode
            } else {
              K.model <- ksvm(Kmatrix, Y, type="C-svc",kernel="matrix",prob.model=prob,C=c,cross=k) # Rular el mètode

            }
            outer.error[o] <- cross(K.model) # La mitjana dels errors és l'error de CV
          }
          v.error <- mean(outer.error)
          print(v.error)
          if (min.error > v.error) {   # < o <= ???
             min.error <- v.error
             best.cost <- c
          }
    }
  }
  best.hyp <- data.frame(cost=best.cost,error= min.error)
  return(best.hyp)
}

## K-fold cross- validation (regression)
#' @keywords internal
#' @importFrom kernlab ksvm cross

kCV.reg <- function(GAMMA, EPS, COST, K, Yresp, k, R) {

  # on Y és el vector resposta, i K.train la submatriu amb els individus de training
  min.error <- Inf
  for (g in GAMMA){
    if (g == 0) {
      Kmatrix <- K
    } else { Kmatrix <- exp(g*K)/exp(g) #Standardized Kernel Matrix. Otherwise exp(g*K)
    }
  for (c in COST){
    for (e in EPS) {
      Kmatrix <- K
      Y <- Yresp
      outer.error <- vector(mode="numeric",length=R)
      for (o in 1:R) {
        unordered <- sample.int(nrow(Kmatrix))
        Kmatrix <- Kmatrix[unordered,unordered]
        Y <- Y[unordered]
        K.model <- ksvm(Kmatrix, Y, type="eps-svr",kernel="matrix",C=c,cross=k) # Rular el mètode
        outer.error[o] <- cross(K.model) # La mitjana dels errors és l'error de CV
      }
      v.error <- mean(outer.error)
      print(v.error)
      if (min.error > v.error) {   # < o <= ???
        min.error <- v.error
        best.cost <- c
        }
      }
    }
  }
  best.hyp <- data.frame(cost=best.cost,error= min.error)
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

## Accuracy
#' @keywords internal
Acc <- function(ct) sum(diag(ct))/sum(ct)

## Precision
#' @keywords internal
Prec <- function(ct) {
  pr <- ct[2,2]/sum(ct[,2])
  if(is.nan(pr)) pr <- 0
  return(pr)
}

## Recall
#' @keywords internal
Rec <-  function(ct) {
  rc <- ct[2,2]/sum(ct[2,])
  if(is.nan(rc)) rc <- 0
  return(rc)
}

## F1
#' @keywords internal
F1 <-  function(Prec,Rec) { (2*Prec*Rec)/(Prec+Rec) }
######

