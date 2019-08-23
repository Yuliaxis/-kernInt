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

## Harmonic mean, to compute the F1 accuracy:
#' @keywords internal
harm <- function (a,b) { 2/(1/a+1/b) }


## K-fold cross- validation
#' @keywords internal
#' @importFrom kernlab ksvm cross

kCV <- function(COST, K, Yresp, k, R) {

  # on Y és el vector resposta, i K.train la submatriu amb els individus de training
  min.error <- Inf

  for (c in COST){
    Kmatrix <- K
    Y <- Yresp
    outer.error <- vector(mode="numeric",length=R)
      for (o in 1:R) {
        unordered <- sample.int(nrow(Kmatrix))
        Kmatrix <- Kmatrix[unordered,unordered]
        Y <- Y[unordered]
        K.model <- ksvm(Kmatrix, Y, type="C-svc",kernel="matrix",C=c,cross=k) # Rular el mètode
        outer.error[o] <- cross(K.model) # La mitjana dels errors és l'error de CV
      }
      v.error <- mean(outer.error)
      print(v.error)
      if (min.error > v.error) {   # < o <= ???
         min.error <- v.error
         best.cost <- c
       }
    }
  best.hyp <- data.frame(cost=best.cost,error= min.error)
  return(best.hyp)
}
