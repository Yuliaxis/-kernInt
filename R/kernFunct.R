# KERNEL FUNCTIONS

#' Quantitative Jaccard
#'
#' This function delivers the quantitative Jaccard kernel matrix, also known as Ruzicka similarity.
#'
#' @param data A matrix or data.frame containing nonnegative values.
#' @return The quantitative Jaccard kernel matrix
#' @examples
#' example <- matrix(abs(rnorm(12)),nrow=4,ncol=3)
#' qJaccMatrix <- qJacc(data=example)
#' @export


qJacc <- function(data) {

  data <- as.matrix(data)
  n <- nrow(data)
  ids <- expand.grid.mod(1:n,rep = FALSE)

  Comp <- rowSums(pmin(data[ids[,1],],data[ids[,2],]))/rowSums(pmax(data[ids[,1],],data[ids[,2],]))

  Ncomb <- 1:((n^2-n)/2)
  K <- matrix(0,nrow=n,ncol=n)
  for (i in Ncomb)  K[ids[i,1],ids[i,2]] <- Comp[i]
  colnames(K) <- rownames(data)
  rownames(K) <- colnames(K)
  tK <- t(K)
  K <- tK + K # Upper triangular matrix to symmetric matrix
  diag(K) <- 1
  return(K)
}

#' Quantitative Jaccard with weights
#'
#' This function delivers the quantitative Jaccard kernel matrix, also known as Ruzicka similarity,
#' with the option to weight each variable.
#'
#' @param data A matrix or data.frame containing nonnegative values.
#' @param w  a vector of weights as long as the number of variables.
#' Class should be set to "numeric". If empty, weights will be authomathically computed as the RF mean
#' decrease in Gini index (see below)
#' @param y If a vector of weights is not provided, a vector of responses for computing the RF
#' weights should be provided instead.
#' @return The weighted quantitative Jaccard kernel matrix
#' @examples
#' example <- matrix(abs(rnorm(12)),nrow=4,ncol=3)
#' weights <- c(0.1,0.7,0.2)
#' qJaccMatrix <- wqJacc(data=example,w = weights)
#' @importFrom methods hasArg
#' @importFrom catkern rfweight
#' @export


wqJacc <- function(data, w, y) {
  if(!hasArg(w)) {
    if(hasArg(y)) {
      w <- rfweight(x=data,y=y,plot=FALSE)
      w <- as.vector(w)
    } else {
      stop("Either a vector of weights or the response variable should be provided")
    }
  }
  if(length(w) != ncol(data)) stop(paste("Number of weights and number of variables do not coincide."))

  weidata <- t(w * t(data))
  return(qJacc(data=weidata))
}
