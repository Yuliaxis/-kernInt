# KERNEL FUNCTIONS

#' Quantitative Jaccard
#'
#' This function delivers the quantitative Jaccard kernel matrix, also known as Ruzicka similarity.
#'
#' @param data A matrix or data.frame containing nonnegative values.
#' @param h An hyperparametr
#' @return The quantitative Jaccard kernel matrix
#' @examples
#' example <- matrix(abs(rnorm(12)),nrow=4,ncol=3)
#' qJaccMatrix <- qJacc(data=example)
#' @export


qJacc <- function(data,h) {

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
  if(hasArg(h))  K <- exp(h*K)/exp(h)
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
#' @param h An hyperparametr
#' @return The weighted quantitative Jaccard kernel matrix
#' @examples
#' example <- matrix(abs(rnorm(12)),nrow=4,ncol=3)
#' weights <- c(0.1,0.7,0.2)
#' qJaccMatrix <- wqJacc(data=example,w = weights)
#' @importFrom methods hasArg
#' @importFrom catkern rfweight
#' @export


wqJacc <- function(data, w, y, h) {
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

  return(qJacc(data=weidata,h=h))
}

#' Aitchison distance kernel
#'
#' Aitchison distance kernel
#'
#' @param data A matrix or data.frame containing only positive values.
#' @return Aitchison distance matrix
#' @examples
#' example <- matrix(abs(rnorm(12)),nrow=4,ncol=3)
#' aitch.dist(data=example)
#' @importFrom robCompositions cenLR
#' @importFrom stats dist
#' @export
#'

aitch.dist <- function(data) {
  minv <- min(  data[data!=min(data)] )
  minv <- minv/10
  clrEucl <- cenLR(data+minv)
  clrPROK <- clrEucl$x.clr
  aitch <- dist(clrPROK, method = "euclidean",diag=TRUE,upper = TRUE)
  return(as.matrix(aitch))
}


#' clr transf + RBF kernel
#'
#' @param data A matrix or data.frame containing only positive values.
#' @param h Kernel hyperparameter (gamma)
#' @return RBF kernel over a clr transformated data
#' @examples
#' example <- matrix(abs(rnorm(12)),nrow=4,ncol=3)
#' kernelMatrix <- clrRBF(data=example)
#' @importFrom robCompositions cenLR
#' @importFrom stats dist
#' @export
#'

clrRBF <- function(data,h) {
  K <- aitch.dist(data)
  if(hasArg(h)) K <- exp(-h*(K^2)) #RBF
  return(K)
}

# Covariance matrix
#' @keywords internal
covmat <- function(ind, visits, h) {
  totalrows <- visits*ind
  K <- matrix(0,nrow=totalrows,ncol=totalrows)
  rownames(K) <- 1:totalrows
  colnames(K) <- rownames(K)
  if(!hasArg(h)) h <- 0.5
  for(i in 1:ind) {
    index <- ((i-1)*visits+1):(i*visits)
    K[index,index] <- h
  }
  diag(K) <- 1
  return(K)
}


#Kernel selection
#' @keywords internal
kernelSelect <- function(kernel,data,y,h) { #h és un hiperparàmetre
  if(kernel == "qJac") {
    cat("quantJaccard kernel\n")
    return(qJacc(data,h))
  } else if(kernel == "wqJac") {
    cat("quantJaccard kernel + weights \n")
    return(wqJacc(data,y=y,h))
  }  else if(kernel == "cRBF") {
    cat("clr + RBF \n")
    return(aitch.dist(data))
  } else if(kernel == "time") {
    cat("Time matrix \n")
    return(TimeK(data,h))
  }else if(kernel == "cov") {
    cat("Covariance matrix \n")
    return(covmat(ind=data[1],visits=data[2],h))
  } else if(kernel == "matrix" | kernel == "time2") {
    cat("Pre-computed kernel matrix given \n")
    return(as.matrix(data))
  }   else {
    cat("standard RBF \n")
    # Jmatrix <-  kernelMatrix(rbfdot(sigma = G),data)
    stop("Kernel not available.") ##temporalment

  }
}


# Kernel computing for different types of data
#' @keywords internal
seqEval <- function(DATA,kernels,y,h) {
  m <- length(DATA)
  n <- nrow(DATA[[1]])
  K <- array(0,dim=c(n,n,m))
  for(i in 1:m) {
    K[,,i] <- kernelSelect(data = DATA[[i]],  kernel = kernels[i],h)
    print( K[,,i])
  }
  return(K)
}

