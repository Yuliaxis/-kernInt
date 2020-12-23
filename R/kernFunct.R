# KERNEL FUNCTIONS


############# STANDARD FUNCTIONS ###################

## Cosine normalization
#' @keywords internal

cosNorm <- function(K) {
  D <- diag(1/sqrt(diag(K)))
  K <- D %*% K %*% D
  return(K)
}


#' Linear kernel
#'
#' @param data A matrix or data.frame with real numbers
#' @param cos.norm If TRUE the cosine normalization is applied
#' @return The linear kernel matrix
#' @examples
#' example <- matrix(rnorm(12),nrow=4,ncol=3)
#' kmatrix <- Linear(data=example)
#' @export


Linear <- function(data,cos.norm=TRUE) {
  data <- as.matrix(data)
  K <- tcrossprod(data)
  if(cos.norm) K <- cosNorm(K)
  rownames(K) <- rownames(data)
  colnames(K) <- rownames(data)
  return(K)
}


#' RBF kernel
#'
#' @param data A matrix or data.frame with real numbers
#' @param h Gamma hyerparameter. If NULL, the euclidian distances are returned
#' @return The RBF kernel matrix
#' @examples
#' example <- matrix(rnorm(12),nrow=4,ncol=3)
#' kmatrix <- RBF(data=example)
#' @export

RBF <- function(data,h=NULL) {
  data <- as.matrix(data)
  N <- nrow(data)
  kk <- tcrossprod(data)
  dd <- diag(kk)
  K <- 2*kk-matrix(dd,N,N)-t(matrix(dd,N,N))
  if(is.null(h) || h == 0 ) return(-K)
  return(exp(h*K))
}

####### Microbiome  beta diversities #########


# Bray Curtis
#
# This function returns the Bray Curtis or Sorensen-Dice kernel
#
# @param data A matrix or data.frame containing nonnegative values.
# @return The Bray Curtis kernel matrix
# @examples
# example <- matrix(abs(rnorm(12)),nrow=4,ncol=3)
# kmatrix <- BCK(data=example)
# @importFrom vegan vegdist
# @export
#
# BCK <- function(data) {
#   bc <- vegdist(data,method="bray",diag=TRUE,upper=TRUE)
#   bc <- as.matrix(bc)
#   return(1-bc)
# }

#' Quantitative Jaccard
#'
#' This function returns the quantitative Jaccard or min-max kernel, also known as Ruzicka similarity.
#'
#' @param data A matrix or data.frame containing nondnegative values.
#' @return The quantitative Jaccard kernel matrix
#' @examples
#' example <- matrix(abs(rnorm(12)),nrow=4,ncol=3)
#' kmatrix <- qJacc(data=example)
#' @importFrom vegan vegdist
#' @export


qJacc <- function(data) {
  ruzicka <- vegdist(data,method="jaccard",diag=TRUE,upper=TRUE)
  ruzicka <- as.matrix(ruzicka)
  return(1-ruzicka)
}


# Jensen-Shannon kernel
#'
#' This function returns the Jensen Shannon Kernel
#'
#' @param data A matrix or data.frame containing positive values.
#' @return The JS kernel matrix
#' @examples
#' example <- matrix(abs(rnorm(12)),nrow=4,ncol=3)
#' kmatrix <- JSK(data=example)
#' @importFrom philentropy JSD
#' @export

JSK <- function(data)    {
  if(sum(rowSums(data)>1)>0) {
    div <- JSD(data,est.prob = "empirical")
  } else {
    div <- JSD(data)
  }
  return(1-div)
}

#' Compositional Linear kernel
#'
#' @param data A matrix or data.frame containing compositiona data in raw counts
#' @param cos.norm If TRUE the cosine normalization is applied
#' @return The compositional linear kernel matrix
#' @examples
#' kmatrix <- clrLin(data=soil$abund)
#' @export

clrLin <- function(data,cos.norm=TRUE) {
  data <- clr(data)
  return(Linear(data=data,cos.norm = cos.norm))
}

#' Aitchison RBF kernel
#'
#' @param data  A matrix or data.frame containing compositiona data in raw counts
#' @param h Gamma hyerparameter. If NULL, the Aitchison distances are returned

#' @return The Aitchison RBF kernel.
#' @examples
#' kmatrix <- clrRBF(data=soil$abund,h=0.001)
#' @export
#'

clrRBF <- function(data,h=NULL) {
  data <- clr(data)
  return(RBF(data=data,h=h) )#RBF

}

######### Kernels for functions ################


#' Functional linear kernel
#'
#' @param data  Output of lsq()
#' @param domain Domain
#' @param cos.norm If TRUE, perform cosine normalization
#' @return A kernel matrix
#' @examples
#' growth2 <- growth
#' colnames(growth2) <-  c( "time", "height")
#' fitted <- lsq(data=growth2,degree=2)
#' Kfun(fitted,domain=c(1,18))
#' @importFrom polynom integral polynomial
#' @importFrom methods is
#' @export

Kfun <- function(data,domain,cos.norm=FALSE) {

  if(!is(data, "lsq")) stop("data should be of class lsq")

  degree <- data$degree

  inf <- min(domain)
  sup <- max(domain)

  d <- length(data$y.names)

  data <- data$coef

  ids <- expand.grid.mod(1:nrow(data),rep = FALSE)
  jds <- seq(from=1,to=ncol(data),by=degree+1)

  if(degree != 1) { ## a l'espera d'una implementació millor:
    Comp <- matrix(0,nrow=nrow(ids),ncol=d)
    for(i in 1:nrow(ids)) {
      for(j in 1:length(jds)) {
        Comp[i,j] <-  integral((polynomial(data[ids[i,1], j:(j+degree)]) *  polynomial(data[ids[i,2],j:(j+degree)])),limits=domain)

      }
    }
  } else {

  ### helper
  evenID <- function(number) seq(from=2,to=number,by=2)
  oddID <- function(number) seq(from=1,to=number,by=2)


   newF <- data[ids[,1],] * data[ids[,2],]

  A <- (newF[,evenID(ncol(data))])/3
  B <- (newF[,oddID(ncol(data))])

  newF <- data[ids[,1],-1] * data[ids[,2],-ncol(data)] + data[ids[,1],-ncol(data)] * data[ids[,2],-1]

  C <-  (newF[,oddID(ncol(data)-1)])/2

  X <- (sup^3-inf^3)*A + (sup^2-inf^2)*C + (sup-inf)*B
  }

  if(d>1) Comp <- rowSums(X)

  K <- matrix(0,ncol=nrow(data),nrow=nrow(data))
  for(i in 1:nrow(ids)) K[ids[i,1],ids[i,2]] <- K[ids[i,2],ids[i,1]] <- Comp[i]
  if(cos.norm) K <- cosNorm(K)
  rownames(K) <- colnames(K) <- rownames(data)
  return(K)
}


#' RBF kernel for functions
#'
#' @param data  Output of lsq()
#' @param domain Domain
#' @param h Gamma hyerparameter. If NULL, the L2 norms are returned
#' @return The L2 norm (euclidean distance) between to matrices
#' @examples
#' growth2 <- growth
#' colnames(growth2) <-  c( "time", "height")
#' fitted <- lsq(data=growth2,degree=2)
#' RBFun(fitted,domain=c(1,18),h=0.00001)
#' @importFrom polynom integral polynomial
#' @importFrom methods is

#' @export

RBFun <- function(data,h=NULL,domain) {

  if(!is(data,"lsq")) stop("data should be of class lsq")

  degree <- data$degree

  inf <- min(domain)
  sup <- max(domain)

  d <- length(data$y.names)

  data <- data$coef

  ids <- expand.grid.mod(1:nrow(data),rep = FALSE)
  jds <- seq(from=1,to=ncol(data),by=degree+1)

  if(degree != 1) { ## a l'espera d'una implementació millor:
    Comp <- matrix(0,nrow=nrow(ids),ncol=d)
    for(i in 1:nrow(ids)) {
      for(j in 1:length(jds)) {
        Comp[i,j] <-  integral((polynomial(data[ids[i,1], j:(j+degree)]) -  polynomial(data[ids[i,2],j:(j+degree)]))^2,limits=domain)

      }
    }
  } else {

    evenID <- function(number) seq(from=2,to=number,by=2)
    oddID <- function(number) seq(from=1,to=number,by=2)


    newF <- data[ids[,1],] - data[ids[,2],]

    A <- (newF[,evenID(ncol(data))])
    B <- (newF[,oddID(ncol(data))])
    C <- A*B
    A <- 1/3*(A^2)
    B <- B^2

    X <- (sup^3-inf^3)*A + (sup^2-inf^2)*C + (sup-inf)*B
  }

  if(d>1) Comp <- rowSums(X)

  K <- matrix(0,ncol=nrow(data),nrow=nrow(data))
  for(i in 1:nrow(ids)) K[ids[i,1],ids[i,2]] <- K[ids[i,2],ids[i,1]] <- Comp[i]
  rownames(K) <- colnames(K) <- rownames(data)
  if(is.null(h) || h == 0 ) return(K)
  return(exp(-h*K))
}



#Kernel selection
#' @keywords internal
kernelSelect <- function(kernel,data,domain=NULL,h=NULL) { #h és un hiperparàmetre
  if(kernel == "jac") {
    cat("quantitative Jaccard (Ruzicka) kernel\n")
    return(qJacc(data=data))
  # }else  if(kernel == "bck") {
    # cat("Bray Curtis kernel\n")
    # return(BCK(data=data))
  } else if(kernel == "jsk") {
    cat("Jensen Shannon kernel \n")
    return(JSK(data=data))
  } else if (kernel=="lin") {
    cat("Linear kernel \n")
    return(Linear(data=data))
  }  else if (kernel=="rbf") {
    cat("RBF kernel \n")
    return(RBF(data=data,h=h))
  } else if(kernel == "clin") {
    cat("Compositional Linear kernel \n")
    return(clrLin(data=data))
  } else if(kernel == "crbf") {
    cat("Aitchison-RBF kernel\n")
    return(clrRBF(data=data,h=h))
  } else if(kernel == "flin") {
    cat("Functional Linear kernel \n")
    return(Kfun(data=data,domain=domain))
  } else if(kernel == "frbf") {
    cat("Functional RBF kernel\n")
    return(RBFun(data=data,h=h,domain=domain))
  } else if(kernel == "matrix") {
    cat("Pre-computed kernel matrix given \n")
    return(as.matrix(data))
  }  else {
    stop("Kernel not available.") ##temporalment
  }
}


# Kernel computing for different types of data
#' @keywords internal

seqEval <- function(DATA,kernels,domain,h) UseMethod("seqEval",DATA)

seqEval.array <- function(DATA,kernels, domain=NULL, h=NULL) {
  m <- dim(DATA)[3]
  if(is(DATA,"matrix")) return(kernelSelect(data = DATA, domain = domain, kernel = kernels, h = h))
  n <- nrow(DATA[,,1])
  K <- array(0,dim=c(n,n,m))
  for(i in 1:m) K[,,i] <- kernelSelect(data = DATA[,,i],domain=domain,kernel = kernels[i], h=h[i])
  return(K)
}
seqEval.list <- function(DATA,kernels,domain=NULL, h=NULL) {
  m <- length(DATA)
  n <- nrow(DATA[[1]])
  K <- array(0,dim=c(n,n,m))
  for(i in 1:m) K[,,i] <- kernelSelect(data = DATA[[i]],domain=domain, kernel = kernels[i], h=h[i])
  return(K)
}

seqEval.default <- function(DATA,kernels,domain=NULL, h=NULL) return(kernelSelect(data = DATA, domain = domain, kernel = kernels, h = h))

