# KERNEL MATRIX INTEGRATION

#' Kernel Matrix Integration
#'
#' @param data A three-dimensional ({n·d·m}) array containing {m} kernel matrices
#' @param coeff A vector of length {m} with the weight of each kernel matrix
#' @return A kernel matrix
#' @examples
#' DATA <- array(dim=c(4,4,3))
#' DATA[,,1] <- matrix(abs(rnorm(16)),nrow=4,ncol=4)
#' diag(DATA[,,1]) <- 1
#' DATA[,,2] <- matrix(abs(rnorm(16)),nrow=4,ncol=4)
#' diag(DATA[,,2]) <- 1
#' DATA[,,3] <- matrix(abs(rnorm(16)),nrow=4,ncol=4)
#' diag(DATA[,,3]) <- 1
#' KInt(data=DATA)
#' KInt(data=DATA,coeff=c(0.5,0.3,0.2))
#'@importFrom methods is
#' @export


KInt <- function(data,coeff="mean") {

  d <- aperm(data,c(3,1,2))

  if(!is(coeff,"numeric")) return(colMeans(d))

  if(dim(d)[1] != length(coeff)) stop("Length of the coefficients
                                        vector different to the number of matrices")
  if(sum(coeff) != 1) coeff <- coeff / sum(coeff)
  return(colSums(coeff*d))
}

#' Fuse Data
#'
#' @param DATA A list of the *m* data to fuse
#' @param kernels A vector of length *m* with the kernels to use in each data
#' @param h Kernel hyperparameter (gamma)
#' @param coeff A vector of length *m* with the weight of each kernel data,
#' or length(m) - 1 if the last coefficient is 1 - sum(coeff).
#' If absent,  all data is considered equally important.
#' @return A consensus kernel matrix (via calling KInt)
#' @examples
#' d <- list()
#' d[[1]] <- matrix(abs(rnorm(20)),nrow=4,ncol=5)
#' d[[2]] <- matrix(abs(rnorm(20)),nrow=4,ncol=5)
#' fuseData(DATA=d,kernel=c("jac","crbf"))
#' fuseData(DATA=d,kernel=c("jac","crbf"),coeff=c(0.9,0.1))
#' @export
#'

fuseData <- function(DATA,coeff="mean",kernels, h=NULL) {
  Kmatrix <- seqEval(DATA=DATA,kernels=kernels,h=h)
  m <- length(DATA)
  return(KInt(Kmatrix,coeff=coeff))
}

