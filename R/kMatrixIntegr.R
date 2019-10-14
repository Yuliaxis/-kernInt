# KERNEL MATRIX INTEGRATION

#' Kernel Matrix Integration
#'
#' @param data A three-dimensional ({n·d·m}) array containing {m} kernel matrices
#' @param coeff A vector of length {m} with the weight of each kernel matrix
#' If absent, the mean of all matrices is computed instead.
#' @return A consensus kernel matrix
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
#' @export


KInt <- function(data,coeff) {
  data[which(is.na(data))] <- 0
  d <- aperm(data,c(3,1,2))
  if(hasArg(coeff)) {
    if(dim(d)[1] != length(coeff)) stop("Length of the coefficients vector different to the number of matrices")
    if(sum(coeff) != 1) coeff <- coeff / sum(coeff)
    IntMatrix <-  colSums(coeff*d)
  } else {
    IntMatrix <-  colMeans(d)
  }
  return(IntMatrix)
}


#' Fuse Data
#'
#' @param DATA A list of the *m* data to fuse
#' @param kernels A vector of length *m* with the kernels to use in each data
#' @param y Only if "wqJac" is chosen: response variable
#' @param h Kernel hyperparameter (gamma)
#' @param coeff A vector of length *m* with the weight of each kernel data,
#' or length(m) - 1 if the last coefficient is 1 - sum(coeff).
#' If absent,  all data is considered equally important.
#' @return A consensus kernel matrix (via calling KInt)
#' @examples
#' d <- list()
#' d[[1]] <- matrix(abs(rnorm(20)),nrow=4,ncol=5)
#' d[[2]] <- matrix(abs(rnorm(20)),nrow=4,ncol=5)
#' fuseData(DATA=d,kernel=c("qJac","cRBF"))
#' fuseData(DATA=d,kernel=c("qJac","cRBF"),coeff=c(0.9,0.1))
#' @export
#'

fuseData <- function(DATA,coeff,kernels,y,h) {
  Kmatrix <- seqEval(DATA,kernels,y,h)
  m <- length(DATA)

  if(hasArg(coeff)) {
    print("Coeff")
    if(length(coeff) == (m - 1) && coeff < 1) coeff <- c(coeff,1-coeff)
    return(KInt(Kmatrix,coeff=coeff))
  } else {
    return(KInt(Kmatrix))
  }
}
