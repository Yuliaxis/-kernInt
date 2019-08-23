# KERNEL MATRIX INTEGRATION


#' Kernel Matrix Integration (via the mean)
#'
#' The most simple approach: mean of matrices of equal dimensions
#'
#' @param data A three-dimensional ({n·d·m}) array containing {m} kernel matrices
#' @return The consensus kernel matrix
#' @examples
#' DATA <- array(dim=c(4,5,3))
#' DATA[,,1] <- matrix(abs(rnorm(20)),nrow=4,ncol=5)
#' DATA[,,2] <- matrix(abs(rnorm(20)),nrow=4,ncol=5)
#' DATA[,,3] <- matrix(abs(rnorm(20)),nrow=4,ncol=5)
#' intMatrix <- meanInt(data=DATA)
#' @export


meanInt <- function(data) {
  d <- aperm(data,c(3,1,2))
  return(colMeans(d))
}
