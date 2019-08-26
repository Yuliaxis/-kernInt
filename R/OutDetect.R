# OUTLIER / NOVELTY DETECTION

#' SVM outlier detection
#'
#'
#' @param data Input data
#' @param kernel "qJac" for quantitative Jaccard and "wqJacc" for quant Jaccard with weights
#' @param nu Hyperparameter nu
#' @return The indexes of the outliers
#' @examples
#' outliers(data=speMGX[,7:ncol(speMGX)],kernel="qJac",nu=0.9)
#' @importFrom kernlab ksvm predict
#' @export


outliers <- function(data,kernel,nu) {

  if(kernel == "qJac") {
    Jmatrix <- qJacc(data)
  }
  model <- ksvm(Jmatrix,nu=nu, type="one-svc", kernel="matrix")
  get_index <- kernlab::predict(model)
  return(which(get_index[,1]==TRUE))
}
