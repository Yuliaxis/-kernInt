# LSQ for time series


#### lsq helper
#' @keywords internal

timeserieshelp <- function(data) {
   condition <- colnames(data) %in% "time"
   if(sum(condition) != 1) stop("Time column needed")
   condition <- which(condition)
   data <- data[,c(condition,(1:ncol(data))[-condition])]
  return(data)
}

#####

### LSQ coefficients
#' Least Squares Coefficients
#'
#' @param data  A matrix with longitudinal data, with a column named "time" containing the time
#' (in arbitrary units) of the measure. lsq() will consider that the row names that share id are repeated
#' measures at different times from the same individual.
#' @param degree The degree of the polynomial to be fitted.
#' @param orthog If true, use orthogonal polynomials.
#' @param scale If true, divide each coefficient by its variance
#' @return An object of class lsq containing the intercept and coefficients over time of for each variable and individual.
#' @examples
#' growth2 <- growth[,2:3]
#' colnames(growth2) <-  c( "time", "height")
#' lsq(data=growth2,degree=2)
#' @importFrom stats as.formula coefficients lm var
#' @export

lsq <- function(data, degree=1, orthog=FALSE,scale=FALSE) {

  data <- timeserieshelp(data)
  ids <- unique(rownames(data))
  cnames <- colnames(data)[-1]
  COEFF <- matrix(NA,nrow=length(ids),ncol=length(cnames)*(degree+1))
  rownames(COEFF) <- ids
  colnames(COEFF) <- rep(c("Intercept",paste("x",seq(1:degree),sep=".")),length(cnames))
  for(ID in ids) {
    subdata <- data[which(rownames(data) == ID),]
    subdata <- as.data.frame(subdata)
    for(cn in cnames) {
      model <- lm(as.formula(paste(paste("subdata",cn,sep="$")," ~  poly(subdata$time,degree,raw=!orthog)")))
      cols <- which(cn==cnames)
      COEFF[ID, ((degree+1)*(cols-1)+1):(cols*(degree+1))] <- coefficients(model)
    }
  }
  if(scale) {
    std <- sqrt(apply(COEFF,2,var))
    COEFF <- t(t(COEFF)/std)
  }
  ob <- list(coef=COEFF,degree=degree,y.names=cnames,domain=range(data[,1]))
  attr(ob, "class") <- "lsq"
  return(ob)
}

