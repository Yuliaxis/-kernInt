# TIME KERNEL MATRICES

#' Time Kernel Matrix
#'
#' @param data  A n·d data.frame or matrix with n individuals and d ordered visits, stating if the individual came to the
#' visit (TRUE or 1) or no (FALSE or 0). The time must be specified as the column names.
#' @param h If argument, the functions returns exp(-h * |tj-ti|)
#' @return A time matrix
#' @examples
#' DATA <- matrix(c(rep(1,5),rep(0,4),rep(1,4),rep(0,4),1,0,0,1,1,0,1,1,rep(0,4),1),ncol=6,nrow=5)
#' colnames(DATA) <- c(0:3,5,9) # Measured in weeks
#' rownames(DATA) <- paste("Indiv",1:5,sep="")
#' TimeK(data=DATA)
#' @export

TimeK <- function(data,h) {
  totaltime <- as.numeric(colnames(data))
  totalrows <- rowSums(data)
  data[!data] <- NA
  Kprova <- (t(data) * totaltime)
  visits <- Kprova[!is.na(Kprova)]
  sta <- 1
  end <- 0
  index <- matrix(ncol=2,nrow=0) #indexos de K que siran diferents de 0
  for(i in 1:length(totalrows)) {
    end <- end + totalrows[i]
    cosa <- expand.grid.mod(sta:end,rep=FALSE)
    sta <- sta + totalrows[i]
    index <- rbind(index,cosa)
  }
  values <- visits[index[,2]] - visits[index[,1]]
  K <- matrix(0,nrow=sum(totalrows),ncol=sum(totalrows))
  colnames(K) <- visits
  rownames(K) <- colnames(K)
  if(hasArg(h)) {
    for(i in 1:length(values)) K[index[i,1],index[i,2]] <- exp(-h*values[i])
  } else {
      for(i in 1:length(values)) K[index[i,1],index[i,2]] <- abs(values[i])
    }
  K <- K + t(K)
  if(hasArg(h))  { diag(K) <- 1
  } else {
    K[K == 0] <- NA
    diag(K) <- 0 }
  return(K)
}

#' Minumum Separation Samples
#'
#' @param table  A n·d data.frame or matrix with n individuals and d ordered visits, stating if the individual came to the
#' visit (TRUE or 1) or no (FALSE or 0). The time must be specified as the column names.
#' @param n Number of samples
#' @return A vector containing, for each individual, the name of the first sample of the minimum period
#' @examples
#' individuals <- which(rowSums(visitTable) > 3)
#' minSep(table=visitTable[individuals,],n=4) # The four samples with minimum separation
#' @export

minSep <- function(table,n) {
  tMatrix <- TimeK(table)
  t <- n-1
  points <- c(1,which(as.numeric(colnames(tMatrix))[-1] - as.numeric(colnames(tMatrix))[-ncol(tMatrix)]<0),ncol(tMatrix))
  tMatrix <- tMatrix[-1,]
  tMatrix <- tMatrix[,-ncol(tMatrix)]
  INDEX <- c()
  j <- 0
  for(i in 1:(length(points)-1)) {
    indiv <-   (points[i]+j):(points[i+1]-1)
    diagon <- diag( tMatrix[indiv,indiv] )
    print(diagon)
    ind <- length(diagon)
    ind.sum <- Inf
    for(k in length(diagon):t) {
      if(sum(diagon[(k-t+1):k]) <= ind.sum) {
        ind.sum <- sum(diagon[(k-t+1):k])
        ind <- k-t+1
      }
    }
    print(ind)
    print(colnames(tMatrix)[indiv[ind]])
    INDEX <- c(INDEX,colnames(tMatrix)[indiv[ind]])
    j <- 1
  }
  names(INDEX) <-  rownames(table)
  return(INDEX)

}
