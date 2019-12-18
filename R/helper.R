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


#Final tr and test indices
#' @keywords internal
ids <- function(x) UseMethod("ids",x)

ids.list <- function(x) {
  is <- rownames(x[[1]])
  if(is.null(is)) is <- 1:nrow(x[[1]])
  return(as.factor(is))
}
ids.array <- function(x) {
  is <- dimnames(x)[[1]]
  if(is.null(is)) is <- 1:nrow(x)
  return(as.factor(is))
}
ids.data.frame <- function(x) {
  is <-  rownames(x)
  if(is.null(is)) is <- 1:nrow(x)
  return(as.factor(is))
}
ids.matrix <- function(x) {
  is <- rownames(x)
  if(is.null(is))  is <- 1:nrow(x)
  return(as.factor(is))
}

#Final tr and test indices
#' @keywords internal

finalTRTE  <- function(data,p) {

  id <-ids(data)
  N <-  nlevels(id)
  # N <- nrow(data)
  all.indexes <- 1:N

  learn.indexes <- trainIndx(n=N,ptrain=p)
  test.indexes <- all.indexes[-learn.indexes]

  ##Mostres vinculades
  if(length(id) > nlevels(id)) {
    trNames <- levels(id)[learn.indexes]
    teNames <-  levels(id)[test.indexes]
    learn.indexes <- which(id %in% trNames)
    test.indexes <- which(id %in% teNames)
    ## Unique test - si es vol llevar, comentar aquestes 4 línies.
    names(test.indexes) <- id[which(id %in% teNames)]
    test.indexes <- sample(test.indexes)[teNames]
    names(test.indexes) <- NULL
    test.indexes <- sort(test.indexes)
    ## Per forçar una visita concreta:
    # test.indexes <- test.indexes[seq(from=8,to=length(test.indexes),by=8)]
  }
  return(list(li=learn.indexes,ti=test.indexes))
}

#' @keywords internal
longTRTE <- function(data,plong) {
  id <-ids(data)
  total <- length(id)
  id <- as.numeric(summary(id,maxsum=length(id)))
  help1 <- cumsum(id)-id[1]
  if(plong=="random") {
    spl <- sapply(id,function(x)sample(x,1))
  } else {
    if(length(plong) == 1) {
      spl <- rep(plong,total)
      }
    else {
      spl <- plong
      }
  }
  test.indexes <- help1 + spl
  learn.indexes <- (1:total)[-test.indexes]
  return(list(li=learn.indexes,ti=test.indexes))
}

## Kernel matrices with weights
#' @keywords internal
#' @importFrom catkern rfweight

compuKerWei <- function(data, train, y, kernel) UseMethod("compuKerWei",data)

compuKerWei.array <- function(data, train, y, kernel) {
  indw <- grep("w",kernel)
  if(sum(indw)>0) {
    w <- vector("list", dim(data)[3])
    for(i in indw) {
      w <- rfweight(x=data[train,,i],y=y,plot=FALSE)
      w <- as.vector(w)
      # names(w) <- colnames(data[,,i])
      w[[i]] <- w
    }
    return(w)
  } else {
    return(NULL)
  }
}

compuKerWei.list <- function(data, train, y, kernel) {

  indw <- grep("w",kernel)
  if(sum(indw)>0) {
    w <- vector("list", length(data))
    for(i in indw) {
      w[[i]] <- rfweight(x=data[[i]][train,],y=y,plot=FALSE)
      w[[i]] <- as.vector(w[[i]])
      # names(w) <- colnames(data[[i]])
      w[[i]] <- w}
    return(w)
  } else {
    return(NULL)
  }
}

compuKerWei.default <- function(data, train, y, kernel) {
  data <- data[train,]
  indw <- grep("w",kernel)
  if(sum(indw)==1) {
    w <- rfweight(x=data,y=y,plot=FALSE)
    w <- as.vector(w)
    # names(w) <- colnames(data)
    return(w)
  } else {
    return(NULL)
  }
}

# Class imbalance: data approach
#' @keywords internal
smote <- function(data=data, diagn=diagn, nlearn=nlearn, N=N, learn.indexes,test.indexes) {
  dades <- data[c(learn.indexes,test.indexes),]
  SobrDadesTr <- ubBalance(as.data.frame(dades[1:nlearn,]), diagn[1:nlearn], type="ubSMOTE", positive=2)

  data <- rbind(SobrDadesTr$X,dades[(nlearn+1):N,])
  nlearn <- length(SobrDadesTr$Y)
  diagn <- c(SobrDadesTr$Y, diagn[test.indexes])
  diagn <- as.factor(diagn)
  return(list(data=data,diagn=diagn,nlearn=nlearn))
}
#
smoteSample <- function(data, diagn, learn.indexes, m, test.indexes, kernel) {
  nlearn <- length(learn.indexes)
  ntest <- length(test.indexes)
  N <- nlearn + ntest
  diagn <- diagn[c(learn.indexes,test.indexes)]
  if(m>1) {
    stop("Option not available.")
  } else {
    Sample <- smote(data=data, diagn=diagn, nlearn=nlearn, N=N, learn.indexes,test.indexes)
  }
  data <- Sample$data
  diagn <- Sample$diagn
  nlearn <- Sample$nlearn
  N <- nrow(data)
  learn.indexes <- 1:nlearn
  test.indexes <- (nlearn+1):N
  return(list(data=data,y=diagn,li=learn.indexes,ti=test.indexes))
}

#' @keywords internal
dataSampl <- function(data, diagn, tedata, kernel, type)  UseMethod("dataSampl",data)

dataSampl.array <- function(data, tedata, diagn, kernel, type) {

    dadespr <- data[,,1]
    rownames(dadespr) <- 1:nrow(dadespr)

    SobrDadesTr <- ubBalance(as.data.frame(dadespr), diagn, type=type, positive=2, k=0)

    ii <- as.numeric(rownames(SobrDadesTr$X))
    data <- data[ii,ii,]
    tedata <- tedata[,ii,]
    diagn <- diagn[ii]
  return(list(data=data,tedata=tedata,diagn=SobrDadesTr$Y))
}

dataSampl.default <- function(data, tedata, diagn, kernel, type) {

  rownames(data) <- 1:nrow(data)
  SobrDadesTr <- ubBalance(as.data.frame(data), diagn, type=type, positive=2, k=0)

  ii <- as.numeric(rownames(SobrDadesTr$X))
  data <- data[ii,ii]
  tedata <- tedata[,ii]
  return(list(data=data,tedata=tedata,diagn=SobrDadesTr$Y))
}


## Hyperparameters depending on kernel selected
#' @keywords internal

hyperkSelection <- function(K, h, kernel) {
  if (is.null(h) || h==0) {
    return(K)
  }
  if(length(h)>1) {
    paste("More kernel hyperparameters than kernel functions - Only the first element will be used")
    h <- h[1]
  }
  if(kernel == "qJac" | kernel == "wqJac") {
    Kmatrix <- exp(h*K)/exp(h) #Standardized Kernel Matrix. Otherwise exp(g*K)
  } else if(kernel == "cRBF" | kernel == "time" | kernel == "time2") {
    Kmatrix <- exp(-h*K)
    Kmatrix[is.na(Kmatrix)] <- 0
  } else if(kernel == "rbf") {
    Kmatrix <- exp(h*K)
  } else if(kernel == "cov") {
    Kmatrix <- K
    Kmatrix[Kmatrix!=0] <- h
    diag(Kmatrix) <- 1
  } else if(kernel == "matrix" | kernel == "linear") {
    Kmatrix <- K
  } else {
    stop("Kernel not available.")
  }
  return(Kmatrix)
}
