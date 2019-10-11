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

ids.list <- function(x) return(as.factor(rownames(x[[1]])))
ids.array <- function(x) return(as.factor(dimnames(x)[[1]]))
ids.data.frame <- function(x) return( as.factor(rownames(x)))
ids.matrix <- function(x) return(as.factor(rownames(x)))

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
    #Unique test - si es vol llevar, comentar aquestes 4 línies.
    names(test.indexes) <- id[which(id %in% teNames)]
    test.indexes <- sample(test.indexes)[teNames]
    names(test.indexes) <- NULL
    test.indexes <- sort(test.indexes)
  }
  return(list(li=learn.indexes,ti=test.indexes))
}

#Final tr and test indices
#' @keywords internal

sampl <- function(data, diagn, learn.indexes, test.indexes, kernel, type) {
  nlearn <- length(learn.indexes)
  ntest <- length(test.indexes)
  N <- nlearn + ntest
  diagn <- diagn[c(learn.indexes,test.indexes)]
  print(nlearn)
  Sample <- dataSampl(data, diagn, nlearn=nlearn, N=N, learn.indexes,test.indexes, kernel=kernel, type)

  data <- Sample$data
  diagn <- Sample$diagn
  nlearn <- Sample$nlearn
  N <- nrow(data)

  learn.indexes <- 1:nlearn
  test.indexes <- (nlearn+1):N
  return(list(data=data,y=diagn,li=learn.indexes,ti=test.indexes))
}

dataSampl <- function(data, diagn, nlearn, N, learn.indexes,test.indexes, kernel, type)  UseMethod("dataSampl",data)

dataSampl.array <- function(data, diagn, nlearn, N, learn.indexes,test.indexes, kernel, type) {

  if(kernel == "matrix") {
    if(type == "ubSMOTE") stop("Kernel matrix as input is not compatible with SMOTE. Original dataset is required.")

    dades <- data[c(learn.indexes,test.indexes),c(learn.indexes,test.indexes),]
    dadespr <- dades[,,1]
    rownames(dadespr) <- 1:N

    if(type == "ubOver")  SobrDadesTr <- ubBalance(dadespr[1:nlearn,], diagn[1:nlearn], type=type, positive=2, k=0)
    if(type == "ubUnder")  SobrDadesTr <- ubBalance(dadespr[1:nlearn,], diagn[1:nlearn], type=type, positive=2)

    ii <- c(as.numeric(rownames(SobrDadesTr$X)),(nlearn+1):N)
    data <- data[ii,ii,]
    diagn <- diagn[ii]
    # N <- nrow(data)
    nlearn <- length(SobrDadesTr$Y)
  } else {
    stop("Option not available yet.")
  }
  return(list(data=data,diagn=diagn,nlearn=nlearn))
}


dataSampl.list <- function(data, diagn, nlearn, N, learn.indexes,test.indexes, kernel, type) {
  m <- length(data)

    if(type == "ubSMOTE") {
      SobrDadesTr <- list()
      for(i in 1:m) SobrDadesTr[[i]] <- ubBalance(dades[[i]][1:nlearn,], diagn[1:nlearn], type=type, positive=2)
    } else {
      dadespr <- dades[[1]]
      rownames(dadespr) <- 1:N
      if(type == "ubOver")  SobrDadesTr <- ubBalance(dadespr[1:nlearn,], diagn[1:nlearn], type=type, positive=2, k=0)
      if(type == "ubUnder")  SobrDadesTr <- ubBalance(dadespr[1:nlearn,], diagn[1:nlearn], type=type, positive=2)
      ii <- c(as.numeric(rownames(SobrDadesTr$X)),(nlearn+1):N)
      diagn <- diagn[ii]
      nlearn <- length(SobrDadesTr$Y)
      for(i in 1:m) data[[i]] <- data[[i]][ii,]
    }
  return(list(data=data,diagn=diagn,nlearn=nlearn))
}

dataSampl.default <- function(data, diagn, nlearn, N, learn.indexes,test.indexes, kernel, type) {
  if(kernel == "matrix") {
    if(type == "ubSMOTE") stop("Kernel matrix as input is not compatible with SMOTE. Original dataset is required.")

    dades <- data[c(learn.indexes,test.indexes),c(learn.indexes,test.indexes)]
    rownames(dades) <- 1:N

    if(type == "ubOver")  SobrDadesTr <- ubBalance(dades[1:nlearn,], diagn[1:nlearn], type=type, positive=2, k=0)
    if(type == "ubUnder")  SobrDadesTr <- ubBalance(dades[1:nlearn,], diagn[1:nlearn], type=type, positive=2)

    ii <- c(as.numeric(rownames(SobrDadesTr$X)),(nlearn+1):N)
    data <- data[ii,ii]
    diagn <- diagn[ii]
    # N <- nrow(data)
    nlearn <- length(SobrDadesTr$Y)
  } else {
    dades <- data[c(learn.indexes,test.indexes),]

    if(type == "ubUnder") SobrDadesTr <- ubBalance(dades[1:nlearn,], diagn[1:nlearn], type=type, positive=2)
    if(type == "ubOver") SobrDadesTr <- ubBalance(dades[1:nlearn,], diagn[1:nlearn], type=type, positive=2,  k=0)
    if(type == "ubSMOTE") SobrDadesTr <- ubBalance(dades[1:nlearn,], diagn[1:nlearn], type=type, positive=2)
    data <- rbind(SobrDadesTr$X,dades[(nlearn+1):N,])
    nlearn <- length(SobrDadesTr$Y)
    # N <- nrow(data)
    diagn <- c(SobrDadesTr$Y, diagn[test.indexes])
    diagn <- as.factor(diagn)
  }
  return(list(data=data,diagn=diagn,nlearn=nlearn))

}


## Hyperparameters depending on kernel selected
#' @keywords internal

hyperkSelection <- function(K, h, kernel) {
  if(kernel == "qJac" | kernel == "wqJac") {
    if (h == 0) {
      Kmatrix <- K
    } else { Kmatrix <- exp(h*K)/exp(h) #Standardized Kernel Matrix. Otherwise exp(g*K)
    }
  } else if(kernel == "cRBF" | kernel == "time" | kernel == "time2") {
    Kmatrix <- exp(-h*K)
  } else if(kernel == "cov") {
    Kmatrix <- K
    Kmatrix[Kmatrix!=0] <- h
    diag(Kmatrix) <- 1
  } else if(kernel == "matrix") {
    Kmatrix <- K
  } else {
    stop("Kernel not available.")
  }
  return(Kmatrix)
}

##  NMSE (regression)
#' @keywords internal
error.norm <- function(target,prediction) {
  N <- length(target)
  error <- sum((target-prediction)^2)/((N-1)*var(target)) ##(norm.mse <- model$deviance/((N-1)*var(target)))
  return(error)
}


## Harmonic mean, to compute the F1 accuracy:
#' @keywords internal
harm <- function (a,b) { 2/(1/a+1/b) }
