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


# Check input data
#' @keywords internal
checkinput <- function(data,y, kernel) {

  if(class(data) == "list") {
    m <- length(data)
    if(m < 2) {
      data <- unlist(data)
    } else {
      if(unique(sapply(data,class)) != "lsq") {
      ## Comprova tots els elements tenen el mateix nombre de files
      n <- unique(sapply(data,nrow))
      } else {
        n <- rep(0,m)
        for(x in 1:m) {
          n[x] <- data[[x]]$coeff
        }
        n <- unique(n)
      }
      if(length(n) != 1) stop("Elements of the list have different number of rows")
    }

  } else if(class(data) == "array") {
    n <- dim(data)[1]
    m <- dim(data)[3]
    if(m < 2) data <- matrix(data[,,1],ncol=dim(data)[2],nrow=n)
  } else if(class(data) == "data.frame" | class(data) == "matrix") {
    m <- 1
    n <- nrow(data)
  } else if(class(data) == "lsq") {
    m <- 1
    n <- nrow(data$coef)
  } else {
    stop("Wrong input data class.")
  }

  if(length(y) != n) stop("Length of the target variable do not match with the row number of predictors")

  if(length(kernel)>m) stop("kernel argument is longer than number of different datasets")

  if(length(kernel)<m) kernel <- rep(kernel,ceiling(m/length(kernel)))[1:(m)]

  return(list(data=data,m=m,kernel=kernel))
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
    # names(test.indexes) <- id[which(id %in% teNames)]
    # test.indexes <- sample(test.indexes)[teNames]
    # names(test.indexes) <- NULL
    # test.indexes <- sort(test.indexes)
    ## Per forçar una visita concreta:
    # test.indexes <- test.indexes[seq(from=8,to=length(test.indexes),by=8)]
  }
  return(list(li=learn.indexes,ti=test.indexes))
}

#' @keywords internal
longTRTE <- function(data,plong) {
  id <-ids(data)
  total <- length(id)
  id <- as.numeric(summary(id,maxsum=length(id))[unique(id)])
  help1 <- cumsum(id)-id
  if(plong=="random") {
    spl <- sapply(id,function(x)sample(x,1)) # test a l'atzar
  } else {
    if(length(plong) == 1) {
      spl <- rep(plong,length(id)) # la mateixa mostra per a tots els individus
      }
    else {
      spl <- plong  ## si volem indicar una mostra específica per individu per usar de test
      }
  }
  test.indexes <- help1 + spl
  learn.indexes <- (1:total)[-test.indexes]
  return(list(li=learn.indexes,ti=test.indexes))
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
  if (kernel == "rbf" | kernel == "crbf" | kernel == "frbf") {
    Kmatrix <- exp(-h*K)
    Kmatrix[is.na(Kmatrix)] <- 0
  }  else {
    Kmatrix <- K
  }
  return(Kmatrix)
}
