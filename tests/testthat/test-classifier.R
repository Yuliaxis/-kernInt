context("Classification procedure")

test_that("Classification output is formally correct (i.e. confusion table)", {
  formal_out <- function(result) {
    expect_is(result,"table")
    expect_equal(dim(result),c(2,2))
    expect_true(!is.null(rownames(result)))
    expect_equal(rownames(result),colnames(result))
  }

  # KERNEL = MATRIX

  ## Cas bàsic
  formal_out(classify(data=grKern[[2]], y=growth[,1], kernel="matrix"))

  ## Classification with cross-val (Cost)
  formal_out(classify(data=grKern[[2]], y=growth[,1],  kernel="matrix", C=c(0.05,0.01,0.1), k=5))

  # Unbalanced Data treatment -undersample
  formal_out(classify(data=grKern[[2]],y=growth[,1], kernel="matrix",classimb="data",type="ubUnder",C=c(0.001,0.01),k=5))

  # Unbalanced Data treatment - oversample
  formal_out(classify(data=grKern[[2]],y=growth[,1], kernel="matrix",classimb="data",type="ubOver",C=c(0.001,0.01),k=5))

  # Unbalanced Data treatment - algorithm weighting
  formal_out( classify(data=grKern[[2]],y=growth[,1], kernel="matrix",classimb="weights",C=c(0.001,0.01),k=5))

  # Unbalanced Data treatment - soft classification
  formal_out(classify(data=grKern[[2]],y=growth[,1], prob=TRUE, kernel="matrix",CUT=0.4, C=c(0.001,0.01),k=5))

  ## Basic MKL Classification (long type A approach)

  TM <- TimeK(grKern[[5]])
  TM <- exp(-0.5*TM)
  TM[is.na(TM)] <- 0

  INPUT <- list()
  INPUT[[1]] <- grKern[[2]]
  INPUT[[2]] <- covmat(ind=93,visits=8)
  INPUT[[3]] <- TM

  formal_out(classify(data=INPUT,y=growth[,1] ,kernel=c("matrix","matrix","matrix"), C=0.5))

  ## MKL Classification (type A approach) with cross-validation
  formal_out(classify(data=INPUT,y=growth[,1] ,kernel=c("matrix","matrix","matrix"), C=c(1,0.5),k=5))

  ## MKL Classification (type A approach) with different coefficients
  formal_out(classify(data=INPUT,y=growth[,1],coeff=c(0.5,0.25,0.25),kernel=c("matrix","matrix","matrix"),C=0.1))

  ## MKL Classification (type A approach) with different coefficients AND cross-val
  w <- matrix(c(0.8,0.1,0.1,0.4,0.3,0.3),ncol=3,nrow=2,byrow = TRUE)
  formal_out(classify(data=INPUT,y=growth[,1],coeff=w,kernel=c("matrix","matrix","matrix"), C=c(0.01,0.1,0.5),k=5))

  ## Basic MKL Classification (long type A approach) - soft and cross-val
  formal_out(classify(data=INPUT,y=growth[,1], prob=TRUE, CUT=0.4, kernel=c("matrix","matrix","matrix"), C=c(0.1,0.5),k=5))

  ## Basic MKL Classification (long type A approach) - oversampling
  formal_out(classify(data=INPUT,y=growth[,1],kernel=c("matrix","matrix","matrix"), classimb="data",type="ubOver",C=0.1))

  #KERNEL != MATRIX

  ## Cas bàsic
  classify(data=growth[,2:3], y=growth[,1],  kernel="linear", C=0.1)

  #### Classificació simple amb cross-val (Cost)
  classify(data=growth[,2:3], y=growth[,1],  kernel="linear", C=c(0.01,0.001),k=5)

  #### Classificació simple amb cross-val (kernel hyperparameter)
  classify(data=growth[,2:3], y=growth[,1],  kernel="rbf", C=c(0.01,0.1), H=c(0.05,0.01), k=5)

  ## Classificació amb component longitudinal ("type A approach")

  ####### mescla matrius precomputades i de computació interna
  INPUT <- list()
  INPUT[[1]] <- grKern[[2]]
  INPUT[[2]] <- c(93,8)
  INPUT[[3]] <- grKern[[5]]
  classify(data=INPUT,y=growth[,1],kernel=c("matrix","cov","matrix"), C=0.5,H=list(a=0,b=0.25,c=0))

  #### Classificació tipus A amb cross-val (Cost)
  classify(data=INPUT,y=growth[,1] ,kernel=c("matrix","cov","matrix"), C=c(0.1,0.05),H=list(a=0,b=c(0.25,0.5),c=0),k=5)

  ####### tot amb matrius de computació interna
  INPUT2 <- list()
  INPUT2[[1]] <- growth[,2:3]
  INPUT2[[2]] <- grKern[[4]]
  INPUT2[[3]] <- grKern[[5]]
  classify(data=INPUT2,y=growth[,1] ,kernel=c("rbf","cov","time"), C=0.5,H=list(a=0.01,b=0.25,c=0.5))

  #### Amb cross-val
  classify(data=INPUT2,y=growth[,1] ,kernel=c("rbf","cov","time"), C=0.5, H=list(a=c(0.1,0.01),b=c(0.75,0.5,0.25),c=0.5),k=5)
})
