context("Regression procedure")

test_that("Regression output is formally correct (i.e. NMSE)", {
  formal_out <- function(result) {
    expect_is(result,"numeric")
    expect_equal(length(result),1)
    expect_true(result >=0 )
  }

  #KERNEL = MATRIX

  km <- (grKern[[1]] + grKern[[3]])/2

  ## Basic
  formal_out(regress(data=km, y=growth[,2],  kernel="matrix", C=0.1))

  ## Regressió simple amb cross-val (Cost i epsilon)
  formal_out(regress(data=km, y=growth[,2], kernel="matrix", C=c(0.1,1),k=5))
  formal_out(regress(data=km, y=growth[,2], kernel="matrix", C=c(0.1,1), E=c(0.001,0.01), k=5))

  ## Type A approach: tenim en compte el component longitudinal però usam el model per a predir un nou individu.
  formal_out(regress(data=grKern[c(1,3,4,5)], y=growth[,2], kernel=c("matrix","matrix","matrix","matrix"), C=1))

  ## Amb Cros-Val
  formal_out(regress(data=grKern[c(1,3,4,5)], y=growth[,2], coeff=c(0.33,0.33,0.17,0.17), kernel=c("matrix","matrix","matrix","matrix"), C=c(0.1,1),E=c(0.01,0.1),k=5))
  w <- matrix(c(0.7,rep(0.1,3),0.4,0.2,0.25,0.15),ncol=4,nrow=2,byrow = TRUE)
  formal_out(regress(data=grKern[c(1,3,4,5)], y=growth[,2], coeff=w, kernel=c("matrix","matrix","matrix","matrix"), C=c(0.1,1),k=5))

  ## Type B:  predim l'altària d'un individu ja conegut en una visita X desconeguda.
  formal_out(regress(data=grKern[c(1,3,4,5)], y=growth[,2], coeff=c(0.33,0.33,0.17,0.17), kernel=c("matrix","matrix","matrix","time2"), C=1))

  #KERNEL != MATRIX

  TM <- TimeK(grKern[[5]])
  TM <- exp(-0.5*TM)
  TM[is.na(TM)] <- 0

  INPUT <- list()
  INPUT[[1]] <- grKern[[1]]
  INPUT[[2]] <-  grKern[[3]]
  INPUT[[3]] <-  grKern[[4]]
  INPUT[[4]] <- TM

  regress(data=INPUT, y=growth[,2], coeff=c(0.3,0.3,0.2,0.2), kernel=c("matrix","matrix","cov","matrix"),C=0.1)

  ######### Amb Cros-Val

  w <- matrix(c(0.4,0.4,0.1,0.1,0.33,0.33,0.17,0.17),ncol=4,nrow=2,byrow = TRUE)
  regress(data=INPUT, y=growth[,2], coeff=w, kernel=c("matrix","matrix","cov","matrix"),C=0.1,k=5)
  regress(data=INPUT, y=growth[,2],  kernel=c("matrix","matrix","cov","matrix"), C=0.1, H=list(a=0,b=0,c=c(0.25,0.5),d=0), k=5)
  regress(data=INPUT, y=growth[,2],  coeff=w,  kernel=c("matrix","matrix","cov","matrix"),C=0.1,H=list(a=0,b=0,c=c(0.25,0.5),d=0), k=5)
  regress(data=INPUT, y=growth[,2],  coeff=w,  kernel=c("matrix","matrix","cov","matrix"),H=list(a=0,b=0,c=c(0.25,0.5),d=0),C=c(0.1,0.01), k=5)


  ##### Type B:  predim l'altària d'un individu ja conegut en una visita X desconeguda.

  regress(data=grKern[c(1,3:5)], y=growth[,2], coeff=c(0.33,0.33,0.17,0.17), kernel=c("matrix","matrix","cov","time"), C=0.1)

  })
