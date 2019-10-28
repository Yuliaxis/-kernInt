context("Kernel functions")

test_that("Ruzicka kernel works", {
  m <- matrix(c(0.333,0.2,0.1,0.333,0.5,0.1,0.333,0.3,0.8),nrow=3,ncol=3)
  res <- matrix(c(1,0.71441,0.36357,0.71441,1,0.33333,0.36357,0.33333,1),nrow=3,ncol=3)
  resg1 <-  round(exp(res)/exp(1),digits=5)
  expect_equal(round(qJacc(m),digits=5), res)
  expect_equal(round(qJacc(m,h=0),digits=5), res)
  expect_equal(round(qJacc(m,h=1),digits=5), resg1)
})


test_that("General Evaluation kernel procedure works", {
  m <- matrix(c(0.333,0.2,0.1,0.333,0.5,0.1,0.333,0.3,0.8),nrow=3,ncol=3)
  res <- matrix(c(1,0.71441,0.36357,0.71441,1,0.33333,0.36357,0.33333,1),nrow=3,ncol=3)
  resg1 <-  round(exp(res)/exp(1),digits=5)
  resl <- matrix(c(1,0.93659,0.71067,0.93659,1,0.61901,0.71067,0.61901,1),nrow=3,ncol=3)
  expect_equal(round(kernelSelect(kernel="qJac",data=m),digits=5), res)
  expect_equal(round(kernelSelect(kernel="qJac",data=m,h=1),digits=5), resg1)
  expect_equal(round(kernelSelect(kernel="linear",data=m),digits=5),resl)
  expect_error(kernelSelect(kernel="none",data=m),"Kernel not available.")
})


test_that("Sequential Evaluation kernel procedure works", {
  M <-  array(dim=c(3,3,2))
  R <- M
  M[,,1] <- matrix(c(0.333,0.2,0.1,0.333,0.5,0.1,0.333,0.3,0.8),nrow=3,ncol=3)
  M[,,2] <- M[,,1]
  R[,,1] <- matrix(c(1,0.71441,0.36357,0.71441,1,0.33333,0.36357,0.33333,1),nrow=3,ncol=3)
  R[,,2] <- matrix(c(1,0.93659,0.71067,0.93659,1,0.61901,0.71067,0.61901,1),nrow=3,ncol=3)
  expect_equal(round(seqEval(DATA=M,kernels = c("qJac","linear")),digits=5), R)
  R[,,1] <-  round(exp(R[,,1])/exp(1),digits=5)
  expect_equal(round(seqEval(DATA=M,kernels = c("qJac","linear"),h=1),digits=5), R)

})
