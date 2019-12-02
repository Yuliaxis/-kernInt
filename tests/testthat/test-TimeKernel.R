context("Least Squares Coefficients procedure")

test_that("lqs() output is formally correct", {
  formal_out <- function(result,degree) {
    expect_is(result,"matrix")
    expect_equal(ncol(result),d+1)
    expect_false(is.null(rownames(result)))
    expect_equal(length(unique(rownames(result))),nrow(result))

  }

  colnames(growth) <-  c("sex", "time", "height")

  d <- 2
  formal_out(lsq(data=growth[,2:3],degree=d),d)

  d <- 4
  formal_out(lsq(data=growth[,2:3],degree=d),d)

  d <- 1
  formal_out(lsq(data=growth[,2:3],degree=d),d)

})

test_that("lsq() output", {
  out <- matrix(c(131.4376, 3.0872, -5.0727, 129.5317, 2.8262, -5.7132, 131.0097, 1.2394, -3.7522),
                nrow=3,ncol=3,byrow=TRUE)
  rownames(out) <- c("boy36","boy20","girl43")
  colnames(out) <- c("Intercept", "x.1", "x.2")

  colnames(growth) <-  c("sex", "time", "height")

  d <- 2
  expect_equal(round(lsq(data=growth[1:24,2:3],degree=2,orthog=TRUE),digits=4),out)
})
