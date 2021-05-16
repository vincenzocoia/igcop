context("Check that the quantile functions are the inverses of the conditional distributions.")

test_that("qcondig21 is the inverse of pcondig21", {
  for (cpar_ in .cpar){
    tau2 <- pcondig(qcondig(.v, .u, cpar = cpar_), .u, cpar = cpar_)
    expect_equal(.v, tau2)
  }
})


test_that("qcondigl21 is the inverse of pcondigl21", {
  for (alpha_ in .alpha){
    tau2 <- pcondigl(qcondigl(.v, .u, cpar = alpha_), .u, cpar = alpha_)
    expect_equal(.v, tau2)
  }
})

test_that("qcondig12 is the inverse of pcondig12", {
  for (cpar_ in .cpar){
    tau2 <- pcondig12(qcondig12(.u, .v, cpar = cpar_), .v, cpar = cpar_)
    expect_equal(.u, tau2)
  }
})


test_that("qcondigl12 is the inverse of pcondigl12", {
  for (alpha_ in .alpha){
    tau2 <- pcondigl12(qcondigl12(.u, .v, cpar = alpha_), .v, cpar = alpha_)
    expect_equal(.u, tau2)
  }
})
