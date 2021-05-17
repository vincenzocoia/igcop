context("Check that the quantile functions are the inverses of the conditional distributions.")

test_that("qcondig21 is the inverse of pcondig21", {
  for (cpar_ in .cpar){
    qu <- qcondig(.v, .u, theta = cpar_[1], alpha = cpar_[2])
    tau2 <- pcondig(qu, .u, theta = cpar_[1], alpha = cpar_[2])
    expect_equal(.v, tau2)
  }
})


test_that("qcondigl21 is the inverse of pcondigl21", {
  for (alpha_ in .alpha){
    tau2 <- pcondigl(qcondigl(.v, .u, alpha = alpha_), .u, alpha = alpha_)
    expect_equal(.v, tau2)
  }
})

test_that("qcondig12 is the inverse of pcondig12", {
  for (cpar_ in .cpar){
    qu <- qcondig12(.u, .v, theta = cpar_[1], alpha = cpar_[2])
    tau2 <- pcondig12(qu, .v, theta = cpar_[1], alpha = cpar_[2])
    expect_equal(.u, tau2)
  }
})


test_that("qcondigl12 is the inverse of pcondigl12", {
  for (alpha_ in .alpha){
    qu <- qcondigl12(.u, .v, alpha = alpha_)
    tau2 <- pcondigl12(qu, .v, alpha = alpha_)
    expect_equal(.u, tau2)
  }
})
