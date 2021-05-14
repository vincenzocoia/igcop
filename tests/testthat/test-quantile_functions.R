context("Check that the quantile functions are the inverses of the conditional distributions.")

test_that("qcondigcop21 is the inverse of pcondigcop21", {
  for (cpar_ in .cpar){
    tau2 <- pcondigcop(qcondigcop(.v, .u, cpar = cpar_), .u, cpar = cpar_)
    expect_equal(.v, tau2)
  }
})


test_that("qcondiglcop21 is the inverse of pcondiglcop21", {
  for (alpha_ in .alpha){
    tau2 <- pcondiglcop(qcondiglcop(.v, .u, cpar = alpha_), .u, cpar = alpha_)
    expect_equal(.v, tau2)
  }
})

test_that("qcondigcop12 is the inverse of pcondigcop12", {
  for (cpar_ in .cpar){
    tau2 <- pcondigcop12(qcondigcop12(.u, .v, cpar = cpar_), .v, cpar = cpar_)
    expect_equal(.u, tau2)
  }
})


test_that("qcondiglcop12 is the inverse of pcondiglcop12", {
  for (alpha_ in .alpha){
    tau2 <- pcondiglcop12(qcondiglcop12(.u, .v, cpar = alpha_), .v, cpar = alpha_)
    expect_equal(.u, tau2)
  }
})
