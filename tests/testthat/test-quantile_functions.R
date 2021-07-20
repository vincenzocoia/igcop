context("Check that the quantile functions are the inverses of the
        conditional distributions.")

test_that("qcondig21 is the inverse of pcondig21", {
  for (cpar_ in .cpar){
    i <- 7L
    print(cpar_)
    qu <- qcondig(.v[i], .u[i], theta = cpar_[1], alpha = cpar_[2])
    print(paste("qu =", qu))
    .v2 <- pcondig(qu, .u[i], theta = cpar_[1], alpha = cpar_[2])
    print(paste(".v2 =", .v2))
    expect_equal(.v[i], .v2)
  }
})


test_that("qcondigl21 is the inverse of pcondigl21", {
  for (alpha_ in .alpha){
    v2 <- pcondigl(qcondigl(.v, .u, alpha = alpha_), .u, alpha = alpha_)
    expect_equal(.v, v2)
  }
})

test_that("qcondig12 is the inverse of pcondig12", {
  for (cpar_ in .cpar){
    qu <- qcondig12(.u, .v, theta = cpar_[1], alpha = cpar_[2])
    u2 <- pcondig12(qu, .v, theta = cpar_[1], alpha = cpar_[2])
    expect_equal(.u, u2)
  }
})


test_that("qcondigl12 is the inverse of pcondigl12", {
  for (alpha_ in .alpha){
    qu <- qcondigl12(.u, .v, alpha = alpha_)
    u2 <- pcondigl12(qu, .v, alpha = alpha_)
    expect_equal(.u, u2)
  }
})
