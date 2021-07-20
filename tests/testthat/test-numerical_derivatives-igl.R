# Density of IGL copula with N(0,1) margins
digl_gaussian <- function(u, v, alpha) {
  digl(u, v, alpha = alpha) *
    dnorm(qnorm(u)) *
    dnorm(qnorm(v))
}

test_that("density matches the numerical density obtained from the cdf", {
  uv <- expand.grid(u = 1:4 / 5, v = 1:4 / 5)
  .u <- uv$u
  .v <- uv$v
  #' Function to calculate numerical derivative
  digl_gaussian_numerical <- function(u, v, alpha, eps = 1.e-5) {
    x <- qnorm(u)
    y <- qnorm(v)
    cdf11 <- pigl(u, v, alpha = alpha)
    cdf22 <- pigl(pnorm(x + eps), pnorm(y + eps), alpha = alpha)
    cdf21 <- pigl(pnorm(x + eps), v, alpha = alpha)
    cdf12 <- pigl(u, pnorm(y + eps), alpha = alpha)
    (cdf22 + cdf11 - cdf12 - cdf21) / eps ^ 2
  }
  for (alpha_ in .alpha) {
    print(alpha_)
    pdf1 <- digl_gaussian_numerical(.u, .v, alpha = alpha_)
    pdf2 <- digl_gaussian(.u, .v, alpha = alpha_)
    expect_equal(pdf1, pdf2, tolerance = 1e-4)
  }
})

test_that("the 2|1 cdf matches the numerically obtained cdf", {
  uv <- expand.grid(u = 1:4 / 5, v = 1:4 / 5)
  .u <- uv$u
  .v <- uv$v
  #' Calculate numerical derivative
  pcondigl21_numerical <- function(v, u, alpha, eps = 1.e-8) {
    cdf11 <- pigl(u, v, alpha = alpha)
    cdf21 <- pigl(u + eps, v, alpha = alpha)
    (cdf21 - cdf11) / eps
  }
  for (alpha_ in .alpha) {
    pcond1 <- pcondigl21_numerical(.v, .u, alpha = alpha_)
    pcond2 <- pcondigl21(.v, .u, alpha = alpha_)
    expect_equal(pcond1, pcond2, tolerance = 1e-5)
  }
})

test_that("the 1|2 cdf matches the numerically obtained cdf", {
  uv <- expand.grid(u = 1:4 / 5, v = 1:4 / 5)
  .u <- uv$u
  .v <- uv$v
  #' Calculate numerical derivative
  pcondigl12_log_numerical <- function(u, v, alpha, eps = 1.e-9) {
    x <- -log(u)
    y <- -log(v)
    cdf11 <- pigl(u, v, alpha = alpha)
    cdf12 <- pigl(u, exp(-y - eps), alpha = alpha)
    (cdf12 - cdf11) / eps
  }
  for (alpha_ in .alpha) {
    pcond1 <- pcondigl12_log_numerical(.u, .v, alpha = alpha_)
    pcond2 <- pcondigl12(.u, .v, alpha = alpha_) * (-.v)
    expect_equal(pcond1, pcond2, tolerance = 1e-5)
  }
})


test_that("density matches the numerical density obtained from the 2|1 conditional distribution", {
  uv <- expand.grid(u = 1:4 / 5, v = 1:4 / 5)
  .u <- uv$u
  .v <- uv$v
  #' Calculate numerical derivative
  digl_numerical_from_2g1 <- function(u, v, alpha, eps = 1.e-10) {
    conda <- pcondigl21(v, u, alpha = alpha)
    condb <- pcondigl21(v + eps, u, alpha = alpha)
    (condb - conda) / eps
  }
  for (alpha_ in .alpha) {
    pdf1 <- digl_numerical_from_2g1(.u, .v, alpha = alpha_)
    pdf2 <- digl(.u, .v, alpha = alpha_)
    expect_equal(pdf1, pdf2, tolerance = 1e-5)
  }
})

rm("digl_gaussian")
