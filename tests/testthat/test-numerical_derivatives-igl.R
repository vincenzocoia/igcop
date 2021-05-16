# Density of IGL copula with N(0,1) margins
digl_gaussian <- function(u, v, cpar) {
  digl(u, v, cpar = cpar) * dnorm(qnorm(u)) * dnorm(qnorm(v))
}

test_that("density matches the numerical density obtained from the cdf", {
  #' Function to calculate numerical derivative
  digl_gaussian_numerical <- function(u, v, cpar, eps = 1.e-5) {
    x <- qnorm(u)
    y <- qnorm(v)
    cdf11 <- pigl(u, v, cpar = cpar)
    cdf22 <- pigl(pnorm(x + eps), pnorm(y + eps), cpar = cpar)
    cdf21 <- pigl(pnorm(x + eps), v, cpar = cpar)
    cdf12 <- pigl(u, pnorm(y + eps), cpar = cpar)
    (cdf22 + cdf11 - cdf12 - cdf21) / eps ^ 2
  }
  for (alpha_ in .alpha) {
    pdf1 <- digl_gaussian_numerical(.u, .v, cpar = alpha_)
    pdf2 <- digl_gaussian(.u, .v, cpar = alpha_)
    expect_equal(pdf1, pdf2, tolerance = 1e-4)
  }
})

test_that("the 2|1 cdf matches the numerically obtained cdf", {
  #' Calculate numerical derivative
  pcondigl21_numerical <- function(v, u, cpar, eps = 1.e-8) {
    cdf11 <- pigl(u, v, cpar = cpar)
    cdf21 <- pigl(u + eps, v, cpar = cpar)
    (cdf21 - cdf11) / eps
  }
  for (alpha_ in .alpha) {
    pcond1 <- pcondigl21_numerical(.v, .u, alpha_)
    pcond2 <- pcondigl21(.v, .u, alpha_)
    expect_equal(pcond1, pcond2, tolerance = 1e-5)
  }
})

test_that("the 1|2 cdf matches the numerically obtained cdf", {
  #' Calculate numerical derivative
  pcondigl12_log_numerical <- function(u, v, cpar, eps = 1.e-9) {
    x <- -log(u)
    y <- -log(v)
    cdf11 <- pigl(u, v, cpar = cpar)
    cdf12 <- pigl(u, exp(-y - eps), cpar = cpar)
    (cdf12 - cdf11) / eps
  }
  for (alpha_ in .alpha) {
    pcond1 <- pcondigl12_log_numerical(.u, .v, cpar = alpha_)
    pcond2 <- pcondigl12(.u, .v, cpar = alpha_) * (-.v)
    expect_equal(pcond1, pcond2, tolerance = 1e-5)
  }
})


test_that("density matches the numerical density obtained from the 2|1 conditional distribution", {
  #' Calculate numerical derivative
  digl_numerical_from_2g1 <- function(u, v, cpar, eps = 1.e-10) {
    conda <- pcondigl21(v, u, cpar = cpar)
    condb <- pcondigl21(v + eps, u, cpar = cpar)
    (condb - conda) / eps
  }
  for (alpha_ in .alpha) {
    pdf1 <- digl_numerical_from_2g1(.u, .v, cpar = alpha_)
    pdf2 <- digl(.u, .v, cpar = alpha_)
    expect_equal(pdf1, pdf2, tolerance = 1e-5)
  }
})

rm("digl_gaussian")
