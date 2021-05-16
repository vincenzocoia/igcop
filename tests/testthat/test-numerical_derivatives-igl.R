# Density of IGL copula with N(0,1) margins
diglcop_gaussian <- function(u, v, cpar) {
  diglcop(u, v, cpar = cpar) * dnorm(qnorm(u)) * dnorm(qnorm(v))
}

test_that("density matches the numerical density obtained from the cdf", {
  #' Function to calculate numerical derivative
  diglcop_gaussian_numerical <- function(u, v, cpar, eps = 1.e-5) {
    x <- qnorm(u)
    y <- qnorm(v)
    cdf11 <- piglcop(u, v, cpar = cpar)
    cdf22 <- piglcop(pnorm(x + eps), pnorm(y + eps), cpar = cpar)
    cdf21 <- piglcop(pnorm(x + eps), v, cpar = cpar)
    cdf12 <- piglcop(u, pnorm(y + eps), cpar = cpar)
    (cdf22 + cdf11 - cdf12 - cdf21) / eps ^ 2
  }
  for (alpha_ in .alpha) {
    pdf1 <- diglcop_gaussian_numerical(.u, .v, cpar = alpha_)
    pdf2 <- diglcop_gaussian(.u, .v, cpar = alpha_)
    expect_equal(pdf1, pdf2, tolerance = 1e-4)
  }
})

test_that("the 2|1 cdf matches the numerically obtained cdf", {
  #' Calculate numerical derivative
  pcondiglcop21_numerical <- function(v, u, cpar, eps = 1.e-8) {
    cdf11 <- piglcop(u, v, cpar = cpar)
    cdf21 <- piglcop(u + eps, v, cpar = cpar)
    (cdf21 - cdf11) / eps
  }
  for (alpha_ in .alpha) {
    pcond1 <- pcondiglcop21_numerical(.v, .u, alpha_)
    pcond2 <- pcondiglcop21(.v, .u, alpha_)
    expect_equal(pcond1, pcond2, tolerance = 1e-5)
  }
})

test_that("the 1|2 cdf matches the numerically obtained cdf", {
  #' Calculate numerical derivative
  pcondiglcop12_log_numerical <- function(u, v, cpar, eps = 1.e-9) {
    x <- -log(u)
    y <- -log(v)
    cdf11 <- piglcop(u, v, cpar = cpar)
    cdf12 <- piglcop(u, exp(-y - eps), cpar = cpar)
    (cdf12 - cdf11) / eps
  }
  for (alpha_ in .alpha) {
    pcond1 <- pcondiglcop12_log_numerical(.u, .v, cpar = alpha_)
    pcond2 <- pcondiglcop12(.u, .v, cpar = alpha_) * (-.v)
    expect_equal(pcond1, pcond2, tolerance = 1e-5)
  }
})


test_that("density matches the numerical density obtained from the 2|1 conditional distribution", {
  #' Calculate numerical derivative
  diglcop_numerical_from_2g1 <- function(u, v, cpar, eps = 1.e-10) {
    conda <- pcondiglcop21(v, u, cpar = cpar)
    condb <- pcondiglcop21(v + eps, u, cpar = cpar)
    (condb - conda) / eps
  }
  for (alpha_ in .alpha) {
    pdf1 <- diglcop_numerical_from_2g1(.u, .v, cpar = alpha_)
    pdf2 <- diglcop(.u, .v, cpar = alpha_)
    expect_equal(pdf1, pdf2, tolerance = 1e-5)
  }
})

rm("diglcop_gaussian")
