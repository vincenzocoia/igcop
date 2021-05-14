context("Check derivatives of IG copula numerically")

# Density of IG copula with N(0,1) margins
digcop_gaussian <- function(u, v, cpar) {
  digcop(u, v, cpar) * dnorm(qnorm(u)) * dnorm(qnorm(v))
}

test_that("density matches the numerical density obtained from the cdf", {
  #' Function to calculate numerical derivative
  digcop_gaussian_numerical <- function(u, v, cpar, eps = 1.e-5) {
    x <- qnorm(u)
    y <- qnorm(v)
    cdf11 <- pigcop(u, v, cpar)
    cdf22 <- pigcop(pnorm(x + eps), pnorm(y + eps), cpar)
    cdf21 <- pigcop(pnorm(x + eps), v, cpar)
    cdf12 <- pigcop(u, pnorm(y + eps), cpar)
    (cdf22 + cdf11 - cdf12 - cdf21) / eps ^ 2
  }
  for (cpar_ in .cpar) {
    pdf1 <- digcop_gaussian_numerical(.u, .v, cpar_)
    pdf2 <- digcop_gaussian(.u, .v, cpar_)
    expect_equal(pdf1, pdf2, tolerance = 1e-4)
  }
})

test_that("the 2|1 cdf matches the numerically obtained cdf", {
  #' Calculate numerical derivative
  pcondigcop21_numerical <- function(v, u, cpar, eps = 1.e-8) {
    cdf11 <- pigcop(u, v, cpar)
    cdf21 <- pigcop(u + eps, v, cpar)
    (cdf21 - cdf11) / eps
  }
  for (cpar_ in .cpar) {
    pcond1 <- pcondigcop21_numerical(.v, .u, cpar_)
    pcond2 <- pcondigcop21(.v, .u, cpar_)
    expect_equal(pcond1, pcond2, tolerance = 1e-6)
  }
})

test_that("the 1|2 cdf matches the numerically obtained cdf", {
  #' Calculate numerical derivative
  pcondigcop12_numerical <- function(u, v, cpar, eps = 1.e-8) {
    cdf11 <- pigcop(u, v, cpar)
    cdf12 <- pigcop(u, v + eps, cpar)
    (cdf12 - cdf11) / eps
  }
  for (cpar_ in .cpar) {
    pcond1 <- pcondigcop12_numerical(.u, .v, cpar_)
    pcond2 <- pcondigcop12(.u, .v, cpar_)
    expect_equal(pcond1, pcond2, tolerance = 1e-6)
  }
})

test_that("density matches the numerical density obtained from the 2|1 conditional distribution", {
  #' Calculate numerical derivative
  digcop_gaussian_numerical_from_2g1 <- function(u, v, cpar, eps = 1.e-8) {
    x <- qnorm(u)
    y <- qnorm(v)
    conda <- pcondigcop21(v, u, cpar)
    condb <- pcondigcop21(pnorm(y + eps), u, cpar)
    (condb - conda) / eps * dnorm(x)
  }
  for (cpar_ in .cpar) {
    pdf1 <- digcop_gaussian_numerical_from_2g1(.u, .v, cpar_)
    pdf2 <- digcop_gaussian(.u, .v, cpar_)
    expect_equal(pdf1, pdf2, tolerance = 1e-7)
  }
})

rm("digcop_gaussian")
