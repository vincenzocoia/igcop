cpar <- c(12.3, 4.3, 2.3, 0.3, 0.1, 3.1, 7.1, 14.1, 34.1, 54.1, 44.1)
u <- c(10^(-(5:1)), 2:8/10, 1 - 10^(-(1:5)))
set.seed(5)
v <- sample(u)

# Density of IGL copula with N(0,1) margins
diglcop_gaussian <- function(u, v, cpar) {
  diglcop(u, v, cpar) * dnorm(qnorm(u)) * dnorm(qnorm(v))
}

test_that("numerical density from cdf", {
  #' Function to calculate numerical derivative
  diglcop_gaussian_numerical <- function(u, v, cpar, eps = 1.e-5) {
    x <- qnorm(u)
    y <- qnorm(v)
    cdf11 <- piglcop(u, v, cpar)
    cdf22 <- piglcop(pnorm(x + eps), pnorm(y + eps), cpar)
    cdf21 <- piglcop(pnorm(x + eps), v, cpar)
    cdf12 <- piglcop(u, pnorm(y + eps), cpar)
    (cdf22 + cdf11 - cdf12 - cdf21) / eps ^ 2
  }
  for (cpar_ in cpar) {
    pdf1 <- diglcop_gaussian_numerical(u, v, cpar_)
    pdf2 <- diglcop_gaussian(u, v, cpar_)
    expect_equal(pdf1, pdf2, tolerance = 1e-4)
  }
})

test_that("numerical 2|1 from cdf", {
  #' Calculate numerical derivative
  pcondiglcop21_numerical <- function(v, u, cpar, eps = 1.e-8) {
    cdf11 <- piglcop(u, v, cpar)
    cdf21 <- piglcop(u + eps, v, cpar)
    (cdf21 - cdf11) / eps
  }
  for (cpar_ in cpar) {
    pcond1 <- pcondiglcop21_numerical(v, u, cpar_)
    pcond2 <- pcondiglcop21(v, u, cpar_)
    expect_equal(pcond1, pcond2, tolerance = 1e-6)
  }
})

test_that("numerical density from 2|1", {
  #' Calculate numerical derivative
  diglcop_gaussian_numerical_from_2g1 <- function(u, v, cpar, eps = 1.e-8) {
    x <- qnorm(u)
    y <- qnorm(v)
    conda <- pcondiglcop21(v, u, cpar)
    condb <- pcondiglcop21(pnorm(y + eps), u, cpar)
    (condb - conda) / eps * dnorm(x)
  }
  for (cpar_ in cpar) {
    pdf1 <- diglcop_gaussian_numerical_from_2g1(u, v, cpar_)
    pdf2 <- diglcop_gaussian(u, v, cpar_)
    plot(pdf1, pdf2)
    expect_equal(pdf1, pdf2, tolerance = 1e-7)
  }
})
