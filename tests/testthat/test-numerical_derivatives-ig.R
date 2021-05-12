cpar <- list(
  c(2.5, 12.3), # nearly independent
  c(2.5, 4.3), # weak dependence spear=.040
  c(2.5, 2.3), # spear=.145
  c(2.5, 0.3), # spear=.244##
  c(2.5, 0.1), # unreliable (for small k)###
  c(5.5,3.1),  # spear=.231
  c(8.5,3.1),  # spear=.330
  c(15.5,3.1), # spear=.462
  c(25.5,3.1), # spear=.553
  c(45.5,3.1), # spear=.629
  c(75.5,3.1), # spear=.674
  c(75.5,7.1), # spear=.661
  c(155.5,3.1), # spear=.714
  c(300.5,3.1), # spear=.733
  c(600.5,3.1), # spear=.743
  c(600.5,7.1), # spear=.829
  c(600.5,14.1), # spear=.853
  c(600.5,34.1), # spear=.811
  c(1200.5,14.1), # spear=.883
  c(2400.5,14.1), # spear=.899
  c(2400.5,34.1), # spear=.918
  c(2400.5,54.1), # spear=.??? unstable for rhoS
  c(2400.5,44.1), # spear=.914
  c(4800.5,34.1), # spear=.937
  c(9600.5,34.1) # spear=.947  # try inference with this parameter vector
)
u <- c(10^(-(5:1)), 2:8/10, 1 - 10^(-(1:5)))
set.seed(5)
v <- sample(u)

# Density of IG copula with N(0,1) margins
digcop_gaussian <- function(u, v, cpar) {
  digcop(u, v, cpar) * dnorm(qnorm(u)) * dnorm(qnorm(v))
}

test_that("numerical density from cdf", {
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
  for (cpar_ in cpar) {
    pdf1 <- digcop_gaussian_numerical(u, v, cpar_)
    pdf2 <- digcop_gaussian(u, v, cpar_)
    expect_equal(pdf1, pdf2, tolerance = 1e-4)
  }
})

test_that("numerical 2|1 from cdf", {
  #' Calculate numerical derivative
  pcondigcop21_numerical <- function(v, u, cpar, eps = 1.e-8) {
    cdf11 <- pigcop(u, v, cpar)
    cdf21 <- pigcop(u + eps, v, cpar)
    (cdf21 - cdf11) / eps
  }
  for (cpar_ in cpar) {
    pcond1 <- pcondigcop21_numerical(v, u, cpar_)
    pcond2 <- pcondigcop21(v, u, cpar_)
    expect_equal(pcond1, pcond2, tolerance = 1e-6)
  }
})

test_that("numerical density from 2|1", {
  #' Calculate numerical derivative
  digcop_gaussian_numerical_from_2g1 <- function(u, v, cpar, eps = 1.e-8) {
    x <- qnorm(u)
    y <- qnorm(v)
    conda <- pcondigcop21(v, u, cpar)
    condb <- pcondigcop21(pnorm(y + eps), u, cpar)
    (condb - conda) / eps * dnorm(x)
  }
  for (cpar_ in cpar) {
    pdf1 <- digcop_gaussian_numerical_from_2g1(u, v, cpar_)
    pdf2 <- digcop_gaussian(u, v, cpar_)
    plot(pdf1, pdf2)
    expect_equal(pdf1, pdf2, tolerance = 1e-7)
  }
})
