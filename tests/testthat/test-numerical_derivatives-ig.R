context("Check derivatives of IG copula numerically")

# Density of IG copula with N(0,1) margins
dig_gaussian <- function(u, v, cpar) {
  dig(u, v, theta = cpar[1], alpha = cpar[2]) *
    dnorm(qnorm(u)) *
    dnorm(qnorm(v))
}

test_that("density matches the numerical density obtained from the cdf", {
  uv <- expand.grid(u = 1:4 / 5, v = 1:4 / 5)
  .u <- uv$u
  .v <- uv$v
  #' Function to calculate numerical derivative
  dig_gaussian_numerical <- function(u, v, cpar, eps = 1.e-5) {
    x <- qnorm(u)
    y <- qnorm(v)
    theta <- cpar[1]
    alpha <- cpar[2]
    cdf11 <- pig(u, v, theta = theta, alpha = alpha)
    cdf22 <- pig(pnorm(x + eps), pnorm(y + eps), theta = theta, alpha = alpha)
    cdf21 <- pig(pnorm(x + eps), v, theta = theta, alpha = alpha)
    cdf12 <- pig(u, pnorm(y + eps), theta = theta, alpha = alpha)
    (cdf22 + cdf11 - cdf12 - cdf21) / eps ^ 2
  }
  for (cpar_ in .cpar) {
    print(cpar_)
    pdf1 <- dig_gaussian_numerical(.u, .v, cpar = cpar_)
    pdf2 <- dig_gaussian(.u, .v, cpar = cpar_)
    expect_equal(pdf1, pdf2, tolerance = 1e-4)
  }
})

test_that("the 2|1 cdf matches the numerically obtained cdf", {
  uv <- expand.grid(u = 1:4 / 5, v = 1:4 / 5)
  .u <- uv$u
  .v <- uv$v
  #' Calculate numerical derivative
  pcondig21_numerical <- function(v, u, cpar, eps = 1.e-8) {
    theta <- cpar[1]
    alpha <- cpar[2]
    cdf11 <- pig(u, v, theta = theta, alpha = alpha)
    cdf21 <- pig(u + eps, v, theta = theta, alpha = alpha)
    (cdf21 - cdf11) / eps
  }
  for (cpar_ in .cpar) {
    pcond1 <- pcondig21_numerical(.v, .u, cpar = cpar_)
    pcond2 <- pcondig21(.v, .u, theta = cpar_[1], alpha = cpar_[2])
    expect_equal(pcond1, pcond2, tolerance = 1e-5)
  }
})

test_that("the 1|2 cdf matches the numerically obtained cdf", {
  uv <- expand.grid(u = 1:4 / 5, v = 1:4 / 5)
  .u <- uv$u
  .v <- uv$v
  #' Calculate numerical derivative
  pcondig12_numerical <- function(u, v, cpar, eps = 1.e-8) {
    theta <- cpar[1]
    alpha <- cpar[2]
    cdf11 <- pig(u, v, theta = theta, alpha = alpha)
    cdf12 <- pig(u, v + eps, theta = theta, alpha = alpha)
    (cdf12 - cdf11) / eps
  }
  for (cpar_ in .cpar) {
    pcond1 <- pcondig12_numerical(.u, .v, cpar = cpar_)
    pcond2 <- pcondig12(.u, .v, theta = cpar_[1], alpha = cpar_[2])
    expect_equal(pcond1, pcond2, tolerance = 1e-6)
  }
})

test_that("density matches the numerical density obtained from the 2|1 conditional distribution", {
  uv <- expand.grid(u = 1:4 / 5, v = 1:4 / 5)
  .u <- uv$u
  .v <- uv$v
  #' Calculate numerical derivative
  dig_gaussian_numerical_from_2g1 <- function(u, v, cpar, eps = 1.e-8) {
    theta <- cpar[1]
    alpha <- cpar[2]
    x <- qnorm(u)
    y <- qnorm(v)
    conda <- pcondig21(v, u, theta = theta, alpha = alpha)
    condb <- pcondig21(pnorm(y + eps), u, theta = theta, alpha = alpha)
    (condb - conda) / eps * dnorm(x)
  }
  for (cpar_ in .cpar) {
    pdf1 <- dig_gaussian_numerical_from_2g1(.u, .v, cpar = cpar_)
    pdf2 <- dig_gaussian(.u, .v, cpar = cpar_)
    expect_equal(pdf1, pdf2, tolerance = 1e-6)
  }
})

rm("dig_gaussian")
