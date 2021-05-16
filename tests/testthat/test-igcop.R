# ## Check that the conditional IG copula distributions work.
# library(plyr)
# library(copsupp)
# library(testthat)
#
# ## Set up parameters
# theta <- c(10, 100) # c(0.1, 1, 10, 100)
# k <- c(2.1, 5, 50, 100) # c(1.1, 2.1, 5, 50, 100)
# par <- expand.grid(theta=theta, k=k)
#
# ## Generate data
# set.seed(123)
# u <- runif(1000)
# tau <- runif(1000)
#
# ## Check IG generator inverse
# d_ply(par, names(par), function(dfrow) {
#     print(dfrow)
#     theta_ <- dfrow$theta
#     k_ <- dfrow$k
#     taunew <- interp_gen(interp_gen_inv(tau, theta_, k_), theta_, k_)
#     expect_equal(tau, taunew)
# })
#
#
# ## Check dcop
# cdfmax <- pmin(u, tau)  # Comonotonicity
# cdfmin <- pmax(u + tau - 1, 0)  # Countermonotonicity
# d_ply(par, names(par), function(dfrow) {
#     print(dfrow)
#     theta_ <- dfrow$theta
#     k_ <- dfrow$k
#     cdf1 <- pig(u, tau, c(theta_, k_))
#     if (any(cdf1>1) | any(cdf1<0)) stop("pig not in [0,1].")
#     if (any(cdfmax - cdf1 < 0)) stop("pig exceeds comonotonicity.")
#     if (any(cdfmin - cdf1 > 0)) stop("pig doesn't always exceed countermonotonicity.")
# })
#
#
# d_ply(par, names(par), function(dfrow) {
#     print(dfrow)
#     theta_ <- dfrow$theta
#     k_ <- dfrow$k
#     v <- qcondig(tau, u, c(theta_, k_))
#     tau2 <- pcondig(v, u, c(theta_, k_))
#     # expect_that(tau, equals(tau2))
#     print(plot(tau, tau2, main=paste(theta_, k_)))
# })
