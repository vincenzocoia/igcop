#' @rdname interpolator
#' @export
interp_kappa <- function(t, eta, k) {
  igl_kappa(1 / (eta * log(t)), k) / t
}



#' @param p Vector of values in [0,1] to evaluate the inverse function at.
#' @rdname interpolator
#' @export
interp_kappa_inv_uniroot <- function(p, eta, k, mxiter = 80, eps = 1.e-12, bd = 5) {
  sol <- numeric()
  at_2 <- interp_kappa(2, eta = eta, k = k)
  expit <- function(x) 1 / (exp(-x) + 1)
  logit <- function(x) log(x / (1 - x))
  ## Get starting values
  xp1 <- 1 / p
  xp2 <- exp(stats::qgamma(1 - p, k - 1) / eta)
  xpm <- pmin(xp1, xp2)
  upper <- pmax(xp1, xp2)
  for (i in seq_along(p)) {
    this_p <- p[i]
    x0 <- xpm[i]
    this_upper <- upper[i]


    # if (this_p > at_2) {
    #   ## Solution is between 1 and 2, but function may be steep here.
    #   ## logit/expit transform. Tolerance is respected since max slope = 0.25
    #   f2 <- function(t) (interp_kappa(expit(t) + 1, eta = eta, k = k) - this_p) ^ 2
    #   cat("L")
    #   res <- nlm(f2, logit(1.5 - 1), steptol = 1e-8)
    #   cat(res$code)
    #   sol[i] <- expit(res$estimate) + 1
    # } else {
    #   f2 <- function(t) (interp_kappa(t, eta = eta, k = k) - this_p) ^ 2
    #   cat("U")
    #   res <- nlm(f2, x0, steptol = 1e-8)
    #   cat(res$code)
    #   sol[i] <- res$estimate
    # }


    f2 <- function(t) {
      if (t <= 1) {
        return(log(this_p) ^ 2 + 1000)
      } else {
        abs(log(interp_kappa(t, eta = eta, k = k)) - log(this_p))
      }
    }
    res <- nlm(f2, x0, steptol = 1e-12)
    cat(res$code)
    sol[i] <- res$estimate
  }
  sol
}


# Instead of finding the t where H(t) = p, find s, where s = log(t - 1). Means t = exp(s) + 1
# Function: H_new(s) = H(exp(s) + 1) = 1 / (exp(s) + 1) * kappa(1 / (eta log(exp(s) + 1))), set = p
# ==> kappa(1 / (eta log(exp(s) + 1))) = p (exp(s) + 1)
# ==> p (exp(s) + 1) - kappa(1 / (eta log(exp(s) + 1))) = 0. Define g(s) as the LHS.
# g'(s) = p * exp(s) - kappa'(same arg) (-1) (same arg) ^ (-2) eta / (exp(s) + 1) * exp(s)
#       = p * exp(s) + kappa'(same arg) 1 / (eta log(exp(s) + 1) ^ 2)  / (exp(s) + 1) * exp(s)
#' @param p Vector of values in [0,1] to evaluate the inverse function at.
#' @rdname interpolator
#' @export
interp_kappa_inv <- function(p, eta, k, mxiter = 80, eps = 1.e-12, bd = 5) {
  sol <- numeric()
  for (i in seq_along(p)) {
    p_ <- p[i]
    ## Algorithm:
    ## Get starting values
    xp1 <- 1 / p_
    xp2 <- exp(stats::qgamma(1 - p_, k - 1) / eta)
    xpm <- pmin(xp1, xp2)
    # print(mean(xpm == xp1))
    tt <- pmax(xpm - eps, 1 + (xpm - 1) / 2) # xpm-eps might overshoot left of 1.
    tt <- pmax(1 + eps, tt)  # tt might be 1.
    # shift to new scale:
    ss <- log(tt - 1)
    iter <- 0
    diff <- 1
    ## Begin Newton-Raphson algorithm
    while(iter < mxiter & max(abs(diff)) > eps) {
      ## Helpful quantities
      es <- exp(ss)
      esp1 <- es + 1
      logesp1 <- log(esp1)
      etalogesp1inv <- 1 / (eta * logesp1)

      logt <- log(tt)
      etalog <- eta * logt
      etaloginv <- 1 / etalog
      ## Evaluate functions
      g <- p_ * esp1 - igl_kappa(etalogesp1inv, k)
      gp <- p_ * es + igl_kappa_D(etalogesp1inv, k) * es * etalogesp1inv / logesp1  / esp1

      diff <- g / gp
      ss <- ss - diff
      while(max(abs(diff)) > bd) {
        diff <- diff / 2
        ss <- ss + diff
      }
      iter <- iter + 1
    }
    # print(iter)
    sol[i] <- tt
  }
  exp(sol) + 1
}



#' @param p Vector of values in [0,1] to evaluate the inverse function at.
#' @rdname interpolator
#' @export
interp_kappa_inv_old <- function(p, eta, k, mxiter = 80, eps = 1.e-12, bd = 5) {
  sol <- numeric()
  for (i in seq_along(p)) {
    p_ <- p[i]
    ## Algorithm:
    ## Get starting values
    xp1 <- 1 / p_
    xp2 <- exp(stats::qgamma(1 - p_, k - 1) / eta)
    xpm <- pmin(xp1, xp2)
    # print(mean(xpm == xp1))
    tt <- pmax(xpm - eps, 1 + (xpm - 1) / 2) # xpm-eps might overshoot left of 1.
    tt <- pmax(1 + eps, tt)  # tt might be 1.
    iter <- 0
    diff <- 1
    ## Begin Newton-Raphson algorithm
    while(iter < mxiter & max(abs(diff)) > eps) {
      ## Helpful quantities
      logt <- log(tt)
      etalog <- eta * logt
      etaloginv <- 1 / etalog
      ## Evaluate functions
      g <- tt * p_ - igl_kappa(etaloginv, k)
      gp <- p_ + etaloginv / tt / logt * igl_kappa_D(etaloginv, k)
      diff <- g / gp
      flag <- diff > tt - 1
      diff[flag] <- (tt[flag] - 1) / 2
      tt <- tt - diff
      while(max(abs(diff)) > bd | any(tt <= 1)) {
        diff <- diff / 2
        tt <- tt + diff
      }
      iter <- iter + 1
      # cat(iter, ",", tt, "\n")
    }
    # print(iter)
    sol[i] <- tt
  }
  sol
}

# 1/t K(1 / (eta * log(t)))
#
#
#
# g1 <- function(x, eta, k, p) { # Support (1, Inf)
#   logt <- log(x)
#   etalog <- eta * logt
#   etaloginv <- 1 / etalog
#   x * p - igl_kappa(etaloginv, k)
# }
# g1_ <- function(x) g1(x, eta = 2, k = 2, p = 0.2)
#
# g2 <- function(x, eta, k, p) { # Support (0, Inf) via log
#   x <- exp(x)
#   logt <- log(x)
#   etalog <- eta * logt
#   etaloginv <- 1 / etalog
#   x * p - igl_kappa(etaloginv, k)
# }
# g2_ <- function(x) g2(x, eta = 2, k = 2, p = 0.2)
#
# g3 <- function(x, eta, k, p) { # Support (-Inf, Inf) via shift then log
#   x <- exp(x) + 1
#   logt <- log(x)
#   etalog <- eta * logt
#   etaloginv <- 1 / etalog
#   x * p - igl_kappa(etaloginv, k)
# }
# g3_ <- function(x) g3(x, eta = 2, k = 2, p = 0.2)
#
# curve(g3_, -5, 2)
