#' @rdname interpolator
#' @export
interp_kappa <- function(t, eta, k) {
  igl_kappa(1 / (eta * log(t)), k) / t
}


#' @rdname interpolator
#' @export
interp_kappa_inv_uniroot_algo <- function(p, eta, k, mxiter = 80, eps = 1.e-12, bd = 5) {
  if (length(p) != 1L) stop("Algorithm requires a single `p`.")
  if (length(eta) != 1L) stop("Algorithm requires a single `eta`.")
  if (length(k) != 1L) stop("Algorithm requires a single `k`.")
  if (p == 0) return(Inf)
  if (p == 1) return(1)
  upper <- 1 / p
  lower <- 1
  f <- function(t) interp_kappa(t, eta = eta, k = k) - p
  fit <- uniroot(f, c(lower, upper))
  cat(fit$code)
  fit$estimate
}

#' @rdname interpolator
#' @export
interp_kappa_inv_uniroot <- function(p, eta, k, mxiter = 80, eps = 1.e-12, bd = 5) {
  lengths <- c(p = length(p), eta = length(eta), k = length(k))
  l <- max(lengths)
  if (lengths[["p"]] == 1) p <- rep(p, l)
  if (lengths[["eta"]] == 1) eta <- rep(eta, l)
  if (lengths[["k"]] == 1) k <- rep(k, l)
  sol <- numeric()
  for (i in 1:l) {
    sol[i] <- interp_kappa_inv_uniroot_algo(
      p[i], eta = eta[i], k = k[i], mxiter = mxiter, eps = eps, bd = bd
    )
  }
  sol
}


#' @rdname interpolator
#' @export
interp_kappa_inv_nr_algo <- function(p, eta, k, mxiter = 80, eps = 1.e-12, bd = 5) {
  if (length(p) != 1L) stop("Algorithm requires a single `p`.")
  if (length(eta) != 1L) stop("Algorithm requires a single `eta`.")
  if (length(k) != 1L) stop("Algorithm requires a single `k`.")
  if (p == 0) return(Inf)
  if (p == 1) return(1)
  ## Starting value: use H_k^{-1}
  x <- interp_gen_inv(p, eta, k, mxiter = mxiter, eps = eps, bd = bd)
  x <- pmax(x - eps, 1 + (x - 1) / 2) # x-eps might overshoot left of 1.
  x <- pmax(1 + eps, x)  # x might *be* 1.
  ## Begin Newton-Raphson algorithm
  iter <- 0
  diff <- 1
  while(iter < mxiter & abs(diff) > eps) {
    eta_recip <- 1 / eta
    logx_recip <- 1 / log(x)
    eta_logx_recip <- eta_recip * logx_recip
    kappa <- igl_kappa(eta_logx_recip, k = k)
    kappap <- igl_kappa_D(eta_logx_recip, k = k)
    g <- kappa - x * p
    gp <- - eta_recip / x * logx_recip ^ 2 * kappap - p
    diff <- g / gp
    if (diff > x - 1) diff <- (x - 1) / 2
    x <- x - diff
    while(abs(diff) > bd || x < 1) {
      diff <- diff / 2
      x <- x + diff
    }
    iter <- iter + 1
  }
  x
}


#' @rdname interpolator
#' @export
interp_kappa_inv_nr <- function(p, eta, k, mxiter = 80, eps = 1.e-12, bd = 5) {
  lengths <- c(p = length(p), eta = length(eta), k = length(k))
  l <- max(lengths)
  if (lengths[["p"]] == 1) p <- rep(p, l)
  if (lengths[["eta"]] == 1) eta <- rep(eta, l)
  if (lengths[["k"]] == 1) k <- rep(k, l)
  sol <- numeric()
  for (i in 1:l) {
    sol[i] <- interp_kappa_inv_nr_algo(
      p[i], eta = eta[i], k = k[i], mxiter = mxiter, eps = eps, bd = bd
    )
  }
  sol
}

#' @rdname interpolator
#' @export
interp_kappa_inv <- interp_kappa_inv_nr
