#' @rdname interpolator
interp_kappa <- function(x, eta, alpha) {
  exp(-x) * igl_kappa(eta * x, alpha)
}

#' @rdname interpolator
interp_kappa_D1 <- function(x, eta, alpha) {
  -exp(-x) * (igl_kappa(eta * x, alpha) - eta * igl_kappa_D(eta * x, alpha))
}


#' @rdname interpolator
interp_kappa_inv_algo <- function(p, eta, alpha, mxiter = 80, eps = 1.e-12, bd = 5) {
  if (length(p) != 1L) stop("Algorithm requires a single `p`.")
  if (length(eta) != 1L) stop("Algorithm requires a single `eta`.")
  if (length(alpha) != 1L) stop("Algorithm requires a single `alpha`.")
  prod <- alpha * eta * p
  if (is.na(prod)) return(prod)
  if (p == 0) return(Inf)
  if (p == 1) return(0)
  x1 <- -log(p)
  x2 <- igl_kappa_inv(p, alpha = alpha) / eta
  p1 <- interp_kappa(x1, eta = eta, alpha = alpha)
  p2 <- interp_kappa(x2, eta = eta, alpha = alpha)
  if (abs(p1 - p) < abs(p2 - p)) {
    x <- x1
  } else {
    x <- x2
  }
  iter <- 0
  diff <- 1
  while(iter < mxiter & abs(diff) > eps) {
    pex <- p * exp(x)
    g <- igl_kappa(eta * x, alpha) - pex
    gp <- eta * igl_kappa_D(eta * x, alpha) - pex
    diff <- g / gp
    if (x - diff < 0) diff <- x / 2
    x <- x - diff
    while(abs(diff) > bd) {
      diff <- diff / 2
      x <- x + diff
    }
    iter <- iter + 1
  }
  x
}

#' @rdname interpolator
interp_kappa_inv <- function(p, eta, alpha, mxiter = 80, eps = 1.e-12, bd = 5) {
  l <- vctrs::vec_size_common(p, eta, alpha)
  if (l == 0L) return(numeric(0L))
  args <- vctrs::vec_recycle_common(p = p, eta = eta, alpha = alpha)
  with(args, {
    x <- numeric(0L)
    for (i in 1:l) {
      x[i] <- interp_kappa_inv_algo(
        p[i], eta = eta[i], alpha = alpha[i], mxiter = mxiter, eps = eps, bd = bd
      )
    }
    x
  })
}
