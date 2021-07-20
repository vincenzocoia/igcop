#' Algorithm: inverse of interpolator of kappa function
#' @inheritParams interp_kappa_inv
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
