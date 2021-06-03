#' Algorithm: inverse of interpolator of generating function
#' @inheritParams interp_gen_inv
interp_gen_inv_algo <- function(p, eta, alpha, mxiter = 40, eps = 1.e-12, bd = 5) {
  if (length(p) != 1L) stop("Algorithm requires a single `p`.")
  if (length(eta) != 1L) stop("Algorithm requires a single `eta`.")
  if (length(alpha) != 1L) stop("Algorithm requires a single `alpha`.")
  prod <- alpha * eta * p
  if (is.na(prod)) return(prod)
  if (p == 0) return(Inf)
  if (p == 1) return(0)
  x1 <- -log(p)
  x2 <- igl_gen_inv_algo(p, alpha) / eta
  p1 <- interp_gen(x1, eta = eta, alpha = alpha)
  p2 <- interp_gen(x2, eta = eta, alpha = alpha)
  if (abs(p - p1) < abs(p - p2)) {
    x <- x1
  } else {
    x <- x2
  }
  x <- log(x)
  iter <- 0
  diff <- 1
  while(iter < mxiter & abs(diff) > eps) {
    ex <- exp(x)
    g <- interp_gen(ex, eta = eta, alpha = alpha) - p
    gp <- interp_gen_D1(ex, eta = eta, alpha = alpha) * ex

    # pex <- p * exp(x)
    # g <- igl_gen(x * eta, alpha) - pex
    # gp <- eta * igl_gen_D(x * eta, alpha) - pex
    diff <- g / gp
    # if (x - diff < 0) diff <- x / 2
    x <- x - diff
    while (abs(diff) > bd) {
      diff <- diff / 2
      x <- x + diff
    }
    iter <- iter + 1
  }
  exp(x)
}
