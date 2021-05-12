#' @rdname interpolator
#' @export
interp_kappa <- function(x, eta, alpha) {
  exp(-x) * igl_kappa(eta * x, alpha)
}

#' @rdname interpolator
#' @export
interp_kappa_D1 <- function(x, eta, alpha) {
  -exp(-x) * (igl_kappa(eta * x, alpha) - eta * igl_kappa_D(eta * x, alpha))
}


#' @rdname interpolator
#' @export
interp_kappa_inv_algo <- function(p, eta, alpha, mxiter = 80, eps = 1.e-12, bd = 5) {
  if (length(p) != 1L) stop("Algorithm requires a single `p`.")
  if (length(eta) != 1L) stop("Algorithm requires a single `eta`.")
  if (length(alpha) != 1L) stop("Algorithm requires a single `alpha`.")
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

# tibble(x = seq(1, 2, length.out = 1000)) %>%
#   mutate(f_interp_kappa = interp_kappa(x, 420.35, 35.1),
#          f_kappa = igl_kappa(1/(420.35 * log(x)), 35.1),
#          f_interp_gen = interp_gen(x, 420.35, 35.1),
#          f_inv = 1/x) %>%
#   pivot_longer(cols = starts_with("f_"), names_to = "fun") %>%
#   ggplot(aes(x, value)) +
#   geom_line(aes(group = fun, colour = fun)) +
#   geom_hline(yintercept = 0.7, linetype = "dotted")


#' @rdname interpolator
#' @export
interp_kappa_inv <- function(p, eta, alpha, mxiter = 80, eps = 1.e-12, bd = 5) {
  lengths <- c(p = length(p), eta = length(eta), alpha = length(alpha))
  l <- max(lengths)
  if (lengths[["p"]] == 1) p <- rep(p, l)
  if (lengths[["eta"]] == 1) eta <- rep(eta, l)
  if (lengths[["alpha"]] == 1) alpha <- rep(alpha, l)
  x <- numeric()
  for (i in 1:l) {
    x[i] <- interp_kappa_inv_algo(
      p[i], eta = eta[i], alpha = alpha[i], mxiter = mxiter, eps = eps, bd = bd
    )
  }
  x
}
