#' @rdname interpolator
#' @export
interp_kappa <- function(x, eta, k) {
  igl_kappa(1 / (eta * log(x)), k) / x
}

#' @rdname interpolator
#' @export
interp_kappa_D1 <- function(x, eta, k) {
  logt <- log(x)
  arg <- 1 / eta / logt
  coeff <- arg / logt
  -x ^ (-2) * (igl_kappa(arg, k) + coeff * igl_kappa_D(arg, k))
}


#' @rdname interpolator
#' @export
interp_kappa_inv_algo <- function(p, eta, k, mxiter = 80, eps = 1.e-12, bd = 5) {
  if (length(p) != 1L) stop("Algorithm requires a single `p`.")
  if (length(eta) != 1L) stop("Algorithm requires a single `eta`.")
  if (length(k) != 1L) stop("Algorithm requires a single `k`.")
  if (p == 0) return(Inf)
  if (p == 1) return(1)
  ## Starting value:
  x1 <- exp(1 / igl_kappa_inv(p, k) / eta)
  x2 <- 1 / p
  p1 <- interp_kappa(x1, eta, k)
  p2 <- interp_kappa(x2, eta, k)
  if (abs(p1 - p) < abs(p2 - p)) {
    x <- x1
  } else {
    x <- x2
  }
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
interp_kappa_inv <- function(p, eta, k, mxiter = 80, eps = 1.e-12, bd = 5) {
  lengths <- c(p = length(p), eta = length(eta), k = length(k))
  l <- max(lengths)
  if (lengths[["p"]] == 1) p <- rep(p, l)
  if (lengths[["eta"]] == 1) eta <- rep(eta, l)
  if (lengths[["k"]] == 1) k <- rep(k, l)
  x <- numeric()
  for (i in 1:l) {
    x[i] <- interp_kappa_inv_algo(
      p[i], eta = eta[i], k = k[i], mxiter = mxiter, eps = eps, bd = bd
    )
  }
  x
}
