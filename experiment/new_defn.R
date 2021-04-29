psi_new <- function(x, k) {
  res <- (k - 1) / x * stats::pgamma(x, k) +
    stats::pgamma(x, k - 1, lower.tail = FALSE)
  res[x == 0] <- 1
  res
}

psi_new_D <- function(x, k) {
  - (k - 1) * stats::pgamma(x, k) / x ^ 2
}

H_new <- function(x, eta, k) {
  exp(-x) * psi_new(eta * x, k)
}

H_kappa_new <- function(x, eta, k) {
  exp(-x) * kappa_new(eta * x, k)
}

kappa_new <- function(x, k) {
  stats::pgamma(x, shape = k - 1, lower.tail = FALSE)
}

kappa_new_D <- function(x, k) {
  -stats::dgamma(x, shape = k - 1)
}


H_new_inv_algo <- function(p, eta, k, mxiter = 40, eps = 1.e-12, bd = 5) {
  if (length(p) != 1L) stop("Algorithm requires a single `p`.")
  if (length(eta) != 1L) stop("Algorithm requires a single `eta`.")
  if (length(k) != 1L) stop("Algorithm requires a single `k`.")
  # if (p == 0) return(Inf)
  # if (p == 1) return(1)
  ## Algorithm:
  ## Get starting values
  x <- -log(p)
  # init2 <- exp(1 / eta / igl_gen_inv(p, k))
  # x <- pmin(init1, init2)
  # x <- pmax(x - eps, 1 + (x - 1) / 2) # x-eps might overshoot left of 1.
  # x <- pmax(1 + eps, x)  # x might be 1.
  iter <- 0
  diff <- 1
  ## Begin Newton-Raphson algorithm
  while(iter < mxiter & abs(diff) > eps) {
    pex <- p * exp(x)
    g <- psi_new(eta * x, k) - pex
    gp <- eta * psi_new_D(eta * x, k) - pex
    diff <- g / gp
    if (x - diff <= 0) diff <- x / 2
    x <- x - diff
    while (abs(diff) > bd) {
      diff <- diff / 2
      x <- x + diff
    }
    iter <- iter + 1
  }
  # print(iter)
  x
}

#' @rdname interpolator
#' @export
H_new_inv <- function(p, eta, k, mxiter = 40, eps = 1.e-12, bd = 5) {
  lengths <- c(p = length(p), eta = length(eta), k = length(k))
  l <- max(lengths)
  if (lengths[["p"]] == 1) p <- rep(p, l)
  if (lengths[["eta"]] == 1) eta <- rep(eta, l)
  if (lengths[["k"]] == 1) k <- rep(k, l)
  x <- numeric()
  for (i in 1:l) {
    x[i] <- H_new_inv_algo(
      p[i], eta = eta[i], k = k[i], mxiter = mxiter, eps = eps, bd = bd
    )
  }
  x
}

digcop_new <- function(u, v, cpar) {
  theta <- cpar[1]
  k <- cpar[2]
  y <- H_new_inv(1 - v, theta, k)
  eta <- theta * (1 - u)
  num <- kappa_new(eta * y, k) - eta * kappa_new_D(eta * y, k)
  den <- psi_new(theta * y, k) - theta * psi_new_D(theta * y, k)
  num / den
}

pigcop_new <- function(u, v, cpar) {
  theta <- cpar[1]
  k <- cpar[2]
  y <- H_new_inv(1 - v, theta, k)
  u + v - 1 + (1 - u) * H_new(y, theta * (1 - u), k)
}

# expand_grid(u = 1:9/10, v = 1:9/10) %>%
#   mutate(p = pigcop_new(u, v, c(3.4, 11.1))) %>%
#   pivot_wider(id_cols = u, names_from = v, values_from = p) %>%
#   select(-u) %>%
#   as.matrix() %>%
#   plotly::plot_ly()
