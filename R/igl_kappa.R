#' Kappa version of generating function
#'
#' The kappa version of a generating function \eqn{\psi} is
#' defined as \eqn{\kappa(x) = \psi(x) + x \psi'(x)}.
#' `igl_kappa` takes \eqn{\psi} to be the IGL generating
#' function, `igl_gen()`.
#'
#' `igl_kappa_D()` is the derivative.
#'
#' @param x Numeric argument of the kappa function. Vectorized.
#' @param p Numeric argument of the inverse function. Vectorized. Between 0 and 1.
#' @param alpha Parameter of the IGL generating function, `igl_gen()`, >0. Vectorized.
#' @rdname kappa
igl_kappa <- function(x, alpha) {
  stats::pgamma(x, shape = alpha, lower.tail = FALSE)
}

#' @rdname kappa
igl_kappa_D <- function(x, alpha) {
  -stats::dgamma(x, shape = alpha)
}

#' @rdname kappa
igl_kappa_inv <- function(p, alpha) {
  stats::qgamma(1 - p, shape = alpha)
}
