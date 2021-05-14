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



## Kappa seems to be calculated correctly. I'd trust the simplified version -- calculation is more direct.
# kappa_v2 <- function(x, k) {
#   igl_gen(x, k) - x * igl_gen_D(x, k)
# }
#

## Kappa primed? No difference.
# kappa_D_v2 <- function(x, k) -x * igl_gen_DD(x, k)
#

## Kappa primed, numeric.
# kappa_D_num <- function(x, k, eps = 1e-8) {
#   (igl_kappa(x + eps, k) - igl_kappa(x, k)) / eps
# }
# kappa_D_num(2, 2)
