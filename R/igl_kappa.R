#' Kappa version of generating function
#'
#' The kappa version of a generating function \eqn{\psi} is
#' defined as \eqn{\kappa(x) = \psi(x) - x \psi'(x)}.
#' `igl_kappa` takes \eqn{\psi} to be the IGL generating
#' function, `igl_gen()`.
#'
#' `igl_kappa_D()` is the derivative.
#'
#' @param x Numeric argument of the kappa function. Vectorized.
#' @param k Parameter of the IGL generating function, `igl_gen()`. Vectorized.
#' @return This function is a nearly direct application of the
#' `stats::pgamma()` function.
#' @rdname kappa
#' @export
igl_kappa <- function(x, k) {
  stats::pgamma(1 / x, shape = k - 1, lower.tail = FALSE)
}

#' @rdname kappa
#' @export
igl_kappa_D <- function(x, k) {
  -x * igl_gen_DD(x, k)
}

#' @rdname kappa
#' @export
igl_kappa_inv <- function(p, k) {
  1 / stats::qgamma(1 - p, shape = k - 1)
}


## Kappa seems to be calculated correctly. I'd trust the simplified version -- calculation is more direct.
# kappa_v2 <- function(x, k) {
#   igl_gen(x, k) - x * igl_gen_D(x, k)
# }
#
# diff <- function(x) kappa_v2(x, 1.1) - igl_kappa(x, 1.1)
# curve(diff, 0, 3)


## Kappa primed? No difference.
# kappa_D_v2 <- function(x, k) -x * igl_gen_DD(x, k)
#
# diff <- function(x) igl_kappa_D(x, 1.1) - kappa_D_v2(x, 1.1)
# curve(diff, 0, 10)

## Kappa primed, numeric.
# kappa_D_num <- function(x, k, eps = 1e-8) {
#   (igl_kappa(x + eps, k) - igl_kappa(x, k)) / eps
# }
# kappa_D_num(2, 2)
