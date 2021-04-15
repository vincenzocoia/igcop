#' Kappa version of generating function
#'
#' The kappa version of a generating function \eqn{\psi} is
#' defined as \eqn{\kappa(t) = \psi(t) - t \psi'(t)}.
#' `igl_kappa` takes \eqn{\psi} to be the IGL generating
#' function, `igl_gen()`.
#'
#' `igl_kappa_D()` is the derivative.
#'
#' @param t Numeric argument of the kappa function. Vectorized.
#' @param k Parameter of the IGL generating function, `igl_gen()`. Vectorized.
#' @return This function is a nearly direct application of the
#' `stats::pgamma()` function.
#' @rdname kappa
#' @export
igl_kappa <- function(t, k) {
  stats::pgamma(1 / t, shape = k - 1, lower.tail = FALSE)
}

#' @rdname kappa
#' @export
igl_kappa_D <- function(t, k) {
  -t * igl_gen_DD(t, k)
}

#' @rdname kappa
#' @export
igl_kappa_inv <- function(p, k) {
  1 / stats::qgamma(1 - p, shape = k - 1)
}
