#' Kappa version of generating function
#'
#' Defined as \eqn{\kappa(t) = \psi(t) - t \psi'(t)}, where
#' \eqn{\psi} is the IGL generating function `igl_gen()`.
#' @param t Numeric argument of the kappa function. Vectorized.
#' @param k Parameter of the IGL generating function, `igl_gen()`. Vectorized.
#' @return This function is a nearly direct application of the
#' `stats::pgamma()` function.
#' @export
igl_kappa <- function(t, k) {
  stats::pgamma(1 / t, shape = k - 1, lower.tail = FALSE)
}
