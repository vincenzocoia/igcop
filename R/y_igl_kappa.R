#' Kappa version of generating function
#'
#' The kappa version of a generating function \eqn{\psi} is
#' defined as \eqn{\kappa(x) = \psi(x) + x \psi'(x)}.
#' `igl_kappa` takes \eqn{\psi} to be the IGL generating
#' function, `igl_gen()`.
#'
#' `igl_kappa_D()` is the derivative; `igl_kappa_inv()` is the inverse.
#'
#' @param x Numeric argument of the kappa function. Vectorized.
#' @param p Numeric argument of the inverse function. Vectorized.
#' Between 0 and 1, inclusive.
#' @param alpha Parameter of the IGL generating function, `igl_gen()`, >0.
#' Vectorized.
#' @rdname kappa
y_igl_kappa <- function(x, alpha) {
  v <- vctrs::vec_recycle_common(x, alpha)
  xvec <- v[[1L]]
  avec <- v[[2L]]
  nn <- length(xvec)
  out <- .C("igl_kappa_vec", as.integer(nn), as.double(xvec), as.double(avec),
            out = as.double(rep(0, nn)), NAOK = TRUE)
  out$out
}

#' @rdname kappa
y_igl_kappa_D <- function(x, alpha) {
  v <- vctrs::vec_recycle_common(x, alpha)
  xvec <- v[[1L]]
  avec <- v[[2L]]
  nn <- length(xvec)
  out <- .C("igl_kappa_D_vec", as.integer(nn), as.double(xvec), as.double(avec),
            out = as.double(rep(0, nn)), NAOK = TRUE)
  out$out
}

#' @rdname kappa
y_igl_kappa_inv <- function(p, alpha) {
  v <- vctrs::vec_recycle_common(p, alpha)
  pvec <- v[[1L]]
  avec <- v[[2L]]
  nn <- length(pvec)
  out <- .C("igl_kappa_inv_vec", as.integer(nn), as.double(pvec),
            as.double(avec), out = as.double(rep(0, nn)), NAOK = TRUE)
  out$out
}
