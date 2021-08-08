#' Generating function for the IGL copula family
#'
#' `y_igl_gen()` is the function itself; `igl_gen_inv()` is
#' its inverse; and `igl_gen_D()` is the
#' derivative. Function arguments and parameters are vectorized.
#'
#' @param x Vector of values >=0 to evaluate the function at.
#' @param p Vector of values to evaluate the inverse function at, between
#' 0 and 1 (inclusive).
#' @param alpha Function parameter, \eqn{\alpha > 0}. Vectorized.
#' @rdname igl_gen
y_igl_gen <- function(x, alpha) {
  v <- vctrs::vec_recycle_common(x, alpha)
  xvec <- v[[1L]]
  avec <- v[[2L]]
  nn <- length(xvec)
  out <- .C("igl_gen_vec", as.integer(nn), as.double(xvec), as.double(avec),
            out = as.double(rep(0, nn)), NAOK = TRUE)
  out$out
}

#' @rdname igl_gen
y_igl_gen_D <- function(x, alpha) {
  v <- vctrs::vec_recycle_common(x, alpha)
  xvec <- v[[1L]]
  avec <- v[[2L]]
  nn <- length(xvec)
  out <- .C("igl_gen_D_vec", as.integer(nn), as.double(xvec), as.double(avec),
            out = as.double(rep(0, nn)), NAOK = TRUE)
  out$out
}

#' @rdname igl_gen
y_igl_gen_inv <- function(p, alpha)
{
  recycled <- vctrs::vec_recycle_common(p, alpha)
  pvec <- recycled[[1L]]
  avec <- recycled[[2L]]
  nn <- length(pvec)
  mxiter <- 20
  eps <- 1.e-12
  bd <- 5;
  out <- .C("igl_gen_inv", as.integer(nn), as.double(pvec),
            as.double(avec), as.integer(mxiter), as.double(eps), as.double(bd),
            inv = as.double(rep(0, nn)), NAOK = TRUE)
  out$inv
}
