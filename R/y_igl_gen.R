#' Generating function for the IGL copula family
#'
#' \code{igl_gen} is the function itself; \code{igl_gen_inv} is
#' its inverse; and \code{igl_gen_D} is the
#' derivative. Function arguments and parameters are vectorized.
#'
#' @param x Vector of values >=0 to evaluate the function at.
#' @param p Vector of values to evaluate the inverse function at, between
#' 0 and 1 (inclusive).
#' @param alpha Parameter of the function, alpha > 0. Vectorized.
#' @rdname igl_gen
#' @export
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
