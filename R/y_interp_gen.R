#' Interpolating Functions
#'
#' These interpolating functions
#' depend on a generating function (of a DJ copula).
#' `y_interp_gen()` uses the IGL generating function \eqn{\Psi_k};
#' `y_interp_kappa()` uses the "kappa transform" of that same function.
#'
#' Appending `_inv` to the function name indicates inverse with
#' respect to the first argument. Appending `_D1` indicates
#' derivative with respect to the first argument. Function arguments
#' and parameters are vectorized.
#'
#' @param x Vector of values >=0 to evaluate the function at.
#' @param p Vector of values between 0 and 1 (inclusive) to
#' evaluate the inverse function at.
#' @param eta Vector of values >0 of the interpolating parameter.
#' @param alpha Vector of values >0 corresponding to the \eqn{alpha} parameter
#' of the IGL generating function.
#' @rdname interpolator
y_interp_gen <- function(x, eta, alpha) {
  v <- vctrs::vec_recycle_common(x, eta, alpha)
  xvec <- v[[1L]]
  evec <- v[[2L]]
  avec <- v[[3L]]
  nn <- length(xvec)
  out <- .C("interp_gen_vec", as.integer(nn), as.double(xvec), as.double(evec),
            as.double(avec), out = as.double(rep(0, nn)), NAOK = TRUE)
  out$out
}

#' @rdname interpolator
y_interp_gen_inv <- function(p, eta, alpha) {
  recycled <- vctrs::vec_recycle_common(p, eta, alpha)
  pvec <- recycled[[1L]]
  thvec <- recycled[[2L]]
  avec <- recycled[[3L]]
  nn <- length(pvec)
  mxiter <- 20
  eps <- 1.e-12
  bd <- 5
  out <- .C("interp_gen_inv", as.integer(nn), as.double(pvec), as.double(thvec),
            as.double(avec), as.integer(mxiter), as.double(eps), as.double(bd),
            inv = as.double(rep(0, nn)), NAOK = TRUE)
  out$inv
}

