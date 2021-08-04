#' Interpolating Functions
#'
#' These interpolating functions
#' depend on a generating function (of a DJ copula).
#' `interp_gen()` uses the IGL generating function \eqn{\Psi_k};
#' `interp_kappa()` uses the "kappa version" of that same function.
#'
#' Appending `_inv` to the function name indicates inverse with
#' respect to the first argument. Appending `_D1` indicates
#' derivative with respect to the first argument. Function arguments
#' and parameters are vectorized.
#'
#' @param x Vector of values >=1 to evaluate the interpolating function at.
#' @param p Vector of values between 0 and 1 to evaluate the inverse
#' function at.
#' @param eta Vector of values >0 of second argument of the
#' interpolating function.
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


#' #' @rdname interpolator
#' y_interp_gen_D1 <- function(x, eta, alpha) {
#'   v <- vctrs::vec_recycle_common(x, eta, alpha)
#'   xvec <- v[[1L]]
#'   evec <- v[[2L]]
#'   avec <- v[[3L]]
#'   nn <- length(xvec)
#'   out <- .C("interp_gen_D1_vec", as.integer(nn), as.double(xvec),
#'             as.double(evec), as.double(avec), out = as.double(rep(0, nn)),
#'             NAOK = TRUE)
#'   out$out
#' }
