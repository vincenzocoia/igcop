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
interp_gen <- function(x, eta, alpha) {
  exp(-x) * igl_gen(eta * x, alpha = alpha)
}


#' @rdname interpolator
interp_gen_D1 <- function(x, eta, alpha) {
  res <- - exp(-x) * (igl_gen(eta * x, alpha = alpha) -
                        eta * igl_gen_D(eta * x, alpha = alpha))
  vctrs::vec_assign(res, x == 0, {
    avec <- vctrs::vec_slice(alpha, x == 0)
    -(1 + eta / 2 * stats::dgamma(0, shape = avec))
  })
}

#' @rdname interpolator
interp_kappa <- function(x, eta, alpha) {
  exp(-x) * igl_kappa(eta * x, alpha)
}

#' @rdname interpolator
interp_kappa_D1 <- function(x, eta, alpha) {
  -exp(-x) * (igl_kappa(eta * x, alpha) - eta * igl_kappa_D(eta * x, alpha))
}

