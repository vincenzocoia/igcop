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
y_interp_gen_inv <- function(p, eta, alpha)
{
  recycled <- vctrs::vec_recycle_common(p, eta, alpha)
  pvec <- recycled[[1]]
  thvec <- recycled[[2]]
  avec <- recycled[[3]]
  nn <- length(pvec)
  mxiter <- 20
  eps <- 1.e-12
  bd <- 5
  out <- .C("interp_gen_inv", as.integer(nn), as.double(pvec), as.double(thvec),
            as.double(avec), as.integer(mxiter), as.double(eps), as.double(bd),
            inv = as.double(rep(0, nn)))
  out$inv
}


#' @rdname interpolator
y_interp_kappa_inv <- function(p, eta, alpha)
{
  recycled <- vctrs::vec_recycle_common(p, eta, alpha)
  pvec <- recycled[[1]]
  thvec <- recycled[[2]]
  avec <- recycled[[3]]
  nn <- length(pvec)
  mxiter <- 20
  eps <- 1.e-12
  bd <- 5
  out <- .C("interp_kappa_inv", as.integer(nn), as.double(pvec),
            as.double(thvec), as.double(avec), as.integer(mxiter),
            as.double(eps), as.double(bd), inv = as.double(rep(0, nn)))
  out$inv
}

