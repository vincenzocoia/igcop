#' #' Interpolating Functions
#' #'
#' #' These interpolating functions
#' #' depend on a generating function (of a DJ copula).
#' #' `interp_gen()` uses the IGL generating function \eqn{\Psi_k};
#' #' `interp_kappa()` uses the "kappa version" of that same function.
#' #'
#' #' Appending `_inv` to the function name indicates inverse with
#' #' respect to the first argument. Appending `_D1` indicates
#' #' derivative with respect to the first argument. Function arguments
#' #' and parameters are vectorized.
#' #'
#' #' @param x Vector of values >=1 to evaluate the interpolating function at.
#' #' @param p Vector of values between 0 and 1 to evaluate the inverse
#' #' function at.
#' #' @param eta Vector of values >0 of second argument of the
#' #' interpolating function.
#' #' @param alpha Vector of values >0 corresponding to the \eqn{alpha} parameter
#' #' of the IGL generating function.
#' #' @rdname interpolator
#' #' @export
#' interp_gen <- function(x, eta, alpha) {
#'     exp(-x) * igl_gen(eta * x, alpha = alpha)
#' }
#'
#'
#' #' @rdname interpolator
#' #' @export
#' interp_gen_D1 <- function(x, eta, alpha) {
#'     - exp(-x) * (igl_gen(eta * x, alpha = alpha) -
#'                      eta * igl_gen_D(eta * x, alpha = alpha))
#'     # res[x == 1 & k > 2] <- -1
#'     # res[x == 1 & k == 2] <- -(1 + eta / 2)
#'     # res[x == 1 & k < 2 & k > 1] <- -Inf
#'     # res
#' }
#'
#' #' @param mxiter Maximum number of iterations to run the Newton-Raphson
#' #' algorithm when computing inverse. Positive integer, default 20.
#' #' @param eps The Newton-Raphson algorithm for computing an inverse will
#' #' stop if the step size is less than this small number.
#' #' @param bd The largest acceptable step size in the Newton-Raphson
#' #' algorithm. Step size is reduced if it reaches this large.
#' #' @rdname interpolator
#' #' @export
#' interp_gen_inv <- function(p, eta, alpha, mxiter = 40, eps = 1.e-12, bd = 5) {
#'     l <- vctrs::vec_size_common(p, eta, alpha)
#'     if (l == 0L) return(numeric(0L))
#'     args <- vctrs::vec_recycle_common(p = p, eta = eta, alpha = alpha)
#'     with(args, {
#'         x <- numeric(0L)
#'         for (i in 1:l) {
#'             x[i] <- interp_gen_inv_algo(
#'                 p[i], eta = eta[i], alpha = alpha[i], mxiter = mxiter, eps = eps, bd = bd
#'             )
#'         }
#'         x
#'     })
#' }
#'
#'
#'
