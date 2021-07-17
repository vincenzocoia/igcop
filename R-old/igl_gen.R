#' #' Generating function for the IGL copula family
#' #'
#' #' \code{igl_gen} is the function itself; \code{igl_gen_inv} is
#' #' its inverse; and \code{igl_gen_D} is the
#' #' derivative. Function arguments and parameters are vectorized.
#' #'
#' #' @param x Vector of values >=0 to evaluate the function at.
#' #' @param p Vector of values to evaluate the inverse function at, between
#' #' 0 and 1 (inclusive).
#' #' @param alpha Parameter of the function, alpha > 0. Vectorized.
#' #' @examples
#' #' arg <- c(0, 0.5, 3, Inf, NA)
#' #' igl_gen(arg, alpha = 1)
#' #' igl_gen(arg, alpha = 0.2)
#' #' igl_gen_D(arg, alpha = 1)
#' #' igl_gen_D(arg, alpha = 2)
#' #' igl_gen_inv(c(0, 0.5, 1), alpha = 0.5)
#' #' @export
#' #' @rdname igl_gen
#' igl_gen <- function(x, alpha) {
#'     res <- stats::pgamma(x, shape = alpha, lower.tail = FALSE) +
#'         alpha * (stats::pgamma(x, shape = alpha + 1) / x)
#'     # Problem w/ res[x==0] <- 1: if length(res)=0 & length(x)=1, then length(res) becomes 1.
#'     vctrs::vec_assign(res, x == 0, value = 1)
#' }
#'
#'
#' #' @rdname igl_gen
#' #' @export
#' igl_gen_D <- function(x, alpha) {
#'     - alpha / x ^ 2 * stats::pgamma(x, alpha + 1)
#' }

#' #' @param mxiter Maximum number of iterations to run the Newton-Raphson
#' #' algorithm for computing the inverse. Positive integer, defaults to 20.
#' #' @param eps The Newton-Raphson algorithm for computing an inverse will
#' #' stop if the step size is less than this small number.
#' #' @param bd The largest acceptable step size in the Newton-Raphson
#' #' algorithm. Step size is reduced if it exceeds this bound.
#' #' @rdname igl_gen
#' #' @export
#' igl_gen_inv <- function(p, alpha, mxiter = 20, eps = 1.e-12, bd = 5){
#'     l <- vctrs::vec_size_common(p, alpha)
#'     if (l == 0L) return(numeric(0L))
#'     args <- vctrs::vec_recycle_common(p = p, alpha = alpha)
#'     with(args, {
#'         x <- numeric(0L)
#'         for (i in 1:l) {
#'             x[i] <- igl_gen_inv_algo(
#'                 p[i], alpha = alpha[i], mxiter = mxiter, eps = eps, bd = bd
#'             )
#'         }
#'         x
#'     })
#' }




