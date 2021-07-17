#' #' @rdname interpolator
#' #' @export
#' interp_kappa <- function(x, eta, alpha) {
#'   exp(-x) * igl_kappa(eta * x, alpha)
#' }
#'
#' #' @rdname interpolator
#' #' @export
#' interp_kappa_D1 <- function(x, eta, alpha) {
#'   -exp(-x) * (igl_kappa(eta * x, alpha) - eta * igl_kappa_D(eta * x, alpha))
#' }
#'
#' #' @rdname interpolator
#' #' @export
#' interp_kappa_inv <- function(p, eta, alpha, mxiter = 80, eps = 1.e-12, bd = 5) {
#'   l <- vctrs::vec_size_common(p, eta, alpha)
#'   if (l == 0L) return(numeric(0L))
#'   args <- vctrs::vec_recycle_common(p = p, eta = eta, alpha = alpha)
#'   with(args, {
#'     x <- numeric(0L)
#'     for (i in 1:l) {
#'       x[i] <- interp_kappa_inv_algo(
#'         p[i], eta = eta[i], alpha = alpha[i], mxiter = mxiter, eps = eps, bd = bd
#'       )
#'     }
#'     x
#'   })
#' }
#'
#'
#'
