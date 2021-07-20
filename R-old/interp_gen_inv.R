#' @param mxiter Maximum number of iterations to run the Newton-Raphson
#' algorithm when computing inverse. Positive integer, default 20.
#' @param eps The Newton-Raphson algorithm for computing an inverse will
#' stop if the step size is less than this small number.
#' @param bd The largest acceptable step size in the Newton-Raphson
#' algorithm. Step size is reduced if it reaches this large.
#' @rdname interpolator
#' @export
interp_gen_inv <- function(p, eta, alpha, mxiter = 40, eps = 1.e-12, bd = 5) {
    l <- vctrs::vec_size_common(p, eta, alpha)
    if (l == 0L) return(numeric(0L))
    args <- vctrs::vec_recycle_common(p = p, eta = eta, alpha = alpha)
    with(args, {
        x <- numeric(0L)
        for (i in 1:l) {
            x[i] <- interp_gen_inv_algo(
                p[i], eta = eta[i], alpha = alpha[i], mxiter = mxiter, eps = eps, bd = bd
            )
        }
        x
    })
}



