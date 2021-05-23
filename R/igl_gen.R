#' Generating function for the IGL copula family
#'
#' \code{igl_gen} is the function itself, and \code{igl_gen_inv} is
#' its inverse; \code{igl_gen_D} is the
#' derivative; and \code{igl_gen_DD} is the second derivative.
#'
#' Function arguments and parameters are vectorized, except
#' for the algorithms (marked by `_algo`).
#'
#' @param x Vector of values >=0 to evaluate the function at, .
#' @param p Vector of values to evaluate the inverse function at, between
#' 0 and 1 (inclusive).
#' @param alpha Parameter of the function, alpha > 0. Vectorized.
#' @examples
#' ## Some examples of evaluating the functions.
#' arg <- c(0, 0.5, 3, Inf, NA)
#' #igl_gen(arg, alpha = 1)
#' #igl_gen_D(arg, alpha = 0.2)
#' #igl_gen_D(arg, alpha = 1)
#' #igl_gen_D(arg, alpha = 2)
#' #igl_gen_inv(c(0, 0.5, 1), alpha = 0.5)
#'
#' ## Visual
#' #foo <- function(u) igl_gen_inv(u, alpha = 0.5)
#' #curve(foo)
#' @rdname igl_gen
igl_gen <- function(x, alpha) {
    res <- stats::pgamma(x, shape = alpha, lower.tail = FALSE) +
        alpha * (stats::pgamma(x, shape = alpha + 1) / x)
    # Problem w/ res[x==0] <- 1: if length(res)=0 & length(x)=1, then length(res) becomes 1.
    vctrs::vec_assign(res, x == 0, value = 1)
}


#' @rdname igl_gen
igl_gen_D <- function(x, alpha) {
    - alpha / x ^ 2 * stats::pgamma(x, alpha + 1)
}

#' @param mxiter Maximum number of iterations to run the Newton-Raphson
#' algorithm when computing inverse. Positive integer, default 20.
#' @param eps The Newton-Raphson algorithm for computing an inverse will
#' stop if the step size is less than this small number.
#' @param bd The largest acceptable step size in the Newton-Raphson
#' algorithm. Step size is reduced if it reaches this large.
#' @rdname igl_gen
igl_gen_inv_algo <- function(p, alpha, mxiter = 20, eps = 1.e-12, bd = 5){
    if (length(p) != 1L) stop("Algorithm requires a single `p`.")
    if (length(alpha) != 1L) stop("Algorithm requires a single `alpha`.")
    prod <- alpha * p
    if (is.na(prod)) return(prod)
    if (p == 0) return(Inf)
    if (p == 1) return(0)
    x_start <- c(
        1 / ((1 - p) ^ (-1 / alpha) - 1),
        alpha / p,
        stats::qgamma(p, shape = alpha + 1)
    )
    p_start <- igl_gen(x_start, alpha = alpha)
    diff_start <- abs(p_start - p)
    best_start <- which(diff_start == min(diff_start))[1]
    x <- x_start[best_start]
    if (x == 0) x <- eps
    x <- log(x)
    iter <- 0
    diff <- 1
    while (iter < mxiter & abs(diff) > eps) {
        ex <- exp(x)
        g <- igl_gen(ex, alpha) - p
        gp <- igl_gen_D(ex, alpha) * ex
        diff <- g / gp
        # if (x - diff < 0) diff <- x / 2
        x <- x - diff
        while (abs(diff) > bd) {
            diff <- diff / 2
            x <- x + diff
        }
        iter <- iter + 1
    }
    exp(x)
}


#' @rdname igl_gen
igl_gen_inv <- function(p, alpha, mxiter = 20, eps = 1.e-12, bd = 5){
    l <- vctrs::vec_size_common(p, alpha)
    if (l == 0L) return(numeric(0L))
    args <- vctrs::vec_recycle_common(p = p, alpha = alpha)
    with(args, {
        x <- numeric(0L)
        for (i in 1:l) {
            x[i] <- igl_gen_inv_algo(
                p[i], alpha = alpha[i], mxiter = mxiter, eps = eps, bd = bd
            )
        }
        x
    })
}
