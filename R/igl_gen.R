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
#' @param k Parameter of the function, k > 1. Vectorized.
#' @examples
#' ## Some examples of evaluating the functions.
#' arg <- c(0, 0.5, 3, Inf, NA)
#' #igl_gen(arg, k = 2)
#' #igl_gen_D(arg, k = 1.2)
#' #igl_gen_D(arg, k = 2)
#' #igl_gen_D(arg, k = 3)
#' #igl_gen_inv(c(0, 0.5, 1), k = 1.5)
#'
#' ## Visual
#' #foo <- function(u) igl_gen_inv(u, k = 1.5)
#' #curve(foo)
#' @rdname igl_gen
#' @export
igl_gen <- function(x, k) {
    tinv <- 1 / x
    res <- (k - 1) * x * stats::pgamma(tinv, k) +
        stats::pgamma(tinv, k - 1, lower.tail = FALSE)
    res[x == Inf] <- 1
    res
}

# igl_gen_v2 <- function(x, k) {
#     x * (gamma(k) - igamma(k, 1 / x)) / gamma(k - 1) +
#         igamma(k - 1, 1 / x) / gamma(k - 1)
# }
# diff <- function(x) igl_gen(x, 1.1) - igl_gen_v2(x, 1.1)
# curve(diff, 0, 100)

#' @rdname igl_gen
#' @export
igl_gen_D <- function(x, k) {
    (k - 1) * stats::pgamma(1 / x, k)
}

# igl_gen_D_v2 <- function(x, k) {
#     (gamma(k) - igamma(k, 1 / x)) / gamma(k - 1)
# }
# diff <- function(x) igl_gen_D(x, 3) - igl_gen_D_v2(x, 3)
# curve(diff, 0, 10)


#' @rdname igl_gen
#' @export
igl_gen_DD <- function(x, k) {
    - (k - 1) / x ^ 2 * stats::dgamma(1 / x, k)
}

# igl_gen_DD_v2 <- function(x, k) {
#     - x ^ (-k - 1) * exp(-1 / x) / gamma(k - 1)
# }
# diff <- function(x) igl_gen_DD(x, 3) - igl_gen_DD_v2(x, 3)
# curve(diff, 0, 10)


#' @rdname igl_gen
#' @export
igl_gen_inv_algo <- function(p, k, mxiter = 20, eps = 1.e-12, bd = 5){
    if (length(p) != 1L) stop("Algorithm requires a single `p`.")
    if (length(k) != 1L) stop("Algorithm requires a single `k`.")
    if (p == 0) return(0)
    if (p == 1) return(Inf)
    ## Compute gamma(k-1) and gamma(k)
    gkm1 <- gamma(k - 1)
    gk <- (k - 1) * gkm1
    ## Algorithm:
    ## Empirically, it looks like this x is a good start:
    x <- (1 - p) ^ (-1 / (k - 1)) - 1
    ## Lower bound (invert igl_gen after removing (k - 1) * x * pgamma(1/x, k) term)
    # t_low <- 1 / qgamma(1 - p, k - 1)
    iter <- 0
    diff <- 1
    ## Begin Newton-Raphson algorithm
    while (iter < mxiter & abs(diff) > eps) {
        ## Helpful quantities
        igam1 <- gkm1 * stats::pgamma(1 / x, k - 1, lower.tail = FALSE)
        igam0 <- gk * stats::pgamma(1 / x, k)
        ## Evaluate functions
        g <- x * igam0 + igam1 - p * gkm1
        gp <- igam0
        diff <- g / gp
        if (diff > x) diff <- x / 2
        x <- x - diff
        while (abs(diff) > bd | x <= 0) {
            diff <- diff / 2
            x <- x + diff
        }
        iter <- iter + 1
        #cat(iter,diff,x,"\n")
    }
    x
}


#' @rdname igl_gen
#' @export
igl_gen_inv <- function(p, k, mxiter = 20, eps = 1.e-12, bd = 5){
    lengths <- c(p = length(p), k = length(k))
    l <- max(lengths)
    if (lengths[["p"]] == 1) p <- rep(p, l)
    if (lengths[["k"]] == 1) k <- rep(k, l)
    sol <- numeric()
    for (i in 1:l) {
        sol[i] <- igl_gen_inv_algo(
            p[i], k = k[i], mxiter = mxiter, eps = eps, bd = bd
        )
    }
    sol
}

