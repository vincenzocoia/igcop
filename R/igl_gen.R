#' Generating function for the IGL copula family
#'
#' \code{igl_gen} is the function itself, and \code{igl_gen_inv} is
#' its inverse; \code{igl_gen_D} is the
#' derivative; and \code{igl_gen_DD} is the second derivative.
#'
#' Function arguments and parameters are vectorized, except
#' for the algorithms (marked by `_algo`).
#'
#' @param t Vector of values >=0 to evaluate the function at, .
#' @param w Vector of values to evaluate the inverse function at, between
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
igl_gen <- function(t, k) {
    tinv <- 1 / t
    res <- (k - 1) * t * stats::pgamma(tinv, k) +
        stats::pgamma(tinv, k - 1, lower.tail = FALSE)
    res[t == Inf] <- 1
    res
}


#' @rdname igl_gen
#' @export
igl_gen_D <- function(t, k) {
    (k - 1) * stats::pgamma(1 / t, k)
}

#' @rdname igl_gen
#' @export
igl_gen_DD <- function(t, k) {
    stats::dgamma(1 / t, k - 1) / t ^ 2
}


#' @rdname igl_gen
#' @export
igl_gen_inv_algo <- function(w, k, mxiter = 20, eps = 1.e-12, bd = 5){
    if (length(w) != 1L) stop("Algorithm requires a single `w`.")
    if (length(k) != 1L) stop("Algorithm requires a single `k`.")
    if (w == 0) return(0)
    if (w == 1) return(Inf)
    ## Compute gamma(k-1) and gamma(k)
    gkm1 <- gamma(k - 1)
    gk <- (k - 1) * gkm1
    ## Algorithm:
    ## Empirically, it looks like this t is a good start:
    t <- (1 - w) ^ (-1 / (k - 1)) - 1
    ## Lower bound (invert igl_gen after removing (k - 1) * t * pgamma(1/t, k) term)
    # t_low <- 1 / qgamma(1 - w, k - 1)
    iter <- 0
    diff <- 1
    ## Begin Newton-Raphson algorithm
    while (iter < mxiter & abs(diff) > eps) {
        ## Helpful quantities
        igam1 <- gkm1 * stats::pgamma(1 / t, k - 1, lower.tail = FALSE)
        igam0 <- gk * stats::pgamma(1 / t, k)
        ## Evaluate functions
        g <- t * igam0 + igam1 - w * gkm1
        gp <- igam0
        diff <- g / gp
        if (diff > t) diff <- t / 2
        t <- t - diff
        while (abs(diff) > bd | t <= 0) {
            diff <- diff / 2
            t <- t + diff
        }
        iter <- iter + 1
        #cat(iter,diff,t,"\n")
    }
    t
}


#' @rdname igl_gen
#' @export
igl_gen_inv <- function(w, k, mxiter = 20, eps = 1.e-12, bd = 5){
    lengths <- c(w = length(w), k = length(k))
    l <- max(lengths)
    if (lengths[["w"]] == 1) w <- rep(w, l)
    if (lengths[["k"]] == 1) k <- rep(k, l)
    sol <- numeric()
    for (i in 1:l) {
        sol[i] <- igl_gen_inv_algo(
            w[i], k = k[i], mxiter = mxiter, eps = eps, bd = bd
        )
    }
    sol
}

