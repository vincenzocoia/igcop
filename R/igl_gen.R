#' Generating function for the IGL copula family
#'
#' \code{igl_gen} is the function itself, and \code{igl_gen_inv} is
#' its inverse; \code{igl_gen_D} is the
#' derivative; and \code{igl_gen_DD} is the second derivative.
#'
#' @param t Vector of values to evaluate the function at, >=0.
#' @param w Vector of values to evaluate the inverse function at, between
#' 0 and 1 (inclusive)
#' @param k Single numeric >1. Parameter of the function.
#' @examples
#' ## Some examples of evaluating the functions.
#' arg <- c(0, 0.5, 3, Inf, NA)
#' igl_gen(arg, k=2)
#' igl_gen_D(arg, k=1.2)
#' igl_gen_D(arg, k=2)
#' igl_gen_D(arg, k=3)
#' igl_gen_inv(c(0, 0.5, 1), k=1.5)
#'
#' ## Visual
#' foo <- function(u) igl_gen_inv(u, k=1.5)
#' curve(foo)
#' @rdname igl_gen
#' @export
igl_gen <- function(t, k) {
    tinv <- 1 / t
    (k - 1) * t * stats::pgamma(tinv, k) +
        stats::pgamma(tinv, k - 1, lower.tail = FALSE)
}


#' @rdname igl_gen
#' @export
igl_gen_D <- function(t, k) {
    (k - 1) * stats::pgamma(1 / t, k)
}

#' @rdname igl_gen
#' @export
igl_gen_DD <- function(t, k) {
    # -t ^ (-k - 1) * exp(-1 / t) / gamma(k - 1)
    stats::dgamma(1 / t, k - 1) / t ^ 2
}

#' @rdname igl_gen
#' @export
igl_gen_inv <- function(w, k, mxiter = 20, eps = 1.e-12, bd = 5){
    ## Compute gamma(k-1) and gamma(k)
    gkm1 <- gamma(k - 1)
    gk <- (k - 1) * gkm1
    ## Algorithm:
    ## Empirically, it looks like this tt is a good start:
    tt <- (1 - w) ^ (-1 / (k - 1)) - 1
    iter <- 0
    diff <- 1
    ## Begin Newton-Raphson algorithm
    while(iter < mxiter & max(abs(diff)) > eps) {
        ## Helpful quantities
        igam1 <- gkm1 * stats::pgamma(1 / tt, k - 1, lower.tail = FALSE)
        igam0 <- gk * stats::pgamma(1 / tt, k)
        ## Evaluate functions
        g <- tt * igam0 + igam1 - w * gkm1
        gp <- igam0
        diff <- g / gp
        flag <- diff > tt
        diff[flag] <- tt[flag] / 2
        tt <- tt - diff
        while(max(abs(diff)) > bd | any(tt <= 0)) {
            diff <- diff / 2
            tt <- tt + diff
        }
        iter <- iter + 1
        #cat(iter,diff,tt,"\n")
    }
    tt
}
