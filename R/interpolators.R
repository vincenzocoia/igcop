#' Interpolating Functions
#'
#' These interpolating functions, denoted \eqn{H} in the vignette,
#' depend on a generating function (of a DJ copula).
#' `interp_gen()` uses the IGL generating function \eqn{\Psi_k};
#' `interp_kappa()` uses the "kappa version" of that same function.
#'
#' Appending `_inv` to the function name indicates inverse with
#' respect to the first argument. Appending `_D1` indicates
#' derivative with respect to the first argument.
#'
#' @param t Vector of values >=1 to evaluate the interpolating function at.
#' @param eta Value of second argument of the interpolating function, >0.
#' This is allowed to be a vector, except in the inverse functions.
#' @param k Single numeric >1 corresponding to the \eqn{k} parameter
#' of the IGL generating function.
#' @rdname interpolator
#' @export
interp_gen <- function(t, eta, k) {
    igl_gen(1 / (eta * log(t)), k) / t
}

#' @rdname interpolator
#' @export
interp_kappa <- function(t, eta, k) {
    igl_kappa(1 / (eta * log(t)), k) / t
}

#' @rdname interpolator
#' @export
interp_kappa_inv <- function(p, eta, k) {
    stop("This function still needs to be programmed.")
}

#' @rdname interpolator
#' @export
interp_gen_D1 <- function(t, eta, k) {
    ## Deal with t=1 separately -- its limit depends on k.
    if (k > 2) {
        replf <- -1
    } else if (k == 2) {
        replf <- -(1 + eta / 2)
    } else {
        replf <- -Inf
    }
    logt <- log(t)
    arg <- 1 / eta / logt
    coeff <- 1 / eta / logt ^ 2
    -t ^ (-2) * (igl_gen(arg, k) + coeff * igl_gen_D(arg, k))
}


#' @param p Vector of values in [0,1] to evaluate the inverse function at.
#' @rdname interpolator
#' @export
interp_gen_inv <- function(p, eta, k, mxiter = 40, eps = 1.e-12, bd = 5) {
    ## Algorithm:
    ## Get starting values
    xp1 <- 1 / p
    xp2 <- exp(1 / eta / igl_gen_inv(p, k))
    xpm <- pmin(xp1, xp2)
    tt <- pmax(xpm - eps, 1 + (xpm - 1) / 2) # xpm-eps might overshoot left of 1.
    tt <- pmax(1 + eps, tt)  # tt might be 1.
    iter <- 0
    diff <- 1
    ## Begin Newton-Raphson algorithm
    while(iter < mxiter & max(abs(diff)) > eps) {
        ## Helpful quantities
        logt <- log(tt)
        etalog <- eta * logt
        etaloginv <- 1 / etalog
        ## Evaluate functions
        g <- tt * p - igl_gen(etaloginv, k)
        gp <- p + etaloginv / tt / logt * igl_gen_D(etaloginv, k)
        diff <- g / gp
        flag <- diff > tt - 1
        diff[flag] <- (tt[flag] - 1) / 2
        tt <- tt - diff
        while(max(abs(diff)) > bd | any(tt <= 1)) {
            diff <- diff / 2
            tt <- tt + diff
        }
        iter <- iter + 1
    }
    tt
}
