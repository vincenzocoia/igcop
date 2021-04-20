#' Interpolating Functions
#'
#' These interpolating functions, denoted \eqn{H} in the vignette,
#' depend on a generating function (of a DJ copula).
#' `interp_gen()` uses the IGL generating function \eqn{\Psi_k};
#' `interp_kappa()` uses the "kappa version" of that same function.
#'
#' Appending `_inv` to the function name indicates inverse with
#' respect to the first argument. Appending `_D1` indicates
#' derivative with respect to the first argument. Function arguments
#' and parameters are vectorized, except for the algorithms (marked by
#' `_algo`).
#'
#' @param t Vector of values >=1 to evaluate the interpolating function at.
#' @param p Vector of values in [0,1] to evaluate the inverse function at.
#' @param eta Vector of values >0 of second argument of the
#' interpolating function.
#' @param k Vector of values >1 corresponding to the \eqn{k} parameter
#' of the IGL generating function.
#' @rdname interpolator
#' @export
interp_gen <- function(t, eta, k) {
    igl_gen(1 / (eta * log(t)), k) / t
}


#' @rdname interpolator
#' @export
interp_gen_D1 <- function(t, eta, k) {
    ## Deal with t=1 separately -- its limit depends on k.
    # if (k > 2) {
    #     replf <- -1
    # } else if (k == 2) {
    #     replf <- -(1 + eta / 2)
    # } else {
    #     replf <- -Inf
    # }
    logt <- log(t)
    arg <- 1 / eta / logt
    coeff <- 1 / eta / logt ^ 2
    -t ^ (-2) * (igl_gen(arg, k) + coeff * igl_gen_D(arg, k))
}


#' @rdname interpolator
#' @export
interp_gen_inv_algo <- function(p, eta, k, mxiter = 40, eps = 1.e-12, bd = 5) {
    if (length(p) != 1L) stop("Algorithm requires a single `p`.")
    if (length(eta) != 1L) stop("Algorithm requires a single `eta`.")
    if (length(k) != 1L) stop("Algorithm requires a single `k`.")
    if (p == 0) return(Inf)
    if (p == 1) return(1)
    ## Algorithm:
    ## Get starting values
    xp1 <- 1 / p
    xp2 <- exp(1 / eta / igl_gen_inv(p, k))
    xpm <- pmin(xp1, xp2)
    t <- pmax(xpm - eps, 1 + (xpm - 1) / 2) # xpm-eps might overshoot left of 1.
    t <- pmax(1 + eps, t)  # t might be 1.
    iter <- 0
    diff <- 1
    ## Begin Newton-Raphson algorithm
    while(iter < mxiter & max(abs(diff)) > eps) {
        ## Helpful quantities
        logt <- log(t)
        etalog <- eta * logt
        etaloginv <- 1 / etalog
        ## Evaluate functions
        g <- t * p - igl_gen(etaloginv, k)
        gp <- p + etaloginv / t / logt * igl_gen_D(etaloginv, k)
        diff <- g / gp
        if (diff > t - 1) diff <- (t - 1) / 2
        t <- t - diff
        while (abs(diff) > bd | t <= 1) {
            diff <- diff / 2
            t <- t + diff
        }
        iter <- iter + 1
    }
    # print(iter)
    t
}

#' @rdname interpolator
#' @export
interp_gen_inv <- function(p, eta, k, mxiter = 40, eps = 1.e-12, bd = 5) {
    vapply(p, function(.p) {
        interp_gen_inv_algo(.p, eta = eta, k = k, mxiter = mxiter, eps = eps, bd = bd)
    },
    FUN.VALUE = numeric(1L))
}
