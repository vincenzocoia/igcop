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
#' @param x Vector of values >=1 to evaluate the interpolating function at.
#' @param p Vector of values in [0,1] to evaluate the inverse function at.
#' @param eta Vector of values >0 of second argument of the
#' interpolating function.
#' @param k Vector of values >1 corresponding to the \eqn{k} parameter
#' of the IGL generating function.
#' @rdname interpolator
#' @export
interp_gen <- function(x, eta, k) {
    igl_gen(1 / (eta * log(x)), k) / x
}


#' @rdname interpolator
#' @export
interp_gen_D1 <- function(x, eta, k) {
    logt <- log(x)
    arg <- 1 / eta / logt
    coeff <- 1 / eta / logt ^ 2
    res <- -x ^ (-2) * (igl_gen(arg, k) + coeff * igl_gen_D(arg, k))
    res[x == 1 & k > 2] <- -1
    res[x == 1 & k == 2] <- -(1 + eta / 2)
    res[x == 1 & k < 2 & k > 1] <- -Inf
    res
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
    init1 <- 1 / p
    init2 <- exp(1 / eta / igl_gen_inv(p, k))
    x <- pmin(init1, init2)
    x <- pmax(x - eps, 1 + (x - 1) / 2) # x-eps might overshoot left of 1.
    x <- pmax(1 + eps, x)  # x might be 1.
    iter <- 0
    diff <- 1
    ## Begin Newton-Raphson algorithm
    while(iter < mxiter & max(abs(diff)) > eps) {
        ## Helpful quantities
        logt <- log(x)
        etalog <- eta * logt
        etaloginv <- 1 / etalog
        ## Evaluate functions
        g <- x * p - igl_gen(etaloginv, k)
        gp <- p + etaloginv / x / logt * igl_gen_D(etaloginv, k)
        diff <- g / gp
        if (diff > x - 1) diff <- (x - 1) / 2
        x <- x - diff
        while (abs(diff) > bd | x <= 1) {
            diff <- diff / 2
            x <- x + diff
        }
        iter <- iter + 1
    }
    # print(iter)
    x
}

#' @rdname interpolator
#' @export
interp_gen_inv <- function(p, eta, k, mxiter = 40, eps = 1.e-12, bd = 5) {
    lengths <- c(p = length(p), eta = length(eta), k = length(k))
    l <- max(lengths)
    if (lengths[["p"]] == 1) p <- rep(p, l)
    if (lengths[["eta"]] == 1) eta <- rep(eta, l)
    if (lengths[["k"]] == 1) k <- rep(k, l)
    x <- numeric()
    for (i in 1:l) {
        x[i] <- interp_gen_inv_algo(
            p[i], eta = eta[i], k = k[i], mxiter = mxiter, eps = eps, bd = bd
        )
    }
    x
}
