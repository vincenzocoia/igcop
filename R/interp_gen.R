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
#' @param alpha Vector of values >0 corresponding to the \eqn{alpha} parameter
#' of the IGL generating function.
#' @rdname interpolator
#' @export
interp_gen <- function(x, eta, alpha) {
    exp(-x) * igl_gen(eta * x, alpha)
}


#' @rdname interpolator
#' @export
interp_gen_D1 <- function(x, eta, alpha) {
    - exp(-x) * (igl_gen(eta * x, alpha) - eta * igl_gen_D(eta * x, alpha))
    # res[x == 1 & k > 2] <- -1
    # res[x == 1 & k == 2] <- -(1 + eta / 2)
    # res[x == 1 & k < 2 & k > 1] <- -Inf
    # res
}


#' @rdname interpolator
#' @export
interp_gen_inv_algo <- function(p, eta, alpha, mxiter = 40, eps = 1.e-12, bd = 5) {
    if (length(p) != 1L) stop("Algorithm requires a single `p`.")
    if (length(eta) != 1L) stop("Algorithm requires a single `eta`.")
    if (length(alpha) != 1L) stop("Algorithm requires a single `alpha`.")
    if (p == 0) return(Inf)
    if (p == 1) return(0)
    x1 <- -log(p)
    x2 <- igl_gen_inv_algo(p, alpha) / eta
    p1 <- interp_gen(x1, eta = eta, alpha = alpha)
    p2 <- interp_gen(x2, eta = eta, alpha = alpha)
    if (abs(p - p1) < abs(p - p2)) {
        x <- x1
    } else {
        x <- x2
    }
    iter <- 0
    diff <- 1
    while(iter < mxiter & abs(diff) > eps) {
        pex <- p * exp(x)
        g <- igl_gen(x * eta, alpha) - pex
        gp <- eta * igl_gen_D(x * eta, alpha) - pex
        diff <- g / gp
        if (x - diff < 0) diff <- x / 2
        x <- x - diff
        while (abs(diff) > bd) {
            diff <- diff / 2
            x <- x + diff
        }
        iter <- iter + 1
    }
    x
}

#' @rdname interpolator
#' @export
interp_gen_inv <- function(p, eta, alpha, mxiter = 40, eps = 1.e-12, bd = 5) {
    lengths <- c(p = length(p), eta = length(eta), alpha = length(alpha))
    l <- max(lengths)
    if (lengths[["p"]] == 1) p <- rep(p, l)
    if (lengths[["eta"]] == 1) eta <- rep(eta, l)
    if (lengths[["alpha"]] == 1) alpha <- rep(alpha, l)
    x <- numeric()
    for (i in 1:l) {
        x[i] <- interp_gen_inv_algo(
            p[i], eta = eta[i], alpha = alpha[i], mxiter = mxiter, eps = eps, bd = bd
        )
    }
    x
}
