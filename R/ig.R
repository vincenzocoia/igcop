#' IG Copula Family Functions
#'
#' Functions related to the IG copula family, denoted  by \code{'ig'}.
#'
#' @param u,v Vectors of values between 0 and 1 representing values of the first
#' and second copula variables.
#' @param tau Vector of quantile levels between 0 and 1
#'  to evaluate a quantile function at.
#' @param cpar Vector of length 2 corresponding to the copula
#' parameters \code{theta>0} and \code{alpha>0}, respectively.
#' @param n Positive integer. Number of observations to randomly draw.
#' @param mxiter Maximum number of iterations to run the Newton-Raphson
#' algorithm when computing inverse. Positive integer, default 20.
#' @param eps The Newton-Raphson algorithm for computing an inverse will
#' stop if the step size is less than this small number.
#' @param bd The largest acceptable step size in the Newton-Raphson
#' algorithm. Step size is reduced if it reaches this large.
#' @note Inputting two vectors greater than length 1 is allowed, if they're
#' the same length.
#' Also, \code{qcondig21} and \code{pcondig21} are the same as
#' \code{qcondig} and \code{pcondig} -- they're the distributions of
#' variable 2 given 1.
#' @return Numeric vector of length equal to the length of the input vector(s).
#' @rdname ig
#' @examples
#' set.seed(1)
#' u <- runif(10)
#' v <- runif(10)
#' pig(u, v, cpar = c(5, 1))
#' dig(u, v, cpar = c(2, 2))
#' logdig(u, v, cpar = c(2, 2))
#' pcondig21(v, u, cpar = c(3, 6))
#' qcondig21(v, u, cpar = c(3, 6))
#' pcondig12(u, v, cpar = c(3, 6))
#' qcondig12(u, v, cpar = c(3, 6))
#' rig(10, cpar = c(3, 3))
#' @export
pcondig21 <- function(v, u, cpar) {
    theta <- cpar[1]
    alpha <- cpar[2]
    if (theta == Inf) return(pcondigl(v, u, cpar = alpha))
    y <- interp_gen_inv(1 - v, eta = theta, alpha = alpha)
    1 - interp_kappa(y, eta = theta * (1 - u), alpha = alpha)
}

#' @param tau Vector of quantile levels between 0 and 1
#' to evaluate a quantile function at.
#' @rdname ig
#' @export
qcondig21 <- function(tau, u, cpar) {
    theta <- cpar[1]
    alpha <- cpar[2]
    if (theta == Inf) return(qcondigl21(tau, u, cpar = alpha))
    inner <- interp_kappa_inv(1 - tau, eta = theta * (1 - u), alpha = alpha)
    1 - interp_gen(inner, eta = theta, alpha = alpha)
}

#' @rdname ig
#' @export
qcondig <- qcondig21

#' @rdname ig
#' @export
pcondig <- pcondig21


qcondig12_algo <- function(tau, v, cpar, mxiter = 80, eps = 1.e-12, bd = 5) {
    if (length(tau) != 1L) stop("Algorithm requires a single `tau`.")
    if (length(v) != 1L) stop("Algorithm requires a single `v`.")
    if (tau == 0) return(0)
    if (tau == 1) return(1)
    theta <- cpar[1]
    alpha <- cpar[2]
    y <- interp_gen_inv(1 - v, eta = theta, alpha = alpha)
    denom <- igl_gen(theta * y, alpha = alpha) -
        theta * igl_gen_D(theta * y, alpha = alpha)
    x0 <- c(tau, 1:99/100)
    tau0 <- pcondig12(x0, v, cpar = cpar)
    diff0 <- abs(tau - tau0)
    i0 <- which(diff0 == min(diff0))[1]
    x <- x0[i0]
    x <- -log(x)
    iter <- 0
    diff <- 1
    while (iter < mxiter & abs(diff) > eps) {
        ex <- exp(-x)
        g <- pcondig12(ex, v, cpar = cpar) - tau
        gp <- - dig(ex, v, cpar = cpar) * ex
        diff <- g / gp
        if (x - diff < 0) diff <- x / 2
        # if (x - diff > 1) diff <- (1 + x) / 2
        x <- x - diff
        while(abs(diff) > bd) {
            diff <- diff / 2
            x <- x + diff
        }
        iter <- iter + 1
    }
    exp(-x)
}

#' @rdname ig
#' @export
qcondig12 <- function(tau, v, cpar, mxiter = 80, eps = 1.e-12, bd = 5) {
    if (cpar[1] == Inf) return(qcondigl12(tau, v, cpar = cpar[2]))
    lengths <- c(tau = length(tau), v = length(v))
    l <- max(lengths)
    if (lengths[["tau"]] == 1) tau <- rep(tau, l)
    if (lengths[["v"]] == 1) v <- rep(v, l)
    x <- numeric()
    for (i in 1:l) {
        x[i] <- qcondig12_algo(
            tau[i], v[i], cpar = cpar, mxiter = mxiter, eps = eps, bd = bd
        )
    }
    x
}



#' @rdname ig
#' @export
pcondig12 <- function(u, v, cpar) {
    theta <- cpar[1]
    alpha <- cpar[2]
    if (theta == Inf) return(pcondigl12(u, v, cpar = alpha))
    y <- interp_gen_inv(1 - v, eta = theta, alpha = alpha)
    1 - (1 - u) *
        interp_gen_D1(y, eta = theta * (1 - u), alpha = alpha) /
        interp_gen_D1(y, eta = theta, alpha = alpha)
}

#' @rdname ig
#' @export
dig <- function(u, v, cpar) {
    theta <- cpar[1]
    alpha <- cpar[2]
    if (theta == Inf) return(digl(u, v, cpar = alpha))
    y <- interp_gen_inv(1 - v, eta = theta, alpha = alpha)
    interp_kappa_D1(y, eta = (1 - u) * theta, alpha = alpha) /
        interp_gen_D1(y, eta = theta, alpha = alpha)
}

#' @rdname ig
#' @export
logdig <- function(u, v, cpar) {
    theta <- cpar[1]
    alpha <- cpar[2]
    if (theta == Inf) return(logdigl(u, v, cpar = alpha))
    y <- interp_gen_inv(1 - v, eta = theta, alpha = alpha)
    eta2 <- (1 - u) * theta
    log(igl_kappa(eta2 * y, alpha) - eta2 * igl_kappa_D(eta2 * y, alpha)) -
        log(igl_gen(theta * y, alpha) - theta * igl_gen_D(theta * y, alpha))
}

#' @rdname ig
#' @export
pig <- function(u, v, cpar) {
    theta <- cpar[1]
    alpha <- cpar[2]
    if (theta == Inf) return(pigl(u, v, cpar = alpha))
    y <- interp_gen_inv(1 - v, eta = theta, alpha = alpha)
    u + v - 1 + (1 - u) * interp_gen(y, eta = theta * (1 - u), alpha = alpha)
}

#' @rdname ig
#' @export
rig <- function(n, cpar) {
    if (cpar[1] == Inf) return(rigl(n, cpar = cpar[2]))
    u <- stats::runif(n)
    tau <- stats::runif(n)
    v <- qcondig(tau, u, cpar = cpar)
    res <- matrix(c(u, v), ncol = 2)
    colnames(res) <- c("u", "v")
    if (requireNamespace("tibble", quietly = TRUE)) {
        res <- tibble::as_tibble(res)
    }
    res
}
