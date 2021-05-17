#' IG Copula Family Functions
#'
#' Functions related to the IG copula family, denoted  by \code{'ig'}.
#'
#' @param u,v Vectors of values between 0 and 1 representing values of the first
#' and second copula variables.
#' @param tau Vector of quantile levels between 0 and 1
#'  to evaluate a quantile function at.
#' @param theta Parameter of the IG copula family. Vectorized; >0.
#' @param alpha Parameter of the IG copula family. Vectorized; >0.
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
#' u <- runif(10)
#' v <- runif(10)
#' pig(u, v, theta = 5, alpha = 1)
#' dig(u, v, theta = 2, alpha = 2)
#' logdig(u, v, theta = 2, alpha = 2)
#' pcondig21(v, u, theta = 3, alpha = 6)
#' qcondig21(v, u, theta = 3, alpha = 6)
#' pcondig12(u, v, theta = 3, alpha = 6)
#' qcondig12(u, v, theta = 3, alpha = 6)
#' rig(10, theta = 3, alpha = 3)
#' @export
pcondig21 <- function(v, u, theta, alpha) {
    check_theta(theta)
    check_alpha(alpha)
    y <- interp_gen_inv(1 - v, eta = theta, alpha = alpha)
    1 - interp_kappa(y, eta = theta * (1 - u), alpha = alpha)
}

#' @param tau Vector of quantile levels between 0 and 1
#' to evaluate a quantile function at.
#' @rdname ig
#' @export
qcondig21 <- function(tau, u, theta, alpha) {
    check_theta(theta)
    check_alpha(alpha)
    inner <- interp_kappa_inv(1 - tau, eta = theta * (1 - u), alpha = alpha)
    1 - interp_gen(inner, eta = theta, alpha = alpha)
}

#' @rdname ig
#' @export
qcondig <- qcondig21

#' @rdname ig
#' @export
pcondig <- pcondig21


qcondig12_algo <- function(tau, v, theta, alpha, mxiter = 80, eps = 1.e-12, bd = 5) {
    if (length(tau) != 1L) stop("Algorithm requires a single `tau`.")
    if (length(v) != 1L) stop("Algorithm requires a single `v`.")
    if (length(theta) != 1L) stop("Algorithm requires a single `theta`.")
    if (length(alpha) != 1L) stop("Algorithm requires a single `alpha`.")
    prod <- alpha * theta * v * tau
    if (is.na(prod)) return(prod)
    if (tau == 0) return(0)
    if (tau == 1) return(1)
    y <- interp_gen_inv(1 - v, eta = theta, alpha = alpha)
    denom <- igl_gen(theta * y, alpha = alpha) -
        theta * igl_gen_D(theta * y, alpha = alpha)
    x0 <- c(tau, 1:99/100)
    tau0 <- pcondig12(x0, v, theta = theta, alpha = alpha)
    diff0 <- abs(tau - tau0)
    i0 <- which(diff0 == min(diff0))[1]
    x <- x0[i0]
    x <- -log(x)
    iter <- 0
    diff <- 1
    while (iter < mxiter & abs(diff) > eps) {
        ex <- exp(-x)
        g <- pcondig12(ex, v, theta = theta, alpha = alpha) - tau
        gp <- - dig(ex, v, theta = theta, alpha = alpha) * ex
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
qcondig12 <- function(tau, v, theta, alpha, mxiter = 80, eps = 1.e-12, bd = 5) {
    check_theta(theta)
    check_alpha(alpha)
    l <- vctrs::vec_size_common(tau, v, theta, alpha)
    if (l == 0L) return(numeric(0L))
    args <- vctrs::vec_recycle_common(
        tau = tau, v = v, theta = theta, alpha = alpha
    )
    with(args, {
        x <- numeric()
        for (i in 1:l) {
            x[i] <- qcondig12_algo(
                tau[i], v[i], theta = theta[i], alpha = alpha[i],
                mxiter = mxiter, eps = eps, bd = bd
            )
        }
        x
    })
}



#' @rdname ig
#' @export
pcondig12 <- function(u, v, theta, alpha) {
    check_theta(theta)
    check_alpha(alpha)
    y <- interp_gen_inv(1 - v, eta = theta, alpha = alpha)
    1 - (1 - u) *
        interp_gen_D1(y, eta = theta * (1 - u), alpha = alpha) /
        interp_gen_D1(y, eta = theta, alpha = alpha)
}

#' @rdname ig
#' @export
dig <- function(u, v, theta, alpha) {
    check_theta(theta)
    check_alpha(alpha)
    y <- interp_gen_inv(1 - v, eta = theta, alpha = alpha)
    interp_kappa_D1(y, eta = (1 - u) * theta, alpha = alpha) /
        interp_gen_D1(y, eta = theta, alpha = alpha)
}

#' @rdname ig
#' @export
logdig <- function(u, v, theta, alpha) {
    check_theta(theta)
    check_alpha(alpha)
    y <- interp_gen_inv(1 - v, eta = theta, alpha = alpha)
    eta2 <- (1 - u) * theta
    log(igl_kappa(eta2 * y, alpha) - eta2 * igl_kappa_D(eta2 * y, alpha)) -
        log(igl_gen(theta * y, alpha) - theta * igl_gen_D(theta * y, alpha))
}

#' @rdname ig
#' @export
pig <- function(u, v, theta, alpha) {
    check_theta(theta)
    check_alpha(alpha)
    y <- interp_gen_inv(1 - v, eta = theta, alpha = alpha)
    u + v - 1 + (1 - u) * interp_gen(y, eta = theta * (1 - u), alpha = alpha)
}

#' @rdname ig
#' @export
rig <- function(n, theta, alpha) {
    check_theta(theta)
    check_alpha(alpha)
    u <- stats::runif(n)
    tau <- stats::runif(n)
    v <- qcondig(tau, u, theta = theta, alpha = alpha)
    v_na <- vctrs::vec_slice(v, is.na(v))
    u <- vctrs::vec_assign(u, is.na(v), v_na)
    res <- matrix(c(u, v), ncol = 2)
    colnames(res) <- c("u", "v")
    if (requireNamespace("tibble", quietly = TRUE)) {
        res <- tibble::as_tibble(res)
    }
    res
}

check_theta <- function(theta) {
    if (isTRUE(any(theta < 0))) {
        stop("`theta` parameter must be positive.")
    }
}
