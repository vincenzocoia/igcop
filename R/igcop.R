#' IG Copula Family Functions
#'
#' Functions related to the IG copula family, denoted  by \code{'igcop'}.
#'
#' @param u,v Vectors of values between 0 and 1 representing values of the first
#' and second copula variables.
#' @param tau Vector of quantile levels between 0 and 1
#'  to evaluate a quantile function at.
#' @param cpar Vector of length 2 corresponding to the copula
#' parameters \code{theta>0} and \code{alpha>0}, respectively.
#' @param n Positive integer. Number of observations to randomly draw.
#' @note Inputting two vectors greater than length 1 is allowed, if they're
#' the same length.
#' Also, \code{qcondigcop21} and \code{pcondigcop21} are the same as
#' \code{qcondigcop} and \code{pcondigcop} -- they're the distributions of
#' variable 2 given 1.
#' @return Numeric vector of length equal to the length of the input vector(s).
#' @rdname igcop
#' @export
pcondigcop <- function(v, u, cpar) {
    theta <- cpar[1]
    alpha <- cpar[2]
    if (theta == Inf) return(pcondiglcop(v, u, cpar = alpha))
    y <- interp_gen_inv(1 - v, eta = theta, alpha = alpha)
    1 - interp_kappa(y, eta = theta * (1 - u), alpha = alpha)
}

#' @param tau Vector of quantile levels between 0 and 1
#' to evaluate a quantile function at.
#' @rdname igcop
#' @export
qcondigcop <- function(tau, u, cpar) {
    theta <- cpar[1]
    alpha <- cpar[2]
    if (theta == Inf) return(qcondiglcop(tau, u, cpar = alpha))
    inner <- interp_kappa_inv(1 - tau, eta = theta * (1 - u), alpha = alpha)
    1 - interp_gen(inner, eta = theta, alpha = alpha)
}

#' @rdname igcop
#' @export
qcondigcop21 <- qcondigcop

#' @rdname igcop
#' @export
pcondigcop21 <- pcondigcop

#' @rdname igcop
#' @export
pcondigcop12 <- function(u, v, cpar) {
    theta <- cpar[1]
    alpha <- cpar[2]
    if (theta == Inf) return(pcondiglcop12(u, v, cpar = alpha))
    y <- interp_gen_inv(1 - v, eta = theta, alpha = alpha)
    1 - (1 - u) *
        interp_gen_D1(y, eta = theta * (1 - u), alpha = alpha) /
        interp_gen_D1(y, eta = theta, alpha = alpha)
}

#' @rdname igcop
#' @export
digcop <- function(u, v, cpar) {
    theta <- cpar[1]
    alpha <- cpar[2]
    if (theta == Inf) return(diglcop(u, v, cpar = alpha))
    y <- interp_gen_inv(1 - v, eta = theta, alpha = alpha)
    interp_kappa_D1(y, eta = (1 - u) * theta, alpha = alpha) /
        interp_gen_D1(y, eta = theta, alpha = alpha)
}

#' @rdname igcop
#' @export
logdigcop <- function(u, v, cpar) {
    theta <- cpar[1]
    alpha <- cpar[2]
    if (theta == Inf) return(logdiglcop(u, v, cpar = alpha))
    y <- interp_gen_inv(1 - v, eta = theta, alpha = alpha)
    eta2 <- (1 - u) * theta
    log(igl_kappa(eta2 * y, alpha) - eta2 * igl_kappa_D(eta2 * y, alpha)) -
        log(igl_gen(theta * y, alpha) - theta * igl_gen_D(theta * y, alpha))
}

#' @rdname igcop
#' @export
pigcop <- function(u, v, cpar) {
    theta <- cpar[1]
    alpha <- cpar[2]
    if (theta == Inf) return(piglcop(u, v, cpar = alpha))
    y <- interp_gen_inv(1 - v, eta = theta, alpha = alpha)
    u + v - 1 + (1 - u) * interp_gen(y, eta = theta * (1 - u), alpha = alpha)
}

#' @rdname igcop
#' @export
rigcop <- function(n, cpar) {
    u <- stats::runif(n)
    tau <- stats::runif(n)
    v <- qcondigcop(tau, u, cpar = cpar)
    res <- matrix(c(u, v), ncol = 2)
    colnames(res) <- c("u", "v")
    if (requireNamespace("tibble", quietly = TRUE)) {
        res <- tibble::as_tibble(res)
    }
    res
}
