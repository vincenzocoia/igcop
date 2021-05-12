#' IGL Copula Family Functions
#'
#' Functions related to the IGL copula family, denoted  by \code{'iglcop'}.
#'
#' @param u,v Vectors of values between 0 and 1 representing values of the first
#' and second copula variables.
#' @param tau Vector of quantile levels between 0 and 1 to
#' evaluate a quantile function at.
#' @param cpar Single numeric >0; corresponds to parameter \code{alpha} in the
#' IGL copula family.
#' @param n Positive integer. Number of observations to randomly draw.
#' @note Inputting two vectors greater than length 1 is allowed, if they're
#' the same length.
#' Also, \code{qcondiglcop21} and \code{pcondiglcop21} are the same as
#' \code{qcondiglcop} and \code{pcondiglcop} -- they're the distributions of
#' variable 2 given 1.
#' @return Numeric vector of length equal to the length of the input vector(s).
#' @rdname iglcop
#' @export
qcondiglcop <- function(tau, u, cpar) {
    alpha <- cpar
    inner <- igl_kappa_inv(1 - tau, alpha = alpha) / (1 - u)
    1 - igl_gen(inner, alpha = alpha)
}

#' @rdname iglcop
#' @export
pcondiglcop <- function(v, u, cpar) {
    alpha <- cpar
    y <- igl_gen_inv(1 - v, alpha = alpha)
    1 - igl_kappa((1 - u) * y, alpha = alpha)
}

#' @rdname iglcop
#' @export
qcondiglcop21 <- qcondiglcop

#' @rdname iglcop
#' @export
pcondiglcop21 <- pcondiglcop

#' @rdname iglcop
#' @export
pcondiglcop12 <- function(u, v, cpar) {
    alpha <- cpar
    y <- igl_gen_inv(1 - v, alpha)
    1 - (1 - u) ^ 2 *
        igl_gen_D((1 - u) * y, alpha = alpha) /
        igl_gen_D(y, alpha = alpha)
}

#' @rdname iglcop
#' @export
qcondiglcop12 <- function(tau, v, cpar) {
    alpha <- cpar
    y <- igl_gen_inv(1 - v, alpha = alpha)
    inner <- (1 - tau) * stats::pgamma(y, shape = alpha + 1)
    1 - stats::qgamma(inner, shape = alpha + 1) / y
}

#' @rdname iglcop
#' @export
diglcop <- function(u, v, cpar) {
    alpha <- cpar
    y <- igl_gen_inv(1 - v, alpha = alpha)
    (1 - u) *
        igl_kappa_D((1 - u) * y, alpha = alpha) /
        igl_gen_D(y, alpha = alpha)
}

#' @rdname iglcop
#' @export
piglcop <- function(u, v, cpar) {
    alpha <- cpar
    y <- igl_gen_inv(1 - v, alpha = alpha)
    u + v - 1 + (1 - u) * igl_gen((1 - u) * y, alpha = alpha)
}

#' @rdname iglcop
#' @export
riglcop <- function(n, cpar) {
    u <- stats::runif(n)
    tau <- stats::runif(n)
    v <- qcondiglcop(tau, u, cpar = cpar)
    res <- matrix(c(u, v), ncol = 2)
    colnames(res) <- c("u", "v")
    if (requireNamespace("tibble", quietly = TRUE)) {
        res <- tibble::as_tibble(res)
    }
    res
}

#' @rdname iglcop
#' @export
logdiglcop <- function(u, v, cpar) {
    alpha <- cpar
    y <- igl_gen_inv(1 - v, alpha = alpha)
    log(1 - u) +
        log(igl_kappa_D((1 - u) * y, alpha = alpha)) -
        log(igl_gen_D(y, alpha = alpha))
}
