#' IGL Copula Family Functions
#'
#' Functions related to the IGL copula family, denoted  by \code{'igl'}.
#'
#' @param u,v Vectors of values between 0 and 1 representing values of the first
#' and second copula variables.
#' @param tau Vector of quantile levels between 0 and 1 to
#' evaluate a quantile function at.
#' @param alpha Single numeric >0; corresponds to parameter \code{alpha} in the
#' IGL copula family.
#' @param n Positive integer. Number of observations to randomly draw.
#' @note Inputting two vectors greater than length 1 is allowed, if they're
#' the same length.
#' Also, \code{qcondigl21} and \code{pcondigl21} are the same as
#' \code{qcondigl} and \code{pcondigl} -- they're the distributions of
#' variable 2 given 1.
#' @return Numeric vector of length equal to the length of the input vector(s).
#' @rdname igl
#' @examples
#' set.seed(1)
#' u <- runif(10)
#' v <- runif(10)
#' pigl(u, v, alpha = 1)
#' digl(u, v, alpha = 2)
#' logdigl(u, v, alpha = 0.4)
#' pcondigl21(v, u, alpha = 6)
#' qcondigl21(v, u, alpha = 6)
#' pcondigl12(u, v, alpha = 6)
#' qcondigl12(u, v, alpha = 6)
#' rigl(10, alpha = 3)
#' @export
qcondigl <- function(tau, u, alpha) {
    check_alpha(alpha)
    alpha <- alpha
    inner <- igl_kappa_inv(1 - tau, alpha = alpha) / (1 - u)
    1 - igl_gen(inner, alpha = alpha)
}

#' @rdname igl
#' @export
pcondigl <- function(v, u, alpha) {
    check_alpha(alpha)
    alpha <- alpha
    y <- igl_gen_inv(1 - v, alpha = alpha)
    1 - igl_kappa((1 - u) * y, alpha = alpha)
}

#' @rdname igl
#' @export
qcondigl21 <- qcondigl

#' @rdname igl
#' @export
pcondigl21 <- pcondigl

#' @rdname igl
#' @export
pcondigl12 <- function(u, v, alpha) {
    check_alpha(alpha)
    alpha <- alpha
    y <- igl_gen_inv(1 - v, alpha)
    1 - (1 - u) ^ 2 *
        igl_gen_D((1 - u) * y, alpha = alpha) /
        igl_gen_D(y, alpha = alpha)
}

#' @rdname igl
#' @export
qcondigl12 <- function(tau, v, alpha) {
    check_alpha(alpha)
    alpha <- alpha
    y <- igl_gen_inv(1 - v, alpha = alpha)
    inner <- (1 - tau) * stats::pgamma(y, shape = alpha + 1)
    1 - stats::qgamma(inner, shape = alpha + 1) / y
}

#' @rdname igl
#' @export
digl <- function(u, v, alpha) {
    check_alpha(alpha)
    alpha <- alpha
    y <- igl_gen_inv(1 - v, alpha = alpha)
    (1 - u) *
        igl_kappa_D((1 - u) * y, alpha = alpha) /
        igl_gen_D(y, alpha = alpha)
}

#' @rdname igl
#' @export
pigl <- function(u, v, alpha) {
    check_alpha(alpha)
    alpha <- alpha
    y <- igl_gen_inv(1 - v, alpha = alpha)
    u + v - 1 + (1 - u) * igl_gen((1 - u) * y, alpha = alpha)
}

#' @rdname igl
#' @export
rigl <- function(n, alpha) {
    check_alpha(alpha)
    u <- stats::runif(n)
    tau <- stats::runif(n)
    v <- qcondigl(tau, u, alpha = alpha)
    v_na <- vctrs::vec_slice(v, is.na(v))
    u <- vctrs::vec_assign(u, is.na(v), v_na)
    res <- matrix(c(u, v), ncol = 2)
    colnames(res) <- c("u", "v")
    if (requireNamespace("tibble", quietly = TRUE)) {
        res <- tibble::as_tibble(res)
    }
    res
}

#' @rdname igl
#' @export
logdigl <- function(u, v, alpha) {
    check_alpha(alpha)
    alpha <- alpha
    y <- igl_gen_inv(1 - v, alpha = alpha)
    log(1 - u) + 2 * log(y) - log(alpha) +
        stats::dgamma((1 - u) * y, shape = alpha, log = TRUE) -
        stats::pgamma(y, shape = alpha + 1, log.p = TRUE)
}

check_alpha <- function(alpha) {
    if (isTRUE(any(alpha < 0))) {
        stop("`alpha` parameter must be positive.")
    }
}
