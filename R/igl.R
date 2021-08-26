#' IGL Copula Family Functions
#'
#' Functions related to the IGL copula family, denoted  by \code{'igl'}.
#'
#' @param u,v Vectors of values between 0 and 1 representing values of the first
#' and second copula variables.
#' @param p Vector of quantile levels between 0 and 1 to
#' evaluate a quantile function at.
#' @param alpha Single numeric >0; corresponds to parameter \code{alpha} in the
#' IGL copula family.
#' @param n Positive integer. Number of observations to randomly draw.
#' @note Inputting two vectors greater than length 1 is allowed, if they're
#' the same length.
#' Also, \code{qcondigl21} and \code{pcondigl21} are the same as
#' \code{qcondigl} and \code{pcondigl} -- they are the distributions of
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
qcondigl <- function(p, u, alpha) {
    check_alpha(alpha)
    p <- trim_square(p)
    u <- trim_square(u)
    inner <- y_igl_kappa_inv(1 - p, alpha = alpha) / (1 - u)
    1 - y_igl_gen(inner, alpha = alpha)
}


#' @rdname igl
#' @export
pcondigl <- function(v, u, alpha) {
    check_alpha(alpha)
    v <- trim_square(v)
    u <- trim_square(u)
    y <- y_igl_gen_inv(1 - v, alpha = alpha)
    1 - y_igl_kappa((1 - u) * y, alpha = alpha)
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
    v <- trim_square(v)
    u <- trim_square(u)
    y <- y_igl_gen_inv(1 - v, alpha)
    y <- trim_square(y, upper = FALSE, eps = 1e-10)
    num <- y_igl_gen_D((1 - u) * y, alpha = alpha)
    den <- y_igl_gen_D(y, alpha = alpha)
    den <- trim_square(den, upper = FALSE, eps = 1e-10)
    1 - (1 - u) ^ 2 * num / den
}

#' @rdname igl
#' @export
qcondigl12 <- function(p, v, alpha) {
    check_alpha(alpha)
    p <- trim_square(p)
    v <- trim_square(v)
    y <- y_igl_gen_inv(1 - v, alpha = alpha)
    inner <- (1 - p) * pgamma(y, shape = alpha + 1)
    res <- 1 - qgamma(inner, shape = alpha + 1) / y
    trim_square(res)
}

#' @rdname igl
#' @export
digl <- function(u, v, alpha) {
    check_alpha(alpha)
    v <- trim_square(v)
    u <- trim_square(u)
    y <- y_igl_gen_inv(1 - v, alpha = alpha)
    (1 - u) *
        y_igl_kappa_D((1 - u) * y, alpha = alpha) /
        y_igl_gen_D(y, alpha = alpha)
}


#' @rdname igl
#' @export
pigl <- function(u, v, alpha) {
    check_alpha(alpha)
    v <- trim_square(v)
    u <- trim_square(u)
    y <- y_igl_gen_inv(1 - v, alpha = alpha)
    u + v - 1 + (1 - u) * y_igl_gen((1 - u) * y, alpha = alpha)
}

#' @rdname igl
#' @export
rigl <- function(n, alpha) {
    check_alpha(alpha)
    u <- stats::runif(n)
    p <- stats::runif(n)
    p <- trim_square(p)
    u <- trim_square(u)
    v <- qcondigl(p, u, alpha = alpha)
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
    v <- trim_square(v)
    u <- trim_square(u)
    y <- y_igl_gen_inv(1 - v, alpha = alpha)
    log(1 - u) + 2 * log(y) - log(alpha) +
        dgamma((1 - u) * y, shape = alpha, log = TRUE) -
        pgamma(y, shape = alpha + 1, log.p = TRUE)
}


