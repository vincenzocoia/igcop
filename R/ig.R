#' IG Copula Family Functions
#'
#' Functions related to the IG copula family, denoted  by \code{'ig'}.
#'
#' @param u,v Vectors of values between 0 and 1 representing values of the first
#' and second copula variables.
#' @param p Vector of quantile levels between 0 and 1
#'  to evaluate a quantile function at.
#' @param theta Parameter of the IG copula family. Vectorized; >0.
#' @param alpha Parameter of the IG copula family. Vectorized; >0.
#' @param n Positive integer. Number of observations to randomly draw.
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
#'
#' # log density available for extra precision
#' log(dig(0.1, 0.1, 2.5, 12.3)) == logdig(0.1, 0.1, 2.5, 12.3)
#' @export
pcondig21 <- function(v, u, theta, alpha) {
    check_theta(theta)
    check_alpha(alpha)
    y <- interp_gen_inv(1 - v, eta = theta, alpha = alpha)
    eta <- theta * (1 - u)
    1 - interp_kappa(y, eta = eta, alpha = alpha)
}

#' @param p Vector of quantile levels between 0 and 1
#' to evaluate a quantile function at.
#' @rdname ig
#' @export
qcondig21 <- function(p, u, theta, alpha) {
    check_theta(theta)
    check_alpha(alpha)
    inner <- interp_kappa_inv(1 - p, eta = theta * (1 - u), alpha = alpha)
    1 - interp_gen(inner, eta = theta, alpha = alpha)
}

#' @rdname ig
#' @export
qcondig <- qcondig21

#' @rdname ig
#' @export
pcondig <- pcondig21

#' @rdname ig
#' @export
pcondig12 <- function(u, v, theta, alpha) {
  check_theta(theta)
  check_alpha(alpha)
  formals_to("pcondig12_vec")
}

#' @rdname ig
#' @export
qcondig12 <- function(p, v, theta, alpha) {
  check_theta(theta)
  check_alpha(alpha)
  formals_to("qcondig12_vec")
}

#' @rdname ig
#' @export
dig <- function(u, v, theta, alpha) {
  check_theta(theta)
  check_alpha(alpha)
  formals_to("dig_vec")
}

#' @rdname ig
#' @export
logdig <- function(u, v, theta, alpha) {
    check_theta(theta)
    check_alpha(alpha)
    y <- interp_gen_inv(1 - v, eta = theta, alpha = alpha)
    eta2 <- (1 - u) * theta
    k <- igl_kappa(eta2 * y, alpha)
    kD <- igl_kappa_D(eta2 * y, alpha)
    g <- igl_gen(theta * y, alpha)
    gD <- igl_gen_D(theta * y, alpha)
    log(k - eta2 * kD) - log(g - theta * gD)
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
    p <- stats::runif(n)
    v <- qcondig(p, u, theta = theta, alpha = alpha)
    v_na <- vctrs::vec_slice(v, is.na(v))
    u <- vctrs::vec_assign(u, is.na(v), v_na)
    res <- matrix(c(u, v), ncol = 2)
    colnames(res) <- c("u", "v")
    if (requireNamespace("tibble", quietly = TRUE)) {
        res <- tibble::as_tibble(res)
    }
    res
}


