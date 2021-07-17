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
    y <- y_interp_gen_inv(1 - v, eta = theta, alpha = alpha)
    1 - interp_kappa(y, eta = theta * (1 - u), alpha = alpha)
}

#' @param p Vector of quantile levels between 0 and 1
#' to evaluate a quantile function at.
#' @rdname ig
#' @export
qcondig21 <- function(p, u, theta, alpha) {
    check_theta(theta)
    check_alpha(alpha)
    inner <- y_interp_kappa_inv(1 - p, eta = theta * (1 - u), alpha = alpha)
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
    y <- y_interp_gen_inv(1 - v, eta = theta, alpha = alpha)
    1 - (1 - u) *
        interp_gen_D1(y, eta = theta * (1 - u), alpha = alpha) /
        interp_gen_D1(y, eta = theta, alpha = alpha)
}

#' @rdname ig
#' @export
qcondig12 <- function(p, v, theta, alpha, mxiter = 20, eps = 1.e-12, bd = 5)
{
  check_theta(theta)
  check_alpha(alpha)
  recycled <- vctrs::vec_recycle_common(p, v, theta, alpha)
  pvec <- recycled[[1]]
  vvec <- recycled[[2]]
  thvec <- recycled[[3]]
  avec <- recycled[[4]]
  nn <- length(pvec)
  out <- .C("qcondig12", as.integer(nn), as.double(pvec), as.double(vvec),
            as.double(thvec), as.double(avec),
            as.integer(mxiter), as.double(eps), as.double(bd),
            qu = as.double(rep(0, nn)))
  out$qu
}

#' @rdname ig
#' @export
dig <- function(u, v, theta, alpha) {
    check_theta(theta)
    check_alpha(alpha)
    y <- y_interp_gen_inv(1 - v, eta = theta, alpha = alpha)
    interp_kappa_D1(y, eta = (1 - u) * theta, alpha = alpha) /
        interp_gen_D1(y, eta = theta, alpha = alpha)
}

#' @rdname ig
#' @export
logdig <- function(u, v, theta, alpha) {
    check_theta(theta)
    check_alpha(alpha)
    y <- y_interp_gen_inv(1 - v, eta = theta, alpha = alpha)
    eta2 <- (1 - u) * theta
    log(igl_kappa(eta2 * y, alpha) - eta2 * igl_kappa_D(eta2 * y, alpha)) -
        log(igl_gen(theta * y, alpha) - theta * igl_gen_D(theta * y, alpha))
}


#' @rdname ig
#' @export
pig <- function(u, v, theta, alpha) {
    check_theta(theta)
    check_alpha(alpha)
    y <- y_interp_gen_inv(1 - v, eta = theta, alpha = alpha)
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


