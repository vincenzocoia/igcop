# H Joe, May 2021
# Functions to convert from standalone R into links of R to C

# Need interp_gen_inv C routine with pointers
# pig : interp_gen_inv
# dig : interp_gen_inv interp_kappa_D1 interp_gen_D1
# logdig  : interp_gen_inv igl_kappa igl_kappa_D igl_gen igl_gen_D
# pcondig21 : interp_gen_inv interp_kappa
# pcondig12 : interp_gen_inv interp_gen_D1
# qcondig12_algo :  interp_gen_inv igl_gen igl_gen_D pcondig12
# qcondig12 : qcondig12_algo (link to interp_gen_inv C routine)

# Need igl_gen_inv C routine with pointers
# pcondigl21 : igl_gen_inv igl_kappa
# pcondigl12 : igl_gen_inv igl_gen_D
# qcondigl12 : igl_gen_inv pgamma qgamma
# digl : igl_gen_inv igl_kappa_D igl_gen_D
# pigl : igl_gen_inv igl_gen

# Need interp_kappa_inv C routine with pointers
# qcondig21 : interp_kappa_inv  interp_gen

# For fewest code changes, added these three functions:
# y_interp_gen_inv <- function(p, eta, alpha)
# y_igl_gen_inv <- function(p, alpha)
# y_interp_kappa_inv <- function(p, eta, alpha)

# Need qcondig12 C routine with pointers


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
    y <- y_interp_gen_inv(1 - v, eta = theta, alpha = alpha)
    1 - interp_kappa(y, eta = theta * (1 - u), alpha = alpha)
}

#' @param tau Vector of quantile levels between 0 and 1
#' to evaluate a quantile function at.
#' @rdname ig
#' @export
qcondig21 <- function(tau, u, theta, alpha) {
    check_theta(theta)
    check_alpha(alpha)
    inner <- y_interp_kappa_inv(1 - tau, eta = theta * (1 - u), alpha = alpha)
    1 - interp_gen(inner, eta = theta, alpha = alpha)
}

# new function for link to C
# Roxygen2 documentation to be added
y_interp_kappa_inv <- function(p, eta, alpha)
{ # replace with link to C code
  recycled <- vctrs::vec_recycle_common(p, eta, alpha)
  pvec <- recycled[[1]]
  thvec <- recycled[[2]]
  avec <- recycled[[3]]
  nn <- length(pvec)
  mxiter <- 20
  eps <- 1.e-12
  bd <- 5
  out <- .C("interp_kappa_inv", as.integer(nn), as.double(pvec),
            as.double(thvec), as.double(avec), as.integer(mxiter),
            as.double(eps), as.double(bd), inv=as.double(rep(0, nn)))
  out$inv
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


y_interp_gen_inv <- function(p, eta, alpha)
{ # replace with link to C code
  recycled <- vctrs::vec_recycle_common(p, eta, alpha)
  pvec <- recycled[[1]]
  thvec <- recycled[[2]]
  avec <- recycled[[3]]
  nn <- length(pvec)
  mxiter <- 20
  eps <- 1.e-12
  bd <- 5
  out <- .C("interp_gen_inv", as.integer(nn), as.double(pvec), as.double(thvec),
            as.double(avec), as.integer(mxiter), as.double(eps), as.double(bd),
            inv=as.double(rep(0, nn)))
  out$inv
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
