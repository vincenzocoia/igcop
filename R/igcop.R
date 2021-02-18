#' IG Copula Family Functions
#'
#' Functions related to the IG copula family, denoted  by \code{'igcop'}.
#'
#' @param u,v Vectors of values in [0,1] representing values of the first
#' and second copula variables.
#' @param cpar Vector of length 2 corresponding to the copula
#' parameters \code{theta>0} and \code{k>1}, respectively.
#' @note Inputting two vectors greater than length 1 is allowed, if they're
#' the same length.
#' Also, \code{qcondigcop21} and \code{pcondigcop21} are the same as
#' \code{qcondigcop} and \code{pcondigcop} -- their the distributions of
#' variable 2 given 1.
#' @return Numeric vector of length equal to the length of the input vector(s).
#' @rdname igcop
#' @export
pcondigcop <- function(v, u, cpar) {
    theta <- cpar[1]
    k <- cpar[2]
    if (theta == Inf) return(pcondiglcop(v, u, k))
    Hkinv <- interp_gen_inv(1 - v, theta, k)
    1 - interp_kappa(Hkinv, theta * (1 - u), k)
}

#' @param tau Vector of quantile levels in [0,1] to evaluate a quantile function
#' at.
#' @rdname igcop
#' @export
qcondigcop <- function(tau, u, cpar) {
    theta <- cpar[1]
    k <- cpar[2]
    if (theta == Inf) return(qcondiglcop(tau, u, k))
    inv <- interp_kappa_inv(1 - tau, theta * (1 - u), k)
    1 - interp_gen(inv, theta, k)
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
    k <- cpar[2]
    if (theta == Inf) return(pcondiglcop12(u, v, k))
    Hkinv <- interp_gen_inv(1 - v, theta, k)
    1 - (1 - u) *
        interp_gen_D1(Hkinv, theta * (1 - u), k) /
        interp_gen_D1(Hkinv, theta, k)
}

#' @rdname igcop
#' @export
digcop <- function(u, v, cpar) {
    theta <- cpar[1]
    k <- cpar[2]
    if (theta == Inf) return(diglcop(u, v, k))
    t <- interp_gen_inv(1 - v, theta, k)
    logt <- log(t)
    th_logt <- theta * logt
    th_logt_1mu <- th_logt * (1 - u)
    th_logt2 <- th_logt * logt
    th_logt2_1mu <- th_logt2 * (1 - u)
    num1 <- igl_kappa(1 / th_logt_1mu)
    num2 <- igl_kappa_D(1 / th_logt_1mu) / th_logt2_1mu
    den1 <- igl_gen(1 / th_logt)
    den2 <- igl_gen_inv(1 / th_logt) / th_logt2
    (num1 + num2) / (den1 + den2)
}

#' @rdname igcop
#' @export
logdigcop <- function(u, v, cpar) {
    log(digcop(u, v, cpar))
}

#' @rdname igcop
#' @export
pigcop <- function(u, v, cpar) {
    theta <- cpar[1]
    k <- cpar[2]
    Hinv <- interp_gen_inv(1 - v, theta, k)
    u + v - 1 + (1 - u) * interp_gen(Hinv, theta * (1 - u), k)
}

#' @rdname igcop
#' @export
rigcop <- function(n, cpar) {
    u <- runif(n)
    tau <- runif(n)
    v <- qcondigcop(tau, u, cpar)
    res <- matrix(c(u, v), ncol = 2)
    colnames(res) <- c("u", "v")
    if (requireNamespace("tibble", quietly = TRUE)) {
        res <- tibble::as_tibble(res)
    }
    res
}
