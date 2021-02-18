#' #' Functions for IG copula 2|1 computations
#' #'
#' #' @param t Vector of values to evaluate the assistant function at, >=1.
#' #' @param p Vector of values to evaluate the inverse assistant function at,
#' #' in [0,1].
#' #' @param k Single numeric, greater than 1.
#' #' @param eta eta parameter. Vectorized to match \code{t} and \code{p}.
#' #' @param mxiter Maximum number of iterations to run the Newton-Raphson
#' #' algorithm for.
#' #' @param eps The Newton-Raphson
#' #' algorithm is stopped if all step sizes are below this value.
#' #' @param bd Largest acceptable step size in the Newton-Raphson algorithm;
#' #' steps larger than this will be reduced.
#' #' @rdname igcond
#' #' @export
#' igcond <- function(t, k, eta) {
#'   pgamma(eta * log(t), k - 1, lower.tail = FALSE) / t
#' }
#'
#' #' @rdname igcond
#' #' @export
#' igcondinv <- function(p, k, eta, mxiter = 40, eps = 1.e-6, bd = 5) {
#'   ## Algorithm:
#'   ## Get starting values
#'   xp1 <- 1 / p
#'   xp2 <- exp(qgamma(1 - p, k - 1) / eta)
#'   xpm <- pmin(xp1, xp2)
#'   tt <- pmax(xpm - eps, 1 + (xpm - 1) / 2) # xpm-eps might overshoot left of 1.
#'   iter <- 0
#'   diff <- 1
#'   ## Begin Newton-Raphson algorithm
#'   while(iter < mxiter & max(abs(diff)) > eps){
#'     ## Helpful quantities
#'     etalog <- eta * log(tt)
#'     ## Evaluate functions
#'     g <- tt * p - pgamma(etalog, k - 1, lower.tail = FALSE)
#'     ## When eta=0, derivative is NaN when 1<k<2, when should just be p.
#'     gp <- p + dgamma(etalog, k - 1) * eta / tt
#'     diff <- g / gp
#'     tt <- tt - diff
#'     while(max(abs(diff)) > bd | any(tt <= 1)) {
#'       diff <- diff / 2
#'       tt <- tt + diff
#'     }
#'     iter <- iter + 1
#'     # cat(paste0("-----", iter, "-----\n"))
#'     # cat(diff, "\n")
#'     # cat(tt, "\n")
#'   }
#'   tt
#' }
