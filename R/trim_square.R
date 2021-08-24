#' Trim values to (0, 1)
#'
#' Values between 0 and 1 that are to close to either 0 or 1
#' are replaced with epsilon or 1-epsilon, respectively.
#'
#' @param u Vector of values.
#' @param eps Distance from 0 and 1 for which to truncate.
#' @return The vector `u` with its values possibly truncated.
#' @note
#' This function is used to trim copula inputs, and
#' is intended to be a temporary fix until functions like
#' igcop's in-house `qgamma` and inverse functions have their
#' accuracies improved.
trim_square <- function(u, eps = 1e-4) {
  u <- vctrs::vec_assign(u, u < eps, eps)
  u <- vctrs::vec_assign(u, u > 1 - eps, 1 - eps)
  u
}
