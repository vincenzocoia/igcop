#' (Upper) Incomplete Gamma Function
#'
#' Calculates the (upper) incomplete gamma function,
#' based on the Gamma cdf `pgamma`.
#'
#' @param s Scale parameter of the incomplete gamma function.
#' Vectorized, > 0.
#' @param x Lower bound of the incomplete gamma function.
#' Vectorized, > 0.
#'
#' @details
#' The (Upper) Incomplete Gamma Function is defined as
#' \deqn{\Gamma(s, x) = \int_x^{\infty} t^{s - 1} e^{-t} dt}.
#' @examples
#' igamma(1:5, 1:5)
#' igamma(0.4, 1:5)
#' igamma(1:5, 0.4)
#'
#' # Should be identical:
#' igamma(1:5, 0)
#' gamma(1:5)
igamma <- function(s, x) {
  gamma(s) * stats::pgamma(x, shape = s, lower.tail = FALSE)
}
