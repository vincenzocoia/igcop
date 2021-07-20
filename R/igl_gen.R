#' Generating function for the IGL copula family
#'
#' \code{igl_gen} is the function itself; \code{igl_gen_inv} is
#' its inverse; and \code{igl_gen_D} is the
#' derivative. Function arguments and parameters are vectorized.
#'
#' @param x Vector of values >=0 to evaluate the function at.
#' @param p Vector of values to evaluate the inverse function at, between
#' 0 and 1 (inclusive).
#' @param alpha Parameter of the function, alpha > 0. Vectorized.
#' @rdname igl_gen
igl_gen <- function(x, alpha) {
  res <- stats::pgamma(x, shape = alpha, lower.tail = FALSE) +
    alpha * (stats::pgamma(x, shape = alpha + 1) / x)
  # Use vctrs, otherwise problem w/ res[x==0] <- 1:
  #   if length(res)=0 & length(x)=1, then length(res) becomes 1.
  vctrs::vec_assign(res, x == 0, value = 1)
}


#' @rdname igl_gen
igl_gen_D <- function(x, alpha) {
  res <- - alpha / x ^ 2 * stats::pgamma(x, alpha + 1)
  if (length(res) == 0) return(res)
  res[isTRUE(x == 0)] <- - stats::dgamma(0, shape = alpha[isTRUE(x == 0)]) / 2
  # vctrs::vec_assign(res, x == 0, {
  #   avec <- vctrs::vec_slice(alpha, x == 0)
  #   - stats::dgamma(0, shape = avec) / 2
  # })
  res
}
