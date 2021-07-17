#' Check validity of copula parameters
#'
#' Ensures input values are non-negative.
#'
#' @param theta Values of theta to check.
#' @param alpha Values of alpha to check.
#' @return An error if any theta or alpha is negative;
#' an invisible value otherwise. `NA` values do not throw an error.
#' @rdname check_params
check_alpha <- function(alpha) {
  if (isTRUE(any(alpha < 0))) {
    stop("`alpha` parameter must be positive.")
  }
}


#' @rdname check_params
check_theta <- function(theta) {
  if (isTRUE(any(theta < 0))) {
    stop("`theta` parameter must be positive.")
  }
}
