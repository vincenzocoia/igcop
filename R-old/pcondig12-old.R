# This function has been moved to C.

#' @rdname ig
#' @export
pcondig12 <- function(u, v, theta, alpha) {
  check_theta(theta)
  check_alpha(alpha)
  y <- y_interp_gen_inv(1 - v, eta = theta, alpha = alpha)
  1 - (1 - u) *
    y_interp_gen_D1(y, eta = theta * (1 - u), alpha = alpha) /
    y_interp_gen_D1(y, eta = theta, alpha = alpha)
}
