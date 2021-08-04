# This function has been moved to C.

#' @rdname ig
#' @export
dig <- function(u, v, theta, alpha) {
  check_theta(theta)
  check_alpha(alpha)
  y <- y_interp_gen_inv(1 - v, eta = theta, alpha = alpha)
  y_interp_kappa_D1(y, eta = (1 - u) * theta, alpha = alpha) /
    y_interp_gen_D1(y, eta = theta, alpha = alpha)
}
