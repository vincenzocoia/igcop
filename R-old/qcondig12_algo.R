#' #' Algorithm: Inverse of the 1|2 IG distribution function
#' #' @inheritParams qcondig12
#' qcondig12_algo <- function(p, v, theta, alpha, mxiter = 80, eps = 1.e-12, bd = 5) {
#'   if (length(p) != 1L) stop("Algorithm requires a single `p`.")
#'   if (length(v) != 1L) stop("Algorithm requires a single `v`.")
#'   if (length(theta) != 1L) stop("Algorithm requires a single `theta`.")
#'   if (length(alpha) != 1L) stop("Algorithm requires a single `alpha`.")
#'   prod <- alpha * theta * v * p
#'   if (is.na(prod)) return(prod)
#'   if (p == 0) return(0)
#'   if (p == 1) return(1)
#'   y <- interp_gen_inv(1 - v, eta = theta, alpha = alpha)
#'   denom <- igl_gen(theta * y, alpha = alpha) -
#'     theta * igl_gen_D(theta * y, alpha = alpha)
#'   x0 <- c(p, 1:99/100)
#'   tau0 <- pcondig12(x0, v, theta = theta, alpha = alpha)
#'   diff0 <- abs(p - tau0)
#'   i0 <- which(diff0 == min(diff0))[1]
#'   x <- x0[i0]
#'   x <- -log(x)
#'   iter <- 0
#'   diff <- 1
#'   while (iter < mxiter & abs(diff) > eps) {
#'     ex <- exp(-x)
#'     g <- pcondig12(ex, v, theta = theta, alpha = alpha) - p
#'     gp <- - dig(ex, v, theta = theta, alpha = alpha) * ex
#'     diff <- g / gp
#'     if (x - diff < 0) diff <- x / 2
#'     # if (x - diff > 1) diff <- (1 + x) / 2
#'     x <- x - diff
#'     while(abs(diff) > bd) {
#'       diff <- diff / 2
#'       x <- x + diff
#'     }
#'     iter <- iter + 1
#'   }
#'   exp(-x)
#' }
