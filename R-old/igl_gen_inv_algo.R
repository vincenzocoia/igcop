#' #' Algorithm for computing inverse of IGL generator
#' #' @inheritParams igl_gen_inv
#' igl_gen_inv_algo <- function(p, alpha, mxiter = 20, eps = 1.e-12, bd = 5){
#'   if (length(p) != 1L) stop("Algorithm requires a single `p`.")
#'   if (length(alpha) != 1L) stop("Algorithm requires a single `alpha`.")
#'   prod <- alpha * p
#'   if (is.na(prod)) return(prod)
#'   if (p == 0) return(Inf)
#'   if (p == 1) return(0)
#'   x_start <- c(
#'     1 / ((1 - p) ^ (-1 / alpha) - 1),
#'     alpha / p,
#'     stats::qgamma(p, shape = alpha + 1)
#'   )
#'   p_start <- igl_gen(x_start, alpha = alpha)
#'   diff_start <- abs(p_start - p)
#'   best_start <- which(diff_start == min(diff_start))[1]
#'   x <- x_start[best_start]
#'   if (x == 0) x <- eps
#'   x <- log(x)
#'   iter <- 0
#'   diff <- 1
#'   while (iter < mxiter & abs(diff) > eps) {
#'     ex <- exp(x)
#'     g <- igl_gen(ex, alpha) - p
#'     gp <- igl_gen_D(ex, alpha) * ex
#'     diff <- g / gp
#'     # if (x - diff < 0) diff <- x / 2
#'     x <- x - diff
#'     while (abs(diff) > bd) {
#'       diff <- diff / 2
#'       x <- x + diff
#'     }
#'     iter <- iter + 1
#'   }
#'   exp(x)
#' }
