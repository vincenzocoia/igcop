#' @rdname interpolator
y_interp_kappa_inv <- function(p, eta, alpha)
{
  recycled <- vctrs::vec_recycle_common(p, eta, alpha)
  pvec <- recycled[[1L]]
  thvec <- recycled[[2L]]
  avec <- recycled[[3L]]
  nn <- length(pvec)
  mxiter <- 20
  eps <- 1.e-12
  bd <- 5
  out <- .C("interp_kappa_inv", as.integer(nn), as.double(pvec),
            as.double(thvec), as.double(avec), as.integer(mxiter),
            as.double(eps), as.double(bd), inv = as.double(rep(0, nn)),
            NAOK = TRUE)
  out$inv
}
