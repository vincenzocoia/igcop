#' @rdname igl_gen
y_igl_gen_inv <- function(p, alpha)
{
  recycled <- vctrs::vec_recycle_common(p, alpha)
  pvec <- recycled[[1]]
  avec <- recycled[[2]]
  nn <- length(pvec)
  mxiter <- 20
  eps <- 1.e-12
  bd <- 5;
  out <- .C("igl_gen_inv", as.integer(nn), as.double(pvec),
            as.double(avec), as.integer(mxiter), as.double(eps), as.double(bd),
            inv=as.double(rep(0, nn)))
  out$inv
}
