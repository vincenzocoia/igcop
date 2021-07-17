#' @rdname ig
#' @export
qcondig12 <- function(p, v, theta, alpha, mxiter = 20, eps = 1.e-12, bd = 5)
{
  check_theta(theta)
  check_alpha(alpha)
  recycled <- vctrs::vec_recycle_common(p, v, theta, alpha)
  pvec <- recycled[[1]]
  vvec <- recycled[[2]]
  thvec <- recycled[[3]]
  avec <- recycled[[4]]
  nn <- length(pvec)
  out <- .C("qcondig12", as.integer(nn), as.double(pvec), as.double(vvec),
            as.double(thvec), as.double(avec),
            as.integer(mxiter), as.double(eps), as.double(bd),
            qu = as.double(rep(0, nn)))
  out$qu
}


