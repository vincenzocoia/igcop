#' Gamma Distribution Links to C
#'
#' `pgamma()`, `dgamma()`, and `qgamma()` are defined natively in igcop,
#' written in C. The R functions named with a `y_` prefix link to the ones
#' in C.
#' @param x Vector of quantiles.
#' @param p Vector of probability levels.
#' @param shape,scale Vector of shape and scale parameters.
#' @return Vector of the evaluated representation of the Gamma distribution.
#' @rdname gamma
y_pgamma <- function(x, shape, scale = 1) {
  v <- vctrs::vec_recycle_common(x, shape, scale)
  xvec <- v[[1L]]
  shvec <- v[[2L]]
  scvec <- v[[3L]]
  nn <- length(xvec)
  out <- .C("pgamma_void", as.integer(nn), as.double(xvec), as.double(shvec),
            as.double(scvec), out = as.double(rep(0, nn)), NAOK = TRUE)
  out$out
}

#' @rdname gamma
y_dgamma <- function(x, shape, scale = 1) {
  v <- vctrs::vec_recycle_common(x, shape, scale)
  xvec <- v[[1L]]
  shvec <- v[[2L]]
  scvec <- v[[3L]]
  nn <- length(xvec)
  out <- .C("dgamma_void", as.integer(nn), as.double(xvec), as.double(shvec),
            as.double(scvec), out = as.double(rep(0, nn)), NAOK = TRUE)
  out$out
}

#' @rdname gamma
y_qgamma <- function(p, shape, scale = 1) {
  v <- vctrs::vec_recycle_common(p, shape, scale)
  pvec <- v[[1L]]
  shvec <- v[[2L]]
  scvec <- v[[3L]]
  nn <- length(pvec)
  out <- .C("qgamma_void", as.integer(nn), as.double(pvec), as.double(shvec),
            as.double(scvec), out = as.double(rep(0, nn)), NAOK = TRUE)
  out$out
}
