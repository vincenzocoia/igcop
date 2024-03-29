% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{igl_gen_vec}
\alias{igl_gen_vec}
\alias{igl_gen_D_vec}
\alias{igl_gen_inv_vec}
\alias{igl_kappa_vec}
\alias{igl_kappa_D_vec}
\alias{igl_kappa_inv_vec}
\alias{interp_gen_vec}
\alias{interp_gen_inv_vec}
\alias{interp_kappa_vec}
\alias{interp_kappa_inv_vec}
\title{IG/IGL Generators and Related Functions: matching inputs}
\usage{
igl_gen_vec(x, alpha)

igl_gen_D_vec(x, alpha)

igl_gen_inv_vec(p, alpha)

igl_kappa_vec(x, alpha)

igl_kappa_D_vec(x, alpha)

igl_kappa_inv_vec(p, alpha)

interp_gen_vec(x, eta, alpha)

interp_gen_inv_vec(p, eta, alpha)

interp_kappa_vec(x, eta, alpha)

interp_kappa_inv_vec(p, eta, alpha)
}
\arguments{
\item{x}{Function argument. Vector of non-negative values.}

\item{p}{Function inverse argument. Vector of values between 0 and 1.}

\item{eta, alpha}{Function parameters. Vector of positive values.}
}
\description{
These are the psi, H, and kappa functions
of the IG and IGL copula families, but with inputs
needing to be vectors of the same length.
These functions are called by the R functions of the
same name, without the \verb{_vec} suffix.
}
\note{
If calling this function manually, make sure each input
are vectors of a common length.
}
\seealso{
\code{igl_gen()} and family;
\code{dig_vec()}, \code{pcondig12_vec()}, and \code{qcondig12_vec()}.
}
