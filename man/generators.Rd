% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vectorize.R
\name{interp_gen_inv}
\alias{interp_gen_inv}
\alias{interp_kappa}
\alias{interp_kappa_inv}
\alias{interp_gen}
\alias{igl_kappa}
\alias{igl_kappa_D}
\alias{igl_kappa_inv}
\alias{igl_gen}
\alias{igl_gen_D}
\alias{igl_gen_inv}
\title{IG/IGL Generators and Related Functions}
\usage{
interp_gen_inv(p, eta, alpha)

interp_kappa(x, eta, alpha)

interp_kappa_inv(p, eta, alpha)

interp_gen(x, eta, alpha)

igl_kappa(x, alpha)

igl_kappa_D(x, alpha)

igl_kappa_inv(p, alpha)

igl_gen(x, alpha)

igl_gen_D(x, alpha)

igl_gen_inv(p, alpha)
}
\arguments{
\item{p}{Function inverse argument. Vector of values between 0 and 1.}

\item{eta, alpha}{Function parameters. Vector of positive values.}

\item{x}{Function argument. Vector of non-negative values.}
}
\value{
The function values, as a vector.
}
\description{
These are the psi, H, and kappa functions
of the IG and IGL copula families.
}
\details{
Kappa function and its relatives have prefix \code{igl_kappa};
Psi function and its relatives have prefix \code{igl_gen};
Interpolating function H with either kappa or psi has
\code{igl} prefix replaced with \code{interp}. Relatives of these functions:
suffix \code{inv} indicates inverse; suffix \code{D} represents function
derivative, and \code{D1} derivative with respect to the first argument.
. Suffix \verb{_vec} indicates that the entries must be vectors of
the same length; \verb{_single} means entries must be
scalars.
}
\note{
Inputs must be recyclable via \code{vctrs::vec_recycle_common()}.
}
