#' igcop: Computational Tools for the IG and IGL Copula Families
#'
#' Compute distributional quantities for an
#' Integrated Gamma (IG) or IG Limit (IGL) copula, such
#' as a cdf and density, along with conditional quantities
#' such as the cdf, quantiles, and densities. Generate
#' data from a copula.
#'
#' @section Usage:
#' Access copula quantities by starting with the `p`, `d`, `q`, or `r`
#' prefixes, followed by the copula name -- either `ig` or `igl`, or
#' their conditional versions, `condig` or `condigl`.
#' @docType package
#' @name igcop
#' @importFrom Rcpp evalCpp
#' @useDynLib igcop, .registration = TRUE
"_PACKAGE"
