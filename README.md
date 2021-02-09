
<!-- README.md is generated from README.Rmd. Please edit that file -->

# igcop

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/igcop)](https://CRAN.R-project.org/package=igcop)
[![Codecov test
coverage](https://codecov.io/gh/vincenzocoia/igcop/branch/master/graph/badge.svg)](https://codecov.io/gh/vincenzocoia/igcop?branch=master)
[![License:
MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://cran.r-project.org/web/licenses/MIT)
<!-- badges: end -->

The goal of igcop is to provide tools for computing on the Integrated
Gamma (IG) and Integrated Gamma Limit (IGL) copula families.

## Installation

igcop is not yet available on [CRAN](https://CRAN.R-project.org), but
can be downloaded from this repository using devtools. Just execute this
line of code in an R instance, after ensuring you have the devtools R
package installed:

``` r
devtools::install_github("vincenzocoia/igcop")
```

## Definition

The IG copula family is defined by parameters \(\theta > 0\) and
\(k > 1\), although computations are problematic for \(k < 2\), although
mostly just for \(k\) close to 1.

…

What makes the IG copula interesting is in regression analysis, where
the response is the second variable. If the response is heavy-tailed,
and is linked to a predictor via an IG copula, then its conditional
distributions have lighter tails with a non-constant extreme value index
across the predictor space. The IGL copula is interesting in a similar
light, except its conditional distributions are all light-tailed (\!) –
meaning that the predictor is solely responsible for the heavy tail of
the response variable.

## Usage

This package piggybacks on the base R syntax for distributions, whose
functions adopt the convention:

    <prefix><name>

For IG and IGL copulas:

  - `<prefix>` corresponds to one of:
      - `r` for random number generation (currently not supported for
        conditional distributions),
      - `p` for cdf,
      - `d` for density, and
      - `q` for quantile (for conditional distributions only).
  - `<name>` corresponds to the possible names:
      - `igcop` and `iglcop` correspond to an IG copula and IGL copula,
        respectively.
      - `condigcop12` and `condiglcop12` correspond to a conditional
        distribution of the first variable given the second, of an IG
        copula and IGL copula respectively.
      - `condigcop21` and `condiglcop21` correspond to a conditional
        distribution of the second variable given the first, of an IG
        copula and IGL copula respectively (also available as
        `condigcop` and `condiglcop` to match the syntax of the
        CopulaModel package).

All of these functions have a `cpar` argument expecting the value of the
copula parameters. For an IG copula, this is `c(theta, k)`, and just `k`
for an IGL copula.

Here are some examples, starting with evaluating the density of an IG
copula at (0.3, 0.6):

``` r
library(igcop)
digcop(0.3, 0.6, cpar = c(3, 2))
#> [1] 1.047923
```

Computations are vectorized over both `u` and `v` (first and second
variables). Here’s the cdf and density of an IGL copula at different
values:

``` r
u <- seq(0.1, 0.9, length.out = 9)
v <- seq(0.9, 0.5, length.out = 9)
piglcop(u, v, cpar = 4)
#> [1] 0.1000000 0.1999991 0.2998577 0.3973520 0.4830642 0.5410399 0.5597112
#> [8] 0.5411682 0.4994190
diglcop(0.2, v, cpar = 4)
#> [1] 2.609373e-06 1.536488e-03 2.655176e-02 1.202567e-01 2.905683e-01
#> [6] 5.044099e-01 7.261548e-01 9.344919e-01 1.120663e+00
```

It doesn’t make sense to talk about quantiles for a multivariate
distribution, so this is only defined for conditional distributions.
Note that the “2 given 1” distributions swap the `u` and `v` arguments
to better align with the conditioning.

``` r
qcondigcop(v, u, cpar = c(5, 3))
#> [1] 0.7241694 0.7000652 0.6876793 0.6832910 0.6850545 0.6916766 0.7015285
#> [8] 0.7111306 0.7112040
```

Generating 5 values from an IG copula:

``` r
rigcop(5, cpar = c(5, 4))
#> # A tibble: 5 x 2
#>        u      v
#>    <dbl>  <dbl>
#> 1 0.604  0.382 
#> 2 0.733  0.973 
#> 3 0.469  0.0380
#> 4 0.0111 0.508 
#> 5 0.635  0.261
```
