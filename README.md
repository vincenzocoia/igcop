
<!-- README.md is generated from README.Rmd. Please edit that file -->

# igcop <img src="man/figures/igcop-240x278.png" align="right" height="150" />

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

The IG copula family is defined by parameters *θ* &gt; 0 and
*a**l**p**h**a* &gt; 0, with the IGL copula family obtained with *θ* → ∞
(so, only has one parameter *a**l**p**h**a* &gt; 0).

The IG and IGL copula families are unique in that they are not
permutation symmetric. Also, the family of conditional distributions
given the first variable has a non-constant extreme value index (EVI)
across the first variable for the IG family, and is 0 for the IGL
family. This means that, when used in a regression context, the
conditional distribution of the response has an EVI that changes with
the predictor for an IG copula, or reduces a heavy-tailed response to a
light-tailed conditional distribution.

## Usage

This package piggybacks on the base R syntax for distributions, whose
functions adopt the convention:

    <prefix><name>

For IG and IGL copulas:

-   `<prefix>` corresponds to one of:
    -   `r` for random number generation (currently not supported for
        conditional distributions),
    -   `p` for cdf,
    -   `d` for density, and
    -   `q` for quantile (for conditional distributions only).
-   `<name>` corresponds to the possible names:
    -   `ig` and `igl` correspond to an IG copula and IGL copula,
        respectively.
    -   `condig12` and `condigl12` correspond to a conditional
        distribution of the first variable given the second, of an IG
        copula and IGL copula respectively.
    -   `condig21` and `condigl21` correspond to a conditional
        distribution of the second variable given the first, of an IG
        copula and IGL copula respectively (also available as `condig`
        and `condigl` to match the syntax of the
        [CopulaModel](https://github.com/vincenzocoia/CopulaModel)
        package).

All of these functions have a `cpar` argument expecting the value of the
copula parameters. For an IG copula, this is `c(theta, alpha)`, and just
`alpha` for an IGL copula.

Here are some examples, starting with evaluating the density of an IG
copula at (0.3, 0.6):

``` r
library(igcop)
dig(0.3, 0.6, cpar = c(3, 2))
#> [1] 1.096211
```

Computations are vectorized over both `u` and `v` (first and second
variables). Here’s the cdf and density of an IGL copula at different
values:

``` r
u <- seq(0.1, 0.9, length.out = 9)
v <- seq(0.9, 0.5, length.out = 9)
pigl(u, v, cpar = 4)
#> [1] 0.1000000 0.2000000 0.2999711 0.3988536 0.4888134 0.5508382 0.5683229
#> [8] 0.5447653 0.4998090
digl(0.2, v, cpar = 4)
#> [1] 2.213223e-08 1.250384e-04 6.146088e-03 4.943650e-02 1.683705e-01
#> [6] 3.617879e-01 5.964016e-01 8.377819e-01 1.064300e+00
```

It doesn’t make sense to talk about quantiles for a multivariate
distribution, so this is only defined for conditional distributions.

Here is an example of a distribution given the first variable (“2 given
1”). Note that the “2 given 1” distributions swap the `u` and `v`
arguments to better align with the conditioning, and you can either
explicitly include the `21` suffix or not.

``` r
qcondig(v, u, cpar = c(5, 3))
#> [1] 0.7435415 0.7228302 0.7121613 0.7073784 0.7056649 0.7039164 0.6972994
#> [8] 0.6777041 0.6356285
qcondig21(v, u, cpar = c(5, 3))
#> [1] 0.7435415 0.7228302 0.7121613 0.7073784 0.7056649 0.7039164 0.6972994
#> [8] 0.6777041 0.6356285
```

Here is the corresponding “1 given 2” distribution. Since this is less
common in regression scenarios, you have to explicitly add the `12`
prefix for “1 given 2”.

``` r
qcondig12(v, u, cpar = c(5, 3))
#> [1] 0.8896885 0.8114873 0.7297887 0.6598357 0.6097781 0.5811235 0.5749922
#> [8] 0.5976573 0.6689895
```

Generating 5 values from an IG copula:

``` r
rig(5, cpar = c(5, 4))
#> # A tibble: 5 x 2
#>       u     v
#>   <dbl> <dbl>
#> 1 0.682 0.414
#> 2 0.183 0.351
#> 3 0.754 0.138
#> 4 0.216 0.146
#> 5 0.136 0.626
```
