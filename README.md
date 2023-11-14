
<!-- README.md is generated from README.Rmd. Please edit that file -->

# igcop <img src="man/figures/igcop-240x278.png" align="right" height="150" />

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/igcop)](https://CRAN.R-project.org/package=igcop)
[![codecov](https://codecov.io/gh/vincenzocoia/igcop/branch/main/graph/badge.svg?token=RS6IYZFX5U)](https://app.codecov.io/gh/vincenzocoia/igcop)
[![License:
MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://cran.r-project.org/web/licenses/MIT)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![R-CMD-check](https://github.com/vincenzocoia/igcop/workflows/R-CMD-check/badge.svg)](https://github.com/vincenzocoia/igcop/actions)
[![R-CMD-check](https://github.com/vincenzocoia/igcop/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/vincenzocoia/igcop/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of igcop is to provide computational tools for the Integrated
Gamma (IG) and Integrated Gamma Limit (IGL) copula families.

## Installation

igcop is available on [CRAN](https://CRAN.R-project.org), and can be
installed by running

``` r
install.packages("igcop")
```

## Definition

The IG copula family is defined by parameters *θ* \> 0 and *α* \> 0,
with the IGL copula family obtained with *θ* → ∞. See the vignette for a
detailed definition.

Here are some contour plots of some normal scores copula densities.

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

The IG and IGL copula families are unique in that, when used in a
regression context, the conditional distribution of the response (the
2nd copula variable) has an Extreme Value Index that increases with the
predictor for an IG copula, and reduces a heavy-tailed response to a
light-tailed conditional distribution for an IGL copula. Specifically,
the Extreme Value Index of the 2\|1 distribution when Variable 2 has a
Pareto(1) marginal distribution is 0 for an IGL copula, and is
(1+*θ*(1−*u*))<sup>−1</sup> for an IG copula (Coia 2017).

## Usage

``` r
library(igcop)
```

This package piggybacks on the base R syntax for distributions, such as
`dnorm()` and `pexp()`, whose functions adopt the convention:

    <prefix><name>

For IG and IGL copulas:

-   `<prefix>` corresponds to one of:
    -   `p` for cdf,
    -   `d` for density (and `logd` for log density),
    -   `q` for quantile (for conditional distributions only), and
    -   `r` for random number generation (not supported for conditional
        distributions).
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

Here are some examples, starting with the density of an IG copula:

``` r
dig(0.3, 0.6, theta = 3, alpha = 2)
#> [1] 1.096211
```

Computations are vectorized over each argument. Here’s the cdf and
density of an IGL copula at different values:

``` r
u <- seq(0.1, 0.9, length.out = 9)
v <- seq(0.9, 0.5, length.out = 9)
pigl(u, v, alpha = 4)
#> [1] 0.1000000 0.2000000 0.2999711 0.3988536 0.4888134 0.5508382 0.5683229
#> [8] 0.5447653 0.4998090
digl(0.2, v, alpha = u)
#> [1] 0.8522462 0.8230206 0.8471676 0.8915708 0.9458967 1.0058156 1.0691273
#> [8] 1.1345476 1.2012456
```

It doesn’t make sense to talk about quantiles for a multivariate
distribution, so these are only defined for conditional distributions.

Here is an example of a distribution given the first variable (“2 given
1”). Note that the “2 given 1” distributions swap the `u` and `v`
arguments to better align with the conditioning, and you can either
explicitly include the `21` suffix or not.

``` r
qcondig(v, u, theta = 5, alpha = 3)
#> [1] 0.7435415 0.7228302 0.7121613 0.7073784 0.7056649 0.7039164 0.6972994
#> [8] 0.6777041 0.6356285
qcondig21(v, u, theta = 5, alpha = 3)
#> [1] 0.7435415 0.7228302 0.7121613 0.7073784 0.7056649 0.7039164 0.6972994
#> [8] 0.6777041 0.6356285
```

Here is the corresponding “1 given 2” distribution. Since this is less
common in regression scenarios, you have to explicitly add the `12`
prefix for “1 given 2.”

``` r
qcondig12(v, u, theta = 5, alpha = 3)
#> [1] 0.8896885 0.8114873 0.7297887 0.6598357 0.6097781 0.5811235 0.5749922
#> [8] 0.5976573 0.6689895
```

Generating 5 values from an IG copula:

``` r
set.seed(42)
rig(5, theta = 5, alpha = 4)
#> # A tibble: 5 × 2
#>       u     v
#>   <dbl> <dbl>
#> 1 0.915 0.598
#> 2 0.937 0.848
#> 3 0.286 0.134
#> 4 0.830 0.761
#> 5 0.642 0.770
```

## Developers

Besides the copula quantities described above, the generating functions
(as outlined in the vignette) are included in this package as internal
functions, and directly link to C++. The notation is:

-   *ψ*<sub>*α*</sub> is `igl_gen()`;
-   *κ*<sub>*α*</sub> is `igl_kappa()`;
-   *H*<sub>*α*</sub> is `interp_gen()`; and
-   *H*<sub>*κ*<sub>*α*</sub></sub> is `interp_kappa()`.

Related functions have the following suffixes:

-   `_inv`: function inverse.
-   `_D`: function derivative.
-   `_D1`: function derivative with respect to first argument.

There are three functions involved when linking to C:

1.  The R function (such as `igl_gen()`) recycles the arguments by
    passing them through the `formals_to()` function, which uses
    `vctrs::vec_recycle_common()`.
2.  These recycled arguments are passed to the corresponding R function
    with the `_vec` suffix, which passes these functions into C++ (via
    the infrastructure created by running `Rcpp::compileAttributes()`).
3.  The C++ functions that accept vector inputs have the `_vec` suffix.
    These functions loop along each entry, and feeds the scalar values
    into a C++ function for computation (either with the `_single`
    prefix, or the `_algo` prefix when the function contains a
    Newton-Raphson algorithm).

Map of dependencies among functions:

-   `igl_gen` : `pgamma`
-   `igl_gen_D` : `pgamma`
-   `igl_gen_inv_algo` : `qgamma` `igl_gen` `igl_gen_D`
-   `igl_gen_inv` : `igl_gen_inv_algo`
-   `interp_gen` : `igl_gen`
-   `interp_gen_D1` : `igl_gen`
-   `interp_gen_inv_algo` : `igl_gen_inv_algo` `interp_gen`
    `interp_gen_D1`
-   `interp_gen_inv` : `interp_gen_inv_algo`
-   `igl_kappa` : `pgamma`
-   `igl_kappa_D` : `dgamma`
-   `igl_kappa_inv` : `qgamma`
-   `interp_kappa` : `igl_kappa`
-   `interp_kappa_D1` : `igl_kappa` `igl_kappa_D`
-   `interp_kappa_inv_algo` : `igl_kappa_inv` `interp_kappa` `igl_kappa`
    `igl_kappa_D` `interp_kappa_inv`
-   `interp_kappa_inv` : `interp_kappa_inv_algo`
-   `pcondig21` : `interp_gen_inv` `interp_kappa`
-   `qcondig21` : `interp_kappa_inv` `interp_gen`
-   `qcondig12_algo` : `interp_gen_inv` `igl_gen` `igl_gen_D`
    `pcondig12`
-   `qcondig12` : `qcondig12_algo`
-   `pcondig12` : `interp_gen_inv` `interp_gen_D1`
-   `dig` : `interp_gen_inv` `interp_kappa_D1` `interp_gen_D1`
-   `logdig` : `interp_gen_inv` `igl_kappa` `igl_kappa_D` `igl_gen`
    `igl_gen_D`
-   `pig` : `interp_gen_inv`
-   `rig` : `qcondig21`
-   `qcondigl21` : `igl_kappa_inv`
-   `pcondigl21` : `igl_gen_inv` `igl_kappa`
-   `pcondigl12` : `igl_gen_inv` `igl_gen_D`
-   `qcondigl12` : `igl_gen_inv` `pgamma` `qgamma`
-   `digl` : `igl_gen_inv` `igl_kappa_D` `igl_gen_D`
-   `pigl` : `igl_gen_inv` `igl_gen`
-   `rigl` : `qcondigl21`

## Attributions

Package developed and maintained by [Vincenzo
Coia](https://vincenzocoia.com/), with thanks to Harry Joe for his help
converting the Newton Raphson algorithms and related functions to C
(originally coded in R in igcop Version 0.2.0).

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-coia2017" class="csl-entry">

Coia, Vincenzo. 2017. “Forecasting of Nonlinear Extreme Quantiles Using
Copula Models.” PhD Dissertation; The University of British Columbia.

</div>

</div>
