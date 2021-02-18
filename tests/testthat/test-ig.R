library(testthat)
library(tidyverse)
library(copsupp)

# Tests from Harry Joe
# check interp_gen_inv in vectorized mode
# check pigcop, digcop numerical derivatives
# check pcondigcop21 and pcondigcop12
# check Spearman rho via rhoS in CopulaModel

cpar <- list(
    c(2.5, 13.3), # nearly independent
    c(2.5, 5.3), # weak dependence spear=.040
    c(2.5, 3.3), # spear=.145
    #c(2.5, 1.3), # spear=.244
    #c(2.5, 1.1), # unreliable (for small k)
    c(5.5,4.1),  # spear=.231
    c(8.5,4.1),  # spear=.330
    c(15.5,4.1), # spear=.462
    c(25.5,4.1), # spear=.553
    c(45.5,4.1), # spear=.629
    c(75.5,4.1), # spear=.674
    c(75.5,8.1), # spear=.661
    c(155.5,4.1), # spear=.714
    c(300.5,4.1), # spear=.733
    c(600.5,4.1), # spear=.743
    c(600.5,8.1), # spear=.829
    c(600.5,15.1), # spear=.853
    c(600.5,35.1), # spear=.811
    c(1200.5,15.1), # spear=.883
    c(2400.5,15.1), # spear=.899
    c(2400.5,35.1), # spear=.918
    c(2400.5,55.1), # spear=.??? unstable for rhoS
    c(2400.5,45.1), # spear=.914
    c(4800.5,35.1), # spear=.937
    c(9600.5,35.1) # spear=.947  # try inference with this parameter vector
)

u <- c(seq(.1, .9, .2), c(0.2, 0.6, 0.9))
uu <- rep(list(u), length(cpar))
v <- c(seq(.1, .9, .2), c(0.8, 0.2, 0.4))
vv <- rep(list(v), length(cpar))

## ------ TEST THE HELPER FUNCTIONS ------

u2 <- purrr::map(cpar, ~ igl_gen(igl_gen_inv(u, .x[2]), .x[2]))
tibble(cpar = cpar,
       u_original = uu,
       u_through_inverse = u2) %>%
    unnest_wider(cpar, names_sep = "_") %>%
    unnest(starts_with("u_")) %>%
    DT::datatable()
test_that("IGL generator inverse works", {
    expect_equal(u2, uu)
})

u2 <- purrr::map(cpar, ~ {
    interp_gen(interp_gen_inv(u, .x[1], .x[2]), .x[1], .x[2])
})
tibble(cpar = cpar,
       u_original = uu,
       u_through_inverse = u2) %>%
    unnest_wider(cpar, names_sep = "_") %>%
    unnest(starts_with("u_")) %>%
    DT::datatable()
test_that("IG generator inverse works", {
    expect_equal(u2, uu)
})

v2 <- purrr::map(cpar, ~ {
    theta <- .x[1]
    k <- .x[2]
    eta <- theta * (1 - u)
    igcond(igcondinv(v, k, eta), k, eta)
})
tibble(cpar = cpar,
       v_original = vv,
       v_through_inverse = v2) %>%
    unnest_wider(cpar, names_sep = "_") %>%
    unnest(starts_with("v_")) %>%
    DT::datatable()
test_that("igcond inverse works", {
    expect_equal(v2, vv)
})

## ------- COMPATIBILITY-BASED CHECKS ------

pcondigcop_2 <- function(v, u, cpar) {
    theta <- cpar[1]
    k     <- cpar[2]
    Hkinv <- interp_gen_inv(1 - v, theta, k)
    inner <- theta * (1 - u) * log(Hkinv)
    1 - pgamma(inner, k - 1, lower.tail = FALSE) / Hkinv
}
one <- map(cpar, ~ pcondigcop(u, v, .x))
two <- map(cpar, ~ pcondigcop_2(u, v, .x))
tibble(cpar = cpar,
       pcondigcop_original = one,
       pcondigcop_alternative = two) %>%
    unnest_wider(cpar, names_sep = "_") %>%
    unnest(starts_with("pcondigcop")) %>%
    DT::datatable()
test_that("Unsimplified calculation of pcondigcop matches the one used in the function", {
    expect_equal(one, two)
})


v2 <- map(cpar, ~ qcondigcop(pcondigcop(v, u, cpar = .x), u, cpar = .x))
tibble(cpar = cpar,
       v_original = vv,
       v_through_the_inverse = v2) %>%
    unnest_wider(cpar, names_sep = "_") %>%
    unnest(starts_with("v_")) %>%
    DT::datatable()
test_that("Quantile function is inverse of pcond", {
    expect_equal(v2, vv)
})

## ------ VECTORIZATION -------

#' Calculate pigcop looping across u and v
pigcop_looped <- function(u, v, cpar) {
    cdf <- numeric()
    for(i in 1:length(u)) {
        cdf[i] <- pigcop(u[i], v[i], cpar)
    }
    cdf
}
cdf1 <- purrr::map(cpar, ~ pigcop(u, v, .x))
cdf2 <- purrr::map(cpar, ~ pigcop_looped(u, v, .x))
tibble(cdf_vectorized = cdf1,
       cdf_looped = cdf2,
       cpar = cpar) %>%
    unnest(starts_with("cdf")) %>%
    unnest_wider(cpar, names_sep = "_") %>%
    DT::datatable()
test_that("Vectorization works for pigcop", {
    expect_equal(cdf1, cdf2)
})

#' Calculate digcop looping across u and v
digcop_looped <- function(u, v, cpar) {
    pdf <- numeric()
    for(i in 1:length(u)) {
        pdf[i] <- digcop(u[i], v[i], cpar)
    }
    pdf
}
pdf1 <- purrr::map(cpar, ~ digcop(u, v, .x))
pdf2 <- purrr::map(cpar, ~ digcop_looped(u, v, .x))
tibble(pdf_vectorized = pdf1,
       pdf_looped = pdf2,
       cpar = cpar) %>%
    unnest(starts_with("pdf")) %>%
    unnest_wider(cpar, names_sep = "_") %>%
    DT::datatable()
test_that("Vectorization works for digcop", {
    expect_equal(pdf1, pdf2)
})

## ------ NUMERICAL DERIVATIVES ------

#' Calculate numerical derivitive
digcop_numerical <- function(u, v, cpar, eps = 1.e-6) {
    cdf11 <- pigcop(u, v, cpar)
    cdf22 <- pigcop(u + eps, v + eps, cpar)
    cdf21 <- pigcop(u + eps, v, cpar)
    cdf12 <- pigcop(u, v + eps, cpar)
    (cdf22 + cdf11 - cdf12 - cdf21) / eps ^ 2
}
pdf1 <- purrr::map(cpar, ~ digcop_numerical(u, v, .x))
pdf2 <- purrr::map(cpar, ~ digcop(u, v, .x))
tibble(pdf_numeric = pdf1,
       pdf_original = pdf2,
       cpar = cpar) %>%
    unnest(starts_with("pdf")) %>%
    unnest_wider(cpar, names_sep = "_") %>%
    DT::datatable()
test_that("Density is the same as numerical derivitive of cdf", {
    expect_equal(pdf1, pdf2, tolerance = 0.009)
})


#' Calculate numerical derivitive
pcondigcop21_numerical <- function(v, u, cpar, eps = 1.e-8) {
    cdf11 <- pigcop(u, v, cpar)
    cdf21 <- pigcop(u + eps, v, cpar)
    (cdf21 - cdf11) / eps
}
pcond1 <- purrr::map(cpar, ~ pcondigcop21_numerical(v, u, .x))
pcond2 <- purrr::map(cpar, ~ pcondigcop21(v, u, .x))
tibble(pcond_numeric = pcond1,
       pcond_original = pcond2,
       cpar = cpar) %>%
    unnest(starts_with("pcond")) %>%
    unnest_wider(cpar, names_sep = "_") %>%
    DT::datatable()
test_that("2|1 is the same as numerical derivitive of cdf", {
    expect_equal(pcond1, pcond2, tolerance = 1e-7)
})

#' Calculate numerical derivitive
digcop_numerical_from_2g1 <- function(u, v, cpar, eps = 1.e-7) {
    conda <- pcondigcop21(v, u, cpar)
    condb <- pcondigcop21(v + eps, u, cpar)
    (condb - conda) / eps
}
pdf1 <- purrr::map(cpar, ~ digcop_numerical_from_2g1(u, v, .x))
pdf2 <- purrr::map(cpar, ~ digcop(u, v, .x))
tibble(pcond_numeric = pdf1,
       pcond_original = pdf2,
       cpar = cpar) %>%
    unnest(starts_with("pcond")) %>%
    unnest_wider(cpar, names_sep = "_") %>%
    DT::datatable()
test_that("density is the same as the numerical derivitive of 2|1", {
    expect_equal(pdf1, pdf2, tolerance = 1e-5)
})

## ------ SIMULATED DATA -------

set.seed(1)
n <- 1000  # Sample size

dats <- map(cpar, ~ {
    mat <- rigcop(n, cpar = .x)
    colnames(mat) <- c("u", "v")
    as_tibble(mat)
})

#' compute non-parametric cdf
#' - "dat" should be a data frame with names "u" and "v".
#' - "at_u" and "at_v" should be single numerics.
dat_to_cdf <- function(dat, at_u, at_v) {
    dat %>%
        filter(u <= at_u, v <= at_v) %>%
        nrow() / nrow(dat)
}
dat_to_count <- function(dat, at_u, at_v) {
    dat %>%
        filter(u <= at_u, v <= at_v) %>%
        nrow()
}

#cdf_estimate ~ N(cdf, cdf * (1 - cdf) / n )
se <- function(cdf, n) {
    sqrt(cdf * (1 - cdf) / n)
}
eval <- expand_grid(nesting(cpar = cpar, dat = dats),
                    u = 1:4 / 5, v = 1:4 / 5) %>%
    transmute(cdf_simu = pmap_dbl(list(dat, u, v), dat_to_cdf),
              count_simu = pmap_dbl(list(dat, u, v), dat_to_count),
              cdf_true = pmap_dbl(list(u, v, cpar), pigcop),
              se_true = se(cdf_true, n = n),
              norm_score = (cdf_simu - cdf_true) / se_true,
              u_score = pnorm(norm_score),
              u_from_binomial = pbinom(count_simu, size = n, prob = cdf_true))
ggplot(eval, aes(cdf_true, cdf_simu)) +
    geom_point(alpha = 0.2) +
    stat_function(fun = function(x) x + 1.96 * se(x, n = n)) +
    stat_function(fun = function(x) x - 1.96 * se(x, n = n)) +
    labs(x = "actual cdf values",
         y = "nonparametric estimations",
         title = "95% prediction interval from\nnormal-approximated sampling distribution")
ggplot(eval, aes(rank(u_score) / length(u_score), u_score)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0) +
    labs(x = "actual",
         y = "calculated",
         title = "PP plot of estimates using binomial distribution")
ggplot(eval, aes(norm_score)) +
    geom_histogram(aes(y = ..density..)) +
    stat_function(fun = dnorm, colour = "blue") +
    ggtitle("Histogram of normalized estimates.\nShould be N(0, 1).")


#======================================================================

# After fix to IG copula density,
# check with some simulated data.

# If numerical min loglik with nlm is attempted,
# use log(theta) and log(k) as the parameters
# It might be better to implement the negative log-likelihood as a loop
# rather than in vectorized mode.
# Try first on a subset of deseasonailzed data or simulated data.




