test_that("integrated gamma vectorization works", {
  set.seed(1)
  alpha <- rexp(10)
  x <- rexp(10)
  expect_identical(
    purrr::map2_dbl(alpha, x, igamma),
    igamma(alpha, x)
  )
  for (alpha_ in alpha) {
    expect_identical(
      purrr::map_dbl(x, igamma, alpha = alpha_),
      igamma(alpha_, x)
    )
  }
  for (x_ in x) {
    expect_identical(
      purrr::map_dbl(alpha, igamma, x = x_),
      igamma(alpha, x_)
    )
  }
})


test_that("integrated gamma aligns with gamma", {
  set.seed(1)
  alpha <- rexp(10)
  expect_identical(
    gamma(alpha),
    igamma(alpha, 0)
  )
})

test_that("integrated gamma aligns with numeric integral", {
  set.seed(1)
  alpha <- rexp(10)
  x <- rexp(10)
  rel.tol <- 1e-7
  for (alpha_ in alpha) {
    for (x_ in x) {
      f <- function(t) {
        t ^ (alpha_ - 1) * exp(-t)
      }
      int_f <- integrate(f, x_, Inf, rel.tol = rel.tol)
      expect_equal(
        int_f$value,
        igamma(alpha_, x_),
        tolerance = rel.tol
      )
    }
  }
})
