context("Test the incomplete gamma function")

test_that("incomplete gamma vectorization works", {
  set.seed(1)
  s <- rexp(10)
  x <- rexp(10)
  expect_identical(
    mapply(igamma, s, x),
    igamma(s, x)
  )
  for (s_ in s) {
    expect_identical(
      vapply(x, function(x_) igamma(s_, x_), FUN.VALUE = numeric(1L)),
      igamma(s_, x)
    )
  }
  for (x_ in x) {
    expect_identical(
      vapply(s, function(s_) igamma(s_, x_), FUN.VALUE = numeric(1L)),
      igamma(s, x_)
    )
  }
})


test_that("incomplete gamma aligns with gamma", {
  set.seed(1)
  s <- rexp(10)
  expect_identical(
    gamma(s),
    igamma(s, 0)
  )
})

test_that("incomplete gamma aligns with numeric integral", {
  set.seed(1)
  s <- rexp(10)
  x <- rexp(10)
  rel.tol <- 1e-7
  for (s_ in s) {
    for (x_ in x) {
      f <- function(t) {
        t ^ (s_ - 1) * exp(-t)
      }
      int_f <- integrate(f, x_, Inf, rel.tol = rel.tol)
      expect_equal(
        int_f$value,
        igamma(s_, x_),
        tolerance = rel.tol
      )
    }
  }
})
