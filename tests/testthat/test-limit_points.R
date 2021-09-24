context("Check limit points")

test_that("Check boundaries of cdf", {
  u <- c(0, 0.7, 1, 0.7)
  v <- c(0.7, 0, 0.7, 1)
  expect_equal(pig(u, v, theta = 1, alpha = 1), c(0, 0, 0.7, 0.7))
  expect_equal(pigl(u, v, alpha = 1), c(0, 0, 0.7, 0.7))
})

test_that("Limit of igl_gen_D at 0 (from the right) works.", {
  avec <- c(0.2, 0.5, 1, 1, 1, 1.8, 100)
  at_zero <- igl_gen_D(0, alpha = avec)
  near_zero <- igl_gen_D(1e-15, alpha = avec)
  near_zero <- vctrs::vec_assign(near_zero, near_zero < -1e6, -Inf)
  expect_equal(at_zero, near_zero, tolerance = 1e-4)
})

test_that("Limit of interp_gen_D1 at 0 (from the right) works, via pcondig12", {
  cpar <- .cpar[-(4:5)] # alpha < 0 results in Inf/Inf.
  for (cpar_ in cpar) {
    at_zero <- pcondig12(.u, 0, theta = cpar_[1], alpha = cpar_[2])
    near_zero <- pcondig12(.u, 1e-15, theta = cpar_[1], alpha = cpar_[2])
    expect_equal(at_zero, near_zero, tolerance = 1e-4)
  }
})
