test_that("Limit of interp_gen_D1 at 0 (from the right) works.", {
  avec <- c(0.2, 0.5, 1, 1, 1, 1.8, 100)
  evec <- c(1, 0.5, 0.25, 1, 2, 1, 4)
  at_zero <- interp_gen_D1(0, eta = evec, alpha = avec)
  near_zero <- interp_gen_D1(1e-15, eta = evec, alpha = avec)
  near_zero <- vctrs::vec_assign(near_zero, near_zero < -1e6, -Inf)
  expect_equal(at_zero, near_zero, tolerance = 1e-4)
})

test_that("Limit of igl_gen_D at 0 (from the right) works.", {
  avec <- c(0.2, 0.5, 1, 1, 1, 1.8, 100)
  at_zero <- igl_gen_D(0, alpha = avec)
  near_zero <- igl_gen_D(1e-15, alpha = avec)
  near_zero <- vctrs::vec_assign(near_zero, near_zero < -1e6, -Inf)
  expect_equal(at_zero, near_zero, tolerance = 1e-4)
})
