context("Check that the inverses of the helper functions (psi, kappa, H) work.")

test_that("igl_gen_inv is the inverse of igl_gen", {
  for (alpha_ in .alpha){
    u2 <- igl_gen(
      igl_gen_inv(.u, alpha = alpha_),
      alpha = alpha_
    )
    expect_equal(.u, u2)
  }
})

test_that("interp_gen_inv is the inverse of interp_gen", {
  for (cpar_ in .cpar){
    eta <- cpar_[1]
    alpha <- cpar_[2]
    u2 <- interp_gen(
      interp_gen_inv(.u, eta = eta, alpha = alpha),
      eta = eta, alpha = alpha
    )
    expect_equal(.u, u2)
  }
})

test_that("interp_kappa_inv is the inverse of interp_kappa", {
  for (cpar_ in .cpar){
    eta <- cpar_[1]
    alpha <- cpar_[2]
    u2 <- interp_kappa(
      interp_kappa_inv(.u, eta = eta, alpha = alpha),
      eta = eta, alpha = alpha
    )
    expect_equal(.u, u2)
  }
})
