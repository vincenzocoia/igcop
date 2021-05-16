context("Check the simplified log densities against the densities")

test_that("exp(log IG density) returns the IG density", {
  for (cpar_ in .cpar) {
    expect_equal(
      digcop(.u, .v, cpar = cpar_),
      exp(logdigcop(.u, .v, cpar = cpar_))
    )
  }
})

test_that("exp(log IGL density) returns the IGL density", {
  for (alpha_ in .alpha) {
    expect_equal(
      diglcop(.u, .v, cpar = alpha_),
      exp(logdiglcop(.u, .v, cpar = alpha_))
    )
  }
})
