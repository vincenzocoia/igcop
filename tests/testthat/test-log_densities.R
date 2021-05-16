context("Check the simplified log densities against the densities")

test_that("exp(log IG density) returns the IG density", {
  for (cpar_ in .cpar) {
    expect_equal(
      dig(.u, .v, cpar = cpar_),
      exp(logdig(.u, .v, cpar = cpar_))
    )
  }
})

test_that("exp(log IGL density) returns the IGL density", {
  for (alpha_ in .alpha) {
    expect_equal(
      digl(.u, .v, cpar = alpha_),
      exp(logdigl(.u, .v, cpar = alpha_))
    )
  }
})
