context("Check the simplified log densities against the densities")

test_that("exp(log IG density) returns the IG density", {
  for (cpar_ in .cpar) {
    expect_equal(
      dig(.u, .v, theta = cpar_[1], alpha = cpar_[2]),
      exp(logdig(.u, .v, theta = cpar_[1], alpha = cpar_[2]))
    )
  }
})

test_that("exp(log IGL density) returns the IGL density", {
  for (alpha_ in .alpha) {
    expect_equal(
      digl(.u, .v, alpha = alpha_),
      exp(logdigl(.u, .v, alpha = alpha_))
    )
  }
})
