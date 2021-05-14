test_that("log IG density works", {
  for (cpar_ in .cpar) {
    expect_equal(
      digcop(.u, .v, cpar_),
      exp(logdigcop(.u, .v, cpar_))
    )
  }
})

test_that("log IGL density works", {
  for (alpha_ in .alpha) {
    expect_equal(
      diglcop(.u, .v, alpha_),
      exp(logdiglcop(.u, .v, alpha_))
    )
  }
})
