context("Check that negative parameters throw an error.")

test_that("negative IG theta throws an error.", {
  expect_error(dig(0.5, 0.5, theta = -10:10, alpha = 1))
  expect_error(logdig(0.5, 0.5, theta = -10:10, alpha = 1))
  expect_error(pig(0.5, 0.5, theta = -10:10, alpha = 1))
  expect_error(rig(10, theta = -10:10, alpha = 1))
  expect_error(pcondig21(0.5, 0.5, theta = -10:10, alpha = 1))
  expect_error(pcondig12(0.5, 0.5, theta = -10:10, alpha = 1))
  expect_error(qcondig21(0.5, 0.5, theta = -10:10, alpha = 1))
  expect_error(qcondig12(0.5, 0.5, theta = -10:10, alpha = 1))
})

test_that("negative IG alpha throws an error.", {
  expect_error(dig(0.5, 0.5, alpha = -10:10, theta = 1))
  expect_error(logdig(0.5, 0.5, alpha = -10:10, theta = 1))
  expect_error(pig(0.5, 0.5, alpha = -10:10, theta = 1))
  expect_error(rig(10, alpha = -10:10, theta = 1))
  expect_error(pcondig21(0.5, 0.5, alpha = -10:10, theta = 1))
  expect_error(pcondig12(0.5, 0.5, alpha = -10:10, theta = 1))
  expect_error(qcondig21(0.5, 0.5, alpha = -10:10, theta = 1))
  expect_error(qcondig12(0.5, 0.5, alpha = -10:10, theta = 1))
})

test_that("negative IGL alpha throws an error.", {
  expect_error(digl(0.5, 0.5, alpha = -10:10))
  expect_error(logdigl(0.5, 0.5, alpha = -10:10))
  expect_error(pigl(0.5, 0.5, alpha = -10:10))
  expect_error(rigl(10, alpha = -10:10))
  expect_error(pcondigl21(0.5, 0.5, alpha = -10:10))
  expect_error(pcondigl12(0.5, 0.5, alpha = -10:10))
  expect_error(qcondigl21(0.5, 0.5, alpha = -10:10))
  expect_error(qcondigl12(0.5, 0.5, alpha = -10:10))
})
