context("Check the handling of missing and improper inputs: IG copula.")

u_na_nan <- c(0.5, NA_real_, NaN)

test_that("pig() handles missingness appropriately", {
  args <- expand.grid(u = u_na_nan,
                      v = u_na_nan,
                      theta = u_na_nan,
                      alpha = u_na_nan)
  eval <- do.call(pig, args = args)
  expect_true(all(is.na(eval[-1])))
  expect_false(is.na(eval[1]))
  expect_length(
    c(pig(numeric(0L), 0.5, theta = 0.5, alpha = 0.5),
      pig(0.5, numeric(0L), theta = 0.5, alpha = 0.5),
      pig(0.5, 0.5, theta = numeric(0L), alpha = 0.5),
      pig(0.5, 0.5, theta = 0.5, alpha = numeric(0L))),
    0L
  )
})

test_that("dig() handles missingness appropriately", {
  args <- expand.grid(u = u_na_nan,
                      v = u_na_nan,
                      theta = u_na_nan,
                      alpha = u_na_nan)
  eval <- do.call(dig, args = args)
  expect_true(all(is.na(eval[-1])))
  expect_false(is.na(eval[1]))
  expect_length(
    c(dig(numeric(0L), 0.5, theta = 0.5, alpha = 0.5),
      dig(0.5, numeric(0L), theta = 0.5, alpha = 0.5),
      dig(0.5, 0.5, theta = numeric(0L), alpha = 0.5),
      dig(0.5, 0.5, theta = 0.5, alpha = numeric(0L))),
    0L
  )
})

test_that("pcondig21() handles missingness appropriately", {
  args <- expand.grid(u = u_na_nan,
                      v = u_na_nan,
                      theta = u_na_nan,
                      alpha = u_na_nan)
  eval <- do.call(pcondig21, args = args)
  expect_true(all(is.na(eval[-1])))
  expect_false(is.na(eval[1]))
  expect_length(
    c(pcondig21(numeric(0L), 0.5, theta = 0.5, alpha = 0.5),
      pcondig21(0.5, numeric(0L), theta = 0.5, alpha = 0.5),
      pcondig21(0.5, 0.5, theta = numeric(0L), alpha = 0.5),
      pcondig21(0.5, 0.5, theta = 0.5, alpha = numeric(0L))),
    0L
  )
})

test_that("pcondig12() handles missingness appropriately", {
  args <- expand.grid(u = u_na_nan,
                      v = u_na_nan,
                      theta = u_na_nan,
                      alpha = u_na_nan)
  eval <- do.call(pcondig12, args = args)
  expect_true(all(is.na(eval[-1])))
  expect_false(is.na(eval[1]))
  expect_length(
    c(pcondig12(numeric(0L), 0.5, theta = 0.5, alpha = 0.5),
      pcondig12(0.5, numeric(0L), theta = 0.5, alpha = 0.5),
      pcondig12(0.5, 0.5, theta = numeric(0L), alpha = 0.5),
      pcondig12(0.5, 0.5, theta = 0.5, alpha = numeric(0L))),
    0L
  )
})

test_that("qcondig12() handles missingness appropriately", {
  args <- expand.grid(p = u_na_nan,
                      v = u_na_nan,
                      theta = u_na_nan,
                      alpha = u_na_nan)
  eval <- do.call(qcondig12, args = args)
  expect_true(all(is.na(eval[-1])))
  expect_false(is.na(eval[1]))
  expect_length(
    c(qcondig12(numeric(0L), 0.5, theta = 0.5, alpha = 0.5),
      qcondig12(0.5, numeric(0L), theta = 0.5, alpha = 0.5),
      qcondig12(0.5, 0.5, theta = numeric(0L), alpha = 0.5),
      qcondig12(0.5, 0.5, theta = 0.5, alpha = numeric(0L))),
    0L
  )
})

test_that("qcondig21() handles missingness appropriately", {
  args <- expand.grid(p = u_na_nan,
                      u = u_na_nan,
                      theta = u_na_nan,
                      alpha = u_na_nan)
  eval <- do.call(qcondig21, args = args)
  expect_true(all(is.na(eval[-1])))
  expect_false(is.na(eval[1]))
  expect_length(
    c(qcondig21(numeric(0L), 0.5, theta = 0.5, alpha = 0.5),
      qcondig21(0.5, numeric(0L), theta = 0.5, alpha = 0.5),
      qcondig21(0.5, 0.5, theta = numeric(0L), alpha = 0.5),
      qcondig21(0.5, 0.5, theta = 0.5, alpha = numeric(0L))),
    0L
  )
})

test_that("rig() handles missingness appropriately", {
  sample_empty <- rig(numeric(0L), theta = 4, alpha = 5)
  expect_equal(nrow(sample_empty), 0L)
  sample_na <- rig(3, theta = c(1, NA_real_, NaN), alpha = 8)
  expect_true(all(is.na(sample_na)[-1, ]))
  expect_true(all(!is.na(sample_na)[1, ]))
})

rm("u_na_nan")
