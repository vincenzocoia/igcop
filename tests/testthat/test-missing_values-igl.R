context("Check the handling of missing and improper inputs: IGL copula.")

u_na_nan <- c(0.5, NA_real_, NaN)

test_that("pigl() handles missingness appropriately", {
  args <- expand.grid(u = u_na_nan,
                      v = u_na_nan,
                      alpha = u_na_nan)
  eval <- do.call(pigl, args = args)
  expect_true(all(is.na(eval[-1])))
  expect_false(is.na(eval[1]))
  expect_length(
    c(pigl(numeric(0L), 0.5, alpha = 0.5),
      pigl(0.5, numeric(0L), alpha = 0.5),
      pigl(0.5, 0.5, alpha = numeric(0L))),
    0L
  )
})

test_that("digl() handles missingness appropriately", {
  args <- expand.grid(u = u_na_nan,
                      v = u_na_nan,
                      alpha = u_na_nan)
  digl(args$u, args$v, args$alpha)
  eval <- do.call(digl, args = args)
  expect_true(all(is.na(eval[-1])))
  expect_false(is.na(eval[1]))
  expect_length(
    c(digl(numeric(0L), 0.5, alpha = 0.5),
      digl(0.5, numeric(0L), alpha = 0.5),
      digl(0.5, 0.5, alpha = numeric(0L))),
    0L
  )
})

test_that("logdigl() handles missingness appropriately", {
  args <- expand.grid(u = u_na_nan,
                      v = u_na_nan,
                      alpha = u_na_nan)
  eval <- do.call(logdigl, args = args)
  expect_true(all(is.na(eval[-1])))
  expect_false(is.na(eval[1]))
  expect_length(
    c(logdigl(numeric(0L), 0.5, alpha = 0.5),
      logdigl(0.5, numeric(0L), alpha = 0.5),
      logdigl(0.5, 0.5, alpha = numeric(0L))),
    0L
  )
})

test_that("pcondigl21() handles missingness appropriately", {
  args <- expand.grid(u = u_na_nan,
                      v = u_na_nan,
                      alpha = u_na_nan)
  eval <- do.call(pcondigl21, args = args)
  expect_true(all(is.na(eval[-1])))
  expect_false(is.na(eval[1]))
  expect_length(
    c(pcondigl21(numeric(0L), 0.5, alpha = 0.5),
      pcondigl21(0.5, numeric(0L), alpha = 0.5),
      pcondigl21(0.5, 0.5, alpha = numeric(0L))),
    0L
  )
})

test_that("pcondigl12() handles missingness appropriately", {
  args <- expand.grid(u = u_na_nan,
                      v = u_na_nan,
                      alpha = u_na_nan)
  eval <- do.call(pcondigl12, args = args)
  expect_true(all(is.na(eval[-1])))
  expect_false(is.na(eval[1]))
  expect_length(
    c(pcondigl12(numeric(0L), 0.5, alpha = 0.5),
      pcondigl12(0.5, numeric(0L), alpha = 0.5),
      pcondigl12(0.5, 0.5, alpha = numeric(0L))),
    0L
  )
})

test_that("qcondigl12() handles missingness appropriately", {
  args <- expand.grid(p = u_na_nan,
                      v = u_na_nan,
                      alpha = u_na_nan)
  eval <- do.call(qcondigl12, args = args)
  expect_true(all(is.na(eval[-1])))
  expect_false(is.na(eval[1]))
  expect_length(
    c(qcondigl12(numeric(0L), 0.5, alpha = 0.5),
      qcondigl12(0.5, numeric(0L), alpha = 0.5),
      qcondigl12(0.5, 0.5, alpha = numeric(0L))),
    0L
  )
})

test_that("qcondigl21() handles missingness appropriately", {
  args <- expand.grid(p = u_na_nan,
                      u = u_na_nan,
                      alpha = u_na_nan)
  eval <- do.call(qcondigl21, args = args)
  expect_true(all(is.na(eval[-1])))
  expect_false(is.na(eval[1]))
  expect_length(
    c(qcondigl21(numeric(0L), 0.5, alpha = 0.5),
      qcondigl21(0.5, numeric(0L), alpha = 0.5),
      qcondigl21(0.5, 0.5, alpha = numeric(0L))),
    0L
  )
})

test_that("rigl() handles missingness appropriately", {
  sample_empty <- rigl(numeric(0L), alpha = 5)
  expect_equal(nrow(sample_empty), 0L)
  sample_na <- rigl(3, alpha = c(1, NA_real_, NaN))
  expect_true(all(is.na(sample_na)[-1, ]))
  expect_true(all(!is.na(sample_na)[1, ]))
})

rm("u_na_nan")
