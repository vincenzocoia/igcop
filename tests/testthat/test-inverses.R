
tau <- c(10^(-(5:1)), 2:8/10, 1 - 10^(-(1:5)))
set.seed(5)
u <- sample(tau)
v <- sample(tau)


test_that("igl_gen works", {
  for (i in seq_along(alpha)){
    alpha_ <- alpha[i]
    tau2 <- igl_gen(
      igl_gen_inv(tau, alpha = alpha_),
      alpha = alpha_
    )
    expect_equal(tau, tau2)
  }
})

test_that("interp_gen_inv works", {
  for (cpar_ in cpar){
    eta <- cpar_[1]
    alpha <- cpar_[2]
    tau2 <- interp_gen(
      interp_gen_inv(tau, eta = eta, alpha = alpha),
      eta = eta, alpha = alpha
    )
    expect_equal(tau, tau2)
  }
})

test_that("interp_kappa_inv works", {
  for (cpar_ in cpar){
    eta <- cpar_[1]
    alpha <- cpar_[2]
    tau2 <- interp_kappa(
      interp_kappa_inv(tau, eta = eta, alpha = alpha),
      eta = eta, alpha = alpha
    )
    expect_equal(tau, tau2)
  }
})

test_that("qcondigcop works", {
  for (cpar_ in cpar){
    tau2 <- pcondigcop(qcondigcop(tau, u, cpar = cpar_), u, cpar = cpar_)
    test_that("qcondigcop works", {
      expect_equal(tau, tau2)
    })
  }
})


test_that("qcondiglcop works", {
  for (alpha_ in alpha){
    tau2 <- pcondiglcop(qcondiglcop(tau, u, cpar = alpha_), u, cpar = alpha_)
    test_that("qcondiglcop works", {
      expect_equal(tau, tau2)
    })
  }
})
