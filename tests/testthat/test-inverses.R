cpar <- list(
  c(2.5, 12.3), # nearly independent
  c(2.5, 4.3), # weak dependence spear=.040
  c(2.5, 2.3), # spear=.145
  c(2.5, 0.3), # spear=.244##
  c(2.5, 0.1), # unreliable (for small k)###
  c(5.5,3.1),  # spear=.231
  c(8.5,3.1),  # spear=.330
  c(15.5,3.1), # spear=.462
  c(25.5,3.1), # spear=.553
  c(45.5,3.1), # spear=.629
  c(75.5,3.1), # spear=.674
  c(75.5,7.1), # spear=.661
  c(155.5,3.1), # spear=.714
  c(300.5,3.1), # spear=.733
  c(600.5,3.1), # spear=.743
  c(600.5,7.1), # spear=.829
  c(600.5,14.1), # spear=.853
  c(600.5,34.1), # spear=.811
  c(1200.5,14.1), # spear=.883
  c(2400.5,14.1), # spear=.899
  c(2400.5,34.1), # spear=.918
  c(2400.5,54.1), # spear=.??? unstable for rhoS
  c(2400.5,44.1), # spear=.914
  c(4800.5,34.1), # spear=.937
  c(9600.5,34.1) # spear=.947  # try inference with this parameter vector
)
alpha <- unique(purrr::map_dbl(cpar, 2))
theta <- unique(purrr::map_dbl(cpar, 1))
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
