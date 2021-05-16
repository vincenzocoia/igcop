.u <- c(10^(-(5:1)), 2:8/10, 1 - 10^(-(1:5)))
.v <- .u[c(4, 16, 17, 9, 10, 6, 1, 14, 11, 13, 5, 3, 7, 8, 15, 2, 12)]
.cpar <- list(
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
.alpha <- unique(vapply(.cpar, `[`, i = 2L, FUN.VALUE = numeric(1L)))
.theta <- unique(vapply(.cpar, `[`, i = 1L, FUN.VALUE = numeric(1L)))
.tau <- .u
