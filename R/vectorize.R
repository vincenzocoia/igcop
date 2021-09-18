#' Send arguments to a function after vectorizing
#'
#' When used within a (encapsulating) function, `formals_to`
#' recycles the inputs of the encapsulating function so that
#' they are vectors of the same length, and then sends these
#' updated arguments to some specified function.
#'
#' @param .fn The function you want to send the recycled arguments to.
#' @return The function `.fn` evaluated with the arguments given in
#' the encapsulating function.
formals_to <- function(.fn) {
  args_syms <- rlang::fn_fmls_syms(fn = rlang::caller_fn(n = 1))
  args_call <- rlang::call2("vec_recycle_common", !!!args_syms, .ns = "vctrs")
  args <- rlang::eval_tidy(args_call, env = rlang::caller_env(n = 1))
  rlang::exec(.fn, !!!args, .env = rlang::caller_env(n = 1))
}

#' IG/IGL Generators and Related Functions
#'
#' These are the psi, H, and kappa functions
#' of the IG and IGL copula families.
#'
#' @param x Function argument. Vector of non-negative values.
#' @param p Function inverse argument. Vector of values between 0 and 1.
#' @param eta,alpha Function parameters. Vector of positive values.
#' @note Inputs must be recyclable via `vctrs::vec_recycle_common()`.
#' @details Kappa function and its relatives have prefix `igl_kappa`;
#' Psi function and its relatives have prefix `igl_gen`;
#' Interpolating function H with either kappa or psi has
#' `igl` prefix replaced with `interp`. Relatives of these functions:
#' suffix `inv` indicates inverse; suffix `D` represents function
#' derivative, and `D1` derivative with respect to the first argument.
#'. Suffix `_vec` indicates that the entries must be vectors of
#' the same length; `_single` means entries must be
#' scalars.
#' @return The function values, as a vector.
#' @rdname generators
interp_gen_inv <- function(p, eta, alpha) {
  formals_to("interp_gen_inv_vec")
}

#' @rdname generators
interp_kappa <- function(x, eta, alpha) {
  formals_to("interp_kappa_vec")
}

#' @rdname generators
interp_kappa_inv <- function(p, eta, alpha) {
  formals_to("interp_kappa_inv_vec")
}

#' @rdname generators
interp_gen <- function(x, eta, alpha) {
  formals_to("interp_gen_vec")
}

#' @rdname generators
igl_kappa <- function(x, alpha) {
  formals_to("igl_kappa_vec")
}

#' @rdname generators
igl_kappa_D <- function(x, alpha) {
  formals_to("igl_kappa_D_vec")
}

#' @rdname generators
igl_kappa_inv <- function(p, alpha) {
  formals_to("igl_kappa_inv_vec")
}

#' @rdname generators
igl_gen <- function(x, alpha) {
  formals_to("igl_gen_vec")
}

#' @rdname generators
igl_gen_D <- function(x, alpha) {
  formals_to("igl_gen_D_vec")
}

#' @rdname generators
igl_gen_inv <- function(p, alpha) {
  formals_to("igl_gen_inv_vec")
}


