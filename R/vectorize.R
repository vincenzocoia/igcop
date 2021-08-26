y_interp_gen_inv <- function(p, eta, alpha) {
  formals_to("interp_gen_inv_vec")
}

y_interp_kappa <- function(x, eta, alpha) {
  formals_to("interp_kappa_vec")
}

y_interp_kappa_inv <- function(p, eta, alpha) {
  formals_to("interp_kappa_inv_vec")
}

y_interp_gen <- function(x, eta, alpha) {
  formals_to("interp_gen_vec")
}

y_igl_kappa <- function(x, alpha) {
  formals_to("igl_kappa_vec")
}

y_igl_kappa_D <- function(x, alpha) {
  formals_to("igl_kappa_D_vec")
}

y_igl_kappa_inv <- function(p, alpha) {
  formals_to("igl_kappa_inv_vec")
}

y_igl_gen <- function(x, alpha) {
  formals_to("igl_gen_vec")
}

y_igl_gen_D <- function(x, alpha) {
  formals_to("igl_gen_D_vec")
}

y_igl_gen_inv <- function(p, alpha) {
  formals_to("igl_gen_inv_vec")
}

formals_to <- function(.fn) {
  args_syms <- rlang::fn_fmls_syms(fn = caller_fn(n = 1))
  args_call <- rlang::call2("vec_recycle_common", !!!args_syms, .ns = "vctrs")
  args <- rlang::eval_tidy(args_call, env = caller_env(n = 1))
  rlang::exec(.fn, !!!args, .env = caller_env(n = 1))
}
