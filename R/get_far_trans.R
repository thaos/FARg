#' @export
get_far.bc_fit <- function(y_trans, y_fit, pnt0, pnt1, qthreshold=NULL){
  pnt1_bc_fit <- transform_pnt(y_trans, pnt1)
  pnt0_bc_fit <- transform_pnt(y_trans, pnt0)
  far <- get_far(y_fit, pnt0_bc_fit, pnt1_bc_fit)
  far
}

#' @export
get_far.std <- function(y_trans, y_fit, pnt0, pnt1, qthreshold=NULL){
  pnt1_std <- transform_pnt(y_trans, pnt1)
  pnt0_std <- transform_pnt(y_trans, pnt0)
  far <- get_far(y_fit, pnt0_std, pnt1_std)
  far
}

boot_func.bc_fit <- function(y_trans, y_fit, indices, pnt0, pnt1){
  print(indices)
  y <- y_trans$y[indices]
  data <- y_trans$data[indices,]
  y_bc <- bc_fit(y_trans$l_lambda, y, data, mu_mod=y_trans$mu_mod, sig2_mod=y_trans$sig2_mod, time_var=y_trans$time_var)
  fit_method <- eval(parse(text=class(y_fit))) 
  fit_inputs <- as.list(set_fparam(fit_method, y_fit)) 
  fit_inputs$y <- y_bc$y_std 
  fit_inputs$data <- data 
  y_bc_fit  <- do.call(fit_method, fit_inputs)
  get_far.bc_fit(y_bc, y_bc_fit, pnt0, pnt1)
}

boot_func.std <- function(y_trans, y_fit, indices, pnt0, pnt1){
  #   print(indices)
  y <- y_trans$y[indices]
  data <- y_trans$data[indices,]
  y_std <- standardize(y, data, mu_mod=y_trans$mu_mod, sig2_mod=y_trans$sig2_mod )
  fit_method <- eval(parse(text=class(y_fit))) 
  fit_inputs <- as.list(set_fparam(fit_method, y_fit)) 
  fit_inputs$y <- y_std$y_std 
  fit_inputs$data <- data 
  y_std_fit  <- do.call(fit_method, fit_inputs)
  get_far.std(y_std, y_std_fit, pnt0, pnt1)
}

#' @export
boot_ic.bc_fit <- function(xp, t0, t1, y_trans, y_fit, ci_p=0.95, ...){
  pnt0 <- set_pnt(t0, xp, y_fit$time_var, y_fit$data)
  pnt1 <- set_pnt(t1, xp, y_fit$time_var, y_fit$data)
  far_mle <- get_far.bc_fit(y_trans, y_fit, pnt0, pnt1)
  boot_res <- simple_boot(y_fit, boot_func.bc_fit, y_trans=y_trans, pnt0=pnt0, pnt1=pnt1)
  alpha <- 1-ci_p
  far_boot <- boot_res["FAR", ]
  ic_boot <- quantile(far_boot, probs=c(alpha/2, .5, 1-alpha/2))
  ans <- c(ic_boot, far_mle)
  names(ans)[1:3] <- c("IC_inf", "FAR_B", "IC_sup")
  ans
}

#' @export
boot_ic.std <- function(xp, t0, t1, y_trans, y_fit, ci_p=0.95, ...){
  pnt0 <- set_pnt(t0, xp, y_fit$time_var, y_fit$data)
  pnt1 <- set_pnt(t1, xp, y_fit$time_var, y_fit$data)
  far_mle <- get_far.std(y_trans, y_fit, pnt0, pnt1)
  boot_res <- simple_boot(y_fit, boot_func.std, y_trans=y_trans, pnt0=pnt0, pnt1=pnt1)
  alpha <- 1-ci_p
  far_boot <- boot_res["FAR", ]
  ic_boot <- quantile(far_boot, probs=c(alpha/2, .5, 1-alpha/2))
  ans <- c(ic_boot, far_mle)
  names(ans)[1:3] <- c("IC_inf", "FAR_B", "IC_sup")
  ans
}
