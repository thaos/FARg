set_fparam <- function(f, gfit){
  l <- formals(f)
  for( p in names(l)){
    #     print('-------------------') 
    #     print(p)
    #     print(gfit[[p]] )
    l[[p]] <- gfit[[p]]
  }
  l
}

#' @rdname get_far 
#' @export
get_far.trans <- function(object, y_fit, pnt0, pnt1, qthreshold=NULL, ...){
  pnt1_trans <- transform_pnt(object, pnt1)
  pnt0_trans <- transform_pnt(object, pnt0)
  far <- get_far(y_fit, pnt0_trans, pnt1_trans, ...)
  far
}

boot_func.bc_fit <- function(object, y_fit, indices, pnt0, pnt1, use_init=TRUE, ...){
  print(indices)
  y <- object$y[indices]
  data <- object$data[indices,]
  y_bc <- bc_fit(object$l_lambda, y, data, mu_mod=object$mu_mod, sig2_mod=object$sig2_mod, time_var=y_fit$time_var)
  fit_method <- eval(parse(text=class(y_fit))) 
  fit_inputs <- as.list(set_fparam(fit_method, y_fit)) 
  fit_inputs$y <- y_bc$y_std 
  fit_inputs$data <- data 
  if(use_init) fit_inputs$init <- y_fit$par
  y_bc_fit  <- do.call(fit_method, fit_inputs)
  get_far.trans(y_bc, y_bc_fit, pnt0, pnt1, ...)
}

boot_func.std <- function(object, y_fit, indices, pnt0, pnt1, use_init=TRUE, ...){
  #   print(indices)
  y <- object$y[indices]
  data <- object$data[indices,]
  y_std <- standardize(y, data, mu_mod=object$mu_mod, sig2_mod=object$sig2_mod )
  fit_method <- eval(parse(text=class(y_fit))) 
  fit_inputs <- as.list(set_fparam(fit_method, y_fit)) 
  fit_inputs$y <- y_std$y_std 
  fit_inputs$data <- data 
  if(use_init) fit_inputs$init <- y_fit$par
  y_std_fit  <- do.call(fit_method, fit_inputs)
  get_far(y_std, y_std_fit, pnt0, pnt1, ...)
}

#' @rdname boot_ic 
#' @export
boot_ic.trans <- function(object, y_fit, xp, t0, t1, ci_p=0.95, use_init=TRUE, R=250, ...){
  pnt0 <- set_pnt(t0, xp, y_fit$time_var, y_fit$data)
  pnt1 <- set_pnt(t1, xp, y_fit$time_var, y_fit$data)
  far_mle <- get_far.trans(object, y_fit, pnt0, pnt1, ...)
  boot_res <- simple_boot(y_fit, boot_func, object=object, pnt0=pnt0, pnt1=pnt1, use_init=use_init, R=R, ...)
  alpha <- 1-ci_p
  far_boot <- boot_res["FAR", ]
  ic_boot <- quantile(far_boot, probs=c(alpha/2, .5, 1-alpha/2))
  ans <- c(ic_boot, far_mle)
  names(ans)[1:3] <- c("IC_inf", "FAR_B", "IC_sup")
  ans
}
