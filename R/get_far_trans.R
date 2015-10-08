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
#' @param y_fit an object of class gauss_fit, gev_fit or gpd_fit. Only needed if the argument object if of class trans 
#' @export
get_far.trans <- function(object, y_fit, pnt0, pnt1, ...){
  pnt1_trans <- transform_pnt(object, pnt1)
  pnt0_trans <- transform_pnt(object, pnt0)
  far <- get_far(y_fit, pnt0_trans, pnt1_trans, ...)
  far
}

#' @rdname get_p 
#' @param y_fit an object of class gauss_fit, gev_fit or gpd_fit. Only needed if the argument object if of class trans 
#' @export
get_p.trans <- function(object, y_fit, pnt, ...){
  pnt_trans <- transform_pnt(object, pnt)
  p <- get_p(y_fit, pnt_trans, ...)
  p
}

boot_sample.bc_fit <- function(object, y_fit, indices, use_init=TRUE){
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
}

boot_sample.std <- function(object, y_fit, indices, use_init=TRUE){
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
}

boot_fun.trans <- function(object, y_fit, indices, get_fun, use_init=TRUE, ...){
  b_fit <- boot_sample(object, y_fit, indices, use_init)
  get_fun(object, b_fit, ...)
}

#' @rdname boot_far 
#' @export
boot_far.trans <- function(object, y_fit, xp, t0, t1, ci_p=0.95, use_init=TRUE, R=250, ...){
  pnt0 <- set_pnt(t0, xp, y_fit$time_var, y_fit$data)
  pnt1 <- set_pnt(t1, xp, y_fit$time_var, y_fit$data)
  far_mle <- get_far.trans(object, y_fit, pnt0, pnt1, ...)
  boot_res <- simple_boot(y_fit, boot_fun, object=object, R=R, get_fun=get_far.trans, pnt0=pnt0, pnt1=pnt1, use_init=use_init, ...)
  alpha <- 1-ci_p
  far_boot <- boot_res["FAR", ]
  ic_boot <- quantile(far_boot, probs=c(alpha/2, .5, 1-alpha/2))
  ans <- c(ic_boot, far_mle)
  names(ans)[1:3] <- c("IC_inf", "FAR_B", "IC_sup")
  ans
  ans
}

#' @rdname boot_p 
#' @export
boot_p.trans <- function(object, y_fit, xp, t, ci_p=0.95, use_init=TRUE, R=250, ...){
  pnt <- set_pnt(t, xp, y_fit$time_var, y_fit$data)
  p_mle <- get_p.trans(object, y_fit, pnt, ...)
  boot_res <- simple_boot(y_fit, boot_fun, object=object, R=R, get_fun=get_p.trans, pnt=pnt, use_init=use_init, ...)
  alpha <- 1-ci_p
  p_boot <- boot_res["p", ]
  ic_boot <- quantile(p_boot, probs=c(alpha/2, .5, 1-alpha/2))
  ans <- c(ic_boot, p_mle)
  names(ans)[1:3] <- c("IC_inf", "p_B", "IC_sup")
  ans
  ans
}
