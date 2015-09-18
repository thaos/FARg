simple_boot <- function(y_fit, statistic, R=250, seed=1, ...){
  nr <- nrow(y_fit$data)
  #set.seed(1)
  i_matrix <- matrix(sample.int(nr, size=nr*R, replace=TRUE), nrow=nr)
  stat_func <- function(indices,y_fit,...){
    statistic(y_fit, indices, ...)
  }
  apply(i_matrix, 2, stat_func, y_fit, ...)
}

#' @export
boot_ic <- function(object, ...){
  UseMethod("boot_ic")
}

#' @export
boot_ic.default <- function(object, xp,t0, t1, ci_p=0.95, use_init=TRUE, ...){
  pnt0 <- set_pnt(t0, xp, time_var=object$time_var, object$data)
  pnt1 <- set_pnt(t1, xp, time_var=object$time_var, object$data)
  far_mle <- get_far(object, pnt0, pnt1, ...)
  #   log <- capture.output({
    #     boot_res=boot(data=object, statistic=boot_func, R=250,  pnt0=pnt0, pnt1=pnt1, use_init=use_init, parallel="snow")
    boot_res=simple_boot(object, statistic=boot_func, R=250, pnt0=pnt0, pnt1=pnt1, use_init=use_init)
  #   })
  alpha <- 1-ci_p
  far_boot <- boot_res[1,]
  ic_boot <- quantile(far_boot,p=c(alpha/2, 0.5, 1-alpha/2))
  ans <- c(ic_boot, far_mle)
  names(ans)[1:3] <- c("IC_inf", "FAR_B", "IC_sup")
  ans
}

boot_func <- function(object, ...){
  UseMethod("boot_func")
}

boot_func.gpd_fit <- function(object, indices, pnt0, pnt1, use_init=TRUE, ...){
  qthreshold <- object$rq_fitted$tau
  stopifnot(length(qthreshold) == 1)
  init <- NULL
  if(use_init) init <- object$par
  bdat <- object$data[indices,]
  b_fit <- gpd_fit(object$y[indices], bdat, object$mu_mod, object$sig_mod, object$time_var, qthreshold, init)
  get_far(b_fit, pnt0, pnt1, under_threshold=TRUE)
}

boot_func.gev_fit <- function(object, indices, pnt0, pnt1, use_init=TRUE, ...){
  init <- NULL
  if(use_init) init <- object$par
  bdat <- object$data[indices,]
  b_fit <- gev_fit(object$y[indices], bdat, object$mu_mod, object$sig_mod, object$time_var, init)
  get_far(b_fit, pnt0, pnt1)
}

boot_func.gauss_fit <- function(object, indices, pnt0, pnt1, use_init=TRUE, ...){
  init <- NULL
  if(use_init) init <- object$par
  bdat <- object$data[indices,]
  b_fit <- gauss_fit(object$y[indices], bdat, object$mu_mod, object$sig2_mod, object$time_var, init)
  get_far(b_fit, pnt0, pnt1)
}


