simple_boot <- function(y_fit, statistic, R=250, seed=1, ...){
  nr <- nrow(y_fit$data)
  #set.seed(1)
  i_matrix <- matrix(sample.int(nr, size=nr*R, replace=TRUE), nrow=nr)
  stat_func <- function(indices,y_fit,...){
    statistic(y_fit, indices, ...)
  }
  apply(i_matrix, 2, stat_func, y_fit, ...)
}

#' boot_ic generic.
#' 
#' Compute the confidence intervales for the FAR using resampling bootstrap.
#' 
#' Resamplin bootstrap is done by drawing with replacement points (which are at least pairs of date and time serie value) to reconstruct a bootstrap time series. Thus, in the bootstraped sample, there might be dates with no points and dates with several points. 
#' @param object an object of class gauss_fit, gev_fit, gpd_fit or of class trans. If object is of class trans, the argument y_fit must be passed.
#' @param y_fit an object of class gauss_fit, gev_fit or gpd_fit. Only needed if the argument object is of class trans 
#' @param xp the threshold used to compute the probability of exeeding that threshold.
#' @param t0 the time t0 to compute the probability of exeeding xp. If the time t0 is not present in dataset used for fitting the model, the closes time of the dataset is used.
#' @param t1 the time t1 to compute the probability of exeeding xp. If the time t1 is not present in dataset used for fitting the model, the closes time of the dataset is used.
#' @param ci_p the confidence level of the confidence intervals (between 0 and 1).
#' @param use_init whether to use the parameters fitted in object of class gauss_fit, gev_fit, or gpd_fit, to initialize the fir on the bootstrap sample.
#' @param ... Arguments to be passed to methods,
#' @examples
#'data(tas)
#'
#'ge_fit <- gev_fit(eur_tas, data=tas, mu_mod=~avg_gbl_tas, sig_mod=~avg_gbl_tas, time_var="year")
#'gp_fit <- gpd_fit(eur_tas, data=tas, mu_mod=~avg_gbl_tas, sig_mod=~avg_gbl_tas, time_var="year", qthreshold=0.9)
#'ga_fit <- gauss_fit(eur_tas, data=tas, mu_mod=~avg_gbl_tas, sig2_mod=~avg_gbl_tas, time_var="year")
#'
#'t1 <- 2003
#'t0 <- 1990
#'xp <- 1.6
#'pnt1 <- set_pnt(t1, xp, time_var="year", tas)
#'pnt0 <- set_pnt(t0, xp, time_var="year", tas)

#'b_gpd <- boot_ic(gp_fit, xp, t0, t1, ci_p=0.95)
#'b_gev <- boot_ic(ge_fit, xp, t0, t1, ci_p=0.95)
#'b_gauss <- boot_ic(ga_fit, xp, t0, t1, ci_p=0.95)
#'
#' @export
boot_ic <- function(object, ...){
  UseMethod("boot_ic")
}

#' @rdname boot_ic 
#' @export
boot_ic.default <- function(object, xp,t0, t1, ci_p=0.95, use_init=TRUE, ...){
  pnt0 <- set_pnt(t0, xp, time_var=object$time_var, object$data)
  pnt1 <- set_pnt(t1, xp, time_var=object$time_var, object$data)
  far_mle <- get_far(object, pnt0, pnt1, ...)
  #   log <- capture.output({
    #     boot_res=boot(data=object, statistic=boot_func, R=250,  pnt0=pnt0, pnt1=pnt1, use_init=use_init, parallel="snow")
    boot_res=simple_boot(object, statistic=boot_func, R=250, pnt0=pnt0, pnt1=pnt1, use_init=use_init, ...)
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
  get_far(b_fit, pnt0, pnt1, ...)
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


