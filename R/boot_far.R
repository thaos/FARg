simple_boot <- function(y_fit, statistic, R=250, seed_init=NULL, ...){
  nr <- nrow(y_fit$data)
  if(!is.null(seed_init)) set.seed(seed_init)
  i_matrix <- matrix(sample.int(nr, size=nr*R, replace=TRUE), nrow=nr)
  stat_func <- function(indices,y_fit,...){
    statistic(y_fit, indices, ...)
  }
  apply(i_matrix, 2, stat_func, y_fit, ...)
}

#' boot_far generic.
#' 
#' Compute the confidence intervales for the FAR using resampling bootstrap.
#' 
#' Resamplin bootstrap is done by drawing with replacement points (which are at least pairs of date and time serie value) to reconstruct a bootstrap time series. Thus, in the bootstraped sample, there might be dates with no points and dates with several points. 
#' @param object an object of class gauss_fit, gev_fit, gpd_fit or of class trans. If object is of class trans, the argument y_fit must be passed.
#' @param xp the threshold used to compute the probability of exeeding that threshold.
#' @param t0 the time t0 to compute the probability of exeeding xp. If the time t0 is not present in dataset used for fitting the model, the closes time of the dataset is used.
#' @param t1 the time t1 to compute the probability of exeeding xp. If the time t1 is not present in dataset used for fitting the model, the closes time of the dataset is used.
#' @param ci_p the confidence level of the confidence intervals (between 0 and 1).
#' @param use_init whether to use the parameters fitted in object of class gauss_fit, gev_fit, or gpd_fit, to initialize the fir on the bootstrap sample.
#' @param R the number of bootstrap samples.
#' @param ... Arguments to be passed to methods,
#' @examples
#'data(tas)
#'
#'ge_fit <- gev_fit(eur_tas, data=tas, mu_mod=~gbl_tas, sig_mod=~gbl_tas, time_var="year")
#'gp_fit <- gpd_fit(eur_tas, data=tas, mu_mod=~gbl_tas, sig_mod=~gbl_tas, time_var="year", qthreshold=0.9)
#'ga_fit <- gauss_fit(eur_tas, data=tas, mu_mod=~gbl_tas, sig_mod=~gbl_tas, time_var="year")
#'
#'t1 <- 2003
#'t0 <- 1990
#'xp <- 1.6

#'b_gpd <- boot_far(gp_fit, xp, t0, t1, ci_p=0.95, under_threshold=TRUE)
#'b_gev <- boot_far(ge_fit, xp, t0, t1, ci_p=0.95)
#'b_gauss <- boot_far(ga_fit, xp, t0, t1, ci_p=0.95)
#'
#' @export
boot_far <- function(object, ...){
  UseMethod("boot_far")
}

#' @rdname boot_far 
#' @export
boot_far.default <- function(object, xp, t0, t1, ci_p=0.95, use_init=TRUE, R=250, ...){
  pnt0 <- set_pnt(t0, xp, time_var=object$time_var, object$data)
  pnt1 <- set_pnt(t1, xp, time_var=object$time_var, object$data)
  far_mle <- get_far(object, pnt0, pnt1, ...)
  #   log <- capture.output({
    #     boot_res=boot(data=object, statistic=boot_func, R=250,  pnt0=pnt0, pnt1=pnt1, use_init=use_init, parallel="snow")
    boot_res=simple_boot(object, statistic=boot_fun, R=R, get_fun=get_far, pnt0=pnt0, pnt1=pnt1, use_init=use_init, ...)
  #   })
  a=as.data.frame(t(boot_res))
  alpha <- 1-ci_p
  far_boot <- boot_res[1,]
  ic_boot <- quantile(far_boot,p=c(alpha/2, 0.5, 1-alpha/2))
  ans <- c(ic_boot, far_mle)
  names(ans)[1:3] <- c("IC_inf", "FAR_B", "IC_sup")
  ans
}

#' boot_p generic.
#' 
#' Compute the confidence intervales for the probability of exceeding xp at time t using resampling bootstrap.
#' 
#' Resamplin bootstrap is done by drawing with replacement points (which are at least pairs of date and time serie value) to reconstruct a bootstrap time series. Thus, in the bootstraped sample, there might be dates with no points and dates with several points. 
#' @param object an object of class gauss_fit, gev_fit, gpd_fit or of class trans. If object is of class trans, the argument y_fit must be passed.
#' @param xp the threshold used to compute the probability of exeeding that threshold.
#' @param t the time t to compute the probability of exeeding xp. If the time t0 is not present in dataset used for fitting the model, the closes time of the dataset is used.
#' @param ci_p the confidence level of the confidence intervals (between 0 and 1).
#' @param use_init whether to use the parameters fitted in object of class gauss_fit, gev_fit, or gpd_fit, to initialize the fir on the bootstrap sample.
#' @param R the number of bootstrap samples.
#' @param ... Arguments to be passed to methods,
#' @examples
#'data(tas)
#'
#'ge_fit <- gev_fit(eur_tas, data=tas, mu_mod=~gbl_tas, sig_mod=~gbl_tas, time_var="year")
#'gp_fit <- gpd_fit(eur_tas, data=tas, mu_mod=~gbl_tas, sig_mod=~gbl_tas, time_var="year", qthreshold=0.9)
#'ga_fit <- gauss_fit(eur_tas, data=tas, mu_mod=~gbl_tas, sig_mod=~gbl_tas, time_var="year")
#'
#'t1 <- 2003
#'t0 <- 1990
#'xp <- 1.6

#'b_gpd <- boot_p(gp_fit, xp, t0, ci_p=0.95, under_threshold=TRUE)
#'b_gev <- boot_p(ge_fit, xp, t0, ci_p=0.95)
#'b_gauss <- boot_p(ga_fit, xp, t0, ci_p=0.95)
#' @export
boot_p <- function(object, ...){
  UseMethod("boot_p")
}

#' @rdname boot_p 
#' @export
boot_p.default <- function(object, xp, t, ci_p=0.95, use_init=TRUE, R=250, ...){
  pnt <- set_pnt(t, xp, time_var=object$time_var, object$data)
  p_mle <- get_p(object, pnt, ...)
  #   log <- capture.output({
    #     boot_res=boot(data=object, statistic=boot_func, R=250,  pnt0=pnt0, pnt1=pnt1, use_init=use_init, parallel="snow")
    boot_res=simple_boot(object, statistic=boot_fun, R=R, get_fun=get_p, pnt=pnt, use_init=use_init, ...)
  #   })
  alpha <- 1-ci_p
  p_boot <- boot_res[1,]
  ic_boot <- quantile(p_boot,p=c(alpha/2, 0.5, 1-alpha/2))
  ans <- c(ic_boot, p_mle)
  names(ans)[1:3] <- c("IC_inf", "p_B", "IC_sup")
  ans
}

boot_fun <- function(object, indices, get_fun, use_init=TRUE, ...){
  UseMethod("boot_fun")
}

boot_fun.default <- function(object, indices, get_fun, use_init=TRUE, ...){
  b_fit <- boot_sample(object, indices, use_init)
  get_fun(b_fit, ...)
}

boot_sample <- function(object, indices, use_init=TRUE, ...){
  UseMethod("boot_sample")
}

boot_sample.gpd_fit <- function(object, indices, use_init=TRUE, ...){
  qthreshold <- object$rq_fitted$tau
  stopifnot(length(qthreshold) == 1)
  init <- NULL
  if(use_init) init <- object$par
  bdat <- object$data[indices,]
  b_fit <- gpd_fit(object$y[indices], bdat, object$mu_mod, object$sig_mod, object$sig_link, object$time_var, qthreshold, init)
}

boot_sample.gev_fit <- function(object, indices, pnt0, pnt1, use_init=TRUE, ...){
  init <- NULL
  if(use_init) init <- object$par
  bdat <- object$data[indices,]
  b_fit <- gev_fit(object$y[indices], bdat, object$mu_mod, object$sig_mod, object$sig_link, object$time_var, init)
}

boot_sample.gauss_fit <- function(object, indices, pnt0, pnt1, use_init=TRUE, ...){
  init <- NULL
  if(use_init) init <- object$par
  bdat <- object$data[indices,]
  b_fit <- gauss_fit(object$y[indices], bdat, object$mu_mod, object$sig_mod, object$sig_link, object$time_var, init)
}

