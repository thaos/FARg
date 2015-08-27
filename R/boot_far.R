simple_boot <- function(y_fit, statistic, R=250, seed=1, ...){
  nr <- nrow(y_fit$ydat)
  #set.seed(1)
  i_matrix <- matrix(sample.int(nr, size=nr*R, replace=TRUE), nrow=nr)
  stat_func <- function(indices,y_fit,...){
    statistic(y_fit, indices, ...)
  }
  apply(i_matrix, 2, stat_func, y_fit, ...)
}

#' @export
boot_ic <- function(xp,t0, t1, y_fit, ci_p=0.95, use_init=TRUE, ...){
  pnt0 <- set_pnt(t0, xp, y_fit$ydat)
  pnt1 <- set_pnt(t1, xp, y_fit$ydat)
  far_mle <- get_far(y_fit, pnt0, pnt1, ...)
  log <- capture.output({
    #     boot_res=boot(data=y_fit, statistic=boot_func, R=250,  pnt0=pnt0, pnt1=pnt1, use_init=use_init, parallel="snow")
    boot_res=simple_boot(y_fit, statistic=boot_func, R=250,  y_fit=y_fit, pnt0=pnt0, pnt1=pnt1, use_init=use_init)
  })
  alpha <- 1-ci_p
  far_boot <- boot_res[1,]
  ic_boot <- quantile(far_boot,p=c(alpha/2, 0.5, 1-alpha/2))
  ans <- c(ic_boot, far_mle)
  names(ans)[1:3] <- c("IC_inf", "FAR_B", "IC_sup")
  ans
}

#' @export
boot_func <- function(y_fit, ...){
  UseMethod("boot_func")
}

#' @export
boot_func.gpd_fit <- function(y_fit, indices, pnt0, pnt1, use_init=TRUE, ...){
  qthreshold <- y_fit$rq_fitted$tau
  stopifnot(length(qthreshold) == 1)
  init <- NULL
  if(use_init) init <- y_fit$par
  bdat <- y_fit$ydat[indices,]
  print(head(bdat))
  b_fit <- gpd_fit(bdat, qthreshold, init)
  get_far(b_fit, pnt0, pnt1, under_threshold=TRUE)
}

#' @export
boot_func.gev_fit <- function(y_fit, indices, pnt0, pnt1, use_init=TRUE, ...){
  init <- NULL
  if(use_init) init <- y_fit$par
  bdat <- y_fit$ydat[indices,]
  b_fit <- gev_fit(bdat, init)
  get_far(b_fit, pnt0, pnt1)
}

#' @export
boot_func.gauss_fit <- function(y_fit, indices, pnt0, pnt1, use_init=TRUE, ...){
  init <- NULL
  if(use_init) init <- y_fit$par
  bdat <- y_fit$ydat[indices,]
  b_fit <- gauss_fit(bdat, init)
  get_far(b_fit, pnt0, pnt1)
}


