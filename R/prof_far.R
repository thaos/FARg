#' @importFrom parallel mclapply
NULL

#' Compute Confidence Intervals for the FAR using the profile likelihood methodology.
#' 
#' Compute Confidence Intervals for the FAR based on the profile likelihood methodology and constrained optimization.
#' 
#' The profile likelihood for the  FAR is computed using constrained optimization which allow us to minimize the negative likelihood for a given value of the FAR. The method for the optimization was the Augmented Lagrangian and Adaptive Barrier Minimization Algorithm implemeted in the \code{auglag} function of the alabama package.
#' @param y_fit an object of class gauss_fit, gev_fit or gpd_fit.
#' @param xp the threshold used to compute the probability of exeeding that threshold.
#' @param t0 the time t0 to compute the probability of exeeding xp. If the time t0 is not present in dataset used for fitting the model, the closes time of the dataset is used.
#' @param t1 the time t1 to compute the probability of exeeding xp. If the time t1 is not present in dataset used for fitting the model, the closes time of the dataset is used.
#' @param ci_p the confidence level of the confidence intervals (between 0 and 1).
#' @param to_plot if TRUE, the function plots the offset profile likelihood. The profile likelihood plot is made in such a way that the FAR value is considered within the confidences intervals if the offset likelihood is above 0.
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
#'
#'p_gpd <- prof_far(gp_fit, xp, t0, t1, ci_p=0.95 ,to_plot=TRUE)
#'p_gev <- prof_far(ge_fit, xp, t0, t1, ci_p=0.95 ,to_plot=TRUE)
#'p_gauss <- prof_far(ga_fit, xp, t0, t1, ci_p=0.95 ,to_plot=TRUE)
#' @export
prof_far=function(y_fit, xp, t0, t1, ci_p=0.95 ,to_plot=FALSE, ...){
  # 1) fit_gauss/gev/gpd -> init + overallmax likelihood
  # 2) get_far -> p0/p1= 1-FAR dans la fonction contraindre
  # 3) creer une fonction contraindre par type d'ajustement (gev, gpd, gaussien)
  # 4) pareil pour la fonction oprmiser profil
  # 5) on vérifie que l'optimisation sous contrainte du bon far donne les même résultats que l'optimisation sans contraine
  # 6) construire le profil comme précédemment
  # 7) normalement on ne pourrait faire qu'une seule fonction pour construire le profile des 3 types d'ajusteùents
  pnt0 <- set_pnt(t0, xp, time_var=y_fit$time_var, y_fit$data)
  pnt1 <- set_pnt(t1, xp, time_var=y_fit$time_var, y_fit$data)
  init <- y_fit$par
  overall_max <- y_fit$objective
  far <- get_far(y_fit, pnt0, pnt1, ...)
  p0p1 <- 1 - far[1]
  print("--- FAR -----")
  print(far)
  fit_check <- optimize_constr(y_fit, optimize_far_prof, pnt0, pnt1, p0p1)
  if(!isTRUE(all.equal(fit_check$value, overall_max))){
    message("Slightly different results with constrained optimization")
    print(fit_check$value)
    print(overall_max)
    #             stop()
  }
  if(!isTRUE(all.equal(fit_check$par, init))){
    message("Slightly different results with constraint optimization")
    print(fit_check$par)
    print(init)
    #             stop()
  }
  aalpha <- qchisq(ci_p, 1)
  f_roots <- function(ratio){
    print(ratio)
    parmax  <-  -profil_optim(y_fit, optimize_far_prof, pnt0, pnt1, R=ratio)
    parmax + aalpha/2 + overall_max 
  }
  print("BORNE INF-------------------------------------------")
  ic_inf <- binf_time(c(-0.1, p0p1), fun=f_roots, nbdiv=2 ,xmax=p0p1,fmax=aalpha/2)
  print("BORNE SUP-------------------------------------------")
  if(p0p1 <= 1)
    ic_sup <- bsup_time(c(p0p1,1.1),fun=f_roots,nbdiv=2,xmax=p0p1,fmax=aalpha/2)
  else
    ic_sup <- bsup_time(c(p0p1 ,p0p1*2), fun=f_roots, nbdiv=2, xmax=p0p1, fmax=aalpha/2)
  ci <- select_prof_ic(ic_inf, ic_sup, ci_p, to_plot=to_plot)
  ci <- sort(1 - ci)
  out  <- c(ci[1], far[1], ci[2],far[-1])
  names(out)[c(1,3)] <- c("IC_inf", "IC_sup")
  # out  <- c(min(ci),get_far()[1],max(ci),get_far()[-1])
  out
}

#' Compute Confidence Intervals for the the probability of exceeding xp at time t using the profile likelihood methodology.
#' 
#' Compute Confidence Intervals for the probability of exceedance based on the profile likelihood methodology and constrained optimization.
#' 
#' The profile likelihood for the probability of exceedance  is computed using constrained optimization which allow us to minimize the negative likelihood for a given value of the probability. The method for the optimization was the Augmented Lagrangian and Adaptive Barrier Minimization Algorithm implemeted in the \code{auglag} function of the alabama package.
#' @param y_fit an object of class gauss_fit, gev_fit or gpd_fit.
#' @param xp the threshold used to compute the probability of exeeding that threshold.
#' @param t the time t to compute the probability of exeeding xp. If the time t is not present in dataset used for fitting the model, the closes time of the dataset is used.
#' @param ci_p the confidence level of the confidence intervals (between 0 and 1).
#' @param to_plot if TRUE, the function plots the offset profile likelihood. The profile likelihood plot is made in such a way that the FAR value is considered within the confidences intervals if the offset likelihood is above 0.
#' @param ... Arguments to be passed to methods.
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
#'
#'p_gpd <- prof_p(gp_fit, xp, t0, ci_p=0.95, to_plot=TRUE)
#'p_gev <- prof_p(ge_fit, xp, t0, ci_p=0.95, to_plot=TRUE)
#'p_gauss <- prof_p(ga_fit, xp, t0, ci_p=0.95, to_plot=TRUE)
#' @export
prof_p=function(y_fit, xp, t, ci_p=0.95 ,to_plot=FALSE, ...){
  pnt <- set_pnt(t, xp, time_var=y_fit$time_var, y_fit$data)
  init <- y_fit$par
  overall_max <- y_fit$objective
  prob <- get_p(y_fit, pnt, ...)
  p <- prob[1]
  print("--- Prob -----")
  print(prob)
  fit_check <- optimize_constr(y_fit, optimize_p_prof, pnt=pnt, P=p)
  print(fit_check$value)
  print(overall_max)
  if(!isTRUE(all.equal(fit_check$value, overall_max))){
    message("Slightly different results with constrained optimization")
    print(fit_check$value)
    print(overall_max)
    #             stop()
  }
  if(!isTRUE(all.equal(fit_check$par, init))){
    message("Slightly different results with constraint optimization")
    print(fit_check$par)
    print(init)
    #             stop()
  }
  aalpha <- qchisq(ci_p, 1)
  f_roots <- function(p){
    print(p)
    parmax  <-  -profil_optim(y_fit, optimize_p_prof, pnt, P=p)
    parmax + aalpha/2 + overall_max 
  }
  print("BORNE INF-------------------------------------------")
  ic_inf <- binf_time(c(-0.1, p), fun=f_roots, nbdiv=2 ,xmax=p,fmax=aalpha/2)
  print("BORNE SUP-------------------------------------------")
  ic_sup <- bsup_time(c(p,1.1),fun=f_roots,nbdiv=2,xmax=p,fmax=aalpha/2)
  ci <- select_prof_ic(ic_inf, ic_sup, ci_p, to_plot=to_plot)
  out  <- c(ci[1], prob[1], ci[2],prob[-1])
  names(out)[c(1,3)] <- c("IC_inf", "IC_sup")
  # out  <- c(min(ci),get_far()[1],max(ci),get_far()[-1])
  out
}


select_prof_ic <- function(ic_inf, ic_sup, ci_p, to_plot=FALSE){
  ratio_l <- c(ic_inf$x, ic_sup$x)
  parmax <- c(ic_inf$likel ,ic_sup$likel)
  aalpha <- qchisq(ci_p, 1)
  if(abs(max(parmax)-aalpha/2)>0.0001){
    warning("non equal mle")
    print(max(parmax))
    print(aalpha/2)
  }
  crit <- aalpha/2 - qchisq(max(1-(1-ci_p)/20,0.999), 1)/2
  cond <- parmax > crit
  ratio_l <- ratio_l[cond]
  parmax <- parmax[cond]
  if(to_plot){
    smth <- spline(ratio_l, parmax, n = 500,method="natural")
    plot(ratio_l, parmax, type = "l", xlab = "", ylab = "")
    abline(h = 0, lty = 2, col = 2)
    lines(smth$x,smth$y, lty = 2, col = 2)
  }
  ci <- ratio_l[parmax > 0]
  stopifnot(length(ci)>1)
  ci  <- c(min(ci), max(ci))
}
