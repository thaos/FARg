#' refit generic.
#' 
#' .
#' 
#' \code{refit} recompute the statistical fit for the intercepts of the model fit
#'_m on the data used for the model fit_o.
#' fit_o and fit_m should follows that same model ans should be of the same R
#' classe
#' @param fit_o an object of class gauss_fit, gev_fit, gpd_fit. the data from
#' the fit_o model are used to refit the intercept of the fit_m model.
#' @param fit_o an object of class gauss_fit, gev_fit, gpd_fit. the data from
#' the fit_o model are used to refit the intercept of the fit_m model.
#' @return the fit_o object with the newly fitted parameters for the intercepts
#' and the other parameters being those of fit_m
#' @examples
#' #data(tas)
#' #Example with the same covariate for the mean and variance parameter
#' ga_fit <- gauss_fit(eur_tas, data=tas, mu_mod=~gbl_tas, sig_mod=~gbl_tas, time_var="year")
#' ga_refit <- refit(ga_fit, ga_fit)
#' ge_fit <- gev_fit(eur_tas, data=tas, mu_mod=~gbl_tas, sig_mod=~gbl_tas, time_var="year")
#' ge_refit <- refit(ge_fit, ge_fit)
#' gp_fit <- gpd_fit(eur_tas, data=tas, mu=~gbl_tas, sig_mod=~gbl_tas, time_var="year", qthreshold=0.9)
#' gp_refit <- refit(gp_fit, gp_fit)
#' @export
refit <- function(fit_o, fit_m){
	UseMethod("refit")
}

#' @describeIn refit for a gaussian fit.
#'@export
refit.gauss_fit <- function(fit_o, fit_m){
   y <- fit_o$y
   data <- fit_o$data
   mu_mod <-fit_o$mu_mod
   sig_mod <-fit_o$sig_mod
   mu_mat <-fit_o$mu_mat
   sig_mat <-fit_o$sig_mat
   mu_terms <-fit_o$mu_terms
   sig_terms <-fit_o$sig_terms
   time_var <- fit_o$time_var
   sig_link <- fit_o$sig_link
   link <- make.link(sig_link)
   init_m <- format_init.gauss_fit(fit_m$par, mu_terms, sig_terms)
   init_o <- format_init.gauss_fit(fit_o$par, mu_terms, sig_terms)
   gauss_lik <- function(init){
     init_m$mu[1] <- init[1]
     init_m$sig[1] <- init[2]
     mu <- get_param(init_m$mu, mu_mat)
     sig <- link$linkinv(get_param(init_m$sig, sig_mat))
     ans <- gauss_negll(y, mu, sig) 
     print(ans)
     ans
   }
   init <- c(init_o$mu[1] ,init_o$sig[1])
   y_fit <- nlminb(start=init, gauss_lik)
   print("--- Parameters Optimization -----")
   print(paste("neg-likelihood =", y_fit$objective))
   print(paste("mle =", do.call(paste, as.list(y_fit$par))))
   init_o$mu[1] <- y_fit$par[[1]]
   init_o$sig[1] <- y_fit$par[[2]]
   fit_o$par <- c(init_o$mu, init_o$sig)
   fit_o
}
#data(tas)
#Example with the same covariate for the mean and variance parameter
#ga_fit <- gauss_fit(eur_tas, data=tas, mu_mod=~gbl_tas, sig_mod=~gbl_tas, time_var="year")
#ga_refit <- refit.gauss_fit(ga_fit, ga_fit)


#' @describeIn refit for a gev fit.
#'@export
refit.gev_fit <- function(fit_o, fit_m){
   y <- fit_o$y
   data <- fit_o$data
   mu_mod <-fit_o$mu_mod
   sig_mod <-fit_o$sig_mod
   mu_mat <-fit_o$mu_mat
   sig_mat <-fit_o$sig_mat
   mu_terms <-fit_o$mu_terms
   sig_terms <-fit_o$sig_terms
   time_var <- fit_o$time_var
   sig_link <- fit_o$sig_link
   link <- make.link(sig_link)
   init_m <- format_init.gev_fit(fit_m$par, mu_terms, sig_terms)
   init_o <- format_init.gev_fit(fit_o$par, mu_terms, sig_terms)
	gev_lik <- function(init){
     init_m$mu[1] <- init[1]
     init_m$sig[1] <- init[2]
     mu <- get_param(init_m$mu, mu_mat)
     sig <- link$linkinv(get_param(init_m$sig, sig_mat))
     ans <- gev_negll(y, mu, sig, init[3]) 
     print(ans)
     ans
	}
   init <- c(init_o$mu[1] ,init_o$sig[1], init_o$xi)
	y_fit=nlminb(start=init, gev_lik)
	print("--- Parameters Optimization -----")
	print(paste("neg-likelihood =", y_fit$objective))
	print(paste("mle =", do.call(paste, as.list(y_fit$par))))
  init_o$mu[1] <- y_fit$par[1]
  init_o$sig[1] <- y_fit$par[2]
  fit_o$par <- c(init_o$mu, init_o$sig, y_fit$par[3])
  fit_o
}
#ge_fit <- gev_fit(eur_tas, data=tas, mu_mod=~gbl_tas, sig_mod=~gbl_tas, time_var="year")
#ge_refit <- refit.gev_fit(ge_fit, ge_fit)


#' @describeIn refit for a gpd fit.
#'@export
refit.gpd_fit <- function(fit_o, fit_m){
  y <- fit_o$y
  data <- fit_o$data
  mu_mod <-fit_o$mu_mod
  sig_mod <-fit_o$sig_mod
  mu_mat <-fit_o$mu_mat
  sig_mat <-fit_o$sig_mat
  mu_terms <-fit_o$mu_terms
  sig_terms <-fit_o$sig_terms
  time_var <- fit_o$time_var
  sig_link <- fit_o$sig_link
  qthreshold <- fit_o$qthreshold
  rq_fitted <- fit_m$rq_fitted 
  link <- make.link(sig_link)
  init_m <- format_init.gpd_fit(fit_m$par, sig_terms)
  init_o <- format_init.gpd_fit(fit_o$par, sig_terms)
  threshold <- predict(rq_fitted, newdata=data)
  reshift_threshold <- function(y, threshold, qthreshold){
    foptim <- function(par){
      ans <- (mean(y <= (threshold + par)) - qthreshold)^2
      # plot(y)
      # lines(threshold + par, col="red")
      # print(mean(y <= (threshold+par)))
      ans
    }
    res <- optimize(foptim, interval=range(scale(y, scale=FALSE)), tol=0.0001)
    res <- res$minimum
  }
  mu0 <- reshift_threshold(y=y, threshold=threshold, qthreshold=qthreshold)
  rq_fitted$coefficients[1] <- rq_fitted$coefficients[1] + mu0
  fit_o$rq_fitted <- rq_fitted
  threshold <- predict(rq_fitted, newdata=data)
  print(mean(y <= (threshold)))
  plot(y)
  lines(threshold , col="red")
  gpd_lik <- function(init){
    init_m$sig[1] <- init[1]
    sig <- link$linkinv(get_param(init_m$sig, sig_mat))
    ans <- gpd_negll(y, threshold, sig, init[2]) 
    print(ans)
    ans
  }
  init <- c(init_o$sig[1], init_o$xi)
  y_fit=nlminb(start=init, gpd_lik) 
  print("--- Parameters Optimization -----")
  print(paste("neg-likelihood =", y_fit$objective))
  print(paste("mle =", do.call(paste, as.list(y_fit$par))))
  init_o$sig[1] <- y_fit$par[1]
  fit_o$par <- c(init_o$sig, y_fit$par[2])
  fit_o
}
#gp_fit <- gpd_fit(eur_tas, data=tas, mu=~gbl_tas, sig_mod=~gbl_tas, time_var="year", qthreshold=0.9)
#gp_refit <- refit.gpd_fit(gp_fit, gp_fit)
#gp_refit <- fit.gpd_fit(gp_fit, gp_fit)
