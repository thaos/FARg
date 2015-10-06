#' @importFrom evmix fgpd 
#' @importFrom MASS  boxcox
NULL

#' Apply the boxcox transformation to an univariate time-serie.
#' 
#' Apply the boxcox transformation to an univariate time-serie.
#' 
#' The univarite time-serie must be positive
#' @param x a vector of values.
#' @param lambda the scalar parameter of the boxcox transformation.
#' @return returns a vector of transformed values.
#' @export
bc <- function(x,lambda){
  y=NA
  if(lambda!=0){y=(x^lambda-1)/lambda}
  if(lambda==0){y=log(x)}
  return(y)
} 

lambda_prof <- function(lambda, y, data, mu_mod=~1, sig2_mod=~1, time_var, init=NULL){
	y_lambda <- bc(y,lambda)
  y_fit <- gauss_fit(y_lambda, data, mu_mod, sig2_mod, time_var, init)
	y_fit$objective <- y_fit$objective+(1-lambda)*sum(log(y))
	y_fit
}

#' Lambda selection and application of the boxcox transformation. 
#' 
#' Select a lambda and apply the boxcox transformation to make rhe data more gaussian.
#' 
#' The selection is lambda is done using the profile likelihood method. The aims is to find the lambda which maximize the likelihood of independant gaussian variables. The mean and the variance of the gaussian variables can be modelled as being linearly dependent on a specified set of covariates. The data are first normalized by \code{exp(mean(log(y)))} to make the optimization step easier.
#' @param l_lambda a vector of lambda value  the profile likelihood will be computed at.
#' @param y a univarite time-serie.
#' @param data a data.frame with the values for the covariates used to model the mean and the variance of the gaussian variables.
#' @param mu_mod a formula giving the list of covariates the mean parameter of the gaussians depends on.
#' @param sig2_mod a formula giving the list of covariates the variabces parameter of the gaussians depends on
#' @param time_var the variable used as time for the time serie.
#' @param ci_p the level of the confidence intervals wanted.
#' @param to_plot whether to plot the profile likelihood.
#' @return returns a list with the arguments of the function as well as :
#' \itemize{
#' \item y_eml, the original time serie normalize by  \code{exp(mean(log(y)))}.
#' \item y_lambda, the original time serie using to the boxcox transformation .
#' \item y_std, the standardized time serie.
#' \item mu_fit, the linear model fitted for the trend.
#' \item sig2_fit, the linear model fitted for the variance.
#'}
#' @examples
#'data(tas)
#'eur_tas_positive <- with(tas, eur_tas + abs(min(eur_tas)) +1)
#'y_bc <- bc_fit(l_lambda=seq(-1, 1, 0.02), y=eur_tas_positive, data=tas, mu_mod=~avg_gbl_tas, sig2_mod=~avg_gbl_tas, time_var="year", to_plot=TRUE)
#'y_bc_fit <- gpd_fit(y_bc$y_std ,  data=tas, time_var="year", qthreshold=0.8)
#'t1 <- 2003
#'t0 <- 1990
#'xp <- 1.6
#'pnt1 <- set_pnt(t1, xp, time_var="year", tas)
#'pnt0 <- set_pnt(t0, xp, time_var="year", tas)
#'get_far(y_bc, y_bc_fit, pnt0, pnt1, under_threshold=TRUE)
#'\donttest{ boot_ic(y_bc,  y_bc_fit, xp, t0, t1, under_threshold=TRUE)}
#' @export
bc_fit=function(l_lambda, y, data, mu_mod=~1, sig2_mod=~1, time_var, ci_p=.95, to_plot=FALSE){
  stopifnot(!is.null(data))
  y_name <- paste(deparse(substitute(y)), collapse="")
  #   if(is.element(y_name, names(data)))
  #     y <- data[, y_name]
  #   y_name <- "y"
  #   if(is.element("y", names(data))){
  #     y_name <- random_name(data=data)
  #     assign(y_name, y)
  #   }
  y_name <- check_y_name(y, y_name, data)
	mu_terms <- terms(mu_mod)
	sig2_terms <- terms(sig2_mod)
  stopifnot(all(y > 0))
  #for the stability of the boxcox trannsformation
  # maybe should not erase the original y
  eml <- exp(mean(log(y)))
  y_eml <- y / eml
  lambda_lik <- function(lambda){
    lambda_prof(lambda, y, data, mu_mod, sig2_mod, time_var)$objective
  }
  aalpha <- qchisq(ci_p, 1)
  l_prof <- -unlist(lapply(l_lambda, lambda_lik))
  overallmax <- max(l_prof)
  crit <- overallmax - qchisq(0.999, 1)/2
  cond <- l_prof > crit
  crit_ci <- overallmax - qchisq(ci_p, 1)/2
  cond_ci <- l_prof > crit_ci
  lambda <- optimize(lambda_lik, range(l_lambda), maximum=FALSE)$minimum
  res <- lambda_prof(lambda, y_eml, data, mu_mod, sig2_mod, time_var)
  res$lambda <- lambda
  res$l_lambda <- l_lambda
  res$ci <- range(l_lambda[cond_ci])
  res$y_lambda <- res$y
  res$y_eml <- y_eml
  res$y <- y
  res$data <- data
  init_f  <-  format_init.gauss_fit(res$par, mu_terms, sig2_terms)
  par_gauss <- compute_par.gauss_fit(res, res$data)
  res$y_std <- (res$y_lambda - par_gauss[, 1]) / sqrt(par_gauss[, 2])
  if(to_plot){
    par(mfrow=c(1,2))
    plot(l_lambda[cond], l_prof[cond],type="l",ylab="")
    abline(v=lambda, lty=2)
    abline(v=res$ci[1], lty=2)
    abline(v=res$ci[2], lty=2)
    abline(h=crit_ci, lty=2)
    boxcox(as.formula(complete_formula(y_name,mu_mod)), data=data, lambda=l_lambda)
  }
  class(res) <- "trans"
  class(res) <- append(class(res), "bc_fit")
  res
}

transform_newdat.bc_fit <- function(y_trans, y, newdata, ...){
  y <- y / exp(mean(log(y_trans$y)))
  y_lambda <- bc(y, y_trans$lambda)
  par_gauss <- compute_par.gauss_fit(y_trans, newdata )
  y_std <- (y_lambda - par_gauss[, 1]) / sqrt(par_gauss[, 2])
  data.frame(y, y_lambda, y_std)
}

#' @rdname transform_pnt
#' @export
transform_pnt.bc_fit <- function(y_trans, pnt, ...){
  # if several y in pnt data.frame, it takes the first one
  pnt$y <- as.numeric(transform_newdat.bc_fit(y_trans, y=pnt$y, newdata=pnt[, -1])$y_std)
  pnt
}


zp2yp <- function(lambda, mu, sigma, zp) 
	(lambda * ( mu + sigma * zp) + 1)^(1 / lambda)

