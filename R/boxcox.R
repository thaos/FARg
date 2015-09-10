#' @importFrom evmix fgpd 
#' @importFrom MASS  boxcox

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

bc_fit=function(l_lambda, y, data, mu_mod=~1, sig2_mod=~1, time_var, ci_p=.95, to_plot=FALSE){
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
  init_f  <-  format_init.gauss_fit(res$par, mu_mod, sig2_mod)
  par_gauss <- compute_par.gauss_fit(res, res$data)
  res$y_std <- (res$y_lambda - par_gauss$mu) / sqrt(par_gauss$sig2)
  if(to_plot){
    par(mfrow=c(1,2))
    plot(l_lambda[cond], l_prof[cond],type="l",ylab="")
    abline(v=lambda, lty=2)
    abline(v=res$ci[1], lty=2)
    abline(v=res$ci[2], lty=2)
    abline(h=crit_ci, lty=2)
    boxcox(as.formula(complete_formula(y,mu_mod)),lambda=l_lambda)
  }
  class(res) <- "bc_fit"
  res
}

transform_newdat.bc_fit <- function(y_trans, y, newdata, ...){
  y <- y / exp(mean(log(y_trans$y)))
  y_lambda <- bc(y, y_trans$lambda)
  par_gauss <- compute_par.gauss_fit(y_trans, newdata )
  y_std <- (y_lambda - par_gauss$mu) / sqrt(par_gauss$sig2)
  data.frame(y, y_lambda, y_std)
}

transform_pnt.bc_fit <- function(y_trans, pnt){
  pnt$y <- as.numeric(transform_newdat.bc_fit(y_trans, y=pnt$y, newdata=pnt)$y_std)
  pnt
}


zp2yp <- function(lambda, mu, sigma, zp) 
	(lambda * ( mu + sigma * zp) + 1)^(1 / lambda)

