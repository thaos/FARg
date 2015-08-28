#' @importFrom extRemes levd 
#' @importFrom extRemes fevd 
#' @importFrom extRemes pevd 
#' @importFrom extRemes plot.fevd 
#' @importFrom extRemes plot.fevd.mle 
#' @importFrom quantreg rq 
#' @importFrom quantreg predict.rq 

get_param <- function(par, mod, data){
  mat <- model.matrix(mod, data=data)
  param <- mat %*% par
}

complete_formula <- function(y, uncomplete_f){
  stopifnot(length(uncomplete_f) == 2)
  if(is.character(y)) 
    return(as.formula(paste(y, "~", uncomplete_f[2])))
  as.formula(paste(deparse(substitute(y)), "~", uncomplete_f[2]))
}


#' @export
gauss_negll <- function(y, mu, sig2){
  n <- length(y)
  stopifnot(length(mu) == n & length(sig2) == n)
	negll <- -0.5*(-n*log(2*pi)-sum(log(sig2)) - sum((y-mu)^2/sig2))
	if(is.na(negll)) negll <- 10^6
	negll
}

#' @export
gpd_negll <- function(y, threshold, sig0, sig1, xi, mu_var, sig_var) {
	stopifnot(length(mu_var) == length(sig_var))
	sig_v  <- sig0 + sig1*sig_var
	levd(x=y, threshold=threshold, scale=sig_v, shape=xi, type="GP",npy=1)
}

#' @export
gev_negll <- function(y, mu0, mu1, sig0, sig1, xi, mu_var, sig_var){
	stopifnot(length(mu_var) == length(sig_var))
	sig_v <- sig0 + sig1*sig_var
	mu_v <- mu0 + mu1*mu_var
	levd(x=y, location=mu_v, scale=sig_v, shape=xi, type="GEV")
}

#' @export
gauss_fit <- function(y, data, mu_mod, sig2_mod, init=NULL){
  env <- environment()
  nb_mup <- length(attr(terms(mu_mod), "term.labels"))+attr(terms(mu_mod),"intercept")
	gauss_lik <- function(init){
    mu_init <- init[1:nb_mup]
    sig2_init <- init[-(1:nb_mup)]
    mu <- get_param(mu_init, mu_mod, data)
    sig2 <- get_param(sig2_init, sig2_mod, data)
		gauss_negll(y, mu, sig2) 
	}
	if(is.null(init)){
		y_fit <- lm.fit(x=model.matrix(mu_mod, data=data), y=y)
		fit_res2 <- residuals(y_fit)^2
		var_fit  <- lm.fit(x=model.matrix(sig2_mod, data=data), y=fit_res2)
		init <- c(coefficients(y_fit), coefficients(var_fit))
		print("--- Parameters Initialization -----")
    print(paste("neg-likelihood =", gauss_lik(init)))
		print(paste("mle =", do.call(paste, as.list(init))))
	}
	y_fit <- nlminb(start=init, gauss_lik)
	print("--- Parameters Optimization -----")
	print(paste("neg-likelihood =", y_fit$objective))
	print(paste("mle =", do.call(paste, as.list(y_fit$par))))
  y_fit$y <- y
	y_fit$data <- data
	y_fit$mu_mod <- mu_mod
	y_fit$sig2_mod <- sig2_mod
	attr(y_fit, "class")= "gauss_fit"
	y_fit
}

#' @export
gpd_fit <- function(ydat, qthreshold, init=NULL){
	if(is.null(ydat$sig_var)) ydat$sig_var <- ydat$mu_var
	rq_fitted <- rq(y~mu_var, data=ydat, tau=qthreshold)
	threshold <- predict(rq_fitted)
	if(is.null(init)){
		print("--- Parameters Initialization -----")
		init <- gpd_fevd(ydat, threshold)$results
		print(paste("neg-likelihood =", init$value))
		print(paste("mle =", do.call(paste, as.list(init$par))))
		init <- init$par
	}
	gpd_lik <- function(y, init, threshold, mu_var, sig_var){
		gpd_negll(y, threshold, init[1], init[2], init[3], mu_var, sig_var) 
	}
	y_fit=nlminb(start=init, gpd_lik, y=ydat$y, threshold=threshold, mu_var=ydat$mu_var, sig_var=ydat$sig_var)
	print("--- Parameters Optimization -----")
	print(paste("neg-likelihood =", y_fit$objective))
	print(paste("mle =", do.call(paste, as.list(y_fit$par))))
	y_fit$ydat <- ydat 
	y_fit$rq_fitted <- rq_fitted 
	y_fit$rate <- mean(ydat$y>threshold)
	attr(y_fit, "class")= "gpd_fit"
	y_fit
}
	
#' @export
gpd_fevd <- function(ydat, threshold, init=NULL){
	if(is.null(init)){
		y_fit <- fevd(ydat$y, ydat, threshold=threshold, scale.fun=~sig_var, type="GP", method="MLE")
	} else{
		y_fit <- fevd(ydat$y, ydat, threshold=threshold, scale.fun=~sig_var, type="GP", method="MLE", initial=init)
	}
	y_fit
}

#' @export
gev_fit <- function(ydat, init=NULL){
	if(is.null(ydat$sig_var)) ydat$sig_var <- ydat$mu_var
	if(is.null(init)){
		print("--- Parameters Initialization -----")
		init <- gev_fevd(ydat)$results
		print(paste("neg-likelihood =", init$value))
		print(paste("mle =", do.call(paste, as.list(init$par))))
		init <- init$par
	}
	gev_lik <- function(y, init, mu_var, sig_var){
		gev_negll(y, init[1], init[2], init[3], init[4], init[5], mu_var, sig_var) 
	}
	y_fit=nlminb(start=init, gev_lik, y=ydat$y, mu_var=ydat$mu_var, sig_var=ydat$sig_var)
	print("--- Parameters Optimization -----")
	print(paste("neg-likelihood =", y_fit$objective))
	print(paste("mle =", do.call(paste, as.list(y_fit$par))))
	attr(y_fit, "class")= "gev_fit"
	y_fit$ydat <- ydat 
	y_fit
}
	
#' @export
gev_fevd <- function(ydat, init=NULL){
	if(is.null(init)){
		y_fit=fevd(ydat$y, ydat, location.fun=~mu_var, scale.fun=~sig_var, method="MLE")
	} else{
		y_fit=fevd(ydat$y, ydat, location.fun=~mu_var, scale.fun=~sig_var, method="MLE", initial=init)
	}
	y_fit
}

#' @export
plot.gev_fit <- function(x, ...){
  res <- fevd(x$ydat$y, x$ydat, location.fun=~mu_var, scale.fun=~sig_var, method="MLE")
  print(all.equal(x$par,res$results$par))
  plot(res)
}

#' @export
plot.gpd_fit <- function(x, ...){
  res <- fevd(x$ydat$y, x$ydat, threshold=predict(x$rq_fitted), scale.fun=~sig_var, type="GP", method="MLE")
  print(all.equal(x$par,res$results$par))
  plot(res)
}

#' @export
compute_par.gauss_fit <- function(object, newdata, ...){
  mu_mod <- object$mu_mod
  sig2_mod <- object$sig2_mod
  nb_mup <- length(attr(terms(mu_mod), "term.labels"))+attr(terms(mu_mod),"intercept")
  nb_sig2p <- length(attr(terms(sig2_mod), "term.labels"))+attr(terms(sig2_mod),"intercept")
  mu_par <- object$par[1:nb_mup]
  sig2_par <- object$par[-(1:nb_mup)]
  data.frame("mu"=get_param(mu_par, mu_mod, newdata), "sig2"=get_param(sig2_par, sig2_mod, newdata))
}

#' @export
plot.gauss_fit <- function(x, ...){
	param <- compute_par.gauss_fit(x, x$data)
	mu <- param$mu
	sig2 <- param$sig2
	res <- x$y - mu
	res_std <- res / sqrt(sig2)
	par(mfrow=c(2,2))
	plot(mu, res_std, main="Standardized Residuals vs Fitted", xlab="Fitted", ylab="Standardized Residuals")
	abline(h=0, ,col="blue")
	extRemes::qqnorm(res_std)
	abline(a=0, b=1, col="blue")
	plot(dens <- density(res_std), main="Transformed Data")
	lines(dens$x, dnorm(x=dens$x), col="blue", lty=2)
	plot(x$y, type="l")
	p=1-c(0.5, 0.05, 0.01)
	qXX=mapply(qnorm, mean=mu, sd=sqrt(sig2), MoreArgs=list(p=p))
	for(i in seq_along(p)){
		    lines(qXX[i,], col=i+1)
	}
	legend("topright", legend=paste("p:", p), lty=1, col= seq_along(p)+1)
}

#' @export
print.gauss_fit  <- print.gev_fit <- print.gpd_fit <- function(x, ...){
	print(x[1:6])
}

