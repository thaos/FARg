#' @importFrom extRemes levd 
#' @importFrom extRemes fevd 
#' @importFrom extRemes pevd 
#' @importFrom extRemes plot.fevd 
#' @importFrom extRemes plot.fevd.mle 
#' @importFrom quantreg rq 
#' @importFrom quantreg predict.rq 

#' @export
gauss_negll <- function(y, mu0, mu1, sig0, sig1, mu_var, sig_var){
	stopifnot(length(mu_var) == length(sig_var))
	n <- length(mu_var)
	mu <- mu1*mu_var+mu0
	sig2 <- sig1*sig_var + sig0
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
gauss_fit <- function(ydat, init=NULL){
	gauss_lik <- function(y, init, mu_var, sig_var){
		gauss_negll(y, init[1], init[2], init[3], init[4], mu_var, sig_var)
	}
	if(is.null(init)){
		y_fit <- lm(y~mu_var, data=ydat)
		fit_res <- residuals(y_fit)
		var_fit <- lm(fit_res^2~sig_var, data=ydat)
		init <- c(coefficients(y_fit), coefficients(var_fit))
		print("--- Parameters Initialization -----")
		print(paste("neg-likelihood =", gauss_lik(y=ydat$y, init=init, mu_var=ydat$mu_var, sig_var=ydat$sig_var)))
		print(paste("mle =", do.call(paste, as.list(init))))
	}
	y_fit <- nlminb(start=init, gauss_lik, y=ydat$y, mu_var=ydat$mu_var, sig_var=ydat$sig_var)
	print("--- Parameters Optimization -----")
	print(paste("neg-likelihood =", y_fit$objective))
	print(paste("mle =", do.call(paste, as.list(y_fit$par))))
	y_fit$ydat <- ydat
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
plot.gev_fit <- function(y_fit){
  res <- fevd(y_fit$ydat$y, y_fit$ydat, location.fun=~mu_var, scale.fun=~sig_var, method="MLE")
  print(all.equal(y_fit$par,res$results$par))
  plot(res)
}

#' @export
plot.gpd_fit <- function(y_fit){
  res <- fevd(y_fit$ydat$y, y_fit$ydat, threshold=predict(y_fit$rq_fitted), scale.fun=~sig_var, type="GP", method="MLE")
  print(all.equal(y_fit$par,res$results$par))
  plot(res)
}


#' @export
plot.gauss_fit <- function(y_fit){
	param <- y_fit$par
	pred <- param[1] + param[2] * y_fit$ydat$mu_var
	res <- y_fit$ydat$y - pred
	var_res <- param[3] +param[4] * y_fit$ydat$sig_var
	res_std <- res / sqrt(var_res)
	par(mfrow=c(2,2))
	plot(pred, res_std, main="Standardized Residuals vs Fitted", xlab="Fitted", ylab="Standardized Residuals")
	abline(h=0, ,col="blue")
	extRemes::qqnorm(res_std)
	abline(a=0, b=1, col="blue")
	plot(dens <- density(res_std), main="Transformed Data")
	lines(dens$x, dnorm(x=dens$x), col="blue", lty=2)
	plot(y_fit$ydat$y, type="l")
	p=1-c(0.5, 0.05, 0.01)
	qXX=mapply(qnorm, mean=pred, sd=sqrt(var_res), MoreArgs=list(p=p))
	for(i in seq_along(p)){
		    lines(qXX[i,], col=i+1)
	}
	legend("topright", legend=paste("p:", p), lty=1, col= seq_along(p)+1)
}

#' @export
print.gauss_fit  <- print.gev_fit <- print.gpd_fit <- function(y_fit){
	print(y_fit[1:6])
}

