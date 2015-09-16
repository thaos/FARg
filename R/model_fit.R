#' @importFrom extRemes levd 
#' @importFrom extRemes fevd 
#' @importFrom extRemes pevd 
#' @importFrom extRemes plot.fevd 
#' @importFrom extRemes plot.fevd.mle 
#' @importFrom quantreg rq 
#' @importFrom quantreg predict.rq 

#' @export
exist_time_var <- function(time_var, data){
  if(is.character(time_var))
    return(exists(time_var) || is.element(time_var, names(data)))
  else
    return(exists(paste(deparse(substitute(time_var)), collapse=" ")) || !inherits(time_var, "simpleError"))
}

#' @export
get_param <- function(par, mod, data){
  mat <- model.matrix(mod, data=data)
  param <- mat %*% par
}

#' @export
complete_formula <- function(y, uncomplete_f){
  stopifnot(length(uncomplete_f) == 2)
  if(is.character(y)) 
    return((paste(y, "~", uncomplete_f[2])))
  paste(deparse(substitute(y)), "~", uncomplete_f[2])
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
gpd_negll <- function(y, threshold, sig, xi) {
  n <- length(y)
  stopifnot(length(threshold) == n & length(sig) == n)
	levd(x=y, threshold=threshold, scale=sig, shape=xi, type="GP",npy=1)
}

#' @export
gev_negll <- function(y, mu, sig, xi){
  n <- length(y)
  stopifnot(length(mu) == n & length(sig) == n)
	levd(x=y, location=mu, scale=sig, shape=xi, type="GEV")
}

format_init.gauss_fit <- function(init, mu_mod=~1, sig2_mod=~1){
  nb_mup <- length(attr(terms(mu_mod), "term.labels"))+attr(terms(mu_mod),"intercept")
  mu <- init[1:nb_mup]
  sig2 <- init[-(1:nb_mup)]
  list("mu"=mu, "sig2"=sig2)
}

#' @export
gauss_fit <- function(y, data, mu_mod, sig2_mod, time_var, init=NULL){
  stopifnot(exist_time_var(time_var, data))
  stopifnot(!is.null(data))
  y_name <- paste(deparse(substitute(y)), collapse="")
  if(is.element(y_name, names(data)))
    y <- data[, y_name]
  nb_mup <- length(attr(terms(mu_mod), "term.labels"))+attr(terms(mu_mod),"intercept")
	gauss_lik <- function(init){
    init_f <- format_init.gauss_fit(init, mu_mod, sig2_mod)
    mu <- get_param(init_f$mu, mu_mod, data)
  sig2 <- get_param(init_f$sig2, sig2_mod, data)
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
  y_fit$time_var <- time_var
	attr(y_fit, "class")= "gauss_fit"
	y_fit
}

format_init.gpd_fit <- function(init,  sig_mod){
  nb_sigp <- length(attr(terms(sig_mod), "term.labels"))+attr(terms(sig_mod),"intercept")
  sig <- init[1:nb_sigp]
  xi <- init[-(1:nb_sigp)]
  list("sig"=sig, "xi"=xi)
}

#'@export
gpd_fit <- function(y, data, mu_mod=~1, sig_mod=~1, time_var, qthreshold, init=NULL){
  stopifnot(exist_time_var(time_var, data))
  stopifnot(!is.null(data))
  y_name <- paste(deparse(substitute(y)), collapse="")
  if(is.element(y_name, names(data)))
    y <- data[, y_name]
  completed_formula <- complete_formula(y_name, mu_mod)
  rq_fitted <- rq(as.formula(completed_formula),data=data, tau=qthreshold)
	threshold <- predict(rq_fitted)
	if(is.null(init)){
		print("--- Parameters Initialization -----")
		init <- fevd(as.formula(paste(y_name,"~ 1")), data, threshold, scale.fun=sig_mod, type="GP", method="MLE")$results
		print(paste("neg-likelihood =", init$value))
		print(paste("mle =", do.call(paste, as.list(init$par))))
		init <- init$par
	}
	gpd_lik <- function(init){
    init_f <- format_init.gpd_fit(init,  sig_mod)
    sig <- get_param(init_f$sig , sig_mod, data)
		gpd_negll(y, threshold, sig, init_f$xi) 
	}
	y_fit=nlminb(start=init, gpd_lik) 
	print("--- Parameters Optimization -----")
	print(paste("neg-likelihood =", y_fit$objective))
	print(paste("mle =", do.call(paste, as.list(y_fit$par))))
	y_fit$y <- y 
	y_fit$data <- data 
	y_fit$mu_mod <- mu_mod
	y_fit$sig_mod <- sig_mod
	y_fit$qthreshold <- qthreshold
	y_fit$rq_fitted <- rq_fitted 
	y_fit$rate <- mean(y>threshold)
  y_fit$time_var <- time_var
	attr(y_fit, "class")= "gpd_fit"
	y_fit
}
	
format_init.gev_fit <- function(init, mu_mod, sig_mod){
  nb_mup <- length(attr(terms(mu_mod), "term.labels"))+attr(terms(mu_mod),"intercept")
  nb_sigp <- length(attr(terms(sig_mod), "term.labels"))+attr(terms(sig_mod),"intercept")
  mu <- init[1:nb_mup]
  sig <- init[nb_mup + (1:nb_sigp)]
  xi <- init[-(1:(nb_sigp + nb_mup))]
  list("mu"=mu, "sig"=sig, "xi"=xi)
}

#' @export
gev_fit <- function(y, data, mu_mod=~1, sig_mod=~1, time_var, init=NULL){
  stopifnot(exist_time_var(time_var, data))
  stopifnot(!is.null(data))
  y_name <- paste(deparse(substitute(y)), collapse=" ")
  if(is.element(y_name, names(data)))
    y <- data[, y_name]
	if(is.null(init)){
		print("--- Parameters Initialization -----")
		init  <- fevd(as.formula(paste(y_name, "~ 1")), data, location.fun=mu_mod, scale.fun=sig_mod, method="MLE")$results
		print(paste("neg-likelihood =", init$value))
		print(paste("mle =", do.call(paste, as.list(init$par))))
		init <- init$par
	}
	gev_lik <- function(init){
    init_f <- format_init.gev_fit(init, mu_mod, sig_mod)
    mu <- get_param(init_f$mu, mu_mod, data)
    sig <- get_param(init_f$sig, sig_mod, data)
		gev_negll(y, mu, sig, init_f$xi) 
	}
	y_fit=nlminb(start=init, gev_lik)
	print("--- Parameters Optimization -----")
	print(paste("neg-likelihood =", y_fit$objective))
	print(paste("mle =", do.call(paste, as.list(y_fit$par))))
	attr(y_fit, "class")= "gev_fit"
	y_fit$y <- y 
	y_fit$data <- data 
	y_fit$mu_mod <- mu_mod
	y_fit$sig_mod <- sig_mod
  y_fit$time_var <- time_var
	y_fit
}
	
#' @export
gev_fevd <- function(y, data, mu_mod=~1, sig_mod=~1, init=NULL){
	if(is.null(init)){
		y_fit=fevd(y, data, location.fun=mu_mod, scale.fun=sig_mod, method="MLE")
	} else{
		y_fit=fevd(y, data, location.fun=mu_mod, scale.fun=sig_mod, method="MLE", initial=init)
	}
	y_fit
}

#' @export
plot.gev_fit <- function(x, ...){
  res <- fevd(x$y, x$data, location.fun=x$mu_mod, scale.fun=x$sig_mod, method="MLE", initial=res$results$par)
  print(all.equal(x$par,res$results$par))
  plot(res)
}

#' @export
plot.gpd_fit <- function(x, ...){
  res <- fevd(x$y, x$data, threshold=predict(x$rq_fitted), scale.fun=x$sig_mod, type="GP", method="MLE", initial=res$results$par)
  print(all.equal(x$par,res$results$par))
  plot(res)
}


#' @export
compute_par <- function(object, newdata, ...){
  UseMethod("compute_par")
}

#' @export
compute_par.gpd_fit <- function(object, newdata, ...){
  threshold <- predict(object$rq_fitted, newdata)
  init_f <- format_init.gpd_fit(object$par,  object$sig_mod)
  data.frame("threshold"=threshold, "sig"=get_param(init_f$sig , object$sig_mod, newdata), "xi"=init_f$xi)
}

#' @export
compute_par.gev_fit <- function(object, newdata, ...){
  init_f <- format_init.gev_fit(object$par, object$mu_mod, object$sig_mod)
  data.frame("mu"=get_param(init_f$mu, object$mu_mod, newdata), "sig"=get_param(init_f$sig, object$sig_mod, newdata), "xi"=init_f$xi)
}


#' @export
compute_par.gauss_fit <- function(object, newdata, ...){
  init_f <- format_init.gauss_fit(object$par,  object$mu_mod, object$sig2_mod)
  data.frame("mu"=get_param(init_f$mu, object$mu_mod, newdata), "sig2"=get_param(init_f$sig2, object$sig2_mod, newdata))
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

