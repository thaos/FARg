#' @importFrom alabama auglag

constrain <- function(y_fit, pnt0, pnt1, R=0.5, ...){
	res=function(par){
		tc=try({
			y_fit$par <- par
			far <- get_far(y_fit, pnt0, pnt1, ...)
			ratio <- 1 - far[1]
		},silent=TRUE)
		if(class(tc)=="try-error" )return(10^6)
		ratio-R	
	}
	res
}

optimize_profile <- function(y_fit, pnt0, pnt1, R){
	UseMethod("optimize_profile")
}

optimize_constr <- function(y_fit, pnt0, pnt1, R){
	tc=try({
		ans <- optimize_profile(y_fit, pnt0, pnt1, R)
	})
	if(class(tc)=="try-error" ){
		print("Cant optimize properly for this value of the ratio p0/p1")
		return(10^6)
	}
	ans
}

optimize_profile.gpd_fit <- function(y_fit, pnt0, pnt1, R){
	gpd_lik <- function(init){
    init_f <- format_init.gpd(init, y_fit$sig_mod)
    sig <- get_param(init_f$sig , y_fit$sig_mod, y_fit$data)
		gpd_negll(y_fit$y, predict(y_fit$rq_fitted), sig, init_f$xi) 
	}
	ans <- auglag(par=y_fit$par, fn=gpd_lik, heq=constrain(y_fit, pnt0, pnt1, R=R), control.outer=list(method="nlminb",trace=FALSE))
	ans
}

optimize_profile.gev_fit <- function(y_fit, pnt0, pnt1, R){
	gev_lik <- function(init){
    init_f <- format_init.gev(init, y_fit$mu_mod, y_fit$sig_mod)
    mu <- get_param(init_f$mu, y_fit$mu_mod, y_fit$data)
    sig <- get_param(init_f$sig, y_fit$sig_mod, y_fit$data)
		gev_negll(y_fit$y, mu, sig, init_f$xi) 
	}
	ans <- auglag(par=y_fit$par, fn=gev_lik, heq=constrain(y_fit, pnt0, pnt1, R=R), control.outer=list(method="nlminb",trace=FALSE))
	ans
}

optimize_profile.gauss_fit <- function(y_fit, pnt0, pnt1, R){
	gauss_lik <- function(init){
    init_f <- format_init.gauss(init, y_fit$mu_mod, y_fit$sig2_mod)
    mu <- get_param(init_f$mu, y_fit$mu_mod, y_fit$data)
  sig2 <- get_param(init_f$sig2, y_fit$sig2_mod, y_fit$data)
		gauss_negll(y_fit$y, mu, sig2) 
	}
	ans <- auglag(par=y_fit$par, fn=gauss_lik, heq=constrain(y_fit, pnt0, pnt1, R=R), control.outer=list(method="nlminb",trace=FALSE))
	ans
}

extract_l<- function(x){
	if (length(x) == 1)
		res = x
	else
		res = x$value
	res
}

profil_optim <- function(y_fit, pnt0, pnt1, R){
	fit <- optimize_constr(y_fit, pnt0, pnt1, R)
	extract_l(fit)
}
