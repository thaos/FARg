#' @importFrom alabama auglag

constrain_far <- function(y_fit, pnt0, pnt1, R=0.5, ...){
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

fit2lik <- function(y_fit){
	UseMethod("fit2lik")
}

fit2lik.gpd_fit <- function(y_fit, pnt0, pnt1, R){
	gpd_lik <- function(init){
    init_f <- format_init.gpd_fit(init, y_fit$sig_terms)
    sig <- get_param(init_f$sig , y_fit$sig_mat)
		gpd_negll(y_fit$y, predict(y_fit$rq_fitted), sig, init_f$xi) 
	}
}

fit2lik.gev_fit <- function(y_fit, pnt0, pnt1, R){
	gev_lik <- function(init){
    init_f <- format_init.gev_fit(init, y_fit$mu_terms, y_fit$sig_terms)
    mu <- get_param(init_f$mu, y_fit$mu_mat)
    sig <- get_param(init_f$sig, y_fit$sig_mat)
		gev_negll(y_fit$y, mu, sig, init_f$xi) 
	}
}

fit2lik.gauss_fit <- function(y_fit){
	gauss_lik <- function(init){
      init_f <- format_init.gauss_fit(init, y_fit$mu_terms, y_fit$sig2_terms)
      mu <- get_param(init_f$mu, y_fit$mu_mat)
      sig2 <- get_param(init_f$sig2, y_fit$sig2_mat)
      gauss_negll(y_fit$y, mu, sig2) 
	}
}


optimize_far_prof <- function(y_fit, pnt0, pnt1, R){
	ans <- auglag(par=y_fit$par, fn=fit2lik(y_fit), heq=constrain_far(y_fit, pnt0, pnt1, R=R), control.outer=list(method="nlminb",trace=FALSE))
}

extract_l<- function(x){
	if (length(x) == 1)
		res = x
	else
		res = x$value
	res
}

profil_optim <- function(y_fit, optim_profil, ...){
	fit <- optimize_constr(y_fit, optim_profil, ...)
	extract_l(fit)
}

optimize_constr <- function(y_fit, optim_profil, ...){
	fit <- tryCatch(expr=optim_profil(y_fit, ...), 
                  error=function(e){
                    print("Cant optimize properly for this value of the ratio p0/p1")
                    return(10^6)
                  })
}

