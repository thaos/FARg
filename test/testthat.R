rm(list=ls())

# library(parallel)
library(devtools)
devtools::load_all()
library(roxygen2)
# library(alabama)
# library(extRemes)
# library(quantreg)

data(tas)
ydat <- with(tas,as.data.frame(cbind(year, eur_tas, avg_gbl_tas, avg_gbl_tas)))
names(ydat) <- c("year", "y", "mu_var", "sig_var")

rq_fitted <- rq(y~mu_var, data=ydat, tau=0.90)
threshold <- predict(rq_fitted)
gpd_fevd(ydat, threshold=threshold)
gev_fevd(ydat)
gpd_fit(ydat, qthreshold=0.9)
plot(gpd_fit(ydat, qthreshold=0.9))
gev_fit(ydat)
plot(gev_fit(ydat))
gauss_fit(ydat)
plot(gauss_fit(ydat))


t1 <- 2003
t0 <- 1990
xp <- 1.6
pnt1 <- set_pnt(t1, xp, ydat)
pnt0 <- set_pnt(t0, xp, ydat)

get_p(gauss_fit(ydat), pnt1)
get_far(gauss_fit(ydat), pnt0, pnt1)

get_p(gpd_fit(ydat, qthreshold=0.9), pnt1)
get_far(gpd_fit(ydat, qthreshold=0.9), pnt0, pnt1)
get_far(gpd_fit(ydat, qthreshold=0.9), pnt0, pnt1, under_threshold=TRUE)

get_p(gev_fit(ydat), pnt1)
get_far(gev_fit(ydat), pnt0, pnt1)

R <- 1-get_far(gpd_fit(ydat, qthreshold=0.9), pnt0, pnt1)[1]
constrain(gpd_fit(ydat, qthreshold=0.9), pnt0, pnt1, R=R)
optimize_profile(gpd_fit(ydat, qthreshold=0.9), pnt0, pnt1, R=R)
optimize_constr(gpd_fit(ydat, qthreshold=0.9), pnt0, pnt1, R=R)

constrain(gev_fit(ydat), pnt0, pnt1, R=R)
optimize_profile(gev_fit(ydat), pnt0, pnt1, R=R)
optimize_constr(gev_fit(ydat), pnt0, pnt1, R=R)

constrain(gauss_fit(ydat), pnt0, pnt1, R=R)
optimize_profile(gauss_fit(ydat), pnt0, pnt1, R=R)
optimize_constr(gauss_fit(ydat), pnt0, pnt1, R=R)

p_gpd <- prof_ic(xp, t0, t1, gpd_fit(ydat, qthreshold=0.9), ci_p=0.95 ,to_plot=TRUE)
p_gev <- prof_ic(xp, t0, t1, gev_fit(ydat), ci_p=0.95 ,to_plot=TRUE)
p_gauss <- prof_ic(xp, t0, t1, gauss_fit(ydat), ci_p=0.95 ,to_plot=TRUE)

b_gpd <- boot_ic(xp, t0, t1, gpd_fit(ydat, qthreshold=0.9), ci_p=0.95) 
b_gev <- boot_ic(xp, t0, t1, gev_fit(ydat), ci_p=0.95) 
b_gauss <- boot_ic(xp, t0, t1, gauss_fit(ydat), ci_p=0.95)



t1 <- 1970
t0 <- 2003
xp <- 1.6


set_env_gauss <- function(){
		this_env <- environment()
		repet <- 4
		changement <- 100
		years <- 1850:2200
		n <- length(unique(years))
		t <- seq(1,n)
		mu <- (t>changement)*.020*(t-changement)+100
		sigma <- (t>changement)*.001*(t-changement)+1
		sigma <- sqrt(sigma)
		mu <- rep(mu,repet)
		sigma <- rep(sigma,repet)
		t <- rep(t,repet)
		year <- rep(years,repet)
		covariate <- (mapply(qnorm,mean=mu,sd=sigma,MoreArgs=list("p"=0.5)))
		ydat <- data.frame("year"=year, "y"=NA, "mu_var"=covariate, "sig_var"=covariate)
		xp <- 102
		t0 <- 1900
		t1 <- 1950
		this_env
}
env <- set_env_gauss()

gen_data_gauss <- function(env){
		list2env(as.list(env),envir=environment())
		this_env <- environment()
		ydat$y <- mapply(rnorm, mean=mu, sd=sigma, MoreArgs=list(n=1))
		this_env
}
env <- gen_data_gauss(env)

get_theo_gauss <- function(env){
		list2env(as.list(env),envir=environment())
		get_p_theo <- function(year, ydat){
				i <- min(which.min(abs(year-ydat$year)))
				sig <- sigma[i]
				mu <- mu[i]
				pnorm(xp, mean=mu, sd=sig, lower.tail=FALSE)
		}
		get_far_theo <- function(t0, t1){
				p0 <- get_p_theo(t0, ydat)
				p1 <- get_p_theo(t1, ydat)
				if (p0 == 0 & p1 == p0)
						far <- 1
				else
						far <- 1 - p0/p1
		}
		get_far_theo(t0, t1)
}
theo <- get_theo_gauss(env)

get_ic <- function(fit_func, ic_func){
		function(env, ci_p=0.95, ...){
				list2env(as.list(env),envir=environment())
				ic <- ic_func(xp, t0, t1, fit_func(ydat, ...), ci_p=ci_p) 
		}
}
print(get_ic(gauss_fit, prof_ic)(env))
print(get_ic(gpd_fit, prof_ic)(env, qthreshold=0.90))
print(get_ic(gev_fit, boot_ic)(env))
print(get_ic(gpd_fit, boot_ic)(env, qthreshold=0.90))
print(get_ic(gauss_fit, boot_ic)(env))

cov_pgpd <- coverage(set_env_gauss, gen_data_gauss, get_theo_gauss, get_ic(gpd_fit, prof_ic)) 
cov_pgev <- coverage(set_env_gauss, gen_data_gauss, get_theo_gauss, get_ic(gev_fit, prof_ic)) 
cov_pgauss <- coverage(set_env_gauss, gen_data_gauss, get_theo_gauss, get_ic(gauss_fit, prof_ic)) 
cov_pgpd(2, qthreshold=0.9)
cov_pgev(2)
cov_pgauss(2)


cov_bgpd <- coverage(set_env_gauss, gen_data_gauss, get_theo_gauss, get_ic(gpd_fit, boot_ic)) 
cov_bgev <- coverage(set_env_gauss, gen_data_gauss, get_theo_gauss, get_ic(gev_fit, boot_ic)) 
cov_bgauss <- coverage(set_env_gauss, gen_data_gauss, get_theo_gauss, get_ic(gauss_fit, boot_ic)) 
cov_bgpd(2, qthreshold=0.9)
cov_bgev(2)
cov_bgauss(2)

set_env_gauss_param <- function(repet=4, xp=106, t0=1900, t1=2150){
		function(){
				this_env <- environment()
				list2env(as.list(parent.env(this_env)),envir=environment())
				changement <- 100
				years <- 1850:2200
				n <- length(unique(years))
				t <- seq(1,n)
				mu <- (t>changement)*.020*(t-changement)+100
				sigma <- (t>changement)*.001*(t-changement)+1
				sigma <- sqrt(sigma)
				mu <- rep(mu,repet)
				sigma <- rep(sigma,repet)
				t <- rep(t,repet)
				year <- rep(years,repet)
				covariate <- (mapply(qnorm,mean=mu,sd=sigma,MoreArgs=list("p"=0.5)))
				ydat <- data.frame("year"=year, "y"=NA, "mu_var"=covariate, "sig_var"=covariate)
				this_env
		}
}
test_several_conf <- function(l_repet=1:3, l_xp=c(103, 106), l_t0=1900, l_t1=c(2100, 2150), l_thresh=c(0.9, 0.95)){
		confs <- expand.grid(l_repet, l_xp, l_t0, l_t1, l_thresh)
		names(confs)=c("repet","xp","t1","t0","qthreshold")
		treat_conf <- function(repet, xp, t0, t1, thresh){
				cov_func <- coverage(set_env_gauss_param(repet, xp, t0, t1), gen_data_gauss, get_theo_gauss, get_ic(gpd_fit, boot_ic))
				cov_func(2, qthreshold=thresh)[[1]]
		}	
		res=with(confs, mapply(treat_conf, repet=repet, xp=xp, t0=t0, t1=t1, thresh=qthreshold))
		res=cbind(confs, res)
		res
}
tsc <- test_several_conf() 
		

