rm(list=ls())
library("FARg")
library(quantreg)

data(tas)
ydat <- with(tas,as.data.frame(cbind(year, eur_tas, avg_gbl_tas, avg_gbl_tas)))
names(ydat) <- c("year", "y", "mu_var", "sig_var")

gev_fevd(y, ydat, mu_mod=~mu_var, sig_mod=~sig_var)
ge_fit <- gev_fit(tas$eur_tas, data=tas, mu_mod=~avg_gbl_tas, sig_mod=~avg_gbl_tas)
plot(ge_fit)
rq_fitted <- rq(y~mu_var, data=ydat, tau=0.90)
threshold <- predict(rq_fitted)
gpd_fevd(y, ydat, threshold=threshold, sig_mod=~1)
gpd_fevd(y, ydat, threshold=threshold, sig_mod=~sig_var)
gp_fit <- gpd_fit(tas$eur_tas, data=tas, mu_mod=~avg_gbl_tas, sig_mod=~avg_gbl_tas, qthreshold=0.9)
compute_par(gp_fit, tas)
plot(gp_fit)
ga_fit <- gauss_fit(tas$eur_tas, data=tas, mu_mod=~avg_gbl_tas, sig2_mod=~avg_gbl_tas)
compute_par(ga_fit, tas)
plot(ga_fit)

y=tas$avg_gbl_tas


t1 <- 2003
t0 <- 1990
xp <- 1.6
pnt1 <- set_pnt(t0, xp, ydat)
pnt0 <- set_pnt(t0, xp, ydat)

get_p(gauss_fit(ydat), pnt1)
get_far(gauss_fit(ydat), pnt0, pnt1)

get_p(gpd_fit(ydat, qthreshold=0.9), pnt1)
get_far(gpd_fit(ydat, qthreshold=0.9), pnt0, pnt1)
get_far(gpd_fit(ydat, qthreshold=0.9), pnt0, pnt1, under_threshold=TRUE)

get_p(gev_fit(ydat), pnt1)
get_far(gev_fit(ydat), pnt0, pnt1)

p_gpd <- prof_ic(xp, t0, t1, gpd_fit(ydat, qthreshold=0.9), ci_p=0.95 ,to_plot=TRUE)
p_gev <- prof_ic(xp, t0, t1, gev_fit(ydat), ci_p=0.95 ,to_plot=TRUE)
p_gauss <- prof_ic(xp, t0, t1, gauss_fit(ydat), ci_p=0.95 ,to_plot=TRUE)

b_gpd <- boot_ic(xp, t0, t1, gpd_fit(ydat, qthreshold=0.9), ci_p=0.95, under_threshold=TRUE)
b_gev <- boot_ic(xp, t0, t1, gev_fit(ydat), ci_p=0.95)
b_gauss <- boot_ic(xp, t0, t1, gauss_fit(ydat), ci_p=0.95)



