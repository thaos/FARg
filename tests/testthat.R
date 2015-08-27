rm(list=ls())
library("FARg")
library(quantreg)

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

p_gpd <- prof_ic(xp, t0, t1, gpd_fit(ydat, qthreshold=0.9), ci_p=0.95 ,to_plot=TRUE)
p_gev <- prof_ic(xp, t0, t1, gev_fit(ydat), ci_p=0.95 ,to_plot=TRUE)
p_gauss <- prof_ic(xp, t0, t1, gauss_fit(ydat), ci_p=0.95 ,to_plot=TRUE)

b_gpd <- boot_ic(xp, t0, t1, gpd_fit(ydat, qthreshold=0.9), ci_p=0.95, under_threshold=TRUE)
b_gev <- boot_ic(xp, t0, t1, gev_fit(ydat), ci_p=0.95)
b_gauss <- boot_ic(xp, t0, t1, gauss_fit(ydat), ci_p=0.95)



