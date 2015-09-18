rm(list=ls())
library(quantreg)

data(tas)
ydat <- with(tas,as.data.frame(cbind(year, eur_tas, avg_gbl_tas, avg_gbl_tas)))
names(ydat) <- c("year", "y", "mu_var", "sig_var")

ge_fit <- gev_fit(tas$eur_tas, data=tas, mu_mod=~avg_gbl_tas, sig_mod=~avg_gbl_tas, time_var="year")
plot(ge_fit)
rq_fitted <- rq(y~mu_var, data=ydat, tau=0.90)
threshold <- predict(rq_fitted)
gp_fit <- gpd_fit(tas$eur_tas, data=tas, mu_mod=~avg_gbl_tas, sig_mod=~avg_gbl_tas, time_var="year", qthreshold=0.9)
compute_par(gp_fit, tas)
plot(gp_fit)
ga_fit <- gauss_fit(tas$eur_tas, data=tas, mu_mod=~avg_gbl_tas, sig2_mod=~avg_gbl_tas, time_var="year")
compute_par(ga_fit, tas)
plot(ga_fit)

y=tas$avg_gbl_tas
ydat_p=tas[,-2]

t1 <- 2003
t0 <- 1990
xp <- 1.6
pnt1 <- set_pnt(t1, xp, time_var="year", ydat_p)
pnt0 <- set_pnt(t0, xp, time_var="year", ydat_p)

get_p(ga_fit ,pnt1)
get_far(ga_fit, pnt0, pnt1)

get_p(gp_fit, pnt1)
pnt1.1 <- set_pnt(t0, .6, time_var="year", ydat_p)
get_p(gp_fit, pnt1.1, under_threshold=TRUE)
get_far(gp_fit, pnt0, pnt1)
get_far(gp_fit, pnt0, pnt1, under_threshold=TRUE)

get_p(ge_fit, pnt1)
get_far(ge_fit, pnt0, pnt1)

p_gpd <- prof_ic(gp_fit, xp, t0, t1, ci_p=0.95 ,to_plot=TRUE)
p_gev <- prof_ic(ge_fit, xp, t0, t1, ge_fit, ci_p=0.95 ,to_plot=TRUE)
p_gauss <- prof_ic(ga_fit, t0, t1, ga_fit, ci_p=0.95 ,to_plot=TRUE)

b_gpd <- boot_ic(gp_fit, xp, t0, t1, ci_p=0.95, under_threshold=TRUE)
b_gev <- boot_ic(ge_fit, xp, t0, t1, ci_p=0.95)
b_gauss <- boot_ic(ga_fit, xp, t0, t1, ci_p=0.95)
