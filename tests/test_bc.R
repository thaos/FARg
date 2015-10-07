library(FARg)
library(evmix)
# library(devtools)
# library(boot)
# devtools::load_all("..") 

n <- 350
t <- seq(1,n)
year <- 1850 + t - 1
mu <- (t>100)*.01*(t-100)+100
sigma2 <- (t>100)*.01*(t-100)+1
# sigma2 <- 1
y <- mapply(rnorm,1,mu,sqrt(sigma2))
covariate <- c(scale(mu))
ydat <- data.frame(y, year, covariate)
xp <- 102 
t0 <- 1870
t1 <- 2003
pnt1 <- set_pnt(t1, xp, time_var=1:nrow(ydat), ydat)
pnt0 <- set_pnt(t0, xp, time_var=1:nrow(ydat), ydat)
pnt1 <- set_pnt(t1, xp, time_var=1:nrow(ydat), NULL)
pnt0 <- set_pnt(t0, xp, time_var=1:nrow(ydat), NULL)
pnt1 <- set_pnt(t1, xp, time_var="year", ydat)
pnt0 <- set_pnt(t0, xp, time_var="year", ydat)
pnt1 <- set_pnt(t1, xp, time_var=year, NULL)
pnt0 <- set_pnt(t0, xp, time_var=year, NULL)
pnt1 <- set_pnt(t1, xp, time_var=year, ydat)
pnt0 <- set_pnt(t0, xp, time_var=year, ydat)

get_far_theo <- function(pnt0, pnt1){
  get_p_theo <- function(pnt){
    i <- min(which.min(abs(ydat$year -  pnt$year)))
    print(ydat[i,])
    p <- pnorm(as.numeric(pnt[1]), mean=mu[i], sd=sqrt(sigma2[i]), lower.tail=FALSE)
  }
  p0 <- get_p_theo(pnt0)
  p1 <- get_p_theo(pnt1)
  if( p0 == 0 & p1 == 0)
    return(1)
  else
    return(1 - p0/p1)
}
get_far_theo(pnt0, pnt1)

y_bc <- bc_fit(l_lambda=seq(-10,10,0.1), y, data=ydat, mu_mod=~covariate, sig2_mod=~covariate, time_var="year", to_plot=TRUE)
lambda <- y_bc$lambda
y_std <- y_bc$y_std
tc <- tcplot(y_std,c(quantile(y_std,0.5),max(y_std)))
# Easton avec thy_bchold fixe 
threshold <-  select_thresh(tc) 
rate=mean(y_std >= threshold)
y_bc_fit <- gpd_fit(y_bc$y_std,  ydat, mu_mod=~y, time_var="year", qthreshold=1-rate)
y_bc_fit <- gpd_fit(y_bc$y_std,  ydat, time_var="year", qthreshold=1-rate)
y_bc_fit <- gpd_fit(y_bc$y_std,  ydat, time_var=year, qthreshold=1-rate)
pnt1_bc <- transform_pnt(y_bc, pnt1)
pnt0_bc <- transform_pnt(y_bc, pnt0)
get_far(y_bc_fit, pnt0_bc, pnt1_bc, under_threshold=TRUE)

# get_far.trans(y_bc, y_bc_fit, pnt0, pnt1)
get_far(y_bc, y_bc_fit, pnt0, pnt1, under_threshold=TRUE)

# bf_bc <- boot_func.bc_fit(y_bc, y_bc_fit, indices=1:length(y), pnt0, pnt1)
# sb_bc <- simple_boot(y_bc_fit, boot_func.bc_fit, y_trans=y_bc, pnt0=pnt0, pnt1=pnt1)

# bic_bc <- boot_far.trans(y_bc,  y_bc_fit, xp, t0, t1)
bic_bc <- boot_far(y_bc,  y_bc_fit, xp, t0, t1, under_threshold=TRUE)
bic_bc <- boot_p(y_bc,  y_bc_fit, xp, t0, under_threshold=TRUE)


# y_std <- standardize(y, data=NULL, mu_mod=~covariate, sig2_mod=~covariate)
y_std <- standardize(y, data=ydat, mu_mod=~covariate, sig2_mod=~covariate)
tc=tcplot(y_std$y_std, c(quantile(y_std$y_std, 0.5), max(y_std$y_std)))
threshold <-  select_thresh(tc) 
rate=mean(y_std$y_std >= threshold)
pnt1_std <- transform_pnt(y_std, pnt1)
pnt0_std <- transform_pnt(y_std, pnt0)
y_std_fit <- gpd_fit(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22:nrow(ydat)) ,  ydat, time_var="year", qthreshold=1-rate)
try(y_std_fit <- gpd_fit(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22:nrow(ydat)) ,  data = NULL, time_var="year", qthreshold=1-rate))
y_std_fit <- gpd_fit(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22:nrow(ydat)) ,  data = ydat[,"covariate", drop=FALSE], time_var="year", qthreshold=1-rate)
y_std_fit <- gpd_fit(y_std$y_std ,  ydat, time_var="year", qthreshold=1-rate)
get_far(y_std_fit, pnt0_std, pnt1_std, under_threshold=TRUE)

# get_far.trans(y_std, y_std_fit, pnt0, pnt1)
get_far(y_std, y_std_fit, pnt0, pnt1, under_threshold=TRUE)

# bf_std <- boot_func.std(y_std, y_std_fit, indices=1:length(y), pnt0, pnt1)
# sb_std <- simple_boot(y_std_fit, boot_func.std, y_trans=y_std, pnt0=pnt0, pnt1=pnt1)

print(system.time({
  #   bic_std <- boot_far.trans(y_std,  y_std_fit, xp, t0, t1)
  bic_std <- boot_far(y_std,  y_std_fit, xp, t0, t1, under_threshold=TRUE)
  bic_std <- boot_p(y_std,  y_std_fit, xp, t0, under_threshold=TRUE)
}))

# muscsh=findpars(y_trans.fit)
# muscsh$mu=rep(threshold,length(muscsh$scale))
# z95f=evmix::qgpd(rep(0.95,length(muscsh$mu)),u=muscsh$mu,sigmau=muscsh$scale,xi=muscsh$shape,phiu=rate)
# y95f=zp2yp(lambda,mu=mu.pred,sigma=sigma.pred,zp=z95f)*sf
# q.theo95=mapply(qnorm,mean=mu,sd=sigma,MoreArgs=list("p"=0.95))
# 
# 
# print(shapiro.test(y_std))
# print(Box.test(y_std))
# print(PP.test(c(y_std)))
# irintiadf.test(y_std))
