library(MASS)
library(extRemes)
library(devtools)
devtools::load_all("..") 

n <- 350
t <- seq(1,n)
year <- 1850 + t - 1
mu <- (t>100)*.01*(t-100)+100
sigma2 <- (t>100)*.01*(t-100)+1
sigma2 <- 1
y <- mapply(rnorm,1,mu,sqrt(sigma2))
covariate <- c(scale(mu))
ydat <- data.frame(y, year, covariate)
xp <- 102 
t0 <- 1870
t1 <- 2003
pnt1 <- set_pnt(t1, xp, ydat)
pnt0 <- set_pnt(t0, xp, ydat)

y_bc <- bc_fit(l_lambda=seq(-10,10,0.1), y, data=ydat, mu_mod=~covariate, sig2_mod=~covariate, to_plot=TRUE)
lambda <- y_bc$lambda
y_std <- y_bc$y_std
tc=tcplot(y_std,c(quantile(y_std,0.5),max(y_std)))
# Easton avec thy_bchold fixe 
threshold <-  select_thresh(tc) 
rate=mean(y_std >= threshold)
y_bc_fit <- gpd_fit(y_bc$y_std,  ydat, qthreshold=1-rate)
pnt1_bc <- transform_pnt(y_bc, pnt1)
pnt0_bc <- transform_pnt(y_bc, pnt0)
get_far(y_bc_fit, pnt0_bc, pnt1_bc)

get_far.bc_fit <- function(y_trans, y_fit, pnt0, pnt1, qthreshold=NULL){
  print(deparse(sys.call()))
  pnt1_bc_fit <- transform_pnt(y_trans, pnt1)
  pnt0_bc_fit <- transform_pnt(y_trans, pnt0)
  far <- get_far(y_fit, pnt0_bc_fit, pnt1_bc_fit)
  far
}
get_far.bc_fit(y_bc, y_bc_fit, pnt0, pnt1)

get_far.bc_fit <- function(y_fit, pnt0, pnt1, qthreshold=NULL){
  lambda <- y_fit$lambda
  y_std <- y_fit$y_std
  if(is.null(qthreshold)){
    tc <- tcplot(y_std)
    threshold <- select_thresh(tc)
    qthreshold <- mean(y_std < threshold)
    print(qthreshold)
  }
  # we suppose the boxcox transform made the time series stationnary
  y_std_fit <- gpd_fit(res$y_std,  ydat, qthreshold=qthreshold)
  plot(y_std_fit)
  pnt1_bc <- transform_pnt(y_fit, pnt1)
  pnt0_bc <- transform_pnt(y_fit, pnt0)
  far <- get_far(y_std_fit, pnt0_bc, pnt1_bc)
  attr(far, "gpd_fit") <- y_std_fit
  far
}
far_bc <- get_far.bc_fit(fit_bc, pnt0, pnt1)
far_bc <- get_far.bc_fit(fit_bc, pnt0, pnt1, qthreshold=.75)

boot_bc <- function(data, indices, y, pnt0, pnt1, qthreshold){
  data_b <- data[indices,]
  y_b <- y[indices]
  fit_bc <- bc_fit(l_lambda=seq(-10,10,0.1), y_b, data=data_b, mu_mod=~covariate, sig2_mod=~covariate, to_plot=FALSE)
  get_far.bc_fit(fit_bc, pnt0, pnt1, qthreshold)
}

bres_bc <- boot(ydat, boot_bc, R=100, y=y, pnt0=pnt0, pnt1=pnt1, qthreshold=attr(far_bc, "gpd_fit")$rq_fitted$tau)
colnames(bres_bc$t) <- names(far_bc)
ic_bc <- apply(bres_bc$t, 2, quantile, probs=c(0.05, 0.5, 0.95)) 


y_std <- standardize(y, data=ydat, mu_mod=~covariate, sig2_mod=~covariate, to_plot=TRUE)
tc=tcplot(y_std$y_std, c(quantile(y_std$y_std, 0.5), max(y_std$y_std)))
threshold <-  select_thresh(tc) 
rate=mean(y_std$y_std >= threshold)
pnt1_std <- transform_pnt(y_std, pnt1)
pnt0_std <- transform_pnt(y_std, pnt0)
y_std_fit <- gpd_fit(y_std$y_std ,  ydat, qthreshold=1-rate)
get_far(y_std_fit, pnt0_std, pnt1_std)

get_far.std <- function(y_trans, y_fit, pnt0, pnt1, qthreshold=NULL){
  print(deparse(sys.call()))
  pnt1_std <- transform_pnt(y_trans, pnt1)
  pnt0_std <- transform_pnt(y_trans, pnt0)
  far <- get_far(y_fit, pnt0_std, pnt1_std)
  far
}
get_far.std(y_std, y_std_fit, pnt0, pnt1)

boot_func.std <- function(y_trans, y_fit, indices, pnt0, pnt1){
  print(indices)
  y <- y_trans$y[indices]
  data <- y_trans$data[indices,]
  y_std <- standardize(y, data, mu_mod=y_trans$mu_mod, sig2_mod=y_trans$sig2_mod )
  fit_method <- eval(parse(text=class(y_fit))) 
  fit_inputs <- as.list(set_fparam(fit_method, y_fit)) 
  fit_inputs$y <- y_std$y_std 
  fit_inputs$data <- data 
  y_std_fit  <- do.call(fit_method, fit_inputs)
  get_far.std(y_std, y_std_fit, pnt0, pnt1)
}
bf.std <- boot_func.std(y_std, y_std_fit, indices=1:length(y), pnt0, pnt1)
sb_std <- simple_boot(y_std_fit, boot_func.std, y_trans=y_std, pnt0=pnt0, pnt1=pnt1)

boot_ic.std <- function(xp, t0, t1, y_trans, y_fit, ci_p=0.95, ...){
  pnt0 <- set_pnt(t0, xp, y_fit$data)
  pnt1 <- set_pnt(t1, xp, y_fit$data)
  far_mle <- get_far.std(y_trans, y_fit, pnt0, pnt1)
  boot_res <- simple_boot(y_std_fit, boot_func.std, y_trans=y_std, pnt0=pnt0, pnt1=pnt1)
  alpha <- 1-ci_p
  far_boot <- boot_res["FAR", ]
  ic_boot <- quantile(far_boot, probs=c(alpha/2, .5, 1-alpha/2))
  ans <- c(ic_boot, far_mle)
  names(ans)[1:3] <- c("IC_inf", "FAR_B", "IC_sup")
  ans
}
bic_std <- boot_ic.std(xp, t0, t1, y_std,  y_std_fit)

set_fparam <- function(f, gfit){
  l <- formals(f)
  for( p in names(l)){
    #     print('-------------------') 
    #     print(p)
    #     print(gfit[[p]] )
    l[[p]] <- gfit[[p]]
  }
  l
}

boot_std <- function(data, indices, y, pnt0, pnt1, qthreshold){
  data_b <- data[indices,]
  y_b <- y[indices]
  y_std <- standardize(y_b, data=data_b, mu_mod=~covariate, sig2_mod=~covariate, to_plot=FALSE)
  get_far.std(y_std, pnt0, pnt1, qthreshold)
}

bres_std <- boot(ydat, boot_std, R=100, y=y, pnt0=pnt0, pnt1=pnt1, qthreshold=attr(far_std, "gpd_fit")$rq_fitted$tau)
colnames(bres_std$t) <- names(far_std)
ic_std <- apply(bres_std$t, 2, quantile, probs=c(0.05, 0.5, 0.95)) 



# muscsh=findpars(y_std.fit)
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
