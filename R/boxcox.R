library(MASS)
library(tseries)
library(evmix)
library(extRemes)

bc <- function(x,lambda){
  y=NA
  if(lambda!=0){y=(x^lambda-1)/lambda}
  if(lambda==0){y=log(x)}
  return(y)
} 

lambda_prof <- function(lambda, y, data, mu_mod=~1, sig2_mod=~1, init=NULL){
	y_lambda <- bc(y,lambda)
  y_fit <- gauss_fit(y_lambda, data, mu_mod, sig2_mod, init)
	y_fit$objective <- y_fit$objective+(1-lambda)*sum(log(y))
  class(y_fit) <- "bc_fit"
	y_fit
}

bc_fit=function(l_lambda, y, data, mu_mod=~1, sig2_mod=~1, ci_p=.95, to_plot=FALSE){
  stopifnot(all(y > 0))
  #for the stability of the boxcox trannsformation
  y <- y / exp(mean(log(y)))
  lambda_lik <- function(lambda){
    lambda_prof(lambda, y, data, mu_mod, sig2_mod)$objective
  }
  aalpha <- qchisq(ci_p, 1)
  l_prof <- -unlist(lapply(l_lambda, lambda_lik))
  overallmax <- max(l_prof)
  crit <- overallmax - qchisq(0.999, 1)/2
  cond <- l_prof > crit
  crit_ci <- overallmax - qchisq(ci_p, 1)/2
  cond_ci <- l_prof > crit_ci
  lambda <- optimize(lambda_lik, range(l_lambda), maximum=FALSE)$minimum
  res <- lambda_prof(lambda, y, data, mu_mod, sig2_mod)
  res$lambda <- lambda
  res$ci <- range(l_lambda[cond_ci])
  if(to_plot){
    par(mfrow=c(1,2))
    plot(l_lambda[cond], l_prof[cond],type="l",ylab="")
    abline(v=lambda, lty=2)
    abline(v=res$ci[1], lty=2)
    abline(v=res$ci[2], lty=2)
    abline(h=crit_ci, lty=2)
    boxcox(as.formula(complete_formula(y,mu_mod)),lambda=l_lambda)
  }
  res
}


zp2yp <- function(lambda,mu,sigma,zp) 
	(lambda*(mu+sigma*zp)+1)^(1/lambda)

getFAR.boxcox <- function(pt0,pt1,xp,ydat,L=seq(-10,10,0.1),to.plot=FALSE){
	print("**********************************************************")
  i0 <- min(which.min(pt0[1]-ydat$year))
  i1 <- min(which.min(pt1[1]-ydat$year))
	if(any(ydat$y<= 0)){
		offs=+ceiling(abs(min(ydat$y)))
		ydat$y=ydat$y+offs
		xp=xp+offs
	}
	covariate=ydat$mua
	sf=exp(mean(log(ydat$y)))
	#         res=with(ydat,bc.fit(mua,y,range=prange,to.plot=TRUE))
	res=with(ydat,bc.fit(mua,y,range=L))
	lambda=res$lambda
	print(lambda)
	#         res=with(ydat,lambda_lik(lambda,covariate,y/sf))
	#         res$lambda=lambda
	ylambda=g(ydat$y/sf,lambda)
	xplambda=g(xp/sf,lambda)
	mu.pred=with(res,par[1]+par[2]*covariate)
	sigma2.pred=with(res,par[3]+par[4]*covariate)
	mu.pred0=with(res,par[1]+par[2]*pt0[2])
	mu.pred1=with(res,par[1]+par[2]*pt1[2])
	residus=(ylambda-mu.pred)/sqrt(sigma2.pred)
	sigma2.pred0=with(res,par[3]+par[4]*pt0[2])
	sigma2.pred1=with(res,par[3]+par[4]*pt1[2])
	r0=(xplambda-mu.pred0)/sqrt(sigma2.pred0)
	r1=(xplambda-mu.pred1)/sqrt(sigma2.pred1)
	#         print(shapiro.test(residus))
	#         print(Box.test(residus))
	#         print(PP.test(c(residus)))
	#         print(adf.test(residus))
	tc=tcplot_sthao(residus)
	# Easton avec threshold fixe 
	threshold=select.mu(tc)
	if(is.finite(threshold)){
		print("Path 1")
		rate=mean(residus>=threshold)
		#                 print(rate)
		residus.fit=fevd(residus,threshold=threshold,type="GP",time.units="years")
		residus.fit=fevd(residus, scale=~covariate, threshold=threshold,type="GP",time.units="years")
		mle=residus.fit$results$par
    pars <- findpars(residus.fit)
    #     p0=pevd(r0, scale=mle[1], shape=mle[2], threshold=threshold, type="GP",lower.tail=FALSE)
    #     p1=pevd(r1, scale=mle[1], shape=mle[2], threshold=threshold, type="GP",lower.tail=FALSE)
		p0=pevd(r0, scale=pars$scale[i0], shape=pars$shape[i0], threshold=threshold, type="GP",lower.tail=FALSE)
		p1=pevd(r1, scale=pars$scale[i1], shape=pars$shape[i1], threshold=threshold, type="GP",lower.tail=FALSE)
		#                 print(p0)
		#                 print(p1)
		p0=p0*rate
		p1=p1*rate
	}
	#         print(c(r0,r1))
	#         if(r0<threshold) p0=mean(r0<residus) # ou pnorm(r0,lower.tail=FALSE
	#         if(r0<threshold) p0=mean(r0<residus) # ou pnorm(r0,lower.tail=FALSE
		print("Direct 2 Path 2")
		#         if(r0<threshold) p0=pnorm(r0,lower.tail=FALSE)
		#         if(r1<threshold) p1=pnorm(r1,lower.tail=FALSE)
	if(r0<threshold) p0=mean(r0<residus) # ou pnorm(r0,lower.tail=FALSE
	if(r1<threshold) p1=mean(r1<residus) # ou pnorm(r1,lower.tail=FALSE
	if(to.plot){
		par(mfrow=c(2,2))
		plot(ydat$year,ydat$y)
		abline(h=xp)
		plot(ydat$year,ylambda)
		abline(h=xplambda)
		plot(ydat$year,residus)
		abline(h=r0,col="violet")
		abline(h=r1,col="blue")
		if(is.finite(threshold))
			abline(h=threshold,col="red")
		hist(residus)
		abline(v=r0,col="violet")
		abline(v=r1,col="blue")
		if(is.finite(threshold))
			abline(v=threshold,col="red")
	}
	far=FAR(p0,p1)
	#         print(far)
	#         if(far < 0)
	#                 print(res$par)
	far
}
