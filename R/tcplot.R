# Based on tcplot from evmix : only removed plot generation.
tcplot <- function(data, tlim=NULL, nt=min(100, length(data)), alpha=0.05) {
  invisible(nt)
  # make sure defaults which result from function evaluations are obtained
  data  <-  sort(data)
  if (is.null(tlim)) {
    tlim <- c(median(data) - 2*.Machine$double.eps, data[length(data) - 11])
  }
  # Check properties of inputs
  thresholds <- seq(tlim[1], tlim[2], length.out = nt)
  n <- length(data)
  data <- data[data > min(thresholds)]
  # Trick to evaluate MRL at all datapoints if there are not too many
  udata <- unique(data)
  if (length(udata) <= nt) {
    warning("less data than number of thresholds requested, so will use unique data as thresholds")
    thresholds <- udata[-length(udata)]
  }
  # Check given thresholds
  nminu <- sum(data > min(thresholds))
  if (nminu <= 10)
    stop("data must have more than 10 exceedances of lowest threshold")
  nmaxu <- sum(data > max(thresholds))
  if (nmaxu <= 5) {
    warning("thresholds above 6th largest input data are dropped")
    thresholds <- thresholds[thresholds < data[length(data) - 5]]
    nmaxu <- sum(data > max(thresholds))
  }
  if (nmaxu <= 10) warning("maximum likelihood estimation is unreliable with less than 10 exceedances")
  nt <- length(thresholds)
  mle.calc <- function(x, u, alpha) {
    gpdmle <- fgpd(x, u)  
    if (is.null(gpdmle$se)) stop("no se estimation : cant build CI") #gpdmle$se = rep(NA, 2)
    if (is.null(gpdmle$cov)) {
      gpdmle$cov12 <- NA
    } else {
      gpdmle$cov12 <- gpdmle$cov[1, 2]
    }
    results <- c(u, sum(x > u), gpdmle$mle, gpdmle$se, gpdmle$sigmau + qnorm(c(alpha/2, 1 - alpha/2)) * gpdmle$se[1],
                  gpdmle$xi + qnorm(c(alpha/2, 1 - alpha/2)) * gpdmle$se[2], gpdmle$cov12)      
    return(results)
  }
  mleresults <- matrix(NA, nrow = nt, ncol = 11)
  for (i in 1:nt) {
    tc<-try({mleresults[i,] <- mle.calc(data, thresholds[i], alpha)}, silent=TRUE)
    if(class(tc)=="try-error" ) {
	break
	}
  }  
  mleresults <- as.data.frame(mleresults)
  names(mleresults) <- c("u", "nu", "sigmau", "xi", "se.sigmau", "se.xi", 
      "cil.sigmau", "ciu.sigmau", "cil.xi", "ciu.xi", "cov12")
  mleresults$mod.sigmau <- mleresults$sigmau - mleresults$xi * mleresults$u
  mleresults$mod.se.sigmau <- sqrt(mleresults$se.sigmau^2 - 2 * mleresults$u * mleresults$cov12 + (mleresults$u * mleresults$se.xi)^2)
  mleresults$mod.cil.sigmau <- mleresults$mod.sigmau + qnorm(alpha/2) * mleresults$mod.se.sigmau
  mleresults$mod.ciu.sigmau <- mleresults$mod.sigmau + qnorm(1 - alpha/2) * mleresults$mod.se.sigmau
  return(mleresults)
}

#' Select threshold from tcplot output.
#' 
#' Select trheshold from threshold stability plot from the tcplot function of the evmix package.
#' 
#' The threhshold is selected such as as the minimum threshold for which the transformed parameters of the GPD model fitted are within the confidence intervals for every higher threshold.
#' @param tc result from the tcplot function of the evmix packages
#' @return returns the value of the threshold selected
#' @export
select_thresh <- function(tc){
	i.na <- which(apply(tc,1,function(x)any(is.na(x))))
	ifelse(length(i.na)!= 0,{tcc <- tc[-c(min(i.na):nrow(tc)),]},{tcc <- tc})
	u2rm <- which(apply(tcc,1,function(x)x["nu"]<11))
	if(length(u2rm)!=0) tcc <- tcc[-u2rm,]
	found <- FALSE
	i <- 1
	while(!found){
		if(i> length(tcc$u)){
			message("no threshold found")
			return(Inf)
		}
		u.cur <- tcc$u[i]
		sig.cur <- tcc$mod.sigmau[i]
		xi.cur <- tcc$xi[i]
		i2r <- seq(1,i,1)
		tsig <- all(sig.cur > tcc$mod.cil.sigmau[-i2r] & sig.cur < tcc$mod.ciu.sigmau[-i2r])
		txi <- all(xi.cur > tcc$cil.xi[-i2r] & xi.cur < tcc$ciu.xi[-i2r])
		found <- tsig & txi
		if(!found) i <- i+1
	}
	u.cur
}
