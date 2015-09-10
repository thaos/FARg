#' @export
#' 
#' @title Parameter Threshold Stability Plots
#'
#' @description Plots the MLE of the GPD parameters against threshold
#'
#' @inheritParams mrlplot
#' @param ylim.xi     y-axis limits for shape parameter or \code{NULL}
#' @param ylim.sigmau y-axis limits for scale parameter or \code{NULL}
#' 
#' @details The MLE of the (modified) GPD scale and shape (xi) parameters are
#'   plotted against a set of possible thresholds. If the GPD is a suitable
#'   model for a threshold \eqn{u} then for all higher thresholds \eqn{v > u} it
#'   will also be suitable, with the shape and modified scale being
#'   constant. Known as the threshold stability plots (Coles, 2001). The modified
#'   scale parameter is \eqn{\sigma_u - u\xi}.
#' 
#' In practice there is sample uncertainty in the parameter estimates, which
#' must be taken into account when choosing a threshold.
#' 
#' The usual asymptotic Wald confidence intervals are shown based on the
#' observed information matrix to measure this uncertainty. The sampling density
#' of the Wald normal approximation is shown by a greyscale image, where lighter
#' greys indicate low density.
#' 
#' A pre-chosen threshold (or more than one) can be given in \code{try.thresh}.
#' The GPD is fitted to the excesses using maximum likelihood estimation. The
#' estimated parameters are shown as a horizontal line which is solid above this
#' threshold, for which they should be the same if the GPD is a good model (upto sample uncertainty).
#' The threshold should always be chosen to be as low as possible to reduce sample uncertainty.
#' Therefore, below the pre-chosen threshold, where the GPD should not be a good model, the line
#' is dashed and the parameter estimates should now deviate from the dashed line
#' (otherwise a lower threshold could be used).
# 
#' If no threshold limits are provided \code{tlim = NULL} then the lowest threshold is set
#' to be just below the median data point and the maximum threshold is set to the 11th
#' largest datapoint. This is a slightly lower order statistic compared to that used in the MRL plot 
#' \code{\link[evmix:mrlplot]{mrlplot}} function to account for the fact the maximum likelihood
#' estimation is likely to be unreliable with 10 or fewer datapoints.
#' 
#' The range of permitted thresholds is just below the minimum datapoint and the
#' second largest value. If there are less unique values of data within the threshold
#' range than the number of threshold evalations requested, then instead of a sequence
#' of thresholds they will be set to each unique datapoint, i.e. MLE will only be applied
#' where there is data.
#' 
#' The missing (\code{NA} and \code{NaN}) and non-finite values are ignored.
#' 
#' The lower x-axis is the threshold and an upper axis either gives the number of 
#' exceedances (\code{p.or.n = FALSE}) or proportion of excess (\code{p.or.n = TRUE}).
#' Note that unlike the \code{gpd} related functions the missing values are ignored, so
#' do not add to the lower tail fraction. But ignoring the missing values is consistent
#' with all the other mixture model functions.
#' 
#' @return \code{\link[evmix:tcplot]{tshapeplot}} and 
#' \code{\link[evmix:tcplot]{tscaleplot}} produces the threshold stability plot for the
#' shape and scale parameter respectively. They also returns a matrix containing columns of
#' the threshold, number of exceedances, MLE shape/scale
#' and their standard devation and \eqn{100(1 - \alpha)\%} Wald confidence interval if requested. Where the
#' observed information matrix is not obtainable the standard deviation and confidence intervals
#' are \code{NA}. For the \code{\link[evmix:tcplot]{tscaleplot}} the modified scale quantities
#' are also provided. \code{\link[evmix:tcplot]{tcplot}} produces both plots on one graph and
#' outputs a merged dataframe of results.
#' 
#' @note If the user specifies the threshold range, the thresholds above the sixth
#' largest are dropped. A warning message is given if any thresholds have at most 10
#' exceedances, in which case the maximum likelihood estimation is unreliable. If there
#' are less than 10 exceedances of the minimum threshold then the function will stop.
#' 
#' By default, no legend is included when using \code{\link[evmix:tcplot]{tcplot}} to get
#' both threshold stability plots.
#' 
#' Error checking of the inputs (e.g. invalid probabilities) is carried out and
#' will either stop or give warning message as appropriate.
#' 
#' @references
#' 
#' Scarrott, C.J. and MacDonald, A. (2012). A review of extreme value
#' threshold estimation and uncertainty quantification. REVSTAT - Statistical
#' Journal 10(1), 33-59. Available from \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}
#' 
#' Coles S.G. (2004). An Introduction to the Statistical Modelling of Extreme Values.
#' Springer-Verlag: London.
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @section Acknowledgments: Based on the threshold stability plot function \code{\link[evd:tcplot]{tcplot}} in the 
#' \code{\link[evd:fpot]{evd}} package for which Stuart Coles' and Alec Stephenson's 
#' contributions are gratefully acknowledged.
#' They are designed to have similar syntax and functionality to simplify the transition for users of these packages.
#'   
#' @seealso \code{\link[evmix:mrlplot]{mrlplot}} and \code{\link[evd:tcplot]{tcplot}} from 
#' \code{\link[evd:mrlplot]{evd}} library
#' @aliases tcplot tshapeplot tscaleplot
#' @family tcplot
#' 
#' @examples
#' \dontrun{
#' x = rnorm(1000)
#' tcplot(x)
#' tshapeplot(x, tlim = c(0, 2))
#' tscaleplot(x, tlim = c(0, 2), try.thresh = c(0.5, 1, 1.5))
#' tcplot(x, tlim = c(0, 2), try.thresh = c(0.5, 1, 1.5))
#' }

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

#' @export
select_thresh <- function(tc){
	i.na <- which(apply(tc,1,function(x)any(is.na(x))))
	ifelse(length(i.na)!= 0,{tcc <- tc[-c(min(i.na):nrow(tc)),]},{tcc <- tc})
	u2rm <- which(apply(tcc,1,function(x)x["nu"]<11))
	if(length(u2rm)!=0) tcc <- tcc[-u2rm,]
	list2env(tcc,env <- environment())
	found <- FALSE
	i <- 1
	while(!found){
		if(i> length(u)){
			message("no threshold found")
			return(Inf)
		}
		u.cur <- u[i]
		sig.cur <- mod.sigmau[i]
		xi.cur <- xi[i]
		i2r <- seq(1,i,1)
		tsig <- all(sig.cur > mod.cil.sigmau[-i2r] & sig.cur < mod.ciu.sigmau[-i2r])
		txi <- all(xi.cur >cil.xi[-i2r] & xi.cur < ciu.xi[-i2r])
		found <- tsig & txi
		if(!found) i <- i+1
	}
	u.cur
}
