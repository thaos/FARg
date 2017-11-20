#' @importFrom quantreg rq
#' @importFrom quantreg predict.rq
NULL

#' get_far generic.
#' 
#' Compute the FAR.
#' 
#' Compute the Fraction of Attributable Risk, which compares the probabilities p0 and p1 of exeeding a threshold xp at two differents time (t0 and t1) of the times series : \eqn{FAR = 1 - p0/p1}. 
#' @param object an object of class gauss_fit, gev_fit, gpd_fit or of class trans. If object is of class trans, the argument y_fit must be passed.
#' @param pnt0 point at time t0. A point is defined by a time, a threshold xp and values of covariates used for the fit.
#' @param pnt1 point at time t1. A point is defined by a time, a threshold xp and values of covariates used for the fit.
#' @param under_threshold used for gpd_fit. It tells whether it provides an estimated when one of points is under the GPD threshold. In which case, gives the empirical probability obtained by quantile regression.
#' @param ... Arguments to be passed to methods,
#' @return a vector with the estimated confidence intervals for the far as well as the probability of exceeding the threshold xp and the models parameters a both time t0 and t1.
#' @examples
#'data(tas)
#'
#'ge_fit <- gev_fit(eur_tas, data=tas, mu_mod=~gbl_tas, sig_mod=~gbl_tas, time_var="year")
#'gp_fit <- gpd_fit(eur_tas, data=tas, mu_mod=~gbl_tas, sig_mod=~gbl_tas, time_var="year", qthreshold=0.9)
#'ga_fit <- gauss_fit(eur_tas, data=tas, mu_mod=~gbl_tas, sig_mod=~gbl_tas, time_var="year")
#'
#'t1 <- 2003
#'t0 <- 1990
#'xp <- 1.6
#'pnt1 <- set_pnt(t1, xp, time_var="year", tas)
#'pnt0 <- set_pnt(t0, xp, time_var="year", tas)
#'
#'get_p(ga_fit ,pnt1)
#'get_far(ga_fit, pnt0, pnt1)
#'
#'get_p(gp_fit, pnt1)
#'get_far(gp_fit, pnt0, pnt1)
#'
#'get_p(ge_fit, pnt1)
#'get_far(ge_fit, pnt0, pnt1)
#' @export
get_far <- function(object, ...){
	UseMethod("get_far")
}

#' get_p generic.
#'
#' Compute the probability of exeeding a given point.
#' @param object an object of class gauss_fit, gev_fit, gpd_fit.
#' @param pnt a point which consists of a line of data.frame containing the same variables used for the fit. It can bet set using the function \code{set_pnt}
#' @param under_threshold used for gpd_fit. It tells whether it provides an estimated when the point is under the GPD threshold. In which case, gives the empirical probability obtained by quantile regression.
#' @param ... Arguments to be passed to methods,
#' @examples
#'data(tas)
#'
#'ge_fit <- gev_fit(eur_tas, data=tas, mu_mod=~gbl_tas, sig_mod=~gbl_tas, time_var="year")
#'gp_fit <- gpd_fit(eur_tas, data=tas, mu_mod=~gbl_tas, sig_mod=~gbl_tas, time_var="year", qthreshold=0.9)
#'ga_fit <- gauss_fit(eur_tas, data=tas, mu_mod=~gbl_tas, sig_mod=~gbl_tas, time_var="year")
#'
#'t1 <- 2003
#'t0 <- 1990
#'xp <- 1.6
#'pnt1 <- set_pnt(t1, xp, time_var="year", tas)
#'pnt0 <- set_pnt(t0, xp, time_var="year", tas)
#'
#'get_p(ga_fit ,pnt1)
#'get_far(ga_fit, pnt0, pnt1)
#'
#'get_p(gp_fit, pnt1)
#'get_far(gp_fit, pnt0, pnt1)
#'
#'get_p(ge_fit, pnt1)
#'get_far(ge_fit, pnt0, pnt1)
#' @export
get_p <- function(object, ...){
	UseMethod("get_p")
}

#' @describeIn get_p for a gaussian fit.
#' @export
get_p.gauss_fit <- function(object, pnt, ...){
  par <- compute_par.gauss_fit(object, newdata=pnt[, -1])
	res <- pnorm(pnt$y, mean=par[, 1], sd=par[, 2], lower.tail=FALSE)
	res <- c(res, par[, 1], par[, 2])
	names(res)=c("p","mu","sigma")
	res
}

#' @describeIn get_p for a GPD fit.
#' @export
get_p.gpd_fit <- function(object, pnt, under_threshold=FALSE, ...){
  y_cord <- pnt$y
	threshold <- predict(object$rq_fitted, newdata=pnt[, -1])
	stopifnot(under_threshold | y_cord >= threshold)
  # A remplacer !!!
  par <- compute_par.gpd_fit(object, newdata=pnt[, -1])
	phi <- object$rate
	if (y_cord > threshold)
		res <- pevd(y_cord, threshold=threshold, scale=par[, 2], shape=par[, 3], type="GP", lower.tail=FALSE) * phi
	else {
		findP <- function(par, pnt, object){
      rq_fitted <- rq(as.formula(complete_formula(object$y, object$mu_mod)), data=object$data, tau=par)
    predicted <- predict.rq(rq_fitted, newdata=pnt[, -1])
			abs(y_cord-predicted)
		}
		res <- optimize(findP, interval=c(0,1-phi), pnt=pnt, object=object, tol=0.001)
		res <- res$minimum
		res <- 1-res
	}		
	res <- c(res, threshold, par[, 2], par[ ,3])
	names(res) <- c("p","threshold","sigma","shape")
	res
}

#' @describeIn get_p for a GEV fit.
#' @export
get_p.gev_fit <- function(object, pnt, ...){
  par <- compute_par.gev_fit(object, newdata=pnt[, -1])
	res <- pevd(pnt$y, loc=par[, 1], scale=par[, 2], shape=par[, 3] ,lower.tail=FALSE)
	res <- c(res, par[, 1], par[ ,2], par[, 3])
	names(res) <- c("p","mu","sigma","shape")
	res
}

# get_far.gauss_fit <- function(y_fit, pnt0, pnt1){
#   p0 <- get_p.gauss_fit(y_fit, pnt0)
#   p1 <- get_p.gauss_fit(y_fit, pnt1)
#   ifelse(p0[1] == 0 & p1[1] == 0, far <- 1, far <- 1-(p0[1]/p1[1]))
#   res <- c(far,p0,p1)
#   names(res) <- c("far","p0","mu0","sc0","p1","mu1","sc1")
#   res
# }

#' @describeIn get_far for a gpd fit.
#' @export
####### normaly not useful, the simple get_far handle the three cases
get_far.gpd_fit <- function(object, pnt0, pnt1, under_threshold=FALSE, ...){
	p0 <- get_p.gpd_fit(object, pnt0, under_threshold=under_threshold)
	p1 <- get_p.gpd_fit(object, pnt1, under_threshold=under_threshold)
	ifelse(p0[1] == 0 & p1[1] == 0, FAR <- 0, FAR <- 1-(p0[1]/p1[1]))
	res <- c(FAR,p0,p1)
	names(res) <- c("FAR", "p0", "thresh0", "sc0", "sh0", "p1","thresh1", "sc1", "sh1")
	res
}

# get_far.gev_fit <- function(y_fit, pnt0, pnt1){
#   p0 <- get_p.gev_fit(y_fit, pnt0)
#   p1 <- get_p.gev_fit(y_fit, pnt1)
#   ifelse(p0[1] == 0 & p1[1] == 0, FAR <- 1, FAR <- 1-(p0[1]/p1[1]))
#   res <- c(FAR,p0,p1)
#   names(res) <- c("FAR", "p0", "mu0", "sc0", "sh0", "p1","mu1", "sc1", "sh1")
#   res
# }

#' @describeIn get_far for a gaussian fit or a gev fit.
#' @export
get_far.default <- function(object, pnt0, pnt1, ...){
	p0 <- get_p(object, pnt0, ...)
	p1 <- get_p(object, pnt1, ...)
	name0 <- paste(names(p0),"0",sep="")
	name1 <- paste(names(p1),"1",sep="")
	ifelse(p0[1] == 0 & p1[1] == 0, FAR <- 0, FAR <- 1-(p0[1]/p1[1]))
	res <- c(FAR,p0,p1)
	names(res) <- c("FAR", name0, name1)
	res
}

#' Create a point. 
#' 
#' Create a point at a given time from a given data. 
#' 
#' The function extract the point of the data set, time of which is the closest to the one specified.
#' @param time the time we want the point to be closest to.
#' @param y the value we want to set the point at the given time.
#' @param time_var the variable used as time in data
#' @param data the data.frame containing the relevent information to create the point
#' @examples
#'data(tas)
#'
#'ge_fit <- gev_fit(eur_tas, data=tas, mu_mod=~gbl_tas, sig_mod=~gbl_tas, time_var="year")
#'gp_fit <- gpd_fit(eur_tas, data=tas, mu_mod=~gbl_tas, sig_mod=~gbl_tas, time_var="year", qthreshold=0.9)
#'ga_fit <- gauss_fit(eur_tas, data=tas, mu_mod=~gbl_tas, sig_mod=~gbl_tas, time_var="year")
#'
#'t1 <- 2003
#'t0 <- 1990
#'xp <- 1.6
#'pnt1 <- set_pnt(t1, xp, time_var="year", tas)
#'pnt0 <- set_pnt(t0, xp, time_var="year", tas)
#'
#'get_p(ga_fit ,pnt1)
#'get_far(ga_fit, pnt0, pnt1)
#'
#'get_p(gp_fit, pnt1)
#'get_far(gp_fit, pnt0, pnt1)
#'
#'get_p(ge_fit, pnt1)
#'get_far(ge_fit, pnt0, pnt1)
#' @export
set_pnt <- function(time, y, time_var="", data=NULL){
  try({
      time_name <- time_var
      if(!is.character(time_name))
        time_name <- paste(deparse(substitute(time_name)), collapse=" ")
      time_var <- data[,time_name]
  })
  i <- min(which.min(abs(time - time_var)))
	cbind(y,data[i,])
}
