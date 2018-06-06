#' @importFrom quantreg rq
#' @importFrom quantreg predict.rq
NULL

#' get_q generic.
#'
#' Compute the quantile corresponding to the probability of exceedance p.
#' @param object an object of class gauss_fit, gev_fit, gpd_fit.
#' @param pnt a point which consists of a line of data.frame containing the same variables used for the fit. It can bet set using the function \code{set_pnt}
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
#'p <- 0.01
#'pnt1 <- set_pnt(t1, p, time_var="year", tas)
#'pnt0 <- set_pnt(t0, p, time_var="year", tas)
#'
#'get_q(ga_fit ,pnt1)
#'
#'get_q(gp_fit, pnt1)
#'
#'get_q(ge_fit, pnt1)
#' @export
get_q <- function(object, ...){
	UseMethod("get_q")
}

#' @describeIn get_q for a Gaussian fit.
#' @export
get_q.gauss_fit <- function(object, pnt, ...){
  par <- compute_par.gauss_fit(object, newdata=pnt[, -1])
	res <- qnorm(p=pnt[, 1], mean=par[, 1], sd=par[, 2], lower.tail=FALSE)
	res <- c(res, par[, 1], par[, 2])
	names(res)=c("q","mu","sigma")
	res
}

#' @describeIn get_q for a GPD fit.
#' @export
get_q.gpd_fit <- function(object, pnt, ...){
  p <- pnt[, 1]
  threshold <- predict(object$rq_fitted, newdata=pnt[, -1])
  par <- compute_par.gpd_fit(object, newdata=pnt[, -1])
	phi <- object$rate
	stopifnot(p <= phi)
  res <- qevd(p/phi, threshold=threshold, scale=par[, 2], shape=par[, 3], type="GP", lower.tail=FALSE)
	res <- c(res, threshold, par[, 2], par[ ,3])
	names(res) <- c("q","threshold","sigma","shape")
	res
}

#' @describeIn get_q for a GEV fit.
#' @export
get_q.gev_fit <- function(object, pnt, ...){
  par <- compute_par.gev_fit(object, newdata=pnt[, -1])
	res <- qevd(pnt[, 1], loc=par[, 1], scale=par[, 2], shape=par[, 3] ,lower.tail=FALSE)
	res <- c(res, par[, 1], par[ ,2], par[, 3])
	names(res) <- c("q","mu","sigma","shape")
	res
}

