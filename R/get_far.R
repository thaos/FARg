#' @importFrom quantreg rq
#' @importFrom quantreg predict.rq

#' @export
get_far <- function(y_fit, ...){
	UseMethod("get_far")
}

#' @export
get_p <- function(y_fit, ...){
	UseMethod("get_p")
}

#' @export
get_p.gauss_fit <- function(y_fit, pnt, ...){
  par <- compute_par.gauss_fit(y_fit, newdata=pnt[, -1])
	res <- pnorm(pnt$y, mean=par$mu, sd=sqrt(par$sig2), lower.tail=FALSE)
	res <- c(res, par$mu, par$sig)
	names(res)=c("p","mu","sigma2")
	res
}

#' @export
get_p.gpd_fit <- function(y_fit, pnt, under_threshold=FALSE, ...){
  y_cord <- pnt$y
	threshold <- predict(y_fit$rq_fitted, newdata=pnt[, -1])
	stopifnot(under_threshold | y_cord >= threshold)
  # A remplacer !!!
  par <- compute_par.gpd_fit(y_fit, newdata=pnt[, -1])
	phi <- y_fit$rate
	if (y_cord > threshold)
		res <- pevd(y_cord, threshold=threshold, scale=par$sig, shape=par$xi, type="GP", lower.tail=FALSE) * phi
	else {
		findP <- function(par, pnt, y_fit){
      rq_fitted <- rq(as.formula(complete_formula(y_fit$y, y_fit$mu_mod)), data=y_fit$data, tau=par)
    predicted <- predict.rq(rq_fitted, newdata=pnt[, -1])
			abs(y_cord-predicted)
		}
		res <- optimize(findP, interval=c(0,1-phi), pnt=pnt, y_fit=y_fit, tol=0.001)
		res <- res$minimum
		res <- 1-res
	}		
	res <- c(res, threshold, par$sig, par$xi)
	names(res) <- c("p","treshold","sigma","shape")
	res
}

#' @export
get_p.gev_fit <- function(y_fit, pnt, ...){
  par <- compute_par.gev_fit(y_fit, newdata=pnt[, -1])
	res <- pevd(pnt$y, loc=par$mu, scale=par$sig, shape=par$xi ,lower.tail=FALSE)
	res <- c(res, par$mu, par$sig, par$xi)
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

#' @export
####### normaly not useful, the simple get_far handle the three cases
get_far.gpd_fit <- function(y_fit, pnt0, pnt1, under_threshold=FALSE, ...){
	p0 <- get_p.gpd_fit(y_fit, pnt0, under_threshold=under_threshold)
	p1 <- get_p.gpd_fit(y_fit, pnt1, under_threshold=under_threshold)
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

#' @export
get_far.default <- function(y_fit, pnt0, pnt1, ...){
	p0 <- get_p(y_fit, pnt0, ...)
	p1 <- get_p(y_fit, pnt1, ...)
	name0 <- paste(names(p0),"0",sep="")
	name1 <- paste(names(p1),"1",sep="")
	ifelse(p0[1] == 0 & p1[1] == 0, FAR <- 0, FAR <- 1-(p0[1]/p1[1]))
	res <- c(FAR,p0,p1)
	names(res) <- c("FAR", name0, name1)
	res
}

#' @export
set_pnt <- function(time, y, time_var="", data=NULL){
  try({
      time_name <- time_var
      if(!is.character(time_name))
        time_name <- paste(deparse(substitute(time_name)), collapse=" ")
      time_formula <- paste(time_name, "~ 1", sep="")
      time_var <- as.numeric(model.frame(as.formula(time_formula), data=data)[[1]])
  })
  i <- min(which.min(abs(time - time_var)))
	cbind(y,data[i,])
}
