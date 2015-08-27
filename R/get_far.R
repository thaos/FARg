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
	m_cord=pnt[1]
	s_cord=pnt[2]
	y_cord=pnt[3]
	mle <- y_fit$par
	mu <- mle[1] + m_cord * mle[2]
	sig <- mle[3] + s_cord * mle[4]
	res=pnorm(y_cord, mean=mu, sd=sqrt(sig), lower.tail=FALSE)
	res=c(res, mu, sig)
	names(res)=c("p","mu","sigma2")
	res
}

#' @export
get_p.gpd_fit <- function(y_fit, pnt, under_threshold=FALSE, ...){
	m_cord <- pnt[1]
	s_cord <- pnt[2]
	y_cord <- pnt[3]
	newdat <- data.frame("mu_var"=m_cord)
	threshold <- predict(y_fit$rq_fitted, newdata=newdat)
	stopifnot(under_threshold | y_cord >= threshold)
	mle <- y_fit$par
	sig <- mle[1] + s_cord * mle[2]
	sha <- mle[3]
	phi <- y_fit$rate
	if (y_cord > threshold)
		res <- pevd(y_cord, threshold=threshold, scale=sig, shape=sha, type="GP", lower.tail=FALSE) * phi
	else {
		findP <- function(par, m_cord, y_cord){
			rq_fitted <- rq(y~mu_var, data=y_fit$ydat, tau=par)
			ndata <- data.frame("mu_var"=m_cord)
			predicted <- predict.rq(rq_fitted, newdata=ndata)
			abs(y_cord-predicted)
		}
		res <- optimize(findP, interval=c(0,1-phi), m_cord=m_cord, y_cord=y_cord, tol=0.001)
		res <- res$minimum
		res <- 1-res
	}		
	res <- c(res, threshold, sig, sha)
	names(res) <- c("p","treshold","sigma","shape")
	res
}

#' @export
get_p.gev_fit <- function(y_fit, pnt, ...){
	m_cord <- pnt[1]
	s_cord <- pnt[2]
	y_cord <- pnt[3]
	mle <- y_fit$par
	mu <- mle[1] + m_cord * mle[2]
	sig <- mle[3] + s_cord * mle[4]
	sha <- mle[5]
	res <- pevd(y_cord, loc=mu, scale=sig, shape=sha ,lower.tail=FALSE)
	res <- c(res, mu, sig, sha)
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
	ifelse(p0[1] == 0 & p1[1] == 0, FAR <- 1, FAR <- 1-(p0[1]/p1[1]))
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
	ifelse(p0[1] == 0 & p1[1] == 0, FAR <- 1, FAR <- 1-(p0[1]/p1[1]))
	res <- c(FAR,p0,p1)
	names(res) <- c("FAR", name0, name1)
	res
}

#' @export
set_pnt <- function(year, y, ydat){
	i <- min(which.min(abs(year -ydat$year)))
	c(ydat$mu_var[i], ydat$sig_var[i],y)
}
