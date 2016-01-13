#' Standardize the times series.
#' 
#' Trys to remove the trend and the heteroscedaticity of the time series estimated by linears models.
#' 
#' Compute the linear trend of the time series as well as  the evolution of the variance residuals with respect to the specified covariates and then standardize the time series accordingly.
#' @param y the time series values we want to remove the average trend and the heteroscedasticity from.
#' @param data a data.frame with the covariates needed to estimates the time series trend and variance.
#' @param mu_mod a formula giving the covariates the time-serie trend depends linearly on.
#' @param sig2_mod a formula giving the covariates the time-serie variance depends linearly on.
#' @return returns a list with the arguments of the function as well as :
#' \itemize{
#' \item y_std, the standardized time serie.
#' \item mu_fit, the linear model fitted for the trend.
#' \item sig2_fit, the linear model fitted for the variance.
#'}
#' @examples
#'data(tas)
#'y_std <- standardize(eur_tas, data=tas, mu_mod=~gbl_tas, sig2_mod=~gbl_tas)
#'y_std_fit <- gpd_fit(y_std$y_std ,  data=tas, time_var="year", qthreshold=0.8)
#'t1 <- 2003
#'t0 <- 1990
#'xp <- 1.6
#'pnt1 <- set_pnt(t1, xp, time_var="year", tas)
#'pnt0 <- set_pnt(t0, xp, time_var="year", tas)
#'get_far(y_std, y_std_fit, pnt0, pnt1, under_threshold=TRUE)
#'boot_far(y_std,  y_std_fit, xp, t0, t1, under_threshold=TRUE)
standardize <- function(y, data, mu_mod=~1, sig2_mod=~1){
  stopifnot(!is.null(data))
  y_name <- paste(deparse(substitute(y)), collapse="")
  #   if(is.element(y_name, names(data)))
  #     y <- data[, y_name]
  #   y_name <- "y"
  #   if(is.element("y", names(data))){
  #     y_name <- random_name(data=data)
  #     assign(y_name, y)
  #   }
  y_name <- check_y_name(y, y_name, data)
	mu_terms <- terms(mu_mod)
	sig2_terms <- terms(sig2_mod)
  y_fit <- lm.fit(x=model.matrix(mu_mod, data=data), y=get(y_name))
  r2 <- residuals(y_fit)^2
  r2_fit  <- lm.fit(x=model.matrix(sig2_mod, data=data), y=r2)
  mup <- fitted(y_fit)
  s2p <- fitted(r2_fit)
  y_std <- (y-mup)/sqrt(s2p)
  res <- list("y"=y,
              "y_std"=y_std,
              "data"=data,
              "mu_mod"=mu_mod, 
              "sig2_mod"=sig2_mod,
              "mu_terms"=mu_terms, 
              "sig2_terms"=sig2_terms,
              "mu_fit"=y_fit,
              "sig2_fit"=r2_fit)
  class(res) <- "trans"
  class(res) <- append(class(res), "std")
  res
}

transform_newdat <- function(y_trans, ...){
 UseMethod("transform_newdat")
}

transform_newdat.std <- function(y_trans, y, newdata, ...){
  mu_mod <- y_trans$mu_mod
  sig2_mod <- y_trans$sig2_mod
  mu_fit <- y_trans$mu_fit
  sig2_fit <- y_trans$sig2_fit
  x_mu <- get_mod_mat(y_trans$mu_terms, data=newdata)
  x_sig2 <- get_mod_mat(y_trans$sig2_terms, data=newdata)
  mu <- x_mu %*% coefficients(mu_fit)
  sig2 <- x_sig2 %*% coefficients(sig2_fit)
  y_std <- (y - mu) / sqrt(sig2)
  y_std
}

#' Apply a transformation on a new point.
#' 
#' Apply a transformation on a new point.
#' 
#' Apply to a point the same transformation as in the y_trans object. 
#' @param y_trans an object of class trans.
#' @param pnt a line from a data.frame giving the necessary information to perfom the transformation.
#' @param ... Arguments to be passed to methods,
#' @return returns a vector where the point value has been transformed.
#' @examples
#'data(tas)
#'t1 <- 2003
#'t0 <- 1990
#'xp <- 1.6
#'pnt1 <- set_pnt(t1, xp, time_var="year", tas)
#'pnt0 <- set_pnt(t0, xp, time_var="year", tas)
#'y_std <- standardize(eur_tas, data=tas, mu_mod=~gbl_tas, sig2_mod=~gbl_tas)
#'y_std_fit <- gpd_fit(y_std$y_std ,  data=tas, time_var="year", qthreshold=0.8)
#'pnt1_std <- transform_pnt(y_std, pnt1)
#'pnt0_std <- transform_pnt(y_std, pnt0)
#'get_far(y_std_fit, pnt0_std, pnt1_std, under_threshold=TRUE)
#'
#'eur_tas_positive <- with(tas, eur_tas + abs(min(eur_tas)) +1)
#'y_bc <- bc_fit(l_lambda=seq(-1, 1, 0.02), y=eur_tas_positive, data=tas, mu_mod=~gbl_tas, sig2_mod=~gbl_tas, time_var="year", to_plot=TRUE)
#'y_bc_fit <- gpd_fit(y_bc$y_std ,  data=tas, time_var="year", qthreshold=0.8)
#'pnt1_bc <- transform_pnt(y_bc, pnt1)
#'pnt0_bc <- transform_pnt(y_bc, pnt0)
#'get_far(y_bc_fit, pnt0_bc, pnt1_bc, under_threshold=TRUE)
transform_pnt <- function(y_trans, pnt, ...){
 UseMethod("transform_pnt")
}

#' @rdname transform_pnt
transform_pnt.std <- function(y_trans, pnt, ...){
  pnt$y <- as.numeric(transform_newdat.std(y_trans, y=pnt$y, newdata=pnt[, -1]))
  pnt
}

