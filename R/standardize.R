#' @export
standardize <- function(y, data, mu_mod=~1, sig2_mod=~1,to_plot=FALSE){
  y_fit <- lm.fit(x=model.matrix(mu_mod, data=data), y=y)
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
              "mu_fit"=y_fit,
              "sig2_fit"=r2_fit)
  class(res) <- "trans"
  class(res) <- append(class(res), "std")
  res
}

#' @export
transform_newdat <- function(y_trans, ...){
 UseMethod("transform_newdat")
}

#' @export
transform_newdat.std <- function(y_trans, y, newdata, ...){
  mu_mod <- y_trans$mu_mod
  sig2_mod <- y_trans$sig2_mod
  mu_fit <- y_trans$mu_fit
  sig2_fit <- y_trans$sig2_fit
  x_mu <- model.matrix(mu_mod, data=newdata)
  x_sig2 <- model.matrix(sig2_mod, data=newdata)
  mu <- x_mu %*% coefficients(mu_fit)
  sig2 <- x_sig2 %*% coefficients(sig2_fit)
  y_std <- (y - mu) / sqrt(sig2)
  y_std
}

#' @export
transform_pnt <- function(y_trans, ...){
 UseMethod("transform_pnt")
}

#' @export
transform_pnt.std <- function(y_trans, pnt, ...){
  pnt$y <- as.numeric(transform_newdat.std(y_trans, y=pnt$y, newdata=pnt[, -1]))
  pnt
}

