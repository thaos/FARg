standardize <- function(y, data, mu_mod=~1, sig2_mod=~1,to_plot=FALSE){
  y_fit <- lm.fit(x=model.matrix(mu_mod, data=data), y=y)
  r2 <- residuals(y_fit)^2
  r2_fit  <- lm.fit(x=model.matrix(sig2_mod, data=data), y=fit_res2)
  mup <- fitted(fit_y)
  s2p <- fitted(fit_r2)
  y_std <- (ydat$y-mup)/sqrt(s2p)
  attr(y_std,"mu_mod"=mu_mod)
  attr(y_std, "sig2_mod"=sig2_mod)
  attr(y_std,"mu_fit"=y_fit)
  attr(y_std, "sig2_fit"=r2_fit)
  class(y_std, "std")
  y_std
}

transform_newdat.std <- function(y_trans, y, newdata){
  mu_mod <- attr(y_trans, "mu_mod")
  sig2_mod <- attr(y_trans, "sig2_mod")
  mu_fit <- attr(y_trans, "mu_fit")
  sig2_fit <- attr(y_trans, "sig2_fit")
  x_mu <- model.matrix(mu_mod, data=newdata)
  x_sig2 <- model.matrix(sig2_mod, data=newdata)
  mu <- x_mu %*% coefficients(mu_fit)
  sig2 <- x_sig2 %*% coefficients(sig2_fit)
  y_std <- (y - mu) / sqrt(sig2)
  y_std
}

pnt2bc <- function(bc_fit, pnt){
  pnt$y <- data2bc(bc_fit, y=pnt$y, pnt)$y_std
  pnt
}


zp2yp <- function(lambda, mu, sigma, zp) 
	(lambda * ( mu + sigma * zp) + 1)^(1 / lambda)

