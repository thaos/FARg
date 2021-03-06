% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/refit.R
\name{refit}
\alias{refit}
\alias{refit.gauss_fit}
\alias{refit.gev_fit}
\alias{refit.gpd_fit}
\title{refit generic.}
\usage{
refit(fit_o, fit_m)

\method{refit}{gauss_fit}(fit_o, fit_m)

\method{refit}{gev_fit}(fit_o, fit_m)

\method{refit}{gpd_fit}(fit_o, fit_m)
}
\arguments{
\item{fit_o}{an object of class gauss_fit, gev_fit, gpd_fit. the data from
the fit_o model are used to refit the intercept of the fit_m model.}

\item{fit_o}{an object of class gauss_fit, gev_fit, gpd_fit. the data from
the fit_o model are used to refit the intercept of the fit_m model.}
}
\value{
the fit_o object with the newly fitted parameters for the intercepts
and the other parameters being those of fit_m
}
\description{
.
}
\details{
\code{refit} recompute the statistical fit for the intercepts of the model fit
_m on the data used for the model fit_o.
fit_o and fit_m should follows that same model ans should be of the same R
classe
}
\section{Methods (by class)}{
\itemize{
\item \code{gauss_fit}: for a gaussian fit.

\item \code{gev_fit}: for a gev fit.

\item \code{gpd_fit}: for a gpd fit.
}}

\examples{
#data(tas)
#Example with the same covariate for the mean and variance parameter
ga_fit <- gauss_fit(eur_tas, data=tas, mu_mod=~gbl_tas, sig_mod=~gbl_tas, time_var="year")
ga_refit <- refit(ga_fit, ga_fit)
ge_fit <- gev_fit(eur_tas, data=tas, mu_mod=~gbl_tas, sig_mod=~gbl_tas, time_var="year")
ge_refit <- refit(ge_fit, ge_fit)
gp_fit <- gpd_fit(eur_tas, data=tas, mu=~gbl_tas, sig_mod=~gbl_tas, time_var="year", qthreshold=0.9)
gp_refit <- refit(gp_fit, gp_fit)
}
