% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geo_fit.r
\name{geo.fit}
\alias{geo.fit}
\title{Function to fit a  either a spatial or spatio-temporal model to geo-statistical data with an intercept and covariates}
\usage{
geo.fit(mesh = NULL, locs = NULL, response = NULL, temp = NULL,
  covariates = NULL, family = "gaussian", control.time = list(model =
  "ar1", param = list(theta = list(prior = "pccor1", param = c(0, 0.9)))),
  control.inla = list(strategy = "gaussian", int.strategy = "eb"),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
  non.linear = NULL, prediction = NULL, sig0 = 1, Psig = 0.5,
  rho0 = 0.3, Prho = 0.5, verbose = FALSE, ...)
}
\arguments{
\item{mesh}{a ``mesh'' object i.e. delauney triangulation of the domain, an
object returned by \link{make.mesh}.}

\item{locs}{a matrix of observation locations, where each row corresponds to the observation.}

\item{response}{a vector of response variable, each corresponds to the spatial locations
in \code{locs}.}

\item{temp}{(optional) a numeric vector specifying a temporal index for each observation (starting at 1.....T).}

\item{covariates}{(optional) a named data.frame of covariates.}

\item{family}{a character vector specifying the assumed likelihood of the response, by default is "gaussian".}

\item{control.time}{(optional) supplied if the \code{temp} argumet is given to fit a spatio-temporal model. This argument
controls the model and prior put on the hyperparameters of the model for the temporal component of the spatio-temporal
model. By default this is \code{list(model = 'ar1', param = list(theta = list(prior='pccor1', param = c(0, 0.9))))}
which is a pc.prior put on the rho coefficient of a AR(1) model with P(rho>0)=0.9.
Refer to Simpson, martins, and rue for further details *****put in proper refs*****}

\item{control.inla}{a list which controls the fitting procedures INLA uses see Rue et al. ***ref book***
by default this is \code{list(strategy='gaussian',int.strategy = 'eb')} for quick and dirty fitting.}

\item{control.compute}{a list of fit statistics the user wants INLA to return. By default this
is \code{list(dic = TRUE, waic = TRUE,cpo = TRUE, config = TRUE)}.}

\item{non.linear}{(optional) should be used if the user requires a non-linear covariate to be included in the model
Must be supplied as a named list with elements \code{random.effect} a numeric vector of the random effect indecies,
and \code{model} the random effect model the user wishes to use for \code{random.effect}}

\item{prediction}{(optional) should be used if the uses wnts to run a prediction model. Must be supplied as a
named list with \code{pred.locs} the locations where predictionas are to be made, and only if a spatio-tempral
model is to be fitted \code{pred.temp} the temporal indecies for the predictions (the same length as \code{pred.locs}).}

\item{sig0}{by default = 1, typical standard deviation to use pc priors for hyperparams of spde model}

\item{Psig}{by default = 0.5 prob for sigma of pc prior}

\item{rho0}{by default = 0.3, typical range to use pc priors for hyperparams of spde model}

\item{Prho}{by default = 0.5 prob for rho of pc prior}

\item{verbose}{Logical if \code{TRUE} model fit is output to screen.}

\item{...}{add inla options to speed up computation i.e., by giving starting values from a previos model}
}
\value{
A \code{inla} result object
}
\description{
Function to fit a  either a spatial or spatio-temporal model to geo-statistical data with an intercept and covariates
}

