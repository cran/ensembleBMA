\name{quantileForecastBMA}
\alias{quantileForecastBMA}
\alias{quantileForecastBMA.ensembleBMAgamma0}
\alias{quantileForecastBMA.ensembleBMAnormal}
\alias{quantileForecastBMA.fitBMAgamma0}
\alias{quantileForecastBMA.fitBMAnormal}
\alias{gamma0BMAquant}
\alias{gamma0BMAcdf}
\alias{normalBMAquant}
\alias{normalBMAcdf}
\title{
   Quantile forecasts for ensemble BMA models.
}
\description{
   Computes quantiles for the probability distribution function
   (PDF) of an ensemble BMA mixture model. 
}
\usage{
quantileForecastBMA( object, ensembleData, quantiles = 0.5, 
                     popData = NULL, \dots)
}
\arguments{
  \item{object}{
     An ensemble BMA model fit.
  }
  \item{ensembleData}{
     An \code{ensembleData} object giving including ensemble
     forecasts, observations and dates corresponding to \code{object}.
  }
  \item{quantiles}{
     The vector of desired quantiles for the PDF of the BMA mixture model.
  }
 \item{popData}{
     For \code{ensembleBMAgamma0} and \code{fitBMAgamma0} objects, there 
     is an additional \code{popData} argument for providing predictors
     in the logistic regression for probability of precipitation.
     The value of \code{popData} (if any) used to create \code{object}
     should be supplied here.
 }
 \item{\dots}{
     Included for generic function compatibility.
 }
}
\value{
A vector of forecasts corresponding to the desired quantiles.
}
\details{
   This method is generic, and can be applied to any ensemble BMA forecasting
   model. \cr
   Note the model may have been applied to a transformation of the data,
   but that information is included in the input \code{object}, and
   the output is transformed appropriately.\cr
   This can be used to compute confidence intervals for the PDF. 
}
\seealso{
  \code{\link{ensembleBMA}}, 
  \code{\link{forecastBMA}}, 
  \code{\link{gridForecastBMA}}, 
  \code{\link{crpsANDmae}}
}
\references{
A. E. Raftery, T. Gneiting, F. Balabdaoui and M. Polakowski,
Using Bayesian Model Averaging to Calibrate Forecast Ensembles,
\emph{Monthly Weather Review 133:1155-1174, 2005}.
}

\examples{

  data(slp)

  slpData <- ensembleData(forecasts = slp[c("AVN","GEM","ETA","NGM","NOGAPS")],
                          observations = slp$obs, dates = slp$date)

  slpFit <- ensembleBMAnormal(slpData)

  slpForc <- quantileForecastBMA( slpFit, slpData)
}
\keyword{models}
% docclass is function
