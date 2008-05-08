\name{cdf}
\alias{cdf}
\alias{cdf.ensembleBMAgamma0}
\alias{cdf.ensembleBMAnormal}
\alias{cdf.fitBMAgamma0}
\alias{cdf.fitBMAnormal}
\alias{cdfBMAgamma0}
\alias{cdfBMAnormal}
\title{
   Cummulative Distribution Function for ensemble forcasting models
}
\description{
   Computes the cumulative distribution function (CDF)
   of an ensemble forecasting model at observation locations.
}
\usage{
cdf( fit, ensembleData, values, dates = NULL, \dots)
}
\arguments{
   \item{fit}{
     A model fit to ensemble forecasting data.
  }
  \item{ensembleData}{
     An \code{ensembleData} object that includes ensemble forecasts,
     verification observations and dates. 
     Missing values (indicated by \code{NA}) are allowed. \\
     This need not be the data used for the model \code{fit}, 
     although it must include the same ensemble members. \\
     If \code{ensembleData} includes dates,
     they must be consistent with \code{fit} and \code{dates}. 
     If \code{ensembleData} does not include dates, they will
     be inferred from \code{fit} and \code{dates}. 
  }
  \item{values}{
    The vector of desired values at which the CDF of the ensemble
    forecasting model is to be evaluated.
  }
 \item{dates}{
     The dates for which the CDF will be computed.
     These dates must be consistent with \code{fit} and \code{ensembleData}.
     The default is to use all of the dates in \code{fit}.
     The dates are ignored if \code{fit} originates from \code{fitBMA}, 
     which also ignores date information.
  }
 \item{\dots}{
     Included for generic function compatibility.
 }
}
\value{
  A vector of probabilities corresponding to the CDF at the desired values.
  Useful for determining propability of freezing, precipitation, etc.
}
\details{
   This method is generic, and can be applied to any ensemble forecasting
   model. \cr
   Note the model may have been applied to a transformation of the data,
   but that information is included in the input \code{fit}, and
   the output is transformed appropriately.
}
\seealso{
  \code{\link{ensembleBMA}}, 
  \code{\link{fitBMA}}, 
  \code{\link{quantileForecast}}
}
\references{
A. E. Raftery, T. Gneiting, F. Balabdaoui and M. Polakowski,
Using Bayesian model averaging to calibrate forecast ensembles,
\emph{Monthly Weather Review 133:1155--1174, 2005}.

J. M. Sloughter, A. E. Raftery, T. Gneiting and C. Fraley,
Probabilistic quantitative precipitation forecasting
using Bayesian model averaging,
\emph{Monthly Weather Review 135:3209--3220, 2007}.

C. Fraley, A. E. Raftery, T. Gneiting and J. M. Sloughter,
\code{ensembleBMA}: An \code{R} Package for Probabilistic Forecasting
using Ensemble and Bayesian Model Averaging,
Technical Report No. 516R, Department of Statistics, University of
Washington, May 2008.

C. Fraley, A. E. Raftery, T. Gneiting,
BMA Forecasting with Missing and Exchangeable Ensemble Members,
in preparation.
}

\examples{

  data(slpTest)

  memberLabels <- c("AVN","GEM","ETA","NGM","NOGAPS")

  slpTestData <- ensembleData(forecasts = slpTest[ ,memberLabels],
                         observations = slpTest$obs, dates = slpTest$date)

  slpTestFit <- ensembleBMAnormal(slpTestData,
                   trainingRule = list(length=30, lag = 2))

  slpTestForc <- quantileForecast( slpTestFit, slpTestData)
  range(slpTestForc)

  slpTestCDF <- cdf( slpTestFit, slpTestData, 
                        values = seq(from=1005, to=1025, by = 5))
}
\keyword{models}
% docclass is function
