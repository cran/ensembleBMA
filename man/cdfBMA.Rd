\name{cdfBMA}
\alias{cdfBMA}
\alias{cdfBMA.ensembleBMAgamma0}
\alias{cdfBMA.ensembleBMAnormal}
\alias{cdfBMA.fitBMAgamma0}
\alias{cdfBMA.fitBMAnormal}
\alias{gamma0BMAcdf}
\alias{normalBMAcdf}
\title{
   Cummulative Distribution Function for ensemble BMA models
}
\description{
   Computes the cumulative distribution function (CDF)
   of an ensemble BMA mixture model at observation locations.
}
\usage{
cdfBMA( fit, ensembleData, values, dates = NULL, popData = NULL, \dots)
}
\arguments{
  \item{fit}{
     An ensemble BMA model fit.
  }
  \item{ensembleData}{
     An \code{ensembleData} object including ensemble
     forecasts and observations. It need not be the object
     used to form \code{fit}, although it must include the same
     ensemble members. If \code{ensembleData} includes dates,
     they must be consistent with \code{fit} and \code{dates}. 
     If \code{ensembleData} does not include dates, they will
     be inferred from the \code{fit} and the \code{dates} argument.
  }
  \item{values}{
     The vector of desired values at which the CDF of the BMA mixture model
     is to be evaluated.
  }
 \item{dates}{
     The dates for which the CDF will be computed.
     These dates must be consistent with \code{fit} and \code{ensembleData}.
     The default is to use all of the dates in \code{fit}.
  }
 \item{popData}{ 
    For \emph{gamma0} model fits, there is an additional \code{popData}
    argument for providing predictors in the logistic regression for 
    probability of zero precipitation. If \code{popData} was supplied to
    obtain in the modeling for \code{fit}, then \code{popData} of the
    same kind must be supplied here.
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
   This method is generic, and can be applied to any ensemble BMA forecasting
   model. \cr
   Note the model may have been applied to a transformation of the data,
   but that information is included in the input \code{fit}, and
   the output is transformed appropriately.
}
\seealso{
  \code{\link{ensembleBMA}}, 
  \code{\link{fitBMA}}, 
  \code{\link{quantileForecastBMA}}
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
Technical Report No. 516, Department of Statistics, University of
Washington, August 2007.
}

\examples{

  data(slpTest)

  memberLabels <- c("AVN","GEM","ETA","NGM","NOGAPS")

  slpTestData <- ensembleData(forecasts = slpTest[ ,memberLabels],
                         observations = slpTest$obs, dates = slpTest$date)

  slpTestFit <- ensembleBMAnormal(slpTestData)

  slpTestForc <- quantileForecastBMA( slpTestFit, slpTestData)
  range(slpTestForc)

  slpTestCDF <- cdfBMA( slpTestFit, slpTestData, 
                        values = seq(from=1005, to=1025, by = 5))
}
\keyword{models}
% docclass is function
