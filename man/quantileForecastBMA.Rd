\name{quantileForecastBMA}
\alias{quantileForecastBMA}
\alias{quantileForecastBMA.ensembleBMAgamma0}
\alias{quantileForecastBMA.ensembleBMAnormal}
\alias{quantileForecastBMA.fitBMAgamma0}
\alias{quantileForecastBMA.fitBMAnormal}
\alias{gamma0BMAquant}
\alias{normalBMAquant}
\title{
   Quantile forecasts at observation locations
}
\description{
   Computes quantiles for the probability distribution function
   (PDF) for an ensemble BMA mixture model at observation locations.
}
\usage{
quantileForecastBMA( fit, ensembleData, quantiles = 0.5, 
                     dates = NULL, popData = NULL, \dots)
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
  \item{quantiles}{
     The vector of desired quantiles for the PDF of the BMA mixture model.
  }
 \item{dates}{
     The dates for which the quantile forecasts will be computed.
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
A vector of forecasts corresponding to the desired quantiles.
}
\details{
   This method is generic, and can be applied to any ensemble BMA forecasting
   model. \cr
   Note the model may have been applied to a transformation of the data,
   but that information is included in the input \code{fit}, and
   the output is transformed appropriately. \cr
   This can be used to compute prediction intervals for the PDF. 
}
\seealso{
  \code{\link{ensembleBMA}}, 
  \code{\link{fitBMA}}, 
  \code{\link{cdfBMA}}
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
using Ensembles and Bayesian Model Averaging,
Technical Report No. 516, Department of Statistics, University of
Washington, August 2007.
}

\examples{

  data(slpTest)

  labels <- c("AVN","GEM","ETA","NGM","NOGAPS")
  slpTestData <- ensembleData( forecasts = slpTest[ ,labels],
                         observations = slpTest$obs, dates = slpTest$date)

  slpTestFit <- ensembleBMAnormal(slpTestData)

  slpTestForc <- quantileForecastBMA( slpTestFit, slpTestData)

\dontrun{
  data(srft)

  labels <- c("CMCG","ETA","GASP","GFS","JMA","NGPS","TCWB","UKMO")
  srftData <- ensembleData( forecasts = srft[ ,labels],
                            dates = srft$date, observations = srft$obs,
                            latitude = srft$lat, longitude = srft$lon)

  srftFit <- ensembleBMAnormal(srftData, date = "2004012900")

  data(srftGrid)

  srftGridData <- ensembleData(forecasts = srftGrid[ ,labels],
                           latitude = srftGrid$lat, longitude = srftGrid$lon)

  srftGridForc <- quantileForecastBMA( srftFit, srftGridData, 
                     date = "2004012900")
}
}
\keyword{models}
% docclass is function
