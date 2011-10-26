\name{verifRankHist}
\alias{verifRankHist}
\title{
   Verification Rank and Histogram
}
\description{
  Computes the rank of verifying observations
  relative to the corresponding ensemble forecasts
   and plots its histogram.
}
\usage{
  verifRankHist( forecasts, observations)
}
\arguments{
   \item{forecasts}{
     A matrix of ensemble forecasts, in which the rows corresponds to
     locations and times and the columns correspond to the individual
     ensemble members.
  }
  \item{observations}{
    A vector of observations corresponding to the locations and times
    of the forecasts.
  }
}
\value{
  A vector giving the rank of verifying observations relative
  to the corresponding ensemble forecasts.
  The verification rank historgram is plotted. 
 }
\details{
  The verification rank is used to assess calibration of a forecast
  ensemble. A more uniform verification rank histogram
  indicates better calibartion. 
}
\seealso{
  \code{\link{ensembleData}}, 
  \code{\link{pit}} 
}
\references{
A. E. Raftery, T. Gneiting, F. Balabdaoui and M. Polakowski,
Using Bayesian model averaging to calibrate forecast ensembles,
\emph{Monthly Weather Review 133:1155-1174, 2005}.

T. Gneiting, F. Balabdaoui and A. Raftery,
Probabilistic forecasts, calibration and sharpness.
\emph{Journal of the Royal Statistical Society, Series B 69:243--268, 2007}.

J. M. Sloughter, A. E. Raftery, T. Gneiting and C. Fraley,
Probabilistic quantitative precipitation forecasting
using Bayesian model averaging,
\emph{Monthly Weather Review 135:3209--3220, 2007}.

C. Fraley, A. E. Raftery, T. Gneiting and J. M. Sloughter,
\code{ensembleBMA}: An \code{R} Package for Probabilistic Forecasting
using Ensemble and Bayesian Model Averaging,
Technical Report No. 516R, Department of Statistics, University of
Washington, 2007 (revised 2010).

C. Fraley, A. E. Raftery, T. Gneiting,
Calibrating Multi-Model Forecast Ensembles
with Exchangeable and Missing Members using Bayesian Model Averaging,
\emph{Monthly Weather Review 138:190--202, 2010}.

J. M. Sloughter, T. Gneiting and A. E. Raftery, 
Probabilistic wind speed forecasting 
using ensembles and Bayesian model averaging,
\emph{Journal of the American Statistical Association, 105:25--35, 2010}.
}

\examples{
  data(srft)

  labels <- c("CMCG","ETA","GASP","GFS","JMA","NGPS","TCWB","UKMO")

  srftData <- ensembleData( forecasts = srft[ ,labels],
                            dates = srft$date,
                            observations = srft$obs,
                            latitude = srft$lat,
                            longitude = srft$lon,
                            forecastHour = 48,
                            initializationTime = "00")

  use <- ensembleValidDates(srftData) >= "2004013000"

  verifRankHist( ensembleForecasts(srftData[use,]),
                 dataVerifObs(srftData[use,]))
}
\keyword{models}
% docclass is function
