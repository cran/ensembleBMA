\name{pitHist}
\alias{pitHist}
\title{
   PIT Histogram
}
\description{
  Computes the probability integral transform of the obervations relative to
  the BMA forecast, and  plots its histogram.
}
\usage{
pitHist( fit, ensembleData, dates=NULL)
}
\arguments{
   \item{fit}{
     A model fit to ensemble forecasting data.
  }
  \item{ensembleData}{
     An \code{ensembleData} object that includes ensemble forecasts,
     verification observations and possibly dates.
     Missing values (indicated by \code{NA}) are allowed. \\
     This need not be the data used for the model \code{fit},
     although it must include the same ensemble members.
  }
  \item{dates}{
     The dates for which the CDF will be computed.
     These dates must be consistent with \code{fit} and
     \code{ensembleData}.
     The default is to use all of the dates in \code{fit}.
     The dates are ignored if \code{fit} originates from \code{fitBMA},
     which also ignores date information.
  }
}
\value{
  The value of the BMA cumulative distribution function CDF
  corresponding to the fit at the observed values in \code{ensembleData}.
  The corresponding histogram is also plotted.
  }
\details{
  PIT histograms are used to assess calibration of
  forecasts, in which case the observations in \code{ensembleData} would
  be those used in modeling \code{fit}. \cr
  Instances in \code{ensembleData} without verifying observations
  are ignored. \cr
  In the case of the \code{gamma0} model for precipitation,
  observations of zero precipitation are randomized within their
  probabilistics range to avoid a false 
  impression of bias.\cr
  Note the model may have been applied to a power transformation of the data,
  but that information is included in the input \code{fit}, and
  the output is transformed appropriately. \cr
  The PIT is a continuous analog of the verification rank.
}
\seealso{
  \code{\link{ensembleData}}, 
  \code{\link{pit}},   
  \code{\link{verifRankHist}}. 
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


\dontrun{
# this takes time
# the PIT should be evaluated over relatively long periods 

  srftFITall <- ensembleBMA( srftData, model = "normal", trainingDays = 25)

  srftPIT <- pitHist( srftFITall, srftData)
}
}
\keyword{models}
% docclass is function
