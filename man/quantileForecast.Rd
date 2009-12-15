\name{quantileForecast}
\alias{quantileForecast}
\alias{quantileForecast.ensembleBMAgamma}
\alias{quantileForecast.ensembleBMAgamma0}
\alias{quantileForecast.ensembleBMAnormal}
\alias{quantileForecast.ensembleBMAbvnormal}
\alias{quantileForecast.fitBMAgamma}
\alias{quantileForecast.fitBMAgamma0}
\alias{quantileForecast.fitBMAnormal}
\alias{quantileForecast.fitBMAbvnormal}
\alias{quantBMAgamma}
\alias{quantBMAgamma0}
\alias{quantBMAnormal}
\title{
   Quantile forecasts at observation locations
}
\description{
   Computes quantiles for the probability distribution function
   (PDF) for ensemble forecasting models.
}
\usage{
quantileForecast( fit, ensembleData, quantiles = 0.5, dates=NULL, \dots)
}
\arguments{
  \item{fit}{
    A model fit to ensemble forecasting data.
  }
  \item{ensembleData}{
     An \code{ensembleData} object that includes ensemble
     forecasts, verification observations and possibly dates.
     Missing values (indicated by \code{NA}) are allowed. \\
     This need not be the data used for the model \code{fit}, 
     although it must include the same ensemble members. 
  }
  \item{quantiles}{
     The vector of desired quantiles for the PDF of the BMA mixture model.
  }
 \item{dates}{
     The dates for which the quantile forecasts will be computed.
     These dates must be consistent with \code{fit} and \code{ensembleData}.
     The default is to use all of the dates in \code{fit}.
     If \code{ensembleData} does not include dates, they will
     be inferred from \code{fit} and \code{dates}. 
  }
 \item{\dots}{
     Included for generic function compatibility.
 }
}
\value{
A vector of forecasts corresponding to the desired quantiles.
}
\details{
   This method is generic, and can be applied to any ensemble forecasting
   model. \cr
   Note the model may have been applied to a power transformation of the data,
   but that information is included in the input \code{fit}, and
   the output is transformed appropriately. \cr
   This can be used to compute prediction intervals for the PDF.  \cr
   For the bivariate normal model for wind speed and direction, the
   CRPS is computed for the marginal wind speed distribution.
}
\seealso{
  \code{\link{ensembleBMA}}, 
  \code{\link{fitBMA}}, 
  \code{\link{cdf}}
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
  data(ensBMAtest)

  ensMemNames <- c("gfs","cmcg","eta","gasp","jma","ngps","tcwb","ukmo")

  obs <- paste("T2","obs", sep = ".")
  ens <- paste("T2", ensMemNames, sep = ".")

  tempTestData <- ensembleData( forecasts = ensBMAtest[,ens],
                                dates = ensBMAtest[,"vdate"],
                                observations = ensBMAtest[,obs],
                                station = ensBMAtest[,"station"],
                                forecastHour = 48,
                                initializationTime = "00")

\dontrun{ # R check
  tempTestFit <- ensembleBMAnormal( tempTestData, trainingDays = 30)
}

  tempTestForc <- quantileForecast( tempTestFit, tempTestData)


\dontrun{ # R check

  data(srft)

  labels <- c("CMCG","ETA","GASP","GFS","JMA","NGPS","TCWB","UKMO")

  srftData <- ensembleData( forecasts = srft[ ,labels],
                            dates = srft$date,
                            observations = srft$obs,
                            latitude = srft$lat,
                            longitude = srft$lon,
                            forecastHour = 48,
                            initializationTime = "00")

  srftFit <- ensembleBMAnormal(srftData, date = "2004013100",
                               trainingDays = 25)

  data(srftGrid)

  srftGridData <- ensembleData(forecasts = srftGrid[ ,labels],
                               latitude = srftGrid$lat,
                               longitude = srftGrid$lon,
                               forecastHour = 48,
                               initializationTime = "00")

  srftGridForc <- quantileForecast( srftFit, srftGridData, 
                     date = "2004013100")
}}
\keyword{models}
% docclass is function
