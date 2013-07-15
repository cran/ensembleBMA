\name{cdf}
\alias{cdf}
\alias{cdf.ensembleBMAgamma}
\alias{cdf.ensembleBMAgamma0}
\alias{cdf.ensembleBMAnormal}
\alias{cdf.fitBMAgamma}
\alias{cdf.fitBMAgamma0}
\alias{cdf.fitBMAnormal}
\alias{cdfBMAgamma}
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
     verification observations and possibly dates. 
     Missing values (indicated by \code{NA}) are allowed. \\
     This need not be the data used for the model \code{fit}, 
     although it must include the same ensemble members. 
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
   Note the model may have been applied to a power transformation of the data,
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


\dontrun{# R check
  tempTestFit <- ensembleBMAnormal( tempTestData, trainingDays = 30)
}

# for quick run only; use more training days for forecasting
  tempTestFit <- ensembleBMAnormal( tempTestData[1:20,], trainingDays = 8)

  tempTestForc <- quantileForecast( tempTestFit, tempTestData)
  range(tempTestForc)

  tempTestCDF <- cdf( tempTestFit, tempTestData, 
                      values = seq(from=277, to=282, by = 1))

  tempTestCDF
}
\keyword{models}
% docclass is function
