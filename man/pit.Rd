\name{pit}
\alias{pit}
\alias{pit.ensembleBMAgamma}
\alias{pit.ensembleBMAgamma0}
\alias{pit.ensembleBMAnormal}
\alias{pit.fitBMAgamma}
\alias{pit.fitBMAgamma0}
\alias{pit.fitBMAnormal}
\title{
   Probability Integral Transform for ensemble forcasting models
}
\description{
   Computes the probabilty integral transform (PIT)
   of a BMA ensemble forecasting model at observation locations.
}
\usage{
pit( fit, ensembleData, dates = NULL, \dots)
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
  The value of the BMA cumulative distribution function CDF
  corresponding to the fit at the observed values in \code{ensembleData}. 
  }
\details{
  Most often used for computing PIT histograms to assess calibration of
  forecasts, in which case the observations in \code{ensembleData} awould
  be those used in modeling \code{fit}. \cr
  Instances in \code{ensembleData} without verifying observations 
   are ignored. \cr
   Note the model may have been applied to a power transformation of the data,
   but that information is included in the input \code{fit}, and
   the output is transformed appropriately. \cr
   The PIT is a continuous analog of the verification rank.
}
\seealso{
  \code{\link{ensembleBMA}}, 
  \code{\link{fitBMA}}, 
  \code{\link{quantileForecast}},
  \code{\link{verifRank}},
}
\references{
A. E. Raftery, T. Gneiting, F. Balabdaoui and M. Polakowski,
Using Bayesian model averaging to calibrate forecast ensembles,
\emph{Monthly Weather Review 133:1155--1174, 2005}.

T. Gneiting, F. Balabdaoui and A. Raftery,
Probabilistic forecasts, calibration and sharpness.
\emph{Journal of the Royal Statistical Society, Series B 69:243--268, 2007}.

C. Fraley, A. E. Raftery, T. Gneiting and J. M. Sloughter,
\code{ensembleBMA}: An \code{R} Package for Probabilistic Forecasting
using Ensemble and Bayesian Model Averaging,
Technical Report No. 516R, Department of Statistics, University of
Washington, revised 2009.
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
  range(tempTestForc)

  tempTestPIT <- pit( tempTestFit, tempTestData, 
                        values = seq(from=277, to=282.5, by = .1))

  hist(tempTestPIT, breaks = 7)
}
\keyword{models}
% docclass is function
