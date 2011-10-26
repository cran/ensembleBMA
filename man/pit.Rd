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
pit( fit, ensembleData, dates = NULL, randomizeATzero=FALSE, \dots)
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
 \item{randomizeATzero}{
     For the \code{gamma0} model for precipitation, observations should be
     randomized at 0 for assessing the calibration. Has no effect for the
     other models.
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
  forecasts, in which case the observations in \code{ensembleData} would
  be those used in modeling \code{fit}. \cr
  Instances in \code{ensembleData} without verifying observations 
   are ignored. \cr
   Note the model may have been applied to a power transformation of the data,
   but that information is included in the input \code{fit}, and
   the output is transformed appropriately. \cr
   The PIT is a continuous analog of the verification rank.
}
\seealso{
  \code{\link{pitHist}},
  \code{\link{verifRankHist}},
  \code{\link{ensembleBMA}}, 
  \code{\link{fitBMA}}, 
  \code{\link{quantileForecast}}
}
\references{
A. E. Raftery, T. Gneiting, F. Balabdaoui and M. Polakowski,
Using Bayesian model averaging to calibrate forecast ensembles,
\emph{Monthly Weather Review 133:1155--1174, 2005}.

T. Gneiting, F. Balabdaoui and A. Raftery,
Probabilistic forecasts, calibration and sharpness.
\emph{Journal of the Royal Statistical Society, Series B 69:243--268, 2007}.

J. M. Sloughter, A. E. Raftery, T Gneiting and C. Fraley, 
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

\dontrun{ # R check
  tempTestFit <- ensembleBMAnormal( tempTestData, trainingDays = 30)
}

  tempTestForc <- quantileForecast( tempTestFit, tempTestData)
  range(tempTestForc)

  tempTestPIT <- pit( tempTestFit, tempTestData) 

}
\keyword{models}
% docclass is function
