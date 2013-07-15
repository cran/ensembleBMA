\name{crps}
\alias{crps}
\alias{crps.ensembleBMAgamma}
\alias{crps.ensembleBMAgamma0}
\alias{crps.ensembleBMAnormal}
\alias{crps.ensembleBMAbvnormal}
\alias{crps.fitBMAgamma}
\alias{crps.fitBMAgamma0}
\alias{crps.fitBMAnormal}
\alias{crps.fitBMAbvnormal}
\alias{CRPS}
\alias{CRPS.default}
\title{
   Continuous Ranked Probability Score
}
\description{
   Computes the continuous ranked probability score (CRPS)
   for univariate ensemble forecasting models. 
}
\usage{
crps( fit, ensembleData, dates=NULL, nSamples=NULL, seed=NULL, \dots)
CRPS( fit, ensembleData, dates=NULL, nSamples=NULL, seed=NULL, \dots)
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
  \item{nSamples}{
     The number of simulation samples for CRPS via simulation.
     For the normal model, the default is analytic computation of the CRPS.
     For the gamma model with a point mass at 0 (precipitation), 
     the CRPS is always computed by simulation,
     with default \code{nSamples = 10000}.
  }
  \item{seed}{
     Argument to \code{set.seed} for random number generation in simulation.
  }
  \item{dates}{
     The dates for which the CRPS will be computed.
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
  \code{crps} is a matrix giving the CRPS for each instance in the data
  for both the raw ensemble and the median probabilistic forecast. \cr
  \code{CRPS} is a vector giving the mean of the CRPS over all
  instances for the raw ensemble and the median probabilistic forecast.
}
\details{
   These methods are generic, and can be applied to all ensemble forecasting
   models. \cr
   For \code{gamma0} model for precipitation and the \code{gamma} model
   for wind speed the CRPS is only available through simulation.
   The default number of simulation samples is 10,000. \cr
   Note that the \code{gamma0} model for precipitation and the
   \code{gamma} model for wind speed may have been applied to a power
   transformation of the data. \cr
   For normal models for temperature and pressure, analytic computation 
   of the CRPS is the default. CRPS will be computed via simulation for
   normal models only if \code{nSamples} is set to a positive value. \cr
   For the bivariate normal model for wind speed and direction, the
   CRPS is computed for the marginal wind speed distribution.
}
\seealso{
  \code{\link{ensembleBMA}},
  \code{\link{fitBMA}}
}
\references{
T. Gneiting and A. E. Raftery,
Strictly proper scoring rules, prediction and estimation,
\emph{Journal of the American Statistical Association 102:359--378, 2007}.

C. Fraley, A. E. Raftery, T. Gneiting and J. M. Sloughter,
\code{ensembleBMA}: An \code{R} Package for Probabilistic Forecasting
using Ensembles and Bayesian Model Averaging,
Technical Report No. 516R, Department of Statistics, University of
Washington, 2007 (revised 2010).
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

  crpsValues <- crps( tempTestFit, tempTestData)
  colMeans(crpsValues)

  CRPS( tempTestFit, tempTestData)
}
\keyword{models}
% docclass is function
