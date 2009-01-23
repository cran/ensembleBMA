\name{ensembleData}
\alias{ensembleData}
\alias{ensembleForecasts}
\alias{ensembleFhour}
\alias{ensembleGroups}
\alias{ensembleItime}
\alias{ensembleMemberLabels}
\alias{ensembleNobs}
\alias{ensembleObsLabels}
\alias{ensembleSize}
\alias{ensembleValidDates}
\alias{ensembleVerifObs}
\alias{ensembleMemberLabels.ensembleData}
\alias{ensembleFhour.ensembleData}
\alias{ensembleForecasts.ensembleData}
\alias{ensembleGroups.ensembleData}
\alias{ensembleItime.ensembleData}
\alias{ensembleObsLabels.ensembleData}
\alias{ensembleSize.ensembleData}
\alias{ensembleNobs.ensembleData}
\alias{ensembleValidDates.ensembleData}
\alias{ensembleVerifObs.ensembleData}
\alias{[.ensembleData}
\alias{matchEnsembleMembers}
\alias{matchEnsembleMembers.ensembleBMAgamma}
\alias{matchEnsembleMembers.ensembleBMAgamma0}
\alias{matchEnsembleMembers.ensembleBMAnormal}
\alias{matchEnsembleMembers.fitBMAgamma}
\alias{matchEnsembleMembers.fitBMAgamma0}
\alias{matchEnsembleMembers.fitBMAnormal}
\alias{matchITandFH}
\alias{matchITandFH.default}
\title{
  Create an ensembleData object
}
\description{
   Creates an \code{ensembleData} object including ensemble forecasts along
   with dates and (optionally) observations. Other descriptive information
   such as latitude, longitude, and station type may be included as well.
}
\usage{
ensembleData( forecasts, dates = NULL, observations = NULL, ...,
              forecastHour, initializationTime, exchangeable = NULL)
}
\arguments{
  \item{forecasts}{
     A matrix with columns corresponding to forecasts from individual members 
     of an ensemble and rows corresponding to forecasts for the same date.
  }
  \item{dates}{
    A numeric or character vector or factor specifying the valid dates
    for the forecasts. If numeric, it is interpreted as a Julian date
    if it has an \code{origin} attribute specifying
    the month, day, and year, e.g. \code{c(month = 1, day = 1, year = 2000)}.
    Otherwise the character form of each date must be a string with format
    format YYYYMMDDHH or YYYYMMDD, where YYYY is the year,
    MM the monthn, DD the day, and (optionally) HH the hour.
   }
  \item{observations}{
    Optional vector of observed weather conditions corresponding 
    to the forecasts. Must be supplied if the data is to be
    used for BMA modeling. 
  }
  \item{...}{
    A named list of additional attributes such as latitude
    and longitude.
  }
  \item{forecastHour}{
     A numeric vector giving the forecast hour.
  }
  \item{initializationTime}{
     A numeric or character vector giving the initialization time.
  }
  \item{exchangeable}{
     A numeric or character vector or factor indicating groups of 
     ensemble members that are exchangeable (indistinguishable).
     The models fit will have equal weights and parameters
     within each group.
     The same names/labels should be used as for the forecasts.     
     The default assumes that none of the ensemble members are
     exhangeable.
  }
}
\value{
  An \code{ensembleData} object, incorporating forecasts and
  (optionally) observations with the associated valid dates.
}
\details{
    All instances in an \code{ensembleData} object are assumed to have
    the same forecast hour and initialization time, which should be
    specified as part of the object. \cr  
    Methods for \code{ensembleData} objects include \code{ensembleSize},
    \code{ensembleForecasts}, \code{ensembleValidDates}, and
    \code{ensembleVerifObs}. \cr
    Subsetting is possible, but in the case of columns it applies only to
    the ensemble forecasts.
}

\seealso{
  \code{\link{ensembleBMA}}, 
  \code{\link{ensembleBMAgamma}}, 
  \code{\link{ensembleBMAgamma0}}, 
  \code{\link{ensembleBMAnormal}}
}

\references{
C. Fraley, A. E. Raftery, T. Gneiting and J. M. Sloughter,
\code{ensembleBMA}: An \code{R} Package for Probabilistic Forecasting
using Ensembles and Bayesian Model Averaging,
Technical Report No. 516R, Department of Statistics, University of
Washington, 2009.
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

  obs <- paste("PCP24","obs", sep = ".")
  ens <- paste("PCP24", ensMemNames, sep = ".")

  prcpTestData <- ensembleData( forecasts = ensBMAtest[,ens],
                                dates = ensBMAtest[,"vdate"], 
                                observations = ensBMAtest[,obs], 
                                station = ensBMAtest[,"station"],
                                forecastHour = 48,
                                initializationTime = "00")
\dontrun{ # R check

  prcpTestFit <- ensembleBMAgamma0( prcpTestData, trainingDays = 30)

}

  obs <- paste("MAXWSP10","obs", sep = ".")
  ens <- paste("MAXWSP10", ensMemNames, sep = ".")

  winsTestData <- ensembleData( forecasts = ensBMAtest[,ens],
                                dates = ensBMAtest[,"vdate"], 
                                observations = ensBMAtest[,obs],
                                station = ensBMAtest[,"station"],
                                forecastHour = 48,
                                initializationTime = "00")

\dontrun{ # R check

 winsTestFit <- ensembleBMAgamma(winsTestData, trainingDays = 30)

}

}
\keyword{data}
% docclass is function

