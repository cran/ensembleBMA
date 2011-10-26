\name{ensembleData}
\alias{ensembleData}
\alias{dataNA}
\alias{dataNobs}
\alias{dataObsLabels}
\alias{dataStartupSpeed}
\alias{dataVerifObs}
\alias{ensembleForecasts}
\alias{ensembleFhour}
\alias{ensembleGroups}
\alias{ensembleItime}
\alias{ensembleMembers}
\alias{ensembleSize}
\alias{ensembleValidDates}
\alias{dataObsLabels.ensembleData}
\alias{dataNA.ensembleData}
\alias{dataNobs.ensembleData}
\alias{dataStartupSpeed.ensembleData}
\alias{dataVerifObs.ensembleData}
\alias{ensembleMembers.ensembleData}
\alias{ensembleFhour.ensembleData}
\alias{ensembleForecasts.ensembleData}
\alias{ensembleGroups.ensembleData}
\alias{ensembleItime.ensembleData}
\alias{ensembleSize.ensembleData}
\alias{ensembleValidDates.ensembleData}
\alias{[.ensembleData}
\alias{matchEnsembleMembers.default}
\alias{matchITandFH}
\alias{matchITandFH.default}
\alias{startupSpeed}
\alias{startupSpeed.ensembleData}
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
              forecastHour, initializationTime,
              startupSpeed = NULL, exchangeable = NULL)
}
\arguments{
  \item{forecasts}{
     A matrix or array (for vector qunatities) 
     with columns corresponding to forecasts from individual members 
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
    Optional vector (or matrix for vector quantities) of observed weather 
    conditions corresponding 
    to the forecasts. Must be supplied if the data is to be
    used for BMA modeling. 
  }
  \item{...}{
    A named list of additional attributes such as latitude,
    longitude, and startupSpeed for wind speed.
  }
  \item{forecastHour}{
     A number giving the \emph{forecast hour}, the time interval 
     between the initialization and forecast times, in units of hours.
  }
  \item{initializationTime}{
     A number or character string giving the initialization time.
  }
  \item{startupSpeed}{
    A numeric value specifying a value  below which the anemometer
    readings for wind speed will be recorded as zero. This value
    is used for all stations when the startup speed is not explicity
    specified as part of the data.
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
    For use with  batch processing modeling functions (\code{ensembleBMA}
    etc), instances \code{ensembleData} object are assumed
    the same forecast hour and initialization time, which should be
    specified as part of the object. \cr  
    Methods for \code{ensembleData} objects include \code{ensembleSize},
    \code{ensembleForecasts}, \code{ensembleValidDates}. \cr
    Subsetting is possible, but in the case of columns it applies only to
    the ensemble forecasts. \cr
    For vector wind computations, the velocity should be in the first 
    column and the direction in the second.
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

