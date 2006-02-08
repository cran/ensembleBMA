\name{ensembleData}
\alias{ensembleData}
\alias{precipitationData}
\alias{temperatureData}
\alias{pressureData}
\alias{ensembleObsNames}
\alias{ensembleMemberNames}
\alias{ensembleSize}
\alias{ensembleNobs}
\alias{ensembleForecasts}
\alias{ensembleDates}
\alias{ensembleVerifObs}
\alias{ensembleObsNames.ensembleData}
\alias{ensembleMemberNames.ensembleData}
\alias{ensembleSize.ensembleData}
\alias{ensembleNobs.ensembleData}
\alias{ensembleForecasts.ensembleData}
\alias{ensembleDates.ensembleData}
\alias{ensembleVerifObs.ensembleData}
\alias{[.ensembleData}
\title{
  Creates an ensembleData object.
}
\description{
   Creates a \code{ensembleData} object. 
}
\usage{
ensembleData( forecasts, observations, dates,  ..., na.action = "remove", 
              missingValues = NULL,  labels = NULL)
}
\arguments{
  \item{forecasts}{
     A matrix with columns corresponding to forecasts from individual members 
     of an ensemble and rows corresponding to forecasts for the same date.
  }
  \item{observations}{
    Vector of observed weather conditions corresponding to the forecast
    values.
  }
  \item{dates}{
    Vector of dates corresponding to the forecasts and observations.
    This must be in integer or ordered factor form. If a factor,
    it must sort in chronological order.
   }
  \item{...}{
     A named list of saditional attributes such as latitude,
    longitude, station name, station type.
    }
  \item{na.action}{
     A character string indicating how missing values handled.
     The default is \code{"remove"} means that data items with missing values
     in either forecasts, observations, or dates are removed from the 
     data. If set to \code{NULL}, no missing values are removed.
   }
  \item{missingValues}{
     One or more numeric codes used to indicate missing values in the data.
     These will be replaced by \code{NA} in the output if 
     \code{na.action = NULL}.
     If not supplied, it is assumed that missing values are designated
     by \code{NA}.
   }
  \item{labels}{
     A character vector labeling the observations, forecasts, and dates.
     The default is to take the names from the observations, forecasts,
     or dates in that order.
   }
}
\value{
   An \code{ensembleData} object, incorporating forecasts and observations
   with the associated integer dates.
}
\details{
    Although the input dates may not be consecutive, the dates of an
    \code{ensembleData} object are consecutive integers giving the 
    order of the levels of the \code{dates} input if it were represented
    as a factor.\cr
    Methods for \code{ensembleData} objects include \code{ensembleSize},
    \code{ensembleForecasts}, \code{ensembleVerifObs}, and \code{ensembleDates}. \cr
    Subsetting is possible, but in the case of columns it applies on to
    the ensemble forecasts.
}

\seealso{
  \code{\link{ensembleBMA}}, 
  \code{\link{ensembleBMAgamma0}}, 
  \code{\link{ensembleBMAnormal}}
}
\examples{

 data(slp)

 slpData <- ensembleData(forecasts = slp[c("AVN","GEM","ETA","NGM","NOGAPS")],
                         observations = slp$obs, dates = slp$date)

\dontrun{
 slpFit <- ensembleBMAnormal( slpData)
}
  data(prcp)

  prcpData <- ensembleData( dates = prcp$date, observations = prcp$obs,
                          forecasts = prcp[,c("CENT","AVN","CMCG","ETA",
                                      "GASP","JMA","NGPS","TCWB","UKMO")])

\dontrun{
  prcpFit <- ensembleBMAgamma0(prcpData)
}

}
\keyword{data}
% docclass is function
