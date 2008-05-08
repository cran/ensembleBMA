\name{ensembleData}
\alias{ensembleData}
\alias{ensembleObsLabels}
\alias{ensembleMemberLabels}
\alias{ensembleSize}
\alias{ensembleNobs}
\alias{ensembleForecasts}
\alias{ensembleDates}
\alias{ensembleVerifObs}
\alias{ensembleGroups}
\alias{ensembleObsLabels.ensembleData}
\alias{ensembleMemberLabels.ensembleData}
\alias{ensembleSize.ensembleData}
\alias{ensembleNobs.ensembleData}
\alias{ensembleForecasts.ensembleData}
\alias{ensembleDates.ensembleData}
\alias{ensembleVerifObs.ensembleData}
\alias{ensembleGroups.ensembleData}
\alias{[.ensembleData}
\alias{matchEnsembleMembers}
\alias{matchEnsembleMembers.ensembleBMAgamma0}
\alias{matchEnsembleMembers.ensembleBMAnormal}
\alias{matchEnsembleMembers.fitBMAgamma0}
\alias{matchEnsembleMembers.fitBMAnormal}
\title{
  Create an ensembleData object
}
\description{
   Creates an \code{ensembleData} object including ensemble forecasts along
   with dates and (optionally) observations. Other descriptive information
   such as latitude, longitude, and station type may be included as well.
}
\usage{
ensembleData( forecasts, caseLabels = NULL, memberLabels = NULL, 
              exchangeable = NULL, dates = NULL, observations = NULL, 
              latitude = NULL, longitude = NULL, ...)
}
\arguments{
  \item{forecasts}{
     A matrix with columns corresponding to forecasts from individual members 
     of an ensemble and rows corresponding to forecasts for the same date.
  }
  \item{caseLabels}{
     An optional set of labels assigned to the forecast instances.
     The default is to use the row names corresponding to \code{forecasts}.
  }
  \item{memberLabels}{
     An optional set of labels assigned to the ensemble members.
     The default is to use the column names corresponding to \code{forecasts}.
     These are important for BMA modeling
     because they are used to match ensemble members
     with BMA  modeling weights and parameters.
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
  \item{dates}{
    A character vector or its factor equivalent specifying the dates corresponding
    to the forecasts.
    The character form of each date must be a string with format
    format YYYYMMDDHH or YYYYMMDD, where YYYY is the year,
    MM the monthn, DD the day, and (optionally) HH the hour.
   }
  \item{observations}{
    Optional vector of observed weather conditions corresponding 
    to the forecasts. Must be supplied if the data is to be
    used for BMA modeling. 
  }
  \item{latitude}{
    Optional vector of latitudes of the forecast locations.
   }
  \item{longitude}{
    Optional vector of longitudes of the forecast locations.
   }
  \item{...}{
    A named list of additional attributes such as station name and station type.
   }
}
\value{
   An \code{ensembleData} object, incorporating forecasts and (optionally) observations
   with the associated dates.
}
\details{
    Methods for \code{ensembleData} objects include \code{ensembleSize},
    \code{ensembleForecasts}, \code{ensembleVerifObs}, and \code{ensembleDates}. \cr
    Subsetting is possible, but in the case of columns it applies only to
    the ensemble forecasts.
}

\seealso{
  \code{\link{ensembleBMA}}, 
  \code{\link{ensembleBMAgamma0}}, 
  \code{\link{ensembleBMAnormal}}
}

\references{
C. Fraley, A. E. Raftery, T. Gneiting and J. M. Sloughter,
\code{ensembleBMA}: An \code{R} Package for Probabilistic Forecasting
using Ensembles and Bayesian Model Averaging,
Technical Report No. 516, Department of Statistics, University of
Washington, August 2007.
}

\examples{

 data(slpTest)

 labels <- c("AVN","GEM","ETA","NGM","NOGAPS")
 slpTestData <- ensembleData(forecasts = slpTest[ ,labels],
                        observations = slpTest$obs, dates = slpTest$date)

 slpTestDataX <- ensembleData(forecasts = slpTest[ ,labels],
                        observations = slpTest$obs, dates = slpTest$date,
         exchangeable = c(AVN=1, GEM=2, ETA=1, NGM=3, NOGAPS=4))

\dontrun{
 slpTestFit <- ensembleBMAnormal( slpTestData,
                                  trainingRule = list(length=25,lag=2))
 slpTestFitX <- ensembleBMAnormal( slpTestDataX,
                                   trainingRule = list(length=25,lag=2))
}
  data(prcpTest)

  labels <- c("CENT","AVN","CMCG","ETA","GASP","JMA","NGPS","TCWB","UKMO")
  prcpTestData <- ensembleData( forecasts = prcpTest[ ,labels],
                          dates = prcpTest$date, observations = prcpTest$obs)
                            

\dontrun{
  prcpTestFit <- ensembleBMAgamma0(prcpTestData,
                                   trainingRule = list(length=30,lag=2))
}

}
\keyword{data}
% docclass is function
