\name{trainingData}
\alias{trainingData}
\title{
    Extract Training Data
}
\description{
  Extracts a subset of an \code{ensembleData} object corresponding
  to a given date and number of  training days.
}
\usage{
trainingData( ensembleData, trainingDays, date) 
}
\arguments{
  \item{ensembleData}{
     An \code{ensembleData} object that includes, ensemble
     forecasts, observations and dates.
    }
  \item{trainingDays}{
    An integer specifying the number of days in the training period.
    
  }
  \item{date}{
     The date for which the training data is desired. 
  }
}
\value{
An \code{ensembleData} object corresponding to the training data for
the given date relative to \code{ensembleData}.
}
\details{
  The most recent days are used for training
  regardless of whether or not they are consecutive.
}
\seealso{
  \code{\link{ensembleBMA}}, 
  \code{\link{fitBMA}}
}

\references{
A. E. Raftery, T. Gneiting, F. Balabdaoui and M. Polakowski,
Using Bayesian model averaging to calibrate forecast ensembles,
\emph{Monthly Weather Review 133:1155-1174, 2005}.

J. M. Sloughter, A. E. Raftery, T. Gneiting and C. Fraley, 
Probabilistic quantitative precipitation forecasting 
using Bayesian model averaging,
\emph{Monthly Weather Review 135:3309--3320, 2007}.

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

  ensNames <- c("gfs","cmcg","eta","gasp","jma","ngps","tcwb","ukmo")

  obs <- paste("T2","obs", sep = ".")
  ens <- paste("T2", ensNames, sep = ".")


  tempTestData <- ensembleData( forecasts = ensBMAtest[,ens],
                                observations = ensBMAtest[,obs],
                                station = ensBMAtest[,"station"],
                                dates = ensBMAtest[,"vdate"],
                                forecastHour = 48,
                                initializationTime = "00")

  tempTrain <- trainingData( tempTestData, trainingDays = 30,
                             date  = "2008010100")
 
  tempTrainFit <- fitBMAnormal( tempTrain)

}
\keyword{models}
% docclass is function
