\name{modelParameters}
\alias{modelParameters}
\alias{modelParameters.ensembleBMAgamma}
\alias{modelParameters.ensembleBMAgamma0}
\alias{modelParameters.ensembleBMAnormal}
\alias{modelParameters.fitBMAgamma}
\alias{modelParameters.fitBMAgamma0}
\alias{modelParameters.fitBMAnormal}
\title{
    Extract model parameters
}
\description{
  Extracts model parameters for ensemble forecasting models.
}
\usage{
modelParameters( fit, \dots) 
}
\arguments{
\item{fit}{
     A model fit to ensemble forecasting data.
  }
\item{\dots}{
  For ensemble fits involving dates, there is an additional \code{dates}
  argument, giving a character representation  of the dates for which
  model parameters are desired.  In this case \code{dates} must
  correspond to the models in the fit and the default is to give
  the model parameters for all available dates.
 }
}
\value{
  A list of parameters (including weights) corresponding to the
  ensemble forecasting model for the specified dates. 
  The list may also include a power by which the forecasts were
  transformed to obtain the model parameters.
}
\seealso{
  \code{\link{ensembleBMAgamma}},
  \code{\link{ensembleBMAgamma0}},
  \code{\link{ensembleBMAnormal}}, 
  \code{\link{fitBMAgamma}},
  \code{\link{fitBMAgamma0}},
  \code{\link{fitBMAnormal}} 
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

  modelParameters( tempTestFit, date = "2008010100")

  tempTrain <- trainingData( tempTestData, date = "2008010100",
                             trainingDays = tempTestFit$training$days)

  tempTrainFit <- fitBMAnormal( tempTrain)

  modelParameters( tempTrainFit)
}
\keyword{models}
% docclass is function
