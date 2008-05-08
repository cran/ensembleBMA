\name{modelParameters}
\alias{modelParameters}
\alias{modelParameters.ensembleBMAgamma0}
\alias{modelParameters.ensembleBMAnormal}
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
The list may also include a transformation and its inverse
if the data was transformed to obtain the model parameters.
}
\seealso{
  \code{\link{ensembleBMAgamma0}},
  \code{\link{ensembleBMAnormal}}, 
  \code{\link{fitBMAgamma0}},
  \code{\link{fitBMAnormal}} 
}
\examples{
  data(slpTest)

  memberLabels <- c("AVN","GEM","ETA","NGM","NOGAPS")
  slpTestData <- ensembleData(forecasts = slpTest[ ,memberLabels],
                         observations = slpTest$obs, dates = slpTest$date)

  slpTestFit1 <- ensembleBMAnormal(slpTestData, date = "2000063000",
                            trainingRule = list(length=30, lag=2))

  modelParameters( slpTestFit1, date = "2000063000")

  trainDat <- trainingData( slpTestData, date = "2000063000",
                            trainingRule = slpTestFit1$training)
 
  slpTestFit2 <- fitBMAnormal(trainDat)
 
  modelParameters( slpTestFit2)
}
\keyword{models}
% docclass is function
