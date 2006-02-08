\name{bmaModelParameters}
\alias{bmaModelParameters}
\alias{bmaModelParameters.ensembleBMAgamma0}
\alias{bmaModelParameters.ensembleBMAnormal}
\alias{bmaModelParameters.fitBMAgamma0}
\alias{bmaModelParameters.fitBMAnormal}
\title{
    Extract BMA model parameters
}
\description{
   Extracts model parameters from an ensemble BMA fit.
}
\usage{
bmaModelParameters( fit, \dots) 
}
\arguments{
\item{fit}{
     An ensemble BMA model fit.
  }
\item{\dots}{
     For \code{ensembleBMAgamma0} and \code{ensembleBMAnormal} fits, there
     is an additional \code{dates} argument, giving a character representation 
     of the dates for which BMA model parameters are desired. 
    \code{dates} must be chosen from \code{names(fit$dateTable)}, and the
     default is to give the BMA model parameters for all of the available 
     dates.
 }
}
\value{
A list of parameters (including weights) corresponding to the BMA model 
for the specified dates. 
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

  slpTestFit1 <- ensembleBMAnormal(slpTestData, date = "2000063000")

  bmaModelParameters( slpTestFit1, date = "2000063000")

  trainDat <- trainingData( slpTestData, date = "2000063000",
                            trainingRule = slpTestFit1$training)
 
  slpTestFit2 <- fitBMAnormal(trainDat)
 
  bmaModelParameters( slpTestFit2)
}
\keyword{models}
% docclass is function
