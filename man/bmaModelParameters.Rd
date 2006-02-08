\name{bmaModelParameters}
\alias{bmaModelParameters}
\alias{bmaModelParameters.ensembleBMAgamma0}
\alias{bmaModelParameters.ensembleBMAnormal}
\alias{bmaModelParameters.fitBMAgamma0}
\alias{bmaModelParameters.fitBMAnormal}
\title{
    Extract ensemble BMA model parameters.
}
\description{
   Extracts model parameters from an ensemble BMA fit.
}
\usage{
bmaModelParameters( object, \dots) 
}
\arguments{
\item{object}{
     An ensemble BMA model fit.
  }
\item{\dots}{
     For \code{ensembleBMAgamma0} and \code{ensembleBMAnormal} objects, there 
     is an additional \code{dates} argument, giving a character representation 
     of the dates for which BMA model parameters are desired. 
    \code{dates} must be chosen from \code{names(object$dateTable)}, and the
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
  data(slp)

  slpData <- ensembleData(forecasts = slp[c("AVN","GEM","ETA","NGM","NOGAPS")],
                          observations = slp$obs, dates = slp$date)

  slpFit1 <- ensembleBMAnormal(slpData, date = "2000063000")

  bmaModelParameters( slpFit1, date = "2000063000")

  trainDat <- trainingData( slpData, date = "2000063000",
                            trainingRule = slpFit1$training)
 
  slpFit2 <- fitBMAnormal(trainDat)
 
  bmaModelParameters( slpFit2)
}
\keyword{models}
% docclass is function
