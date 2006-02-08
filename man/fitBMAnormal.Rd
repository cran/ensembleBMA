\name{fitBMAnormal}
\alias{fitBMAnormal}
\title{
    BMA for a mixture of normals.
}
\description{
   Fits a Baysian Model Averaging mixture of normals 
   to ensemble forecasting data.
}
\usage{
fitBMAnormal( ensembleData, control = controlBMAnormal()) 
}
\arguments{
  \item{ensembleData}{
     An \code{ensembleData} object with forecasts, observations and dates
     for precipitation.
  }
  \item{control}{
     A list of control values for the fitting functions. The defaults are
     given by the function \code{controlBMAnormal}.
  }
}
\value{
A list with the following output components:
  \item{biasCoefs}{
    The fitted bias-correction coefficients.
  }
  \item{sd}{
    The fitted standard deviations for the mixture of normals model.
  }
  \item{weights}{
    The fitted weights for the mixture of normals model.
  }
  \item{nIter}{
    The number of EM iterations.
  }
}
\details{
   This function fits a BMA model to a training data set. \cr
   It is called by \code{ensembleBMAnormal}, which produces a sequence
   of fits over a larger precipitation data set.  \cr
   The following methods are available for the output of \code{fitBMA}:
   \code{gridForecastBMA}, \code{quantileForecastBMA}, and 
   \code{bmaModelParameters}. 
}
\seealso{
  \code{\link{ensembleData}}, 
  \code{\link{controlBMAnormal}}, 
  \code{\link{ensembleBMAnormal}},
  \code{\link{gridForecastBMA}},
  \code{\link{quantileForecastBMA}},
  \code{\link{bmaModelParameters}}
}
\references{
A. E. Raftery, T. Gneiting, F. Balabdaoui and M. Polakowski,
Using Bayesian Model Averaging to Calibrate Forecast Ensembles,
\emph{Monthly Weather Review 133:1155-1174, 2005}.
}

\examples{
  data(slp)

  slpData <- ensembleData(forecasts = slp[c("AVN","GEM","ETA","NGM","NOGAPS")],
                          observations = slp$obs, dates = slp$date)

  DATE <- sort(unique(slpData$dates))[27]
  trainDat <- trainingData( slpData, date = DATE, trainingRule=list(length=25,lag=2))
  slpFit25a <- fitBMA(trainDat, model = "normal")

  D <- as.numeric(slpData$dates) <= 25
  slpFit25b <- fitBMA(slpData[D, ], model = "normal")
}
\keyword{models}
% docclass is function
