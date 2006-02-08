\name{crpsNormal}
\alias{crpsNormal}
\title{
   Continuous Ranked Probability Score for Normal Ensemble Data
}
\description{
   Computes the continuous ranked probability score
   for a normal ensemble BMA forecasting model.
}
\usage{
crpsNormal( sd, weights, biasCoefs, ensembleData)
}
\arguments{
  \item{sd}{
     Standard deviation for an ensemble BMA forecasting model. 
     This will be a scalar if the variances are assumed equal 
     across the ensemble, and otherwise a vector.
  }
  \item{weights}{
     Weights for the ensemble BMA forecasting model. 
  }
  \item{biasCoefs}{
     Bias coefficients for the ensemble BMA forecasting model. 
  }
  \item{ensembleData}{
     An \code{ensembleData} object giving including ensemble
     forecasts and observations associated with the ensemble BMA model.
  }
}
\value{
 The mean of the continuous ranked probability scores for
 \code{ensembleData} according to the ensemble BMA model 
  specified by \code{sd}, \code{weights} and \code{biasCoefs}.
}
\details{
 This function is mainly intended for internal use by \code{ensembleBMA}
 and \code{forecastBMA} with the \code{minCRPS} option.
}
\seealso{
  \code{\link{ensembleBMAnormal}},
  \code{\link{forecastBMAnormal}}
}
\references{
D. A. Unger, A method  to estimate the continuous ranked probability score,
\emph{Preprints of the Ninth Conference on Probability and Statistics in
Atmospheric Sciences}, Virginia Beach, VA USA, 206-213,
American Meteorological Society.

H. Hersbach, Decomposition of the continuous ranked probability score
for ensemble prediction systems,
\emph{Weather and Forecasting 15, 559-570}.

T. Gneiting and A. E. Raftery, 
Strictly proper scoring rules, prediction and estimation, 
Technical Report No. 463R, Department of Statistics, University of
Washington, November 2006 (to appear in \emph{Journal of the American 
Statistical Association}).
}

\examples{
  data(slp)

 slpData <- ensembleData(forecasts = slp[c("AVN","GEM","ETA","NGM","NOGAPS")],
                         observations = slp$obs, dates = slp$date)

## default training data for the 32nd date
  trainDat <- trainingData(slpData, date = sort(unique(slpData$dates))[32],
                           trainingRule = list(length=30, lag=2))
  slpFitTD <- fitBMAnormal(trainDat)

  slp32 <- slpData[sort(unique(slpData$dates))[32], ]
  crpsNormal( slpFitTD$sd, slpFitTD$weights, slpFitTD$biasCoefs, slp32)
}
\keyword{models}
% docclass is function
