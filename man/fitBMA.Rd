\name{fitBMA}
\alias{fitBMA}
\title{
    BMA model fitting.
}
\description{
   Fits a Bayesian Modeling Averaging mixture model
   to ensemble forecasting data.
}
\usage{
fitBMA( ensembleData, control = NULL, model = NULL, popData = NULL)
}
\arguments{
  \item{ensembleData}{
     An \code{ensembleData} object with forecasts, observations and dates
     for precipitation.
  }
  \item{control}{
     A list of control values for the fitting functions. 
     The default is \code{controlBMAnormal()} for normal models
     and \code{controlBMAgamma0()} for gamma models with a point mass at 0.
  }
 \item{model}{
      A character string describing the BMA model to be fit. 
      Current choices are \code{"normal"} for temperature or pressure data,
      and \code{"gamma0"} for precipitation data.
  }
 \item{popData}{
      Optional predictors for the logistic regression for probability of 
      precipitation. This option applies to the mixture of gammas model with 
      a point mass at zero that is used for precipitation.
      In this model, the default predictors are an intercept, the
      transformed forecast data, and a logical variable indicating
      whether or not the forecast is equal to 0. 
      If provided, the predictors in \code{popData} would replace the logical 
      variable in the regression. \cr
      To supply one pop predictor per ensemble member,
      \code{popData} can be a matrix or data frame (if the
      predictor is categorical) with number of rows equal to the
      number of observations, and number of columns equal to the
      ensemble size. \cr
      To supply multiple numeric pop predictors per ensemble member,
      \code{popData} can be an array of dimension (number of
      observations) by (ensemble size) by (number of pop predictors). \cr
      To supply multiple pop predictors per ensemble member, some of
      which may be categorical, \code{popData} must be a list of
     (number of observations) by (ensemble size) matrices, one
      for each pop predictor.
 }
}
\value{
A list with the following output components:
  \item{...}{
    One or more components corresponding to the coeffcients of the model.  
  }
  \item{weights}{
    The fitted weights for the BMA mixture.
  }
  \item{nIter}{
    The number of EM iterations.
  }
  \item{transformation}{
    The function corresponding the transformation (if any) of the data used 
    to fit the models for the  point mass at 0 and mean of nonzero 
    observations. 
    The untransformed forecast is used to fit the variance model.
    This is input as part of \code{control}.
   }
  \item{inverseTransformation}{
    The function corresponding to the inverse of \code{transformation} 
   (if any).
    Used in various diagnostic methods for the output.
   }
}
\details{
   This function fits a BMA model to a training data set. \cr
   The following methods are available for the output of \code{fitBMA}:
   \code{gridForecastBMA}, \code{quantileForecastBMA}, and 
   \code{bmaModelParameters}. 
}
\seealso{
  \code{\link{ensembleData}}, 
  \code{\link{ensembleBMA}},
  \code{\link{fitBMAgamma0}}, 
  \code{\link{fitBMAnormal}}, 
  \code{\link{gridForecastBMA}}, 
  \code{\link{quantileForecastBMA}}, 
  \code{\link{bmaModelParameters}}, 
}
\references{
J. M. Sloughter, A. E. Raftery, T. Gneiting and C. Fraley,
Probabilistic Quantitative Precipitation Forecasting
using Bayesian Model Averaging,
Technical Report No. 496R, Department of Statistics, University of
Washington, October 2006 (to appear in \emph{Montly Weather Review}).
}

\examples{
\dontrun{
  data(slp)

  slpData <- ensembleData(forecasts = slp[c("AVN","GEM","ETA","NGM","NOGAPS")],
                          observations = slp$obs, dates = slp$date)

  DATE <- sort(unique(slpData$dates))[27]
  trainDat <- trainingData( slpData, date=DATE, trainingRule=list(length=25,lag=2))
  slpFit25a <- fitBMA(trainDat, model = "normal")

  D <- as.numeric(slpData$dates) <= 25
  slpFit25b <- fitBMA(slpData[D, ], model = "normal")
}
}
\keyword{models}
% docclass is function
