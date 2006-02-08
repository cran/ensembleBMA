\name{fitBMAgamma0}
\alias{fitBMAgamma0}
\title{
    BMA for a mixture of gammas with a point mass at 0.
}
\description{
   Fits a Bayesian Modeling Averaging mixture of gammas 
   with a point mass at 0 to ensemble forecasting data.
   Intended for modeling precipitation.
}
\usage{
fitBMAgamma0( ensembleData, control = controlBMAgamma0(), popData = NULL) 
}
\arguments{
  \item{ensembleData}{
     An \code{ensembleData} object with forecasts, observations and dates
     for precipitation.
  }
  \item{control}{
     A list of control values for the fitting functions. The defaults are
     given by the function \code{controlBMAgamma0}.
  }
  \item{popData}{
      Optional predictors for the logistic regression for probability of 
      precipitation. This option applies to the mixture of gammas model with 
      a point mass at zero that is used for precipitation.
      In this model, the default predictors are an intercept, the
      transformed forecast data, and a logical variable indicating
      whether or not the forecast is equal to 0. 
      In addition, the coefficient of the transformed forecast must
      not be greater than 0 and that of the logical variable
      must not be less than 0 in the default.
      If provided, the predictors in \code{popData} would replace the logical
      variable in the regression, and no constraints are imposed on
      the regression coefficients. \cr    
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
  \item{prob0coefs}{
    The fitted coefficients in the model for the point mass at 0 
   for each member of the ensemble.
  }
  \item{biasCoefs}{
    The fitted coefficients in the model for the mean of nonzero observations
    for each member of the ensemble (used for bias correction).
  }
  \item{varCoefs}{
    The fitted coefficients for the model for the variance of nonzero 
    observations (these are the same for all members of the ensemble).
  }
  \item{weights}{
    The fitted weights for the mixture of gammas model for the
    nonzero observations.
  }
  \item{nIter}{
    The number of EM iterations.
  }
  \item{transformation}{
    The function corresponding the transformation of the data used to fit 
    the models for the  point mass at 0 and mean of nonzero observations. 
    The untransformed forecast is used to fit the variance model.
    This is input as part of \code{control}.
   }
  \item{inverseTransformation}{
    The function corresponding to the inverse of \code{transformation}.
    Used in various diagnostic methods for the output.
   }
}
\details{
   This function fits a BMA model to a training data set. \cr
   It is called by \code{ensembleBMAgamma0}, which produces a sequence
   of fits over a larger precipitation data set.  \cr
 The following methods are available for the output of \code{fitBMA}:
   \code{gridForecastBMA}, \code{quantileForecastBMA}, and 
   \code{bmaModelParameters}. 
}
\seealso{
  \code{\link{ensembleData}}, 
  \code{\link{controlBMAgamma0}}, 
  \code{\link{ensembleBMAgamma0}},
  \code{\link{gridForecastBMA}},
  \code{\link{quantileForecastBMA}},
  \code{\link{bmaModelParameters}}
}
\references{
J. M. Sloughter, A. E. Raftery, T. Gneiting and C. Fraley,
Probabilistic Quantitative Precipitation Forecasting
using Bayesian Model Averaging,
Technical Report No. 496R, Department of Statistics, University of
Washington, October 2006 (to appear in \emph{Montly Weather Review}).
}

\examples{
  data(prcp)

  prcpData <- ensembleData( dates = prcp$date, observations = prcp$obs,
                          forecasts = prcp[,c("CENT","AVN","CMCG","ETA",
                                      "GASP","JMA","NGPS","TCWB","UKMO")])

  DATE <- sort(unique(prcpData$dates))[27]
  trainDat <- trainingData(prcpData, date = DATE,
                           trainingRule = list(length=25,lag=2))
  prcpFit25a <- fitBMAgamma0(trainDat)

  D <- as.numeric(prcpData$dates) <= 25
  prcpFit25b <- fitBMAgamma0(prcpData[D, ])

  prcpFit25c <- fitBMAgamma0(prcpData[D, ], 
                     popData = ensembleForecasts(prcpData[D,]) == 0)
}
\keyword{models}
% docclass is function
