\name{ensembleBMAgamma0}
\alias{ensembleBMAgamma0}
\title{
    BMA mixture of gammas ensemble forecasting models with a point mass at 0.
}
\description{
   Fits a Bayesian Model Averaging mixture of gammas with a point mass at 0 
   to ensemble forecasting data. Intended for modeling precipitation. 
}
\usage{
ensembleBMAgamma0( ensembleData, dates = NULL, 
                   trainingRule = list(length = 30, lag = 2), 
                   control = controlBMAgamma0(), warmStart = FALSE, popData = NULL)
}
\arguments{
  \item{ensembleData}{
     An \code{ensembleData} object including ensemble forecasts, observations 
     and dates of precipitation. 
  }
  \item{dates}{
     The dates for which modeling is desired. By default, this will include all
     dates consistent with the training rule.
  }
  \item{trainingRule}{
     A list giving the \code{length} and \code{lag} for the training period.
     The default is to use a 30 time step training period for a forecast
     2 time steps ahead of the last time step in the training period.
  }
  \item{control}{
     A list of control values for the fitting functions. The defaults are
     given by the function \code{controlBMAgamma0}.
  }
  \item{warmStart}{
      A logical variable indicating whether or not estimation of models for a sequence
      of dates or time steps should be initialized with the weights from the previous
      date or time step. The default is for the initilaizatiin to be independent of
      the result at the previous time step.
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
  \item{dateTable}{
    The table of observations corresponding to the dates in 
    \code{ensembleData}  in chronological order.
  }
  \item{trainingRule}{
    The training rule specified as input.
  }
  \item{prob0coefs}{
    The fitted coefficients in the model for the point mass at 0 
    for each member of the ensemble at each date.
  }
  \item{biasCoefs}{
    The fitted coefficients in the model for the mean of nonzero observations
    for each member of the ensemble at each date (bias correction).
  }
  \item{varCoefs}{
    The fitted coefficients for the model for the variance of nonzero 
    observations at each date.
    The coefficients are the same for all members of the ensemble.
  }
  \item{weights}{
    The fitted weights for the mixture of gammas model for the
    nonzero observations at each date.
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
   The output is for all of the \code{dates} in \code{ensembleBMA}, so there
   will be missing entries denoted by \code{NA} for dates that are too recent 
   to be forecast with  the training rule. \cr
   The following methods are available for \code{ensembleBMAgamma0} objects:
   \code{gridForecastBMA}, \code{quantileForecastBMA}, 
   \code{bmaModelParameters}, \code{brierSkillScores}, and \code{crpsANDmae}.
}
\seealso{
  \code{\link{ensembleData}}, 
  \code{\link{forecastBMAgamma0}}, 
  \code{\link{controlBMAgamma0}}, 
  \code{\link{fitBMAgamma0}},
  \code{\link{gridForecastBMA}},
  \code{\link{quantileForecastBMA}},
  \code{\link{bmaModelParameters}},
  \code{\link{brierSkillScores}},
  \code{\link{crpsANDmae}},
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
  data(prcp)

  prcpData <- ensembleData( dates = prcp$date, observations = prcp$obs,
                          forecasts = prcp[,c("CENT","AVN","CMCG","ETA",
                                      "GASP","JMA","NGPS","TCWB","UKMO")])

  prcpFit <- ensembleBMAgamma0(prcpData)
}
}
\keyword{models}
% docclass is function
