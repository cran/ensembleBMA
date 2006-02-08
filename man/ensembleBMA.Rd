\name{ensembleBMA}
\alias{ensembleBMA}
\alias{getExchangeable}
\title{
    BMA mixture model fit
}
\description{
   Fits a BMA mixture model to ensemble forecasts.
   Allows specification of a model, training rule, and forecasting dates.
}
\usage{
ensembleBMA( ensembleData, dates = NULL, trainingRule = list(length=30,lag=2), 
             control = NULL, model = NULL, warmStart = FALSE, minCRPS = NULL, 
             exchangeable = NULL, popData = NULL)
}
\arguments{
  \item{ensembleData}{
     An \code{ensembleData}, \code{temperatureData}, \code{pressureData}, 
     or \code{precipitationData} object including ensemble forecasts, 
     observations and dates. 
  }
 \item{dates}{
     The dates for which forecasts are desired. 
     By default, this will include all
     dates consistent with the training rule.
  }
  \item{trainingRule}{
     A list giving the \code{length} and \code{lag} for the training period.
     The default is to use a 30 time step training period for a forecast
     2 time steps ahead of the last time step in the training period.
  }
  \item{control}{
     A list of control values for the fitting functions. 
     The default is \code{controlBMAnormal()} for normal models
     and \code{controlBMAgamma0()} for gamma models with a point mass at 0.
  }
  \item{model}{
      A character string describing the BMA model to be fit. 
      Current choices are \code{"normal"}, typically used for temperature 
      or pressure data, and \code{"gamma0"}, 
      typically used for precipitation data.
  }
  \item{warmStart}{
      A logical variable indicating whether or not estimation of models for 
      a sequence of dates or time steps should be initialized with the weights 
      from the previous date or time step. 
      The default is for the initialization to be independent of
      the result at the previous time step.
  }
 \item{minCRPS}{
     A logical variable indicating whether or not to add a postprocessing
     step after a normal BMA fit to choose the standard deviation so as to
     minimize the CRPS for the training data. The default is not to do the CRPS minimization for
     normal data (e.g. temperature or pressure), because it can involve
     considerably more computation, particularly when unequal variances
     are specified in \code{control}. This argument is used only for
     normal models.
  }
  \item{exchangeable}{
     A numeric or character vector or factor indicating groups of
     ensemble members that are exchangeable (indistinguishable).
     The model fit will have equal weights and parameters within each group.
     The default determines exchangeability from \code{ensembleData}.
  }
  \item{popData}{
      Optional predictors for the logistic regression for probability of 
      zero precipitation. This option applies to the mixture of gammas model 
      with a point mass at zero that is used for precipitation.
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
    The table of observations corresponding to the dates in \code{x}   
     in chronological order.
  }
  \item{trainingRule}{
    The training rule specified as input.
  }
  \item{\dots}{
     One or more components corresponding to fitted coefficients for the model.
  }
  \item{weights}{
    The fitted BMA weights for the mixture components for each ensemble member at each date.
  }
  \item{transformation}{
    A function corresponding to the transformation (if any) of the data used to fit 
    the models for the  point mass at 0 and mean of nonzero observations. 
    The untransformed forecast is used to fit the variance model.
    This is input as part of \code{control}.
   }
  \item{inverseTransformation}{
    A function corresponding to the inverse of \code{transformation} (if any).
    Used for quantile forecasts and verification.
    This is input as part of \code{control}.
   }
}
\details{
   If dates are specified in \code{dates} that cannot be forecast 
   with the training rule, the corresponding BMA model parameter outputs will
   be missing (\code{NA}) but not \code{NULL}. \cr
   The training rule uses the number of days corresponding to its 
   \code{length} regardless of whether or not the dates are consecutive.\cr
   The following methods are available for the output of \code{ensembleBMA}:
   \code{cdfBMA}, \code{quantileForecastBMA},  \code{bmaModelParameters},  
   \code{brierScore}, \code{crps} and \code{mae}.
}
\seealso{
  \code{\link{ensembleData}}, 
  \code{\link{ensembleBMAnormal}}, 
  \code{\link{ensembleBMAgamma0}},
  \code{\link{cdfBMA}},
  \code{\link{quantileForecastBMA}},
  \code{\link{bmaModelParameters}}, 
  \code{\link{brierScore}},
  \code{\link{crps}},
  \code{\link{mae}},
  \code{\link{controlBMAnormal}},
  \code{\link{controlBMAgamma0}}
}

\references{
A. E. Raftery, T. Gneiting, F. Balabdaoui and M. Polakowski,
Using Bayesian Model Averaging to calibrate forecast ensembles,
\emph{Monthly Weather Review 133:1155--1174, 2005}.

J. M. Sloughter, A. E. Raftery, T. Gneiting and C. Fraley, 
Probabilistic quantitative precipitation forecasting 
using Bayesian model averaging,
\emph{Monthly Weather Review 135:3209--3220, 2007}.

C. Fraley, A. E. Raftery, T. Gneiting and J. M. Sloughter,
\code{ensembleBMA}: An \code{R} Package for Probabilistic Forecasting
using Ensembles and Bayesian Model Averaging,
Technical Report No. 516, Department of Statistics, University of
Washington, August 2007.
}

\examples{
  data(slpTest)

  labels <- c("AVN","GEM","ETA","NGM","NOGAPS")
  slpTestData <- ensembleData(forecasts = slpTest[ ,labels],
                         observations = slpTest$obs, dates = slpTest$date)

  slpTestFit <- ensembleBMA( slpTestData, model = "normal")

\dontrun{
  data(prcp)

  labels <- c("CENT","AVN","CMCG","ETA","GASP","JMA","NGPS","TCWB","UKMO")
  prcpTestData <- ensembleData( forecasts = prcpTest[ , labels],
                          dates = prcpTest$date, observations = prcpTest$obs)
  prcpTestFit <- ensembleBMA(prcpTestData, model = "gamma0")
}
}
\keyword{models}
% docclass is function
