\name{forecastBMA}
\alias{forecastBMA}
\alias{forecastBMAnormal}
\alias{forecastBMAgamma0}
\title{
    Ensemble BMA forecast for a given date at observation locations.  
}
\description{
   Gives an ensemble BMA forecast for a given date at observation locations.
   The model to be fitted depends on the type of data.
}
\usage{
forecastBMA( date, ensembleData, trainingRule = list(length = 30, lag = 2), 
             control = NULL, quantiles = 0.5, model = NULL, 
             minCRPS = NULL, popData = NULL)
}
\arguments{
  \item{date}{
     An integer or character string givin the date for which the
     forecast is desired.
  }
  \item{ensembleData}{
     An \code{ensembleData} object including ensemble forecasts, observations 
     and dates. Must include the date referenced in the \code{date} argument.
  }
  \item{trainingRule}{
     A list giving the \code{length} and \code{lag} for the training period.
     The default is to use a 30 day training period for a forecast
     2 days ahead of the last day in the training period.
  }
  \item{control}{
     A list of control values for the fitting functions. 
     The default is \code{controlBMAnormal()} for temperature or
     pressure data and \code{controlBMAgamma0()} for precipitation data.
  }
  \item{quantiles}{
     A vector of quantiles for which the forecast is desired.
     The default is 0.5, for the median forecast of the model.
  }
 \item{model}{
      A character string describing the BMA model to be fit. 
      Current choices are \code{"normal"} for temperature or pressure data,
      and \code{"gamma0"} for precipitation data.
      Used only in function \code{forecastBMA}.
  }
  \item{minCRPS}{
     A logical variable indicating whether or not to add a postprocessing
     step after a normal BMA fit to choose the standard deviation so as to
     minimize the CRPS. The default is not to do the CRPS minimization for
     normal data (temperature and pressure), because it can involve
     considerably more computation, particularly when unequal variances
     are specified in \code{control}. This argument applies only to
     normal models.
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
The matrix of quantile forecasts for the given date.
}
\seealso{
  \code{\link{ensembleData}}, 
  \code{\link{controlBMAnormal}}, 
  \code{\link{controlBMAgamma0}},
  \code{\link{gridForecastBMA}},
  \code{\link{quantileForecastBMA}}
}
\references{
A. E. Raftery, T. Gneiting, F. Balabdaoui and M. Polakowski,
Using Bayesian Model Averaging to Calibrate Forecast Ensembles,
\emph{Monthly Weather Review 133:1155-1174, 2005}.

J. M. Sloughter, A. E. Raftery, T. Gneiting and C. Fraley, 
Probabilistic Quantitative Precipitation Forecasting 
using Bayesian Model Averaging,
Technical Report No. 496R, Department of Statistics, University of
Washington, October 2006 (to appear in \emph{Montly Weather Review}).
}

\examples{
  data(slp)

  slpData <- ensembleData(forecasts = slp[c("AVN","GEM","ETA","NGM","NOGAPS")],
                          observations = slp$obs, dates = slp$date)

  unique(slpData$dates)

  slpForc1 <- forecastBMA( "2000063000", slpData, model = "normal",
                               quantiles = c(0.5,0.9))

  slpForc2 <- forecastBMAnormal( "2000063000", slpData, 
                               quantiles = c(0.5,0.9))
}
\keyword{models}
% docclass is function








