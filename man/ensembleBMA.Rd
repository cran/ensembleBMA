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
ensembleBMA( ensembleData, dates = NULL,
             trainingRule = list(length=NA,lag=NA), control = NULL,
             model = NULL, exchangeable = NULL, minCRPS = NULL)
}
\arguments{
  \item{ensembleData}{
     An \code{ensembleData} object including ensemble forecasts, 
     verification observations and dates. 
     Missing values (indicated by \code{NA}) are allowed.
  }
 \item{dates}{
     The dates for which BMA forecasting models are desired. 
     By default, this will include all
     dates consistent with the training rule.
  }
  \item{trainingRule}{
    A list giving the \code{length} and \code{lag} for the training
    period. The length gives the number of time steps (e.g. days)
    in the training period, and the lag gives the number of time steps
    ahead of the most recent date in the training period for which the
    forecast is valid. There is no default.
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
  \item{exchangeable}{
     A numeric or character vector or factor indicating groups of
     ensemble members that are exchangeable (indistinguishable).
     The model fit will have equal weights and parameters within each group.
     The default determines exchangeability from \code{ensembleData}.
  }
 \item{minCRPS}{
     A logical variable indicating whether or not to add a postprocessing
     step after a normal BMA fit to choose the standard deviation so as to
     minimize the CRPS for the training data. This argument is used only
     for normal models, and the default is to not do the CRPS minimization
     for those models because it may require consderably more
     computation time, expecially when there are many ensemble members.
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
   \code{cdf}, \code{quantileForecast},  \code{modelParameters},  
   \code{brierScore}, \code{crps} and \code{mae}.
}
\seealso{
  \code{\link{ensembleData}}, 
  \code{\link{ensembleBMAnormal}}, 
  \code{\link{ensembleBMAgamma0}},
  \code{\link{cdf}},
  \code{\link{quantileForecast}},
  \code{\link{modelParameters}}, 
  \code{\link{brierScore}},
  \code{\link{crps}},
  \code{\link{mae}},
  \code{\link{controlBMAnormal}},
  \code{\link{controlBMAgamma0}},
  \code{\link{trainingControl}}
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
Technical Report No. 516R, Department of Statistics, University of
Washington, May 2008.

C. Fraley, A. E. Raftery, T. Gneiting,
BMA Forecasting with Missing and Exchangeable Ensemble Members,
in preparation.
}

\examples{
  data(slpTest)

  labels <- c("AVN","GEM","ETA","NGM","NOGAPS")
  slpTestData <- ensembleData(forecasts = slpTest[ ,labels],
                         observations = slpTest$obs, dates = slpTest$date)

  slpTestFit <- ensembleBMA( slpTestData, model = "normal",
                             trainingRule = list(length=30,lag=2))

\dontrun{
  data(prcp)

  labels <- c("CENT","AVN","CMCG","ETA","GASP","JMA","NGPS","TCWB","UKMO")
  prcpTestData <- ensembleData( forecasts = prcpTest[ , labels],
                          dates = prcpTest$date, observations = prcpTest$obs)
  prcpTestFit <- ensembleBMA(prcpTestData, model = "gamma0",
                             trainingRule = list(length=30,lag=2))
}
}
\keyword{models}
% docclass is function
