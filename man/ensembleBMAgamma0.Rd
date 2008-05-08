\name{ensembleBMAgamma0}
\alias{ensembleBMAgamma0}
\title{
    BMA precipitation modeling 
}
\description{
   Fits a Bayesian Model Averaging mixture of gammas with a point mass at 0 
   to ensemble forecasts. Intended for predicting precipitation. 
   Allows specification of a training rule and forecasting dates.
}
\usage{
ensembleBMAgamma0( ensembleData, dates = NULL, 
                   trainingRule = list(length = NA, lag = NA), 
                   control = controlBMAgamma0(), warmStart = FALSE, 
                   exchangeable = NULL)
}
\arguments{
  \item{ensembleData}{
      An \code{ensembleData} object including ensemble forecasts, 
      verification observations and dates. 
      Missing values (indicated by \code{NA}) are allowed.
  }
  \item{dates}{
     The dates for which forecasting models are desired. 
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
     A list of control values for the fitting functions. The defaults are
     given by the function \code{controlBMAgamma0}.
  }
  \item{warmStart}{
      A logical variable indicating whether or not estimation of models 
      for a sequence of dates or time steps should be initialized with 
      the weights from the previous date or time step. The default is 
      for the initialization to be independent of the result at the 
      previous time step.
  }
  \item{exchangeable}{
     A numeric or character vector or factor indicating groups of
     ensemble members that are exchangeable (indistinguishable).
     The models fit will have equal weights and parameters within each group.
     The default determines exchangeability from \code{ensembleData}.
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
    (probability of zero precipitaion) for each member of the ensemble 
     at each date. 
  }
  \item{biasCoefs}{
    The fitted coefficients in the model for the mean of the gamma components
    for each member of the ensemble at each date (bias correction).
  }
  \item{varCoefs}{
    The fitted coefficients for the model for the variance of gamma
    components for each date.
    The coefficients are the same for all members of the ensemble.
  }
  \item{weights}{
    The fitted BMA weights for the gamma components for each ensemble member at each date.
  }
  \item{transformation}{
    The function corresponding to the transformation of the data used to fit 
    the models for the  point mass at 0 and the bias model. 
    The untransformed forecast is used to fit the variance model.
    This is input as part of \code{control}.
   }
  \item{inverseTransformation}{
    The function corresponding to the inverse of \code{transformation}.
    Used for quantile forecasts and verification.
    This is input as part of \code{control}.
   }
}
\details{
   The output is for all of the \code{dates} in \code{ensembleBMA}, so there
   will be missing entries denoted by \code{NA} for dates that are too recent 
   to be forecast with  the training rule. \cr
   The following methods are available for \code{ensembleBMAgamma0} objects:
   \code{cdf}, \code{quantileForecast},   \code{modelParameters}, 
   \code{brierScore}, \code{crps} and \code{mae}.
}
\seealso{
  \code{\link{ensembleData}}, 
  \code{\link{controlBMAgamma0}}, 
  \code{\link{trainingControl}}, 
  \code{\link{fitBMAgamma0}},
  \code{\link{cdf}},
  \code{\link{quantileForecast}},
  \code{\link{modelParameters}},
  \code{\link{brierScore}},
  \code{\link{crps}},
  \code{\link{mae}}
}
\references{
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
\dontrun{
  data(prcpTest)

  labels <- c("CENT","AVN","CMCG","ETA","GASP","JMA","NGPS","TCWB","UKMO")

  prcpTestData <- ensembleData( forecasts = prcpTest[ ,labels],
                          dates = prcpTest$date, observations = prcpTest$obs)

\dontrun{
  prcpTestFit <- ensembleBMA( prcpTestData, model = "gamma0",
                              trainingRule = list(length=30,lag=2))
}
  prcpTestFit <- ensembleBMAgamma0( prcpTestData,
                              trainingRule = list(length=30,lag=2))
}
}
\keyword{models}
% docclass is function
