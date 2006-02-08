\name{ensembleBMAnormal}
\alias{ensembleBMAnormal}
\alias{crpsNormal}
\title{
    BMA mixture of normals modeling
}
\description{
   Fits a Bayesian Model Averaging mixture of normals to ensemble forecasts.
   Allows specification of a training rule and forecasting dates.
}
\usage{
ensembleBMAnormal(ensembleData, dates = NULL, 
                  trainingRule = list(length=30, lag=2), 
                  control = controlBMAnormal(), warmStart = FALSE, 
                  minCRPS = FALSE, exchangeable = NULL)
}
\arguments{
  \item{ensembleData}{
     An \code{ensembleData} object including ensemble forecasts, observations 
     and dates of precipitation. 
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
     A list of control values for the fitting functions. The defaults are
     given by the function \code{controlBMAnormal}.
  }
  \item{warmStart}{
      A logical variable indicating whether or not estimation of models 
      for a sequence of dates or time steps should be initialized with 
      the weights from the previous date or time step. The default is 
      for the initialization to be independent of the result at the 
      previous time step.
  }
  \item{minCRPS}{
     A logical variable indicating whether or not to add a postprocessing
     step after the BMA fit to choose the standard deviation so as to
     minimize the CRPS for the training data. The default is not to do the CRPS minimization,
     which can add considerable extra cost to the computation, especially
     when unequal variances are specified in \code{control}.
  }
 \item{exchangeable}{
     A numeric or character vector or factor indicating groups of
     ensemble members that are exchangeable (indistinguishable).
     The modeling will have equal weights and parameters within each group.
     The default determines exchangeability from \code{ensembleData}.
  }
}
\value{
A list with the following output components:
  \item{dateTable}{
    The table of observations corresponding to the dates in 
    \code{ensembleData} in chronological order.
  }
  \item{trainingRule}{
    The training rule specified as input.
  }
 \item{biasCoefs}{
    The fitted bias-correction coefficients for each ensemble member at each date.
  }
  \item{sd}{
    The fitted standard deviations for the mixture of normals model
    at each date.
  }
  \item{weights}{
    The fitted BMA weights for the normal components for each ensemble member at each date.
  }
}
\details{
   The output is for all of the \code{dates} in \code{ensembleData}, so there
   will be missing entries denoted by \code{NA} for dates that are too recent 
   to be forecast with  the training rule. \cr
   The following methods are available for \code{ensembleBMAnormal} objects:
   \code{cdfBMA}, \code{quantileForecastBMA}, \code{bmaModelParameters},  
   \code{brierScore}, \code{crps} and \code{mae}.
}
\seealso{
  \code{\link{ensembleData}}, 
  \code{\link{controlBMAnormal}}, 
  \code{\link{fitBMAnormal}},
  \code{\link{cdfBMA}},
  \code{\link{quantileForecastBMA}},
  \code{\link{bmaModelParameters}},
  \code{\link{brierScore}},
  \code{\link{crps}},
  \code{\link{mae}}
}
\references{
A. E. Raftery, T. Gneiting, F. Balabdaoui and M. Polakowski,
Using Bayesian model averaging to calibrate forecast ensembles,
\emph{Monthly Weather Review 133:1155-1174, 2005}.

C. Fraley, A. E. Raftery, T. Gneiting and J. M. Sloughter,
\code{ensembleBMA}: An \code{R} Package for Probabilistic Forecasting
using Ensembles and Bayesian Model Averaging,
Technical Report No. 516, Department of Statistics, University of
Washington, August 2007.
}

\examples{
\dontrun{
  data(slpTest)

  memberLabels <- c("AVN","GEM","ETA","NGM","NOGAPS")
  slpTestData <- ensembleData(forecasts = slpTest[ ,memberLabels],
                         observations = slpTest$obs, dates = slpTest$date)

\dontrun{
  slpTestFit <- ensembleBMAnormal( slpTestData, model = "normal")
}

  slpTestFit <- ensembleBMAnormal( slpTestData)
}
}
\keyword{models}
% docclass is function
