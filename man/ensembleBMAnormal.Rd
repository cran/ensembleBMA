\name{ensembleBMAnormal}
\alias{ensembleBMAnormal}
\title{
    BMA mixture of normals modeling for ensemble forecasts.  
}
\description{
   Fits a Bayesian Model Averaging mixture of normals to ensemble forecasts.
}
\usage{
ensembleBMAnormal(ensembleData, dates = NULL, 
                  trainingRule = list(length=30, lag=2), 
                  control = controlBMAnormal(), warmStart = FALSE, minCRPS = FALSE)
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
     given by the function \code{controlBMAnormal}.
  }
  \item{warmStart}{
      A logical variable indicating whether or not estimation of models for a sequence
      of dates or time steps should be initialized with the weights from the previous
      date or time step. The default is for the initilaizatiin to be independent of
      the result at the previous time step.
  }
  \item{minCRPS}{
     A logical variable indicating whether or not to add a postprocessing
     step after the BMA fit to choose the standard deviation so as to
     minimize the CRPS. The default is not to do the CRPS minimization,
     which can add considerable extra cost to the computation, especially
     when unequal variances are specified in \code{control}.
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
    The fitted bias-correction coefficients at each date.
  }
  \item{sd}{
    The fitted standard deviations for the mixture of normals model
    at each date.
  }
  \item{weights}{
    The fitted weights for the mixture of normals modelat each date.
  }
}
\details{
   The output is for all of the \code{dates} in \code{ensembleData}, so there
   will be missing entries denoted by \code{NA} for dates that are too recent 
   to be forecast with  the training rule. \cr
   The following methods are available for \code{ensembleBMAnormal} objects:
   \code{gridForecastBMA}, \code{quantileForecastBMA}, 
   \code{bmaModelParameters},  \code{brierSkillScores}, and \code{crpsANDmae}.
}
\seealso{
  \code{\link{ensembleData}}, 
  \code{\link{forecastBMAnormal}}, 
  \code{\link{controlBMAnormal}}, 
  \code{\link{fitBMAnormal}},
  \code{\link{gridForecastBMA}},
  \code{\link{quantileForecastBMA}},
  \code{\link{bmaModelParameters}},
  \code{\link{brierSkillScores}},
  \code{\link{crpsANDmae}}
}
\references{
A. E. Raftery, T. Gneiting, F. Balabdaoui and M. Polakowski,
Using Bayesian Model Averaging to Calibrate Forecast Ensembles,
\emph{Monthly Weather Review 133:1155-1174, 2005}.
}

\examples{
\dontrun{
  data(slp)

  slpData <- ensembleData(forecasts = slp[c("AVN","GEM","ETA","NGM","NOGAPS")],
                          observations = slp$obs, dates = slp$date)

  slpFit <- ensembleBMAnormal( slpData, minCRPS = TRUE)
}
}
\keyword{models}
% docclass is function
