\name{mae}
\alias{mae}
\alias{mae.ensembleBMAgamma0}
\alias{mae.ensembleBMAnormal}
\alias{mae.fitBMAgamma0}
\alias{mae.fitBMAnormal}
\title{
   Mean Absolute Error
}
\description{
   Computes the mean absolute error (MAE) for ensemble BMA forecasting models. 
}
\usage{
mae( fit, ensembleData, dates=NULL, nSamples=NULL, seed=NULL, popData=NULL, \dots)
}
\arguments{
  \item{fit}{
     An ensemble BMA model fit for \code{ensembleData}.
  }
  \item{ensembleData}{
     An \code{ensembleData} object that includes both ensemble
     forecasts and observations. It need not be the object
     used to form \code{fit}, although it must include the same
     ensemble members. If \code{ensembleData} includes dates,
     they must be consistent with \code{fit} and \code{dates}. 
     If \code{ensembleData} does not include dates, they will
     be inferred from \code{fit} and \code{dates}.
  }
 \item{dates}{
     The dates for which the CRPS and MAE will be computed.
     These dates must be consistent with \code{fit} and \code{ensembleData}.
     The default is to use all of the dates in \code{fit}.
  }
  \item{nSamples}{
     The number of simulation samples for MAE via simulation.
     The default is analytic computation of the MAE.
  }
  \item{seed}{
     Argument to \code{set.seed} for random number generation in simulation.
  }
  \item{popData}{ 
    For \emph{gamma0} model fits, there is an additional \code{popData}
    argument for providing predictors in the logistic regression for 
    probability of zero precipitation. If \code{popData} was supplied to
    obtain in the modeling for \code{fit}, then \code{popData} of the
    same kind must be supplied here.
  }
  \item{\dots}{
     Included for generic function compatibility.
  }
}
\value{
A vector giving the MAE for the deterministic forecasts associated with
the raw ensemble and the BMA model. For the raw ensemble, this is the
mean absolute difference of the the ensemble means and the observations.
For the BMA  mixture of normals model,
it is the mean absolute difference of the BMA predictive mean (Raftery et al. 2005)
and the observations.
For the BMA mixture of gammas with a point mass at 0 model used for
precipitation, it is the mean absolute difference of the BMA median forecast
and the observations (Sloughter et al. 2007).
}
\details{
   This method is generic, and can be applied to all ensemble BMA forecasting
   models. \cr
   Note the model may have been applied to a transformation of the data,
   but that information is included in the input \code{fit}, and
   the output is transformed appropriately.
}
\seealso{
  \code{\link{ensembleBMA}} 
}
\references{
C. Fraley, A. E. Raftery, T. Gneiting and J. M. Sloughter,
\code{ensembleBMA}: An \code{R} Package for Probabilistic Forecasting
using Ensembles and Bayesian Model Averaging,
Technical Report No. 516, Department of Statistics, University of
Washington, August 2007.
}

\examples{
  data(slpTest)
  
  memberLabels <- c("AVN","GEM","ETA","NGM","NOGAPS")
  slpTestData <- ensembleData(forecasts = slpTest[ ,memberLabels],
                         observations = slpTest$obs, dates = slpTest$date)

  slpTestFit <- ensembleBMAnormal(slpTestData)

  mae( slpTestFit, slpTestData)
}
\keyword{models}
% docclass is function
