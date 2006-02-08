\name{crps}
\alias{crps}
\alias{crps.ensembleBMAgamma0}
\alias{crps.ensembleBMAnormal}
\alias{crps.fitBMAgamma0}
\alias{crps.fitBMAnormal}
\title{
   Continuous Ranked Probability Score
}
\description{
   Computes the continuous ranked probability score (CRPS)
   for ensemble BMA forecasts. 
}
\usage{
crps( fit, ensembleData, dates=NULL, nSamples=NULL, seed=NULL, 
            popData=NULL, \dots)
}
\arguments{
  \item{fit}{
     An ensemble BMA model fit for \code{ensembleData}.
  }
  \item{ensembleData}{
     An \code{ensembleData} object including ensemble
     forecasts and observations. It need not be the object
     used to form \code{fit}, although it must include the same
     ensemble members. If \code{ensembleData} includes dates,
     they must be consistent with \code{fit} and \code{dates}. 
     If \code{ensembleData} does not include dates, they will
     be inferred from the \code{fit} and the \code{dates} argument.
  }
 \item{dates}{
     The dates for which the CRPS will be computed.
     These dates must be consistent with \code{fit} and \code{ensembleData}.
     The default is to use all of the dates in \code{fit}.
  }
  \item{nSamples}{
     The number of simulation samples for CRPS via simulation.
     For the normal model, the default is analytic computation of the CRPS.
     For the gamma model with a point mass at 0 (precipitation), 
     the CRPS is always computed by simulation,
     with default \code{nSamples = 10000}.
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
A vector giving the CRPS for the raw ensemble
and the BMA probabilistic forecast.
}
\details{
   This method is generic, and can be applied to all ensemble BMA forecasting
   models. \cr
   Note the model may have been applied to a transformation of the data,
   but that information is included in the input \code{fit}, and
   the output is transformed appropriately.\cr
   For the gamma0 model for precipitation, the CRPS is only available
   through simulation,  because a transformation of the forecast values 
   is being modeled  rather than the forecast values themselves (the 
   default transformation is the cube root). The default number of
   simulation samples in this case is 10,000. \cr
   For normal models for temperature and pressure, analytic computation 
   of the CRPS is the default. CRPS will be computed via simulation for
   normal models only if \code{nSamples} is set to a positive value.
}
\seealso{
  \code{\link{ensembleBMA}} 
}
\references{
E. P. Grimit, T. Gneiting, V. J. Berrocal and N. A. Johnson, 
The continuous ranked probability score for circular variables 
and its application to mesoscale forecast ensemble verification,
\emph{Quarterly Journal of the Royal Meteorological Society 
      132:2925-2642 (2006)}.

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

  crps( slpTestFit, slpTestData)
}
\keyword{models}
% docclass is function
