\name{crpsANDmae}
\alias{crpsANDmae}
\alias{crpsANDmae.ensembleBMAgamma0}
\alias{crpsANDmae.ensembleBMAnormal}
\title{
   CRPS and MAE via sampling
}
\description{
   Computes the continuous ranked probability score and median absolute error
   for ensemble BMA forecasting models via sampling. Can be used for
   distributions for which the CDF is not available.
}
\usage{
crpsANDmae( object, ensembleData, nSamples = 10000, seed = NULL, 
            popData = NULL, \dots)
}
\arguments{
  \item{object}{
     An ensemble BMA model fit for \code{ensembleData}.
  }
  \item{ensembleData}{
     An \code{ensembleData} object giving including ensemble
     forecasts, observations and dates corresponding to \code{object}.
  }
  \item{nSamples}{
     The number of simulation samples (gamma0 for precipitation model only)
  }
  \item{seed}{
     Argument to \code{set.seed} for random number generation in simulation.
  }
  \item{popData}{
   For \code{ensembleBMAgamma0} and \code{fitBMAgamma0} objects, there 
     is an additional \code{popData} argument for providing predictors
     in the logistic regression for probability of precipitation.
     The value of \code{popData} (if any) used to create \code{object}
     should be supplied here.
  }
  \item{\dots}{
     Included for generic function compatibility.
  }
}
\value{
A matrix whose rows correspond to the \code{crps} and
\code{mae}, respectively, and whose columns correpond
to these values for climatology,  the ensemble,
the BMA model, and simulations from the BMA model.  
}
\details{
   This method is generic, and can be applied to all ensemble BMA forecasting
   models. \cr
   Note the model may have been applied to a transformation of the data,
   but that information is included in the input \code{x}, and
   the output is transformed appropriately.\cr
   For the gamma0 model for precipitation, the CRPS is only available
   through simulation,  because a transformation of the forecast values 
   is being modeled  rather than the forecast values themselves (the 
   default transformation is the cube root). For normal models for 
   temperature and pressure,  simulation is not need and the 
   \code{nSamples} and \code{seed}  arguments are ignored.
}
\seealso{
  \code{\link{ensembleBMA}} 
}
\references{
D. A. Unger, A method  to estimate the continuous ranked probability score,
\emph{Preprints of the Ninth Conference on Probability and Statistics in
Atmospheric Sciences}, Virginia Beach, VA USA, 206-213,
American Meteorological Society.

H. Hersbach, Decomposition of the continuous ranked probability score
for ensemble prediction systems,
\emph{Weather and Forecasting 15, 559-570}.

T. Gneiting and A. E. Raftery, 
Strictly proper scoring rules, prediction and estimation, 
Technical Report No. 463R, Department of Statistics, University of
Washington, November 2006 
(to appear in \emph{Journal of the American Statistical Association}).
}

\examples{
  data(slp)
  
  slpData <- ensembleData(forecasts = slp[c("AVN","GEM","ETA","NGM","NOGAPS")],
                          observations = slp$obs, dates = slp$date)

  slpFit <- ensembleBMAnormal(slpData)

  crpsANDmae( slpFit, slpData)
}
\keyword{models}
% docclass is function
