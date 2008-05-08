\name{controlBMAnormal}
\alias{controlBMAnormal}
\title{
  Control parameters for BMA mixtures of normals
}
\description{
   Specifies a list of values controling the Bayesian Model Averaging fit 
   of a mixture of normals to ensemble forecasts.
}
\usage{
controlBMAnormal(maxIter, eps, equalVariance, biasCorrection, start) 
}
\arguments{
  \item{maxIter}{
     An integer specifying an upper limit on the number of iterations
     for fitting the BMA mixture via EM. The default is \code{Inf}, 
     which sets no upper limit on the number of iterations, 
     so that the convergence criterion based on \code{eps} is used.
  }
  \item{eps}{
     A numeric convergence tolerance. The EM fit for the mixture 
     model is terminated when the relative error in successive 
     objective values in the M-step falls below \code{eps}.
     The default is \code{sqrt(.Machine$double.eps)}, 
     which is approximately \code{1.e-8} on IEEE compliant machines.
  }
  \item{equalVariance}{
     A logical value indicating whether or not the variances for
     the mixture components should to be equal.
     The default is to constrain them to be equal.
  }
  \item{biasCorrection}{
   A character string describing the type of bias correction to be used.
  \describe{
   \item{\code{"regression"}}{
    The bias correction term is formed by regression on the forecast values 
    (including an intercept).
   }
   \item{\code{"additive"}}{
      The mean of the difference between observations 
      and forecasts is used for bias correction.
   }
   \item{\code{"none"}}{
      No bias correction.
   }
  }
  }
 \item{start}{
     An optional list of starting values for standard deviations 
     and weights. The default is to start with all standard deviations
     equal to 1, and with equal weights for each member of the ensemble.
  }
}
\value{
  A list whose components are the input arguments and their assigned values. 
}
\seealso{
  \code{\link{ensembleBMAnormal}}, 
  \code{\link{fitBMAnormal}}
}
\references{
A. E. Raftery, T. Gneiting, F. Balabdaoui and M. Polakowski,
Using Bayesian model averaging to calibrate forecast ensembles,
\emph{Monthly Weather Review 133:1155-1174, 2005}.

C. Fraley, A. E. Raftery, T. Gneiting and J. M. Sloughter,
\code{ensembleBMA}: An \code{R} Package for Probabilistic Forecasting
using Ensemble and Bayesian Model Averaging,
Technical Report No. 516, Department of Statistics, University of
Washington, August 2007.
}

\examples{
  data(slpTest)

  memberLabels <- c("AVN","GEM","ETA","NGM","NOGAPS")
  slpTestData <- ensembleData(forecasts = slpTest[ ,memberLabels],
                         observations = slpTest$obs, dates = slpTest$date)

  slpTestFit <- ensembleBMAnormal(slpTestData,
    trainingRule = list(length =30, lag = 2),
    control = controlBMAnormal(maxIter = 100, equalVariance = FALSE,
                               biasCorrection = "additive"))
}
\keyword{models}
% docclass is function
