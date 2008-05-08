

\name{controlBMAgamma0}
\alias{controlBMAgamma0}
\title{
  Control parameters for BMA precipitation modeling
}
\description{
   Specifies a list of values controling the Bayesian Model Averaging fit 
   of a mixture gammas with a point mass at 0 to ensemble forecasts for
   precipitation.
}
\usage{
controlBMAgamma0(maxIter, eps, nEsteps, transformation, inverseTransformation, 
                start)
}
\arguments{
  \item{maxIter}{
     An integer specifying an upper limit on the number of iterations
     for fitting  the BMA mixture via EM. The default is 
    \code{Inf}, which sets no upper limit on the number of iterations, 
     so that the convergence criterion based on \code{eps} is used.
  }
  \item{eps}{
     A numeric convergence tolerance. The EM fit for the mixture of
     gammas is terminated when the relative error in successive 
     objective values in the M-step falls below \code{eps}.
     The default is \code{sqrt(.Machine$double.eps)}, 
     which is approximately \code{1.e-8} on IEEE compliant machines.
  }
  \item{nEsteps}{
     The number of E-steps per M-step. The default is 1, but using
     more E-steps per M-step will cause the convergence criterion
     to be satisfied in fewer iterations. This is a heuristic for
     a gain in speed since the results for multiple M-steps
     often are  not much different from the results with one 
     E-step per M-step.
  }
  \item{transformation}{
    The function corresponding to the transformation of the data used to fit 
    the models for the  point mass at 0 and mean of nonzero observations. 
    The untransformed forecast is used to fit the variance model.
   }
  \item{inverseTransformation}{
    The function corresponding to the inverse of \code{transformation}.
    Used for quantile forecasts and verification.
   }
  \item{start}{
     An optional list of starting values for variance coefficients
     and weights. The default is to start with the variance coefficients
     equal to 1, and with equal weights for each member of the ensemble.
  }
}

\value{
  A list whose components are the input arguments and their assigned values. 
}

\seealso{
  \code{\link{ensembleBMAgamma0}}, 
  \code{\link{fitBMAgamma0}}
}

\references{
J. M. Sloughter, A. E. Raftery, T Gneiting and C. Fraley, 
Probabilistic quantitative precipitation forecasting 
using Bayesian model averaging,
\emph{Monthly Weather Review 135:3209--3220, 2007}.

C. Fraley, A. E. Raftery, T. Gneiting and J. M. Sloughter,
\code{ensembleBMA}: An \code{R} Package for Probabilistic Ensemble Forecasting
using Bayesian Model Averaging,
Technical Report No. 516R, Department of Statistics, University of
Washington, May 2008.
}

\examples{
  data(prcpTest)

  labels <- c("CENT","AVN","CMCG","ETA","GASP","JMA","NGPS","TCWB","UKMO")
  prcpTestData <- ensembleData( forecasts = prcpTest[ , labels],
                          dates = prcpTest$date, observations = prcpTest$obs)

  prcpTestFitAlt <- ensembleBMAgamma0(prcpTestData, 
                         trainingRule=list(length=25,lag=2),
       control = controlBMAgamma0(transformation = function(x) x^(1/4),
                              inverseTransformation = function(x) x^4))
}
\keyword{models}
% docclass is function
