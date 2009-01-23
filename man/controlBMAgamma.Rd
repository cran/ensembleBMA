\name{controlBMAgamma}
\alias{controlBMAgamma}
\title{
  Control parameters for BMA wind speed modeling
}
\description{
   Specifies a list of values controling the Bayesian Model Averaging fit 
   of a mixture of gammas to ensemble forecasts for wind speed.
}
\usage{
controlBMAgamma(maxIter, tol, nEsteps, power = 1, start)
}
\arguments{
  \item{maxIter}{
     An integer specifying an upper limit on the number of iterations`
     for fitting  the BMA mixture via EM. The default is 
    \code{Inf}, which sets no upper limit on the number of iterations, 
     so that the convergence criterion based on \code{eps} is used.
  }
  \item{tol}{
     A numeric convergence tolerance. The EM fit for the mixture of
     gammas is terminated when the relative error in successive 
     objective values in the M-step falls below \code{tol}.
     The default is \code{sqrt(.Machine$double.eps)}, 
     which is approximately \code{1.e-8} on IEEE compliant machines.
  }
  \item{nEsteps}{
     The number of E-steps per M-step, after five iterations of standard
     EM. The default is 1, but using more E-steps per M-step causes
     the convergence criterion to be satisfied in fewer iterations.
     This is a heuristic to speed things up since the results for multiple
     E-steps often are not much different from the results with one 
     E-step per M-step. 
  }
  \item{power}{
    A scalar value giving the power by which the data will be
    transformed to fit the model for mean of the observations.
    The default is not to transform the data.
    The untransformed forecast is used to fit the variance model.
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
  \code{\link{ensembleBMAgamma}}, 
  \code{\link{fitBMAgamma}}
}

\references{
J. M. Sloughter, T.Gneiting and A. E. Raftery, 
Probabilistic wind speed forecasting 
using ensembles and Bayesian model averaging,
Technical Report No. 544, Department of Statistics,
University of Washington, October 2008.

C. Fraley, A. E. Raftery, T. Gneiting and J. M. Sloughter,
\code{ensembleBMA}: An \code{R} Package for Probabilistic Ensemble Forecasting
using Bayesian Model Averaging,
Technical Report No. 516R, Department of Statistics, University of
Washington, 2009.
}

\examples{
  data(ensBMAtest)

  ensMemNames <- c("gfs","cmcg","eta","gasp","jma","ngps","tcwb","ukmo")
  
  obs <- paste("MAXWSP10","obs", sep = ".")
  ens <- paste("MAXWSP10", ensMemNames, sep = ".")

  winsTestData <- ensembleData( forecasts = ensBMAtest[,ens],
                                dates = ensBMAtest[,"vdate"], 
                                observations = ensBMAtest[,obs],
                                station = ensBMAtest[,"station"],
                                forecastHour = 48,
                                initializationTime = "00")

   winsTestFit1 <- ensembleBMAgamma(winsTestData, trainingDays = 30,
       control = controlBMAgamma(maxIter = 100, tol = 1.e-6))
}
\keyword{models}
% docclass is function

