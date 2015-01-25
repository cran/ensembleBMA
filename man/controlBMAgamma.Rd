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
controlBMAgamma(maxIter, tol, power = 1, startupSpeed = NULL, init, 
           optim.control = list(ndeps = rep( sqrt(.Machine$double.eps), 2)))
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
  \item{power}{
    A scalar value giving the power by which the data will be
    transformed to fit the model for mean of the observations.
    The default is not to transform the data.
    The untransformed forecast is used to fit the variance model.
   }
  \item{startupSpeed}{
    A scalar value giving a global value for the anemometer startup speed,
    or the threshold below which a value of 0 is recorded. As this can
    vary from station to station and network to network, it may be
    preferable to include \code{startupSpeed} as part of the 
   \code{ensembleData} object.
   }
  \item{init}{
     An optional list of initial values for variance coefficients
     and weights. The default is to start with the variance coefficients
     equal to 1, and with equal weights for each member of the ensemble.
  }
  \item{optim.control}{
     Control parameters for the optim function used in the M-step of EM.
     The default here is list(ndeps = rep( sqrt(.Machine$double.eps), 2)),
     which assigns a smaller finite-difference step size than the
    \code{optim} default of 1.e-3. To use the default control parameters for 
    \code{optim}, set \code{optim.control=NULL}.
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
J. M. Sloughter, T. Gneiting and A. E. Raftery, 
Probabilistic wind speed forecasting 
using ensembles and Bayesian model averaging,
\emph{Journal of the American Statistical Association, 105:25--35, 2010}.

C. Fraley, A. E. Raftery, T. Gneiting and J. M. Sloughter,
\code{ensembleBMA}: An \code{R} Package for Probabilistic Ensemble Forecasting
using Bayesian Model Averaging,
Technical Report No. 516R, Department of Statistics, University of
Washington, 2007 (revised 2010).
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

\dontrun{ # R check
   winsTestFit1 <- ensembleBMAgamma(winsTestData, trainingDays = 30,
       control = controlBMAgamma(maxIter = 100, tol = 1.e-6,
   startupSpeed =1))
}

# for quick run only; use more training days for forecasting
   winsTestFit1 <- ensembleBMAgamma(winsTestData[1:14,], trainingDays = 5,
      control = controlBMAgamma(maxIter = 100, tol = 1.e-6, startupSpeed = 1))
}
\keyword{models}


% docclass is function

