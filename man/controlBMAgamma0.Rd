\name{controlBMAgamma0}
\alias{controlBMAgamma0}
\title{
  Control parameters for BMA precipitation modeling
}
\description{
   Specifies a list of values controling the Bayesian Model Averaging fit 
   of a mixture of gammas with a point mass at 0 to ensemble forecasts for
   precipitation.
}
\usage{
controlBMAgamma0(maxIter = Inf, tol = sqrt(.Machine$double.eps), 
                 power = (1/3), rainobs = 10, 
                 init = list(varCoefs = NULL, weights = NULL),
             optim.control = list(ndeps = rep( sqrt(.Machine$double.eps), 2)))
}
\arguments{
  \item{maxIter}{
     An integer specifying an upper limit on the number of iterations
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
    transformed to fit the models for the point mass at 0
    and mean of nonzero observations. The default is to use the 1/3
    power of the data.
    The untransformed forecast is used to fit the variance model.
   }
  \item{rainobs}{
    An integer specifying a minimum number of observations with
    nonzero precipitation in the training data. When necessary and
    possible, the training period will be extended backward in
    increments of days to meet the minimum requirement.
    It is not possible to fit the BMA model for precipitation without
    sufficient nonzero observations. The default minimum number is 10.
    It many instances fewer nonzero observations may suffice, but it
    could also be that more are needed to model precipitation in some
    datasets.
  }
  \item{init}{
     An optional list of initial values for variance coefficients
     and weights. The default is to start with the variance coefficients
     equal to 1, and with equal weights for each member of the ensemble.
  }
  \item{optim.control}{
     Control parameters	for the	optim function used in the M-step of EM.
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
Washington, 2007 (revised 2010).
}

\examples{
  data(ensBMAtest)

  ensMemNames <- c("gfs","cmcg","eta","gasp","jma","ngps","tcwb","ukmo")
  
  obs <- paste("PCP24","obs", sep = ".")
  ens <- paste("PCP24", ensMemNames, sep = ".")

  prcpTestData <- ensembleData( forecasts = ensBMAtest[,ens],
                                dates = ensBMAtest[,"vdate"], 
                                observations = ensBMAtest[,obs], 
                                station = ensBMAtest[,"station"],
                                forecastHour = 48,
                                initializationTime = "00")

\dontrun{ # R check
  prcpTestFit1 <- ensembleBMAgamma0( prcpTestData, trainingDays = 30,
       control = controlBMAgamma0(power = (1/4)))
}

# for quick run only; use more training days for forecasting
  prcpTestFit1 <- ensembleBMAgamma0( prcpTestData[1:14,], trainingDays = 6,
       control = controlBMAgamma0(power = (1/4)))
}
\keyword{models}
% docclass is function

