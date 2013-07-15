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
controlBMAnormal(maxIter, tol, equalVariance, biasCorrection, init) 
}
\arguments{
  \item{maxIter}{
     An integer specifying an upper limit on the number of iterations
     for fitting the BMA mixture via EM. The default is \code{Inf}, 
     which sets no upper limit on the number of iterations, 
     so that the convergence criterion based on \code{eps} is used.
  }
  \item{tol}{
     A numeric convergence tolerance. The EM fit for the mixture 
     model is terminated when the relative error in successive 
     objective values in the M-step falls below \code{tol}.
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
 \item{init}{
     An optional list of initial values for standard deviations 
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
\emph{Monthly Weather Review 133:1155--1174, 2005}.

C. Fraley, A. E. Raftery, T. Gneiting and J. M. Sloughter,
\code{ensembleBMA}: An \code{R} Package for Probabilistic Forecasting
using Ensemble and Bayesian Model Averaging,
Technical Report No. 516R, Department of Statistics, University of
Washington, 2007 (revised 2010).
}

\examples{
  data(ensBMAtest)

  ensMemNames <- c("gfs","cmcg","eta","gasp","jma","ngps","tcwb","ukmo")

  obs <- paste("T2","obs", sep = ".")
  ens <- paste("T2", ensMemNames, sep = ".")

  tempTestData <- ensembleData( forecasts = ensBMAtest[,ens],
                                dates = ensBMAtest[,"vdate"],
                                observations = ensBMAtest[,obs],
                                station = ensBMAtest[,"station"],
                                forecastHour = 48,
                                initializationTime = "00")

\dontrun{ # R check
  tempTestFit1 <- ensembleBMAnormal(tempTestData, trainingDays = 30,
    control = controlBMAnormal(maxIter = 100, biasCorrection = "additive"))
}

# for quick run only; use more training days for forecasting
  tempTestFit1 <- ensembleBMAnormal(tempTestData[1:20,], trainingDays = 5,
    control = controlBMAnormal(maxIter = 100, biasCorrection = "additive"))

}
\keyword{models}
% docclass is function
