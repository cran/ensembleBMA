\name{fitBMAgamma}
\alias{fitBMAgamma}
\alias{inverseLogit}
\title{
    BMA wind speed model fit to a training set
}
\description{
   Fits a Bayesian Modeling Averaging mixture of gammas.
   Intended for wind speed forecasts.
}
\usage{
fitBMAgamma( ensembleData, control = controlBMAgamma(), exchangeable = NULL) 
}
\arguments{
  \item{ensembleData}{
     An \code{ensembleData} object including ensemble forecasts and 
     verification observations.
     Missing values (indicated by \code{NA}) are allowed. Dates are ignored
     if they are included. This is the training set for the model.
  }
  \item{control}{
     A list of control values for the fitting functions. The defaults are
     given by the function \code{controlBMAgamma}.
  }
  \item{exchangeable}{
     An optional numeric or character vector or factor indicating groups of
     ensemble members that are exchangeable (indistinguishable).
     The model fit will have equal weights and parameters within each group.
     If supplied, this argument will override any specification of 
     exchangeability in \code{ensembleData}.
  }
}
\value{
A list with the following output components:
  \item{biasCoefs}{
    The fitted coefficients in the model for the mean of nonzero observations
    for each member of the ensemble (used for bias correction).
  }
  \item{varCoefs}{
    The fitted coefficients for the model for the variance of nonzero 
    observations (these are the same for all members of the ensemble).
  }
  \item{weights}{
   The fitted BMA weights for the gamma components for each ensemble member.
  }
  \item{nIter}{
    The number of EM iterations.
  }
  \item{power}{
    A scalar value giving to the power by which the data was transformed
    to fit the models for the  point mass at 0 and the bias model. 
    The untransformed forecast is used to fit the variance model.
    This is input as part of \code{control}.
   }
}
\details{
   This function fits a BMA model to a training data set. \cr
   It is called by \code{ensembleBMAgamma}, which can produce a sequence
   of fits over a larger precipitation data set.  \cr
   Methods available for the output of \code{fitBMA} include:
   \code{cdf}, \code{quantileForecast}, and 
   \code{modelParameters}. 
}
\seealso{
  \code{\link{ensembleData}}, 
  \code{\link{controlBMAgamma}}, 
  \code{\link{ensembleBMAgamma}},
  \code{\link{cdf}},
  \code{\link{quantileForecast}},
  \code{\link{modelParameters}}
}
\references{
J. M. Sloughter, T. Gneiting and A. E. Raftery, 
Probabilistic wind speed forecasting 
using ensembles and Bayesian model averaging,
\emph{Journal of the American Statistical Association, 105:25--35, 2010}.

C. Fraley, A. E. Raftery, T. Gneiting and J. M. Sloughter,
\code{ensembleBMA}: An \code{R} Package for Probabilistic Forecasting
using Ensembles and Bayesian Model Averaging,
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
                                startupSpeed = 1,
                                forecastHour = 48,
                                initializationTime = "00")

\dontrun{ # R check
  winsTrain <- trainingData( winsTestData, trainingDays = 30,
                             date = "2008010100")
}

# for quick run only; use more training days for forecasting
  winsTrain <- trainingData( winsTestData, trainingDays = 10,
                             date = "2008010100")

  winsTrainFit <- fitBMAgamma( winsTrain)

## equivalent to
##    winsTrainFit <- fitBMA( winsTrain, model = "gamma")
}
\keyword{models}
% docclass is function
