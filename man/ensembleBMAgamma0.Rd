\name{ensembleBMAgamma0}
\alias{ensembleBMAgamma0}
\alias{[.ensembleBMAgamma0}
\title{
    BMA precipitation modeling 
}
\description{
   Fits a Bayesian Model Averaging mixture of gammas with a point mass at 0 
   to ensemble forecasts. Intended for predicting precipitation. 
   Allows specification of a training rule and forecasting dates.
}
\usage{
ensembleBMAgamma0( ensembleData, trainingDays, dates = NULL, 
                   control = controlBMAgamma0(), exchangeable = NULL)
}
\arguments{
  \item{ensembleData}{
      An \code{ensembleData} object including ensemble forecasts with the 
      corresponding verifying observations and their dates. 
      Missing values (indicated by \code{NA}) are allowed.
  }
  \item{trainingDays}{
    An integer giving the number of time steps (e.g. days)
    in the training period. There is no default.
  }
  \item{dates}{
     The dates for which forecasting models are desired. 
     By default, this will be all dates in \code{ensembleData} for which
     modeling is allowed given the training rule.
  }
  \item{control}{
     A list of control values for the fitting functions. The defaults are
     given by the function \code{controlBMAgamma0}.
  }
  \item{exchangeable}{
     A numeric or character vector or factor indicating groups of
     ensemble members that are exchangeable (indistinguishable).
     The models fit will have equal weights and parameters within each group.
     The default determines exchangeability from \code{ensembleData}.
  }
}
\value{
  A list with the following output components:
  \item{training}{
    A list containing information on the training length and lag and
    the number of instances used for training for each modeling day.
  }
  \item{prob0coefs}{
    The fitted coefficients in the model for the point mass at 0 
    (probability of zero precipitaion) for each member of the ensemble 
     at each date. 
  }
  \item{biasCoefs}{
    The fitted coefficients in the model for the mean of the gamma components
    for each member of the ensemble at each date (bias correction).
  }
  \item{varCoefs}{
    The fitted coefficients for the model for the variance of gamma
    components for each date.
    The coefficients are the same for all members of the ensemble.
  }
  \item{weights}{
    The fitted BMA weights for the gamma components for each ensemble member at each date.
  }
  \item{power}{
    A scalar value giving to the power by which the data was transformed
    to fit the models for the  point mass at 0 and the bias model. 
    The untransformed forecast is used to fit the variance model.
    This is input as part of \code{control}.
   }
}
\details{
   The output is for all of the \code{dates} in \code{ensembleBMA}, so there
   will be missing entries denoted by \code{NA} for dates that are too recent 
   to be forecast with  the training rule. \cr
   The following methods are available for \code{ensembleBMAgamma0} objects:
   \code{cdf}, \code{quantileForecast},   \code{modelParameters}, 
   \code{brierScore}, \code{crps}, \code{CRPS} and \code{MAE}.
}
\seealso{
  \code{\link{ensembleData}}, 
  \code{\link{controlBMAgamma0}}, 
  \code{\link{fitBMAgamma0}},
  \code{\link{cdf}},
  \code{\link{quantileForecast}},
  \code{\link{modelParameters}},
  \code{\link{brierScore}},
  \code{\link{crps}},
  \code{\link{MAE}}
}
\references{
J. M. Sloughter, A. E. Raftery, T. Gneiting and C. Fraley, 
Probabilistic quantitative precipitation forecasting 
using Bayesian model averaging,
\emph{Monthly Weather Review 135:3209--3220, 2007}.

C. Fraley, A. E. Raftery, T. Gneiting and J. M. Sloughter,
\code{ensembleBMA}: An \code{R} Package for Probabilistic Forecasting
using Ensembles and Bayesian Model Averaging,
Technical Report No. 516R, Department of Statistics, University of
Washington, 2007 (revised 2010).

C. Fraley, A. E. Raftery, T. Gneiting,
Calibrating Multi-Model Forecast Ensembles
with Exchangeable and Missing Members using Bayesian Model Averaging,
\emph{Monthly Weather Review 138:190--202, 2010}.
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
  prcpTestFit <- ensembleBMAgamma0( prcpTestData, trainingDays = 30)

## equivalent to
##    prcpTestFit <- ensembleBMA( prcpTestData, trainingDays = 30,
##                                model = "gamma0")
}

# for quick run only; use more training days for forecasting
  prcpTestFit <- ensembleBMAgamma0( prcpTestData[3:16,], trainingDays = 6)
}
\keyword{models}
% docclass is function
