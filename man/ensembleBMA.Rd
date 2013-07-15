\name{ensembleBMA}
\alias{ensembleBMA}
\alias{getDates}
\alias{getExchangeable}
\alias{matchDates}
\alias{match2side}
\alias{modelDates}
\alias{modelDates.default}
\alias{modelMembers}
\alias{modelMembers.default}
\alias{matchEnsembleMembers}
\alias{matchEnsembleMembers.ensembleBMAgamma}
\alias{matchEnsembleMembers.ensembleBMAgamma0}
\alias{matchEnsembleMembers.ensembleBMAnormal}
\alias{matchEnsembleMembers.fitBMAgamma}
\alias{matchEnsembleMembers.fitBMAgamma0}
\alias{matchEnsembleMembers.fitBMAnormal}
\title{
    BMA mixture model fit
}
\description{
   Fits a BMA mixture model to ensemble forecasts.
   Allows specification of a model, training rule, and forecasting dates.
}
\usage{
ensembleBMA( ensembleData, trainingDays, dates = NULL, control = NULL,
             model = NULL, exchangeable = NULL, minCRPS = NULL)
}
\arguments{
  \item{ensembleData}{
     An \code{ensembleData} object including ensemble forecasts with the 
     corresponding verifying observations and their dates. 
     Missing values (indicated by \code{NA}) are allowed.
  }
 \item{dates}{
     The dates for which BMA forecasting models are desired. 
     By default, this will be all dates in \code{ensembleData}
     for which modeling is allowed given training rule.
  }
  \item{trainingDays}{
    An integer giving the number of time steps (e.g. days)
    in the training period. There is no default.
  }
  \item{control}{
     A list of control values for the fitting functions. 
     The default is \code{controlBMAnormal()} for normal models
     and \code{controlBMAgamma0()} for gamma models with a point mass at 0.
  }
  \item{model}{
      A character string describing the BMA model to be fit. 
      Current choices are \code{"normal"}, typically used for temperature 
      or pressure data, and \code{"gamma0"}, 
      typically used for precipitation data.
  }
  \item{exchangeable}{
     A numeric or character vector or factor indicating groups of
     ensemble members that are exchangeable (indistinguishable).
     The model fit will have equal weights and parameters within each group.
     The default determines exchangeability from \code{ensembleData}.
  }
 \item{minCRPS}{
     A logical variable indicating whether or not to add a postprocessing
     step after a normal BMA fit to choose the standard deviation so as to
     minimize the CRPS for the training data. This argument is used only
     for normal models, and the default is to not do the CRPS minimization
     for those models because it may require consderably more
     computation time, expecially when there are many ensemble members.
  }
}
\value{
A list with the following output components:
  \item{dateTable}{
    The table of observations corresponding to the dates in \code{x}   
     in chronological order.
  }
  \item{trainingDays}{
    The number of days in the training period as specified in input.
  }
  \item{\dots}{
     One or more components corresponding to fitted coefficients for the model.
  }
  \item{weights}{
    The fitted BMA weights for the mixture components for each ensemble
    member at each date.
  }
  \item{power}{
    A scalar value giving the power (if any) by which the data was
    transformed for modeling.
    The untransformed forecast is used to fit the variance model.
    This is input as part of \code{control},
    and applies only to certain models.
   }
}
\details{
   If dates are specified in \code{dates} that cannot be forecast 
   with the training rule, the corresponding BMA model parameter outputs will
   be missing (\code{NA}) but not \code{NULL}. \cr
   The training rule uses the number of days corresponding to its 
   \code{length} regardless of whether or not the dates are consecutive.\cr
   The following methods are available for the output of \code{ensembleBMA}:
   \code{cdf}, \code{quantileForecast},  \code{modelParameters},  
   \code{brierScore}, \code{crps}, \code{CRPS}  and \code{MAE}.
}
\seealso{
  \code{\link{ensembleData}}, 
  \code{\link{ensembleBMAnormal}}, 
  \code{\link{ensembleBMAgamma0}},
  \code{\link{ensembleBMAgamma}},
  \code{\link{cdf}},
  \code{\link{quantileForecast}},
  \code{\link{modelParameters}}, 
  \code{\link{brierScore}},
  \code{\link{crps}},
  \code{\link{MAE}},
  \code{\link{controlBMAnormal}},
  \code{\link{controlBMAgamma0}},
  \code{\link{controlBMAgamma}}
}

\references{
A. E. Raftery, T. Gneiting, F. Balabdaoui and M. Polakowski,
Using Bayesian Model Averaging to calibrate forecast ensembles,
\emph{Monthly Weather Review 133:1155--1174, 2005}.

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

J. M. Sloughter, T. Gneiting and A. E. Raftery, 
Probabilistic wind speed forecasting 
using ensembles and Bayesian model averaging,
\emph{Journal of the American Statistical Association, 105:25--35, 2010}.
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

\dontrun{ # R checl
  tempTestFit <- ensembleBMA( tempTestData, trainingDays = 30,
                              model = "normal")

## equivalent to
##    tempTestFit <- ensembleBMAnormal( tempTestData, trainingDays = 30)
}

# for quick run only; use more training days for forecasting
  tempTestFit <- ensembleBMA( tempTestData[1:20,], trainingDays = 8,
                              model = "normal")
}
\keyword{models}
% docclass is function
