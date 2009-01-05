\name{ensembleBMAnormal}
\alias{ensembleBMAnormal}
\alias{[.ensembleBMAnormal}
\alias{crpsNormal}
\title{
    BMA mixture of normals modeling
}
\description{
   Fits a Bayesian Model Averaging mixture of normals to ensemble forecasts.
   Allows specification of a training rule and forecasting dates.
}
\usage{
ensembleBMAnormal(ensembleData, trainingDays, dates = NULL, 
                  control = controlBMAnormal(), exchangeable = NULL,
                  minCRPS = FALSE)
}
\arguments{
  \item{ensembleData}{
     An \code{ensembleData} object including ensemble forecasts, 
     verification observations and dates. 
     Missing values (indicated by \code{NA}) are allowed.
  }
  \item{trainingDays}{
    An integer giving the number of time steps (e.g. days)
    in the training period. There is no default.
  }
  \item{dates}{
     The dates for which BMA forecasting models are desired. 
     By default, this will include all
     dates consistent with the training rule.
  }
  \item{control}{
     A list of control values for the fitting functions. The defaults are
     given by the function \code{controlBMAnormal}.
  }
 \item{exchangeable}{
     A numeric or character vector or factor indicating groups of
     ensemble members that are exchangeable (indistinguishable).
     The modeling will have equal weights and parameters within each group.
     The default determines exchangeability from \code{ensembleData}.
  }
  \item{minCRPS}{
     A logical variable indicating whether or not to add a postprocessing
     step after the BMA fit to choose the standard deviation so as to
     minimize the CRPS for the training data. The default is not to do
     the CRPS minimization, because it can add considerable extra cost
     to the computation, particularly when there are many ensemble members.
  }
}
\value{
 A list with the following output components:
  \item{training}{
    A list containing information on the training length and lag and
    the number of instances used for training for each modeling day.
  }
 \item{biasCoefs}{
    The fitted bias-correction coefficients for each ensemble member at each date.
  }
  \item{sd}{
    The fitted standard deviations for the mixture of normals model
    at each date.
  }
  \item{weights}{
    The fitted BMA weights for the normal components for each ensemble member at each date.
  }
}
\details{
   The output is for all of the \code{dates} in \code{ensembleData}, so there
   will be missing entries denoted by \code{NA} for dates that are too recent 
   to be forecast with  the training rule. \cr
   The following methods are available for \code{ensembleBMAnormal} objects:
   \code{cdf}, \code{quantileForecast}, \code{modelParameters},  
   \code{brierScore}, \code{crps} and \code{mae}.
}
\seealso{
  \code{\link{ensembleData}}, 
  \code{\link{controlBMAnormal}}, 
  \code{\link{fitBMAnormal}},
  \code{\link{cdf}},
  \code{\link{quantileForecast}},
  \code{\link{modelParameters}},
  \code{\link{brierScore}},
  \code{\link{crps}},
  \code{\link{mae}}
}
\references{
A. E. Raftery, T. Gneiting, F. Balabdaoui and M. Polakowski,
Using Bayesian model averaging to calibrate forecast ensembles,
\emph{Monthly Weather Review 133:1155-1174, 2005}.

C. Fraley, A. E. Raftery, T. Gneiting and J. M. Sloughter,
\code{ensembleBMA}: An \code{R} Package for Probabilistic Forecasting
using Ensembles and Bayesian Model Averaging,
Technical Report No. 516R, Department of Statistics, University of
Washington, 2009.

C. Fraley, A. E. Raftery, T. Gneiting,
Using Bayesian Model Averaging to Calibrate Forecast Ensembles
with Missing and Exchangeable Ensemble Members,
(in preparation).
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
                                initializationTime = "00Z")

  tempTestFit <- ensembleBMAnormal( tempTestData, trainingDays = 30)

## equivalent to
##   tempTestFit <- ensembleBMA( tempTestData, trainingDays = 30,
##                                model = "normal")
}
\keyword{models}
% docclass is function
