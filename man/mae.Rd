\name{MAE}
\alias{MAE}
\alias{MAE.default}
\title{
   Mean Absolute Error
}
\description{
   Computes the mean absolute error (MAE) for ensemble forecasting models. 
}
\usage{
MAE( fit, ensembleData, dates=NULL, \dots)
}
\arguments{
  \item{fit}{
     A model fit to ensemble forecasting data.
  }
  \item{ensembleData}{
     An \code{ensembleData} object that includes ensemble
     forecasts, verification observations and dates.
     Missing values (indicated by \code{NA}) are allowed. \\
     This need not be the data used for the model \code{fit}, 
     although it must include the same ensemble members. \\
     If \code{ensembleData} includes dates,
     they must be consistent with \code{fit} and \code{dates}. 
     If \code{ensembleData} does not include dates, they will
     be inferred from \code{fit} and \code{dates}. 
   }
 \item{dates}{
     The dates for which the CRPS and MAE will be computed.
     These dates must be consistent with \code{fit} and \code{ensembleData}.
     The default is to use all of the dates in \code{fit}.
  }
  \item{\dots}{
     Included for generic function compatibility.
  }
}
\value{
A vector giving the MAE for the  deterministic forecasts associated with
the raw ensemble and for the ensemble forecasting model.
This is the mean absolute difference of the raw ensemble medians and the 
observations,
and the mean absolute difference of the median forecast
and the observations (as in Sloughter et al. 2007). \\
Note that Raftery et al. 2005 uses
the mean absolute difference of the raw ensemble means and the 
observations,
and the mean absolute difference of the BMA predictive mean and the
observations.
}
\details{
   This method is generic, and can be applied to all ensemble forecasting
   models. \cr
   Note the model may have been applied to a power transformation of the data,
   but that information is included in the input \code{fit}, and
   the output is transformed appropriately.
}
\seealso{
  \code{\link{ensembleBMA}} 
}
\references{
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
                                initializationTime = "00")
\dontrun{ # R check
  tempTestFit <- ensembleBMAnormal( tempTestData, trainingDays = 30)
}

  MAE( tempTestFit, tempTestData)
}
\keyword{models}
% docclass is function
