\name{fitBMAnormal}
\alias{fitBMAnormal}
\title{
    BMA mixture of normals fit to a training set
}
\description{
   Fits a Bayesian Model Averaging mixture of normals to a given training set.
}
\usage{
fitBMAnormal( ensembleData, control = controlBMAnormal(),
              exchangeable = NULL) 
}
\arguments{
  \item{ensembleData}{
     An \code{ensembleData} object including ensemble forecasts and 
     verification observations.
     Missing values (indicated by \code{NA}) are allowed. Dates are ignored
     if they are included.  This is the training set for the model.
  }
  \item{control}{
     A list of control values for the fitting functions. The defaults are
     given by the function \code{controlBMAnormal}.
  }
  \item{exchangeable}{
     An optional numeric or character vector or factor indicating groups of
     ensemble members that are exchangeable (indistinguishable).
     The models have equal weights and parameters within each group.
     If supplied, this argument will override any specification of 
     exchangeability in \code{ensembleData}.
  }
}
\value{
A list with the following output components:
  \item{biasCoefs}{
    The fitted bias-correction coefficients.
  }
  \item{sd}{
    The fitted standard deviations for the mixture of normals model
    (equal or varying across components according to the \code{equalVariance}
    setting in the \code{control} input).
  }
  \item{weights}{
    The fitted BMA weights for the normal components for each ensemble member.
  }
  \item{nIter}{
    The number of EM iterations.
  }
}
\details{
   This function fits a BMA model to a training data set. \cr
   It is called by \code{ensembleBMAnormal}, which can produce a sequence
   of fits over a larger data set.  \cr
   Methods available for the output of \code{fitBMAnormal} include:
   \code{cdf}, \code{quantileForecast}, and \code{modelParameters}. 
}
\seealso{
  \code{\link{ensembleData}}, 
  \code{\link{controlBMAnormal}}, 
  \code{\link{ensembleBMAnormal}},
  \code{\link{cdf}},
  \code{\link{quantileForecast}},
  \code{\link{modelParameters}}
}
\references{
A. E. Raftery, T. Gneiting, F. Balabdaoui and M. Polakowski,
Using Bayesian Model Averaging to calibrate forecast ensembles,
\emph{Monthly Weather Review 133:1155--1174, 2005}.

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

  ensNames <- c("gfs","cmcg","eta","gasp","jma","ngps","tcwb","ukmo")

  obs <- paste("T2","obs", sep = ".")
  ens <- paste("T2", ensNames, sep = ".")

  tempTestData <- ensembleData( forecasts = ensBMAtest[,ens],
                                observations = ensBMAtest[,obs],
                                station = ensBMAtest[,"station"],
                                dates = ensBMAtest[,"vdate"],
                                forecastHour = 48,
                                initializationTime = "00")

  tempTrain <- trainingData( tempTestData, trainingDays = 30,
                             date  = "2008010100")

  tempTrainFit <- fitBMAnormal( tempTrain)

## equivalent to
##    tempTrainFit <- fitBMA( tempTrain, model = "normal")

}
\keyword{models}
% docclass is function
