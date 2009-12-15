\name{fitBMA}
\alias{fitBMA}
\title{
    BMA model fit to a training set
}
\description{
   Fits a Bayesian Modeling Averaging mixture model to a given training set.

}
\usage{
fitBMA( ensembleData, control = NULL, model = NULL, exchangeable = NULL)
}
\arguments{
  \item{ensembleData}{
     An \code{ensembleData} object including ensemble forecasts and 
     verification observations.
     Missing values (indicated by \code{NA}) are allowed. Dates are ignored
     if they are included.  This is the training set for the model.
  }
  \item{control}{
     A list of control values for the fitting functions. 
     The default is \code{controlBMAnormal()} for normal models
     and \code{controlBMAgamma0()} for gamma models with a point mass at 0.
  }
 \item{model}{
      A character string describing the BMA model to be fit. 
      Current choices are \code{"normal"} for temperature or pressure data,
      and \code{"gamma0"} for precipitation data.
  }
  \item{exchangeable}{
     A numeric or character vector or factor indicating groups of
     ensemble members that are exchangeable (indistinguishable).
     The model fit will have equal weights and parameters 
     within each group.
     The default determines exchangeability from \code{ensembleData}.
  }
}
\value{
A list with the following output components:
  \item{...}{
    One or more components corresponding to the coeffcients of the model.  
  }
  \item{weights}{
    The fitted BMA weights for the mixture components for each ensemble member.
  }
  \item{nIter}{
    The number of EM iterations.
  }
  \item{power}{
    A scalar value giving the power (if any) by which the data was transformed
    for modeling.
    The untransformed forecast is used to fit the variance model.
    This is input as part of \code{control},
    and applies only to certain models.
   }
}
\details{
   This function fits a BMA model to a training data set. \cr
   Methods available for \code{fitBMA} objects (the output of \code{fitBMA})
   include: \code{cdf}, \code{quantileForecast}, and
   \code{modelParameters}.
}
\seealso{
  \code{\link{ensembleData}}, 
  \code{\link{ensembleBMA}},
  \code{\link{fitBMAgamma}}, 
  \code{\link{fitBMAgamma0}}, 
  \code{\link{fitBMAnormal}}, 
  \code{\link{cdf}}, 
  \code{\link{quantileForecast}}, 
  \code{\link{modelParameters}}, 
  \code{\link{controlBMAgamma}},
  \code{\link{controlBMAgamma0}},
  \code{\link{controlBMAnormal}}
}
\references{
A. E. Raftery, T. Gneiting, F. Balabdaoui and M. Polakowski,
Using Bayesian model averaging to calibrate forecast ensembles,
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

  tempTrainFit <- fitBMA( tempTrain, model = "normal")

## equivalent to
##    tempTrainFit <- fitBMAnormal( tempTrain)

}
\keyword{models}
% docclass is function
