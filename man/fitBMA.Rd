\name{fitBMA}
\alias{fitBMA}
\title{
    BMA model fit to a training set
}
\description{
   Fits a Bayesian Modeling Averaging mixture model to a given training set.

}
\usage{
fitBMA( ensembleData, control = NULL, model = NULL, 
        exchangeable = NULL, popData = NULL)
}
\arguments{
  \item{ensembleData}{
     An \code{ensembleData} object with ensemble forecasts, observations 
     and dates. This is the training set for the model.
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
 \item{popData}{
      Optional predictors for the logistic regression for probability of 
      zero precipitation. This option applies to the mixture of gammas model 
      with a point mass at zero that is used for precipitation.
      In this model, the default predictors are an intercept, the
      transformed forecast data, and a logical variable indicating
      whether or not the forecast is equal to 0. 
      If provided, the predictors in \code{popData} would replace the logical 
      variable in the regression. \cr
      To supply one pop predictor per ensemble member,
      \code{popData} can be a matrix or data frame (if the
      predictor is categorical) with number of rows equal to the
      number of observations, and number of columns equal to the
      ensemble size. \cr
      To supply multiple numeric pop predictors per ensemble member,
      \code{popData} can be an array of dimension (number of
      observations) by (ensemble size) by (number of pop predictors). \cr
      To supply multiple pop predictors per ensemble member, some of
      which may be categorical, \code{popData} must be a list of
     (number of observations) by (ensemble size) matrices, one
      for each pop predictor.
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
  \item{transformation}{
    The function corresponding to the transformation (if any) of the data used 
    to fit the models for the  point mass at 0 and mean of nonzero 
    observations. 
    The untransformed forecast is used to fit the variance model.
    This is input as part of \code{control}.
   }
  \item{inverseTransformation}{
    The function corresponding to the inverse of \code{transformation} 
   (if any).
    Used for quantile forecasts and verification.
    This is input as part of \code{control}.
   }
}
\details{
   This function fits a BMA model to a training data set. \cr
   Methods available for \code{fitBMA} objects (the output of \code{fitBMA})
   include:   \code{cdfBMA},  \code{quantileForecastBMA}, and 
   \code{bmaModelParameters}. 
}
\seealso{
  \code{\link{ensembleData}}, 
  \code{\link{ensembleBMA}},
  \code{\link{fitBMAgamma0}}, 
  \code{\link{fitBMAnormal}}, 
  \code{\link{cdfBMA}}, 
  \code{\link{quantileForecastBMA}}, 
  \code{\link{bmaModelParameters}}, 
  \code{\link{controlBMAnormal}},
  \code{\link{controlBMAgamma0}}
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
Technical Report No. 516, Department of Statistics, University of
Washington, August 2007.
}

\examples{
\dontrun{
  data(slpTest)

  memberLabels <- c("AVN","GEM","ETA","NGM","NOGAPS")
  slpTestData <- ensembleData(forecasts = slpTest[ ,memberLabels],
                         observations = slpTest$obs, dates = slpTest$date)

  DATE <- sort(unique(slpTestData$dates))[27]
  trainDat <- trainingData( slpTestData, date=DATE, 
                            trainingRule=list(length=25,lag=2))
  slpFitA <- fitBMA(trainDat, model = "normal")

  D <- as.numeric(slpTestData$dates) <= 25
  slpFitB <- fitBMA(slpTestData[D, ], model = "normal")
}
}
\keyword{models}
% docclass is function
