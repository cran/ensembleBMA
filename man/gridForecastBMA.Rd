\name{gridForecastBMA}
\alias{gridForecastBMA}
\alias{gridForecastBMA.ensembleBMAgamma0}
\alias{gridForecastBMA.ensembleBMAnormal}
\alias{gridForecastBMA.fitBMAgamma0}
\alias{gridForecastBMA.fitBMAnormal}
\title{
   BMA forecasts for a grid of ensemble predictions.
}
\description{
   Computes quantiles for the probability distribution function
   (PDF) of an ensemble BMA mixture model at a grid of 
   locations. 
}
\usage{
gridForecastBMA( object, gridData, quantiles = 0.5, 
                 date = NULL, popData = NULL, \dots)
}
\arguments{
  \item{object}{
     An ensemble BMA model fit.
  }
  \item{gridData}{
     A two- or three dimensional array.
     If two dimensional, the rows represents spatial locations on a grid
     the columns represent the ensemble forecasts at each location.
     If three dimensional, the for two dimension represents spatial locations 
     on a grid and the thrid dimension represents the ensemble forecasts 
     at each location.
  }
  \item{quantiles}{
     The vector of desired quantiles for the PDF of the BMA mixture model.
  }
 \item{date}{
     For \code{ensembleBMAgamma0} and \code{ensembleBMAnormal} objects, there 
     is an additional \code{date} argument, giving a character representation 
     of the date for which the grid forecast is desired. 
    \code{date} must be chosen from \code{names(object$dateTable)}.
 }
 \item{popData}{
     For \code{ensembleBMAgamma0} and \code{fitBMAgamma0} objects, there 
     is an additional \code{popData} argument for providing predictors
     in the logistic regression for probability of precipitation.
     The value of \code{popData} (if any) used to create \code{object}
     should be supplied here.
 }
 \item{\dots}{
     Included for generic function compatibility.
 }
}
\value{
     A three dimensional grid in which the first two dimensions
     represent spatial locations as in the input and the third dimension 
     gives the desired quantile forecasts at each location. 
}
\details{
   This method is generic, and can be applied to any ensemble BMA forecasting
   model. \cr
   Note the model may have been applied to a transformation of the data,
   but that information is included in the input \code{object}, and
   the output is transformed appropriately. \cr
   This can be used to compute confidence intervals for the PDF. 
}
\seealso{
  \code{\link{ensembleBMA}}, 
  \code{\link{fitBMA}},
  \code{\link{quantileForecastBMA}} 
}
\references{
A. E. Raftery, T. Gneiting, F. Balabdaoui and M. Polakowski,
Using Bayesian Model Averaging to Calibrate Forecast Ensembles,
\emph{Monthly Weather Review 133:1155-1174, 2005}.

V. J. Berrocal, A. E. Raftery and T. Gneiting,
Combining Spatial and Ensemble Information in Probabilistic Weather Forecasts,
Technical Report no. 495, Department of Statistics, University of
Washington, February 2006.
}

\examples{
  data(srft)
  data(srftGrid)

  srftData <- ensembleData(dates = srft$date, observations = srft$observation,
  forecasts = srft[c("CMCG","ETA","GASP","GFS","JMA","NGPS","TCWB","UKMO")])

  srftFeb14fit1 <- ensembleBMAnormal( srftData, date = "2004021400",
                                     trainingRule = list(length=25, lag=2))

  gridForecastBMA(srftFeb14fit1, srftGrid[,1:ensembleSize(srftData)], 
                  date = "2004021400")

  trainingFeb14data <- trainingData( srftData, date = "2004021400",
                                 trainingRule = list(length=25, lag=2))

  srftFeb14fit2 <- fitBMAnormal( trainingFeb14data)

  gridForecastBMA(srftFeb14fit2, srftGrid[,1:ensembleSize(srftData)], 
                  date = "2004021400")
}
\keyword{models}
% docclass is function
