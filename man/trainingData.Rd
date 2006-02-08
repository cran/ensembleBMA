\name{trainingData}
\alias{trainingData}
\title{
    Extract Training Data
}
\description{
  Extracts the training data corresponding to a given date and training rule.
}
\usage{
trainingData( ensembleData, date, trainingRule = list(length=30, lag=2)) 
}
\arguments{
  \item{ensembleData}{
     An \code{ensembleData} object including ensemble forecasts, observations 
     and dates. 
  }
  \item{date}{
     The date for which the training data is desired. 
  }
  \item{trainingRule}{
     A list giving the \code{length} and \code{lag} for the training period.
     The default is to use a 30 day training period for a forecast
     2 days ahead of the last day in the training period.
  }
}
\value{
An \code{ensembleData} object corresponding to the training data for
the given date relative to \code{ensembleData}.
}
\details{
  The training rule uses the most recent days for the given period
  regardless of whether or not they are consecutive.
}
\seealso{
  \code{\link{ensembleBMA}}, 
  \code{\link{fitBMA}} 
}

\references{
A. E. Raftery, T. Gneiting, F. Balabdaoui and M. Polakowski,
Using Bayesian model averaging to calibrate forecast ensembles,
\emph{Monthly Weather Review 133:1155-1174, 2005}.

J. M. Sloughter, A. E. Raftery, T. Gneiting and C. Fraley, 
Probabilistic quantitative precipitation forecasting 
using Bayesian model averaging,
\emph{Monthly Weather Review 135:3309--3320, 2007}.

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
  slpTestData <- ensembleData( forecasts = slpTest[ ,memberLabels],
                         observations = slpTest$obs, dates = slpTest$date)

  trainDat <- trainingData(slpTestData, date = "2000063000",
                           trainingRule = list(length=30,lag =2))
 
  slpTestFitTD <- fitBMAnormal(trainDat)
}
}
\keyword{models}
% docclass is function
