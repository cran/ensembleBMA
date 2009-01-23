\name{verifRank}
\alias{verifRank}
\title{
   Verification Rank
}
\description{
  Computes the rank of verifying observations.
  relative to the corresponding ensemble forecasts.
}
\usage{
verifRank( forecasts, observations)
}
\arguments{
   \item{forecasts}{
     A matrix of ensemble forecasts, in which the rows corresponds to
     locations and times and the columns correspond to the individual
     ensemble members.
  }
  \item{observations}{
    A vector of observations corresponding to the locations and times
    of the forecasts.
  }
}
\value{
  A vector giving the rank of verifying observations relative
  to the correspnding ensemble forecasts.
 }
\details{
  The verification rank is used to assess calibration of a forecast
  ensemble. A more uniform verification rank histogram
  indicates better calibartion.
}
\seealso{
  \code{\link{ensembleData}}, 
  \code{\link{pit}} 
}
\references{
A. E. Raftery, T. Gneiting, F. Balabdaoui and M. Polakowski,
Using Bayesian model averaging to calibrate forecast ensembles,
\emph{Monthly Weather Review 133:1155-1174, 2005}.

J. M. Sloughter, A. E. Raftery, T. Gneiting and C. Fraley,
Probabilistic quantitative precipitation forecasting
using Bayesian model averaging,
\emph{Monthly Weather Review 135:3209--3220, 2007}.

J. M. Sloughter, T.Gneiting and A. E. Raftery, 
Probabilistic wind speed forecasting 
using ensembles and Bayesian model averaging,
Technical Report No. 544, Department of Statistics,
University of Washington, October 2008. 

C. Fraley, A. E. Raftery, T. Gneiting and J. M. Sloughter,
\code{ensembleBMA}: An \code{R} Package for Probabilistic Forecasting
using Ensemble and Bayesian Model Averaging,
Technical Report No. 516R, Department of Statistics, University of
Washington, revised 2009.
}

\examples{
  data(prcpDJdata)

  prcpDJdataVR <- verifRank( ensembleForecasts(prcpDJdata),
                               ensembleVerifObs(prcpDJdata))

  hist(prcpDJdataVR, breaks = ensembleSize(prcpDJdata), xlab = "",
       main = "Verification Rank Histogram")

}
\keyword{models}
% docclass is function
