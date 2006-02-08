\name{brierScore}
\alias{brierScore}
\alias{brierScore.ensembleBMAgamma0}
\alias{brierScore.ensembleBMAnormal}
\alias{brierScore.fitBMAgamma0}
\alias{brierScore.fitBMAnormal}
\title{
   Brier Scores
}
\description{
   Computes climatology, ensemble, logistic, and BMA Brier scores
   given observation thresholds.
}
\usage{
brierScore( fit, ensembleData, thresholds, dates = NULL, popData = NULL, 
                \dots)
}
\arguments{
  \item{fit}{
     An ensemble BMA model fit for \code{ensembleData}.
  }
  \item{ensembleData}{
     An \code{ensembleData} object including ensemble
     forecasts and observations. It need not be the object
     used to form \code{fit}, although it must include the same
     ensemble members. If \code{ensembleData} includes dates,
     they must be consistent with \code{fit} and \code{dates}. 
     If \code{ensembleData} does not include dates, they will
     be inferred from the \code{fit} and the \code{dates} argument.
  }
  \item{thresholds}{
     One or more threshold values for the Brier
     score computations.
  }
  \item{dates}{
     The dates for which the Brier score will be computed.
     These dates must be consistent with \code{fit} and \code{ensembleData}.
     The default is to use all of the dates in \code{fit}.
  } 
  \item{popData}{
    For \emph{gamma0} model fits, there is an additional \code{popData}
    argument for providing predictors in the logistic regression for 
    probability of zero precipitation. If \code{popData} was supplied to
    obtain in the modeling for \code{fit}, then \code{popData} of the
    same kind must be supplied here.
  }
  \item{\dots}{
     Included for generic function compatibility.
  }
 }
}
\value{
 A data frame giving the climatology 
(empirical distribution of the verifying observations),
 ensemble (voting), logistic (coefficients determined by logistic regression
 on the training data), and BMA Brier scores for the specified thresholds.
}
\details{
There can be a lot of warnings due to logistic fitting near the
extremes.
}
\seealso{
  \code{\link{ensembleBMA}} 
}
\references{
G. W. Brier,
Verification of forecasts expressed in terms of probability,
\emph{Monthly Weather Review, 78:1--3 (1950)}.

T. Gneiting and A. E. Raftery, 
Strictly proper scoring rules, prediction and estimation, 
\emph{Journal of the American Statistical Association 102:359--378 (2007).}

C. Fraley, A. E. Raftery, T. Gneiting and J. M. Sloughter,
\code{ensembleBMA}: An \code{R} Package for Probabilistic Forecasting
using Ensembles and Bayesian Model Averaging,
Technical Report No. 516, Department of Statistics, University of
Washington, August 2007.
}

\examples{
  data(prcpTest)
                                      
  labels <- c("CENT","AVN","CMCG","ETA","GASP","JMA","NGPS","TCWB","UKMO")
  prcpTestData <- ensembleData( forecasts = prcpTest[ ,labels],
                          dates = prcpTest$date, observations = prcpTest$obs)
\dontrun{
  prcpTestFit <- ensembleBMAgamma0(prcpTestData)
}
 
  hist(prcpTestData$obs)

  brierScore(prcpTestFit, prcpTestData, thresholds = c(0, 5, 10, 15, 20))
}
\keyword{models}
% docclass is function
