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
  Computes Brier Scores for climatology, raw ensemble, and ensemble
   forecasting models  given observation thresholds.
}
\usage{
brierScore( fit, ensembleData, thresholds, dates = NULL, \dots)
}
\arguments{
  \item{fit}{
     A model fit to ensemble forecasting data.
  }
  \item{ensembleData}{
     An \code{ensembleData} object including ensemble
     forecasts, verification observations and dates. 
     Missing values (indicated by \code{NA}) are allowed. \\
     This need not be the data used for the model \code{fit}, 
     although it must include the same
     ensemble members.\\
     If \code{ensembleData} includes dates,
     they must be consistent with \code{fit} and \code{dates}. 
     If \code{ensembleData} does not include dates, they will
     be inferred from \code{fit} and \code{dates}. 
  }
  \item{thresholds}{
     One or more threshold values for the Brier
     score computations.
  }
  \item{dates}{
     The dates for which the Brier score will be computed.
     These dates must be consistent with \code{fit} and \code{ensembleData}.
     The default is to use all of the dates in \code{fit}.
     The dates are ignored if \code{fit} originates from \code{fitBMA}, 
     which also ignores date information.
  } 
  \item{\dots}{
     Included for generic function compatibility.
  }
 }
}
\value{
 A data frame giving the Brier Scores for climatology 
(empirical distribution of the verifying observations),
ensemble (voting), and ensemble foreacsting models
 for the specified thresholds. \cr
 A logistic Brier score is also given for the \emph{BMAgamma0} model.
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
Technical Report No. 516R, Department of Statistics, University of
Washington, May 2008.
}

\examples{
  data(prcpTest)
                                      
  labels <- c("CENT","AVN","CMCG","ETA","GASP","JMA","NGPS","TCWB","UKMO")
  prcpTestData <- ensembleData( forecasts = prcpTest[ ,labels],
                          dates = prcpTest$date, observations = prcpTest$obs)
\dontrun{
  prcpTestFit <- ensembleBMAgamma0(prcpTestData, trainingRule =
                                   list(length = 30, lag = 2))
}
 
  hist(prcpTestData$obs)

  brierScore(prcpTestFit, prcpTestData, thresholds = c(0, 5, 10, 15, 20))
}
\keyword{models}
% docclass is function