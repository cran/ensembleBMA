\name{brierSkillScores}
\alias{brierSkillScores}
\alias{brierSkillScores.ensembleBMAgamma0}
\alias{brierSkillScores.ensembleBMAnormal}
\title{
   Brier Skill Scores
}
\description{
   Computes ensemble, logistic, and BMA Brier Skill Scores
   given observation thresholds.
}
\usage{
brierSkillScores( object, ensembleData, thresholds, popData = NULL, \dots)
}
\arguments{
  \item{object}{
     An ensemble BMA model fit for \code{ensembleData}.
  }
  \item{ensembleData}{
     An \code{ensembleData} object giving including ensemble
     forecasts, observations, and dates corresponding to \code{object}.
  }
  \item{thresholds}{
     One or more threshold values with respect to which the Brier
     Skill Scores will be computed.
  }
  \item{popData}{
    For \code{ensembleBMAgamma0} objects, there  is an additional \code{popData} 
    argument for providing predictors in the logistic regression for probability 
    of precipitation. The value of \code{popData} (if any) used to create \code{object}
    should be supplied here.
  }
  \item{\dots}{
     Included for generic function compatibility.
  }
 }
}
\value{
 A data frame giving the ensemble (voting), logistic,
 and BMA Brier Skill Scores for the specified thresholds.
}
\details{
There will be a lot of warnings due to logistic fitting near the
extremes.
}
\seealso{
  \code{\link{crpsANDmae}}, 
  \code{\link{ensembleBMA}} 
}
\references{
G. W. Brier,
Verification of forecasts expressed in terms of probability,
\emph{Monthly Weather Review, 78:1-3 (1950)}.

T. Gneiting and A. E. Raftery, 
Strictly proper scoring rules, prediction and estimation, 
Technical Report No. 463R, Department of Statistics, University of
Washington, November 2006 (to appear in \emph{Journal of the 
American Statistical Association}).
}

\examples{
  data(prcp)

  prcpData <- ensembleData( dates = prcp$date, observations = prcp$obs,
                          forecasts = prcp[,c("CENT","AVN","CMCG","ETA",
                                      "GASP","JMA","NGPS","TCWB","UKMO")])
\dontrun{
  prcpFit <- ensembleBMAgamma0(prcpData)
}
 
  hist(prcpData$obs)

  brierSkillScores(prcpFit, prcpData, thresholds = c(0, 10, 50, 100))
}
\keyword{models}
% docclass is function
