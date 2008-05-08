\name{trainingControl}
\alias{trainingControl}
\title{
    Set Training Rule Argument
}
\description{
  Function for use in setting training rules in ensembleBMA.
}
\usage{
trainingControl(length = 30, lag = NA)
}
\arguments{
  \item{length}{
     The number of time steps (e.g. days)
     in the training period.
    }
  \item{lag}{
    The number of time steps ahead of the most recent date
    in the training period for which the
    forecast is valid. There is no default.
  }
}
\value{
  An list with components specifying \code{length} and \code{lag} that
  can be use to set the \code{trainingRule} argument in 
  \code{ensembleData}.
}
\seealso{
  \code{\link{ensembleBMA}}
}

\examples{
\dontrun{
  data(slpTest)

  memberLabels <- c("AVN","GEM","ETA","NGM","NOGAPS")
  slpTestData <- ensembleData( forecasts = slpTest[ ,memberLabels],
                         observations = slpTest$obs, dates = slpTest$date)

  trainDat <- trainingData(slpTestData, date = "2000063000",
                           trainingRule = trainingControl(length=25,lag=2))
 
  slpTestFitTD <- fitBMAnormal(trainDat)
}
}
\keyword{models}
% docclass is function
