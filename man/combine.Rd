\name{combine}
\alias{combine}
\alias{combineCORE}
\alias{combine.ensembleBMAgamma}
\alias{combine.ensembleBMAgamma0}
\alias{combine.ensembleBMAnormal}
\title{
   Combine Compatible BMA Models
}
\description{
  Combines BMA models having the same characteristics for different dates. 
}
\usage{
combine( x, y, ...)
}
\arguments{
   \item{x}{
      An \code{ensembleBMA} model. 
    }
  \item{y}{
    An \code{ensembleBMA} model having the same characteristics as
    \code{x} except for dates.
  }
  \item{\dots}{
      Other \code{ensembleBMA} models compatible with \code{x} and \code{y}.
  }   
}
\value{
  An \code{ensembleBMA} model that merges the models from each input into 
  a single model for the common dates.
 }
\details{
  Input models are checked for compatibility, and entries from different
  inputs having the same dates are eliminated. Dates are ordered
  chronologically and ensemble members are ordered in the order in which
  they occur in inout \code{x}.
}
\seealso{
  \code{\link{ensembleBMA}}
}
\references{
C. Fraley, A. E. Raftery, T. Gneiting and J. M. Sloughter,
\code{ensembleBMA}: An \code{R} Package for Probabilistic Forecasting
using Ensemble and Bayesian Model Averaging,
Technical Report No. 516R, Department of Statistics, University of
Washington, 2007 (revised 2010).
}

\examples{
  data(ensBMAtest)

  ensMemNames <- c("gfs","cmcg","eta","gasp","jma","ngps","tcwb","ukmo")

  obs <- paste("T2","obs", sep = ".")
  ens <- paste("T2", ensMemNames, sep = ".")


  tempTestData <- ensembleData( forecasts = ensBMAtest[,ens],
                                dates = ensBMAtest[,"vdate"],
                                observations = ensBMAtest[,obs],
                                station = ensBMAtest[,"station"],
                                forecastHour = 48,
                                initializationTime = "00")

\dontrun{ # R check
  tempTestFit12 <- ensembleBMAnormal( tempTestData, trainingDays = 30,
                   dates = c("2008010100","2008010200"))

  tempTestFit34 <- ensembleBMAnormal( tempTestData, trainingDays = 30,
                   dates = c("2008010300","2008010400"))
}

# for quick run only; use more training days for forecasting
  tempTestFit12 <- ensembleBMAnormal( tempTestData, trainingDays = 8,
                   dates = c("2008010100","2008010200"))

  tempTestFit34 <- ensembleBMAnormal( tempTestData, trainingDays = 8.,
                   dates = c("2008010300","2008010400"))

  tempTestFit <- combine( tempTestFit12, tempTestFit34)
}
\keyword{models}
% docclass is function
