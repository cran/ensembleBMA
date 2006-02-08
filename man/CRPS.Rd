\name{CRPS}
\alias{CRPS}
\title{
  Calculate the Continuous Ranked Probability Score for an ensemble forecast.
}
\description{
  Calculate the Continuous Ranked Probability Score for an ensemble forecast.
}
\usage{
CRPS(a,b,sigma,w,X,Y)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{a}{
vector of K intercepts in the regression bias correction.  If
no regression is desired, 'a' should be a vector of zeros.
}
  \item{b}{
vector of K slopes in the regression bias correction.  If
no regression is desired, 'b' should be a vector of ones.
}
  \item{sigma}{
vector of K standard deviations from the BMA fit
(a,b,sigma are all outputs of EM.normals or EM.for.date).  If there
is only one variance parameter (constant variance), then this can
be a single number.
}
  \item{w}{
vector of K weights from the BMA fit
}
  \item{X}{
matrix of ensemble forecasts. This is an n by K matrix, where there are n
observations to be used in the fitting, and K ensemble members
}
  \item{Y}{
n-vector of observations.
}

}
\details{
  This function calculates the Continuous Ranked Probability Score (CRPS) for an ensemble BMA forecast with the specified forecasts, observations, and BMA parameters
}
\value{
value of the CRPS
}

\references{
Raftery, A. E., T. Gneiting, F. Balabdaoui, & M. Polakowski, "Using Bayesian Model Averaging to calibrate forecast ensembles." Monthly Weather Review, to appear, 2005.
earlier version available at:
http://www.stat.washington.edu/www/research/reports/2003/tr440.pdf

Gneiting, T, and A. E. Raftery, "Strictly proper scoring rules, prediction, and estimation." University of Washington Technical Report 463.
available at:
http://www.stat.washington.edu/www/research/reports/2004/tr463.pdf

Gneiting, T., A. Westveld, A. E. Raftery, and T. Goldman, " Calibrated Probabilistic Forecasting Using Ensemble Model Output Statistics and Minimum CRPS Estimation." Monthly Weather Review, to appear, 2005.
earlier version available at:
http://www.stat.washington.edu/www/research/reports/2004/tr449.pdf
}

\author{ Adrian E. Raftery, J. McLean Sloughter, Michael Polakowski }


\seealso{
  \code{
  \link{EM.normals}
  },
  \code{
  \link{EM.for.date}
  },
  \code{
  \link{bmaquant}
  },
  \code{
  \link{bmacdf}
  }
}

\examples{
#read in the sea-level pressure data and calculate BMA estimates
#for forecasting on the 35th day in the data set
data(slp)
unique.dates <- unique(slp$date)
date.list <- NULL

for(i in 1:length(unique.dates))
{
  date.list[slp$date==unique.dates[i]] <- i
}

X <- cbind(slp$F1,slp$F2,slp$F3,slp$F4,slp$F5)
Y <- slp$Y

EMresult <- EM.for.date(date = 35,date.list = date.list,X = X,Y = Y )

#now calculate the CRPS over the training period (observations 43 through 161)
CRPS(a = EMresult$a,b = EMresult$b, sigma = EMresult$sigma, w =  EMresult$w, X=X[43:161,], Y=Y[43:161])

#calculate the BMA estimates without CRPS minimization, and compare the new CRPS score
EMresult.without <- EM.for.date(date = 35,date.list = date.list,X = X,Y = Y, min.CRPS=FALSE )
CRPS(a = EMresult.without$a,b = EMresult.without$b, sigma = EMresult.without$sigma, w =  EMresult.without$w, X=X[43:161,], Y=Y[43:161])

}

\keyword{file}



