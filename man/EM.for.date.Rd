\name{EM.for.date}
\alias{EM.for.date}
\title{
  EM.normals wrapper
}
\description{
This is a wrapper function for the EM.normals() function - it uses a list of dates
and a specified training period to extract the appropriate training set,
and to pass that training set to the EM.normals() function.
}
\usage{
EM.for.date(date, date.list, X, Y, eps = 1e-005, maxiter=1000,
start.w=NULL, start.sigma=NULL, const.var = TRUE, num.training.days = 25, lead=2,
reg.adjust = TRUE, min.CRPS = TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{date}{ 
The date for which a forecast is to be made.
}
  \item{date.list}{
The list of dates corresponding to the ensemble predictions and observations.
This vector must have the same length as that of the observations vector.
Each element in date.list corresponds to one observation. If there are ten
obervations taken on a particular day, there will have to be ten entries of
this date in date.list corresponding to that day. In addition, differences between 
entries in the date list must correspond to the number of days between the entries - 
that is, the value for the first day of a month/year must be one higher than the value 
for the last day of the previous month/year.
}
  \item{X}{
Matrix of ensemble members. This is an n by K matrix, where there are n
observations to be used in the fitting, and K ensemble members.
}
  \item{Y}{
n-vector of observations.
}
  \item{eps}{
stopping criterion. 
}
  \item{maxiter}{
The maximum number of EM iterations allowed. 
}
  \item{start.w}{
A vector for the starting values for the EM algorithm for the ensemble weights.  
The vector of weights should sum up to one.  If not specified, the EM algorithm will
start with equal weights for all members.
}
  \item{start.sigma}{
A starting value for the EM algorithm for the ensemble standard deviation.  
If const.var = TRUE, this should be a single number.  If const.var = FALSE, 
it should be a vector with a length equal to the number of ensemble members. 
If not specified, the EM algorithm will start with each standard deviation equal to
the overall standard deviation of the observations.
}
  \item{const.var}{
Either TRUE or FALSE. TRUE if each ensemble member should have the same
variance, and FALSE if each member should have its own variance. Default
is TRUE 
}
  \item{num.training.days}{
The number of training days to be used to train the model. Default is 25

}

  \item{lead}{
The lead time for the forecast, in days. A 48-hour ahead forecast would have a lead time of 2 days. Default is 2.

}
  \item{reg.adjust}{
Either TRUE or FALSE. TRUE if there should be a regression adjustment,
FALSE if there should be no adjustment. The regression used is a simple
linear regression for each ensemble member of the form Y = A + B*X. In
this case Y is the new prediction value and X is the original prediction.

}

  \item{min.CRPS}{
Either TRUE or FALSE. If TRUE, the standard deviation returned will be adjusted to minimize the CRPS score
after the EM algorithm has converged. If FALSE, the MLE for standard deviation obtained from the EM
algorithm will be returned without adjustment.

}

}
\details{
This function extracts the appropriate training set from a data set for a 
given desired forecast day, training period, and lead time, and
invokes EM.normals to construct a probabilistic model for an
ensemble composed of a mixture of normals. See EM.normals for details on the
fitting of the model.
}
\value{
A list of values is returned.
  \item{loglik}{
maximized log likelihood.
}
  \item{a}{
vector of intercept values in the linear regression for the ensemble members.
}
\item{b}{
vector of slope coefficient from the linear regression for the ensemble
members.
}
\item{w}{
vector of weights for the ensemble members.
}
\item{sigma}{
standard deviation of each normal distribution. If const.var = TRUE then one number is
returned. If const.var = FALSE then a vector of standard deviations is returned.
}
\item{z}{
A matrix of latent variables.
}
\item{niter}{
number of iterations of the EM algorithm.
}

}
\references{
Raftery, A. E., T. Gneiting, F. Balabdaoui, & M. Polakowski, "Using Bayesian Model Averaging to calibrate forecast ensembles." Monthlly Weather Review, to appear, 2005.
earlier version available at:
http://www.stat.washington.edu/www/research/reports/2003/tr440.pdf
}

\author{ Adrian E. Raftery, J. McLean Sloughter, Michael Polakowski }


\seealso{
  \code{
  \link{EM.normals}
  },
  \code{
  \link{CRPS}
  },
  \code{
  \link{bmacdf}
  },
  \code{
  \link{bmaquant}
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

}


\keyword{file}




