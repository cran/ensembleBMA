\name{EM.normals}
\alias{EM.normals}
\title{
EM algorithm for a BMA mixture of normals
}
\description{
Implements the EM algorithm for a mixture of normals centered at 
(possibly bias-corrected) ensemble-member forecasts.  
}
\usage{
EM.normals(X, Y, eps = 1e-05, maxiter = 1000, start.w = NULL,
start.sigma = NULL, const.var = TRUE, reg.adjust = TRUE, min.CRPS = TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{

  \item{X}{
matrix of ensemble forecasts. This is an n by K matrix, where there are n
observations to be used in the fitting, and K ensemble members
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
Either TRUE or FALSE.  TRUE if each ensemble member should have the same
variance, and FALSE if each member should have its own variance. 
}
  \item{reg.adjust}{
Either TRUE or FALSE.  TRUE if there should be a bias correction to the forecasts,
in the form of a linear regression adjustment, FALSE if there should be no adjustment.  
The regression used is a simple linear regression for each ensemble member of the form Y = A + B*X.  In
this case Y is the new forecast value and X is the original forecast.
}

  \item{min.CRPS}{
Either TRUE or FALSE. If TRUE, the standard deviation returned will be adjusted to minimize the CRPS score
after the EM algorithm has converged. If FALSE, the MLE for standard deviation obtained from the EM
algorithm will be returned without adjustment.

}


}
\details{
  The EM algorithm is used to implement Bayesian Model Averaging (BMA) to
  create a probabilistic forecast as a weighted mixture of normal distributions.
  The EM algorithm runs until either the specified maximum number of iterations has
  been reached, or when the absolute values of the differences in estimated values drops below 
  the specified epsilon. The differences considered for stopping criteria are the differences in weights,
  latent variables, and logarithms of variances.
}
\value{
A list of values is returned.
  \item{loglik}{
maximized log likelihood.
}
  \item{a}{
vector of intercept values in the linear regression for the ensemble forecasts.
}
\item{b}{
vector of slope coefficient from the linear regression for the ensemble
forecasts.
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

\author{ Adrian E. Raftery, J. McLean Sloughter, Michael Polakowski, Fadoua Balabdaoui }



\seealso{
  \code{
  \link{EM.for.date}
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

#create a simulated dataset with equal weights, no bias,
#and standard deviation of 1 in each component
x <- matrix(rnorm(1000,0,2),nrow = 200, ncol = 5)

y.latent <- floor(runif(200,1,6))
y.means <- NULL
for(i in 1:200)
{
  y.means[i] <- x[i,y.latent[i]]
}
y <- rnorm(200,y.means, sd = 1)

#calculate the BMA estimates of the parameters
EMresult <- EM.normals(x, y, reg.adjust=FALSE, min.CRPS=FALSE)

#read in the sea-level pressure data and calculate BMA estimates
#for forecasting on the 35th day in the data set
#(this requires training on observations 43 through 161)
data(slp)
unique.dates <- unique(slp$date)
date.list <- NULL

for(i in 1:length(unique.dates))
{
  date.list[slp$date==unique.dates[i]] <- i
}

X <- cbind(slp$F1,slp$F2,slp$F3,slp$F4,slp$F5)
Y <- slp$Y

#calculate the BMA estimates of the parameters
EMresult <- EM.normals(X = X[43:161,],Y = Y[43:161] )

}

\keyword{file}





