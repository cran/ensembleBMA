
\name{bmaquant}
\alias{bmaquant}
\title{ Find a specific quantile of a BMA mixture of normal distributions }
\description{
Returns the quantile of a mixture of weighted normals pdf.  This
can be used to find the bounds of a confidence interval for such a pdf.
The bisection method is used to find the desired quantile.
}
\usage{
bmaquant(a, b, sigma, w, alpha, X, niter = 14)
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
  \item{alpha}{
quantile desired (.05, .95, etc.).  Can be a vector of quantiles.
}
  \item{X}{
vector of ensemble predictions.
}
  \item{niter}{
number of iterations for the bisection method.  Default is 14.
}

}

\value{
the desired quantile.
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
  \link{EM.for.date}
  },
  \code{
  \link{CRPS}
  },
  \code{
  \link{bmacdf}
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

# 95th percentile
bmaquant(a = EMresult$a,b = EMresult$b, sigma = EMresult$sigma,
w =  EMresult$w, alpha = .95, x[1,])

# 5th percentile
bmaquant(a = EMresult$a,b = EMresult$b, sigma = EMresult$sigma,
w =  EMresult$w, alpha = .05 ,x[1,])

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

# 5th and 95th percentiles
bmaquant(a = EMresult$a,b = EMresult$b, sigma = EMresult$sigma,
w =  EMresult$w, alpha = c(.05,.95) ,X[1,])

}


\keyword{file}

