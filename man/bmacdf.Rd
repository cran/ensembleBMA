\name{bmacdf}
\alias{bmacdf}
\title{
  Find the BMA cdf at x for the K-component BMA mixture of normals
}
\description{
Find the BMA cdf at a given value Y. 
}
\usage{
bmacdf(a, b, sigma, w, X, Y)
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
vector of K forecasts
}
  \item{Y}{
Observation value at which cdf is required.
}
}
\details{
  The BMA pdf is a mixture of weighted normals.  This function returns the
  cdf at a given quantile x.
}
\value{
value of the BMA cdf evaluated at x
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

#evaluate and plot the BMA cdf for the first observation
index <- seq(-5,5,by=.1)
cdf.at.index=NULL
for(i in 1:length(index))
{
  cdf.at.index[i]=bmacdf(a = EMresult$a,b = EMresult$b, sigma = EMresult$sigma,
  w =  EMresult$w, x[1,],index[i])
}
plot(index, cdf.at.index, type="l")

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

#evaluate and plot the BMA cdf for the first observation
index <- seq(1000,1025,by=.1)
cdf.at.index=NULL
for(i in 1:length(index))
{
  cdf.at.index[i]=bmacdf(a = EMresult$a,b = EMresult$b, sigma = EMresult$sigma,
  w =  EMresult$w, X[1,],index[i])
}
plot(index, cdf.at.index, type="l")

}

\keyword{file}



