\name{prcpFit}
\alias{prcpFit}
\title{
  BMA Model Fit to Precipitation Data
}
\description{
 The \code{ensembleBMAgamma0} model fit with a 30 day training period to the 
 precipitation data set from \cr
 \emph{http://www.stat.washington.edu/MURI}, 
 which gives daily daily 48 hour forecasts of 24 hour accumulated 
 precipitation over the US Pacific Northwest region from December 12, 2002
 through March 31, 2005 on a 9 member version of the University of Washington 
 mesoscale 
 ensemble (Grimit and Mass 2002; Eckel and Mass 2005). 
 Precipitation amounts are quantized to hundredths of an inch.
}
\format{
  A list with the following arguments:
\describe{
\item{\code{dateTable}}{
 A named vector in which the names are the dates and the entries are
 the number of observations for each date.
}
\item{\code{trainingRule}}{
 The training rule used to compute the model fits.
}
\item{\code{prob0coefs}}{
 The coefficients in the logistic regression for probability of zero precipitation.
}
\item{\code{biasCoefs}}{
 The coefficients in the linear regression for bias correction.
}
\item{\code{varCoefs}}{
 The variance coefficients of the models.
}
\item{\code{weights}}{
 The BMA weights for the models.
}
\item{\code{power}}{
  An scalar value giving the power by which the forecasts are transformed
  for the BMA fitting.
}
}
}

\references{
E. P. Grimit and C. F. Mass, 
Initial results of a mesoscale short-range ensemble forecasting system over
the Pacific Northwest, 
\emph{Weather and Forecasting 17:192--205, 2002}.

F. A. Eckel and C. F. Mass, 
Effective mesoscale, short-range ensemble forecasting,
\emph{Weather and Forecasting 20:328--350, 2005}.

J. M. Sloughter, A. E. Raftery, T. Gneiting and C. Fraley, 
Probabilistic quantitative precipitation forecasting 
using Bayesian model averaging,
\emph{Monthly Weather Review 135:3309--3320, 2007}.

C. Fraley, A. E. Raftery, T. Gneiting and J. M. Sloughter,
\code{ensembleBMA}: An \code{R} Package for Probabilistic Forecasting
using Ensembles and Bayesian Model Averaging,
Technical Report No. 516R, Department of Statistics, University of
Washington, 2007 (revised 2010).
}

\examples{
\dontrun{ # R check

  data(prcpFit)

  modelParameters(prcpFit, date = "20030113")

  data(prcpGrid) 

  prcpGridData <- ensembleData(forecasts = prcpGrid[,1:9], 
                               latitude = prcpGrid[,"latitude"],
                               longitude = prcpGrid[,"longitude"],
                               forecsatHour = 48,
                               initializationTime = "00")

# probability of precipitation
  1 - cdf( prcpFit, prcpGridData, value = 0)

# probability of precipitation above 0.25 in
  1 - cdf( prcpFit, prcpGridData, date = "20030115", value = 25)
  
}}
\keyword{data}

% docclass is function
