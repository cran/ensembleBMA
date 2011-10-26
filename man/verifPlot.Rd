\name{verifPlot}
\alias{verifPlot}
\title{
   Plot observations along with median, 10th and 90th percentile forecasts.  
}
\description{
   Computes the median, 10th and 90th percentile forecasts,
   and plots the corresponding observations.
}
\usage{
  verifPlot( fit, ensembleData, dates = NULL)
}
\arguments{
   \item{fit}{
     A model fit to ensemble forecasting data.
  }
  \item{ensembleData}{
     An \code{ensembleData} object that includes ensemble forecasts,
     verification observations and possibly dates.
     Missing values (indicated by \code{NA}) are allowed. \\
     This need not be the data used for the model \code{fit},
     although it must include the same ensemble members.
  }
   \item{dates}{
     The dates for which the CDF will be computed.
     These dates must be consistent with \code{fit} and
     \code{ensembleData}.
     The default is to use all of the dates in \code{fit}.
     The dates are ignored if \code{fit} originates from \code{fitBMA},
     which also ignores date information.
  }
}
\value{
  A matrix giving the median, 10th and 90th percentile forecasts for the 
  ensemble data at the specified dates. 
  If observations are available, they are plotted along with the forecasts
  in order of increasing 90th percentile forecast. 
 }
\seealso{
  \code{\link{ensembleData}}, 
  \code{\link{pit}} 
}
\references{
A. E. Raftery, T. Gneiting, F. Balabdaoui and M. Polakowski,
Using Bayesian model averaging to calibrate forecast ensembles,
\emph{Monthly Weather Review 133:1155-1174, 2005}.

T. Gneiting, F. Balabdaoui and A. Raftery,
Probabilistic forecasts, calibration and sharpness.
\emph{Journal of the Royal Statistical Society, Series B 69:243--268, 2007}.

J. M. Sloughter, A. E. Raftery, T. Gneiting and C. Fraley,
Probabilistic quantitative precipitation forecasting
using Bayesian model averaging,
\emph{Monthly Weather Review 135:3209--3220, 2007}.

C. Fraley, A. E. Raftery, T. Gneiting and J. M. Sloughter,
\code{ensembleBMA}: An \code{R} Package for Probabilistic Forecasting
using Ensemble and Bayesian Model Averaging,
Technical Report No. 516R, Department of Statistics, University of
Washington, 2007 (revised 2010).

C. Fraley, A. E. Raftery, T. Gneiting,
Calibrating Multi-Model Forecast Ensembles
with Exchangeable and Missing Members using Bayesian Model Averaging,
\emph{Monthly Weather Review 138:190--202, 2010}.

J. M. Sloughter, T. Gneiting and A. E. Raftery, 
Probabilistic wind speed forecasting 
using ensembles and Bayesian model averaging,
\emph{Journal of the American Statistical Association, 105:25--35, 2010}.
}

\examples{
  data(prcpFit)
  data(prcpDJdata)

  forc <- verifPlot( prcpFit, prcpDJdata, date = "20030113")
}
\keyword{models}
% docclass is function
