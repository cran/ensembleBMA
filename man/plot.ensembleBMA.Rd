\name{plot.ensembleBMA}
\alias{plot.ensembleBMA}
\alias{plot.ensembleBMAgamma}
\alias{plot.ensembleBMAgamma0}
\alias{plot.ensembleBMAnormal}
\alias{plot.fitBMA}
\alias{plot.fitBMAgamma}
\alias{plot.fitBMAgamma0}
\alias{plot.fitBMAnormal}
\alias{plotBMAgamma}
\alias{plotBMAgamma0}
\alias{plotBMAnormal}
\title{
   Plot the Predictive Distribution Function for ensemble forcasting models
}
\description{
   Plots the Predictive Distribution Function (PDF)
   of an ensemble forecasting model.
}
\usage{
\method{plot}{ensembleBMAgamma}( x, ensembleData, dates=NULL, ask=TRUE, ...)
\method{plot}{ensembleBMAgamma0}( x, ensembleData, dates=NULL, ask=TRUE, ...)
\method{plot}{ensembleBMAnormal}( x, ensembleData, dates=NULL, ask=TRUE, ...)
\method{plot}{fitBMAgamma}( x, ensembleData, dates=NULL, ...)
\method{plot}{fitBMAgamma0}( x, ensembleData, dates=NULL, ...)
\method{plot}{fitBMAnormal}( x, ensembleData, dates=NULL, ...)
}
\arguments{
   \item{x}{
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
     The dates for which the PDF will be computed.
     These dates must be consistent with \code{fit} and \code{ensembleData}.
     The default is to use all of the dates in \code{fit}.
     The dates are ignored if \code{fit} originates from \code{fitBMA}, 
     which also ignores date information.
  }
 \item{ask}{
     A logical value indicating whether or not the user should be prompted 
     for the next plot.
  }
 \item{\dots}{
     Included for generic function compatibility.
 }
}

\details{
   This method is generic, and can be applied to any ensemble forecasting
   model. \cr
   The colored curves are the weighted PDFs of the ensemble members,
   and the bold curve is the overall PDF. The vertical black line represents
   the median forecast, and the dotted back lines represent the .1 and .9
   quartiles. The vertical orange line is the verifying observation (if
   any).\cr
   Exchangeable members are represented in the plots by the weighted
   group sum rather than by the indivdual weighted PDFs of each member.
}

\references{
A. E. Raftery, T. Gneiting, F. Balabdaoui and M. Polakowski,
Using Bayesian model averaging to calibrate forecast ensembles,
\emph{Monthly Weather Review 133:1155--1174, 2005}.

J. M. Sloughter, A. E. Raftery, T. Gneiting and C. Fraley,
Probabilistic quantitative precipitation forecasting
using Bayesian model averaging,
\emph{Monthly Weather Review 135:3209--3220, 2007}.

J. M. Sloughter, T. Gneiting and A. E. Raftery, 
Probabilistic wind speed forecasting 
using ensembles and Bayesian model averaging,
\emph{Journal of the American Statistical Association, 105:25--35, 2010}.

C. Fraley, A. E. Raftery, T. Gneiting,
Calibrating Multi-Model Forecast Ensembles
with Exchangeable and Missing Members using Bayesian Model Averaging,
\emph{Monthly Weather Review 138:190-202, 2010}.

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

\dontrun{# R check
  tempTestFit <- ensembleBMAnormal( tempTestData, trainingDays = 30)
  plot(tempTestFit, tempTestData)
}

}
\keyword{models}
% docclass is function
