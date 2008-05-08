\name{prcpTest}
\alias{prcpTest}
\alias{prcpTestFit}
\title{
  Precipitation Test Data Set
}
\description{
 This is a subset of daily 48 hour forecasts of 24 hour accumulated 
 precipitation over the US Pacific Northwest region in 2003 based 
 on a 9 member version of the University of Washington mesoscale ensemble
(Grimit and Mass 2002; Eckel and Mass 2005). 
Precipitation amounts are quantized to hundredths of an inch. \cr
Note that forecasts are not available for some of the interim dates.
}
\format{
A data frame with 175 rows and 15 columns: \cr
\code{CENT,AVN,CMCG,ETA,GASP,JMA,NGAPS,TCWB,UKMO}
forecasts from the 9 members of the ensemble (numeric).\cr
\code{observation} the observed sea-level pressure (numeric).\cr
\code{date} the date of each forecast/observation,
		format YYYYMMDDHH (categorical). \cr
\code{station} weather station identifier (categorical). \cr
\code{latitude} the latitude of each weather station (numeric).\cr
\code{longitude} the longitude of each weather station (numeric). \cr
\code{elevation} the elevation of each weather station (numeric). 
}

\details{
This dataset is a small subset of the data used in Sloughter et al.(2007),
provided for the purposes of testing. Typically forecasting would be
performed on much larger datasets.
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
\emph{Monthly Weather Review 135:3209--3220, 2007}.
}

\examples{
\dontrun{
  data(prcpTest)

  labels <- c("CENT","AVN","CMCG","ETA","GASP","JMA","NGPS","TCWB","UKMO")
  prcpTestData <- ensembleData( forecasts = prcpTest[,labels],
                          dates = prcpTest$date, observations = prcpTest$obs)

  prcpTestFit <- ensembleBMA(prcpTestData, model = "gamma0",
                             trainingRule = list(length = 30, lag = 2))
}
}
\keyword{data}

% docclass is function
