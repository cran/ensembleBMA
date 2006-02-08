\name{slpTest}
\docType{data}
\alias{slpTest}
\title{Sea-Level Pressure Test Data Set}
\description{
This data set gives 48-hour sea-level pressure forecasts for the
US Pacific Northwest from a 5 member ensemble at the University of Washington 
(UW MM5) and the associated observations for five different stations
for the period from April 16, 2000 to June 30, 2000 
at 0 hours GMT, using a 12-km grid resolution. 
Pressure measurements are given in millibars/hectopascales. \cr
Note that forecasts are not available for some of the interim dates.
}
\format{
A data frame with 205 rows and 11 columns: \cr
\code{AVN,GEN,ETA,NGM,NOGAPS} 
forecasts from the 5 members of the ensemble (numeric). \cr
\code{observation} the observed sea-level pressure (numeric). \cr
\code{date}  the date of each forecast/observation, 
		format YYYYMMDDHH (categorical). \cr
\code{latitude}  the latitude of each forecast/observation (numeric). \cr
\code{longitude}  the longitude of each forecast/observation (numeric). \cr
\code{station} weather station identifier (categorical). \cr
\code{type} weather station type (categorical). 
}
\details{
This dataset is a small subset of the data used in Raftery et al. (2005),
provided for the purposes of testing. Typically forecasting would be
performed on much larger datasets.
}
\references{
A. E. Raftery, T. Gneiting, F. Balabdaoui and M. Polakowski,
Using Bayesian Model Averaging to calibrate forecast ensembles,
\emph{Monthly Weather Review 133:1155--1174, 2005}.

E. P. Grimit and C. F. Mass, 
Initial results of a mesoscale short-range ensemble forecasting system over
the Pacific Northwest, 
\emph{Weather and Forecasting 17:192--205, 2002}.
}

\examples{
\dontrun{
  data(slpTest)
 
  memberLabels <- c("AVN","GEM","ETA","NGM","NOGAPS")
  slpTestData <- ensembleData(forecasts = slpTest[ ,memberLabels],
                   observations = slpTest$obs, dates = slpTest$date)

  slpTestBMA <- ensembleBMA( slpTestData, date="2000063000", model="normal")
  slpTestForc <- quantileForecastBMA( slpTestBMA, slpTestData, 
                                      date="2000063000")
}
}
\keyword{datasets}
